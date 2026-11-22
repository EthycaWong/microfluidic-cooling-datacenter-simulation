#include "system/Server.hpp"
#include "Units.hpp"

#include <algorithm>
#include <cmath>

namespace mf {

Server::Server(const ServerParams& p, std::shared_ptr<const Fluid> fluid, int index)
    : p_(p), fluid_(fluid), index_(index),
      gpu_{ChipGrid([&] { auto c = p.gpu; c.name = "GPU1"; return c; }(), fluid),
           ChipGrid([&] { auto c = p.gpu; c.name = "GPU2"; return c; }(), fluid)},
      cpu_(p.cpu, fluid), air_(p.air) {
    setActivity(0.0);
}

std::string Server::branchName(int b) const { return b == 0 ? "GPU1" : b == 1 ? "GPU2" : "CPU"; }

void Server::attach(HydraulicNetwork& net, int supplyNode, int returnNode, double z) {
    const FluidProps ref = referenceWater(p_.balanceT); // component curves are specified for water
    const std::string tag = "S" + std::to_string(index_ + 1) + ".";
    nodeS_ = net.addNode(tag + "in", z);
    nodeR_ = net.addNode(tag + "out", z);
    qdIn_ = net.addEdge(supplyNode, nodeS_,
                        std::make_unique<RefLoss>(tag + "QD supply", p_.qdDpRefPa, lpmToM3s(p_.qdQRefLpm), 0.0, ref));
    qdOut_ = net.addEdge(nodeR_, returnNode,
                         std::make_unique<RefLoss>(tag + "QD return", p_.qdDpRefPa, lpmToM3s(p_.qdQRefLpm), 0.0, ref));
    for (int b = 0; b < kBranches; ++b) {
        auto s = std::make_unique<SeriesEdge>(tag + branchName(b));
        s->add(std::make_unique<PipeEdge>(tag + branchName(b) + " hose", p_.hoseDiameter, p_.hoseLength, p_.hoseK));
        if (b < 2) {
            s->add(std::make_unique<ChannelEdge>(tag + branchName(b) + " microchannels", p_.gpu.channels));
        } else {
            const ColdPlateParams& c = p_.cpu.coldPlate;
            s->add(std::make_unique<RefLoss>(tag + "CPU cold plate", c.dpRefPa, lpmToM3s(c.qRefLpm), c.laminarFraction, ref));
        }
        auto valve = std::make_unique<RefLoss>(tag + branchName(b) + " balancing valve", valveDp_[b],
                                               lpmToM3s(p_.branchTargetLpm), 0.0, ref);
        valve_[b] = static_cast<RefLoss*>(s->add(std::move(valve)));
        branchEdge_[b] = net.addEdge(nodeS_, nodeR_, std::move(s));
    }
}

void Server::clearBalance() {
    for (int b = 0; b < kBranches; ++b) { valveDp_[b] = 0.0; if (valve_[b]) valve_[b]->setDpRef(0.0); }
}

void Server::autoBalance(HydraulicNetwork& net) {
    clearBalance();
    if (!p_.autoBalance) return;
    const FluidProps ref = fluid_->props(p_.balanceT);
    const double m = ref.rho * lpmToM3s(p_.branchTargetLpm);
    std::array<double, kBranches> dp{};
    double mx = 0.0;
    for (int b = 0; b < kBranches; ++b) {
        auto& e = net.edge(branchEdge_[b]);
        e.setTemperature(p_.balanceT);
        dp[b] = e.dp(m, net.fluid());
        mx = std::max(mx, dp[b]);
    }
    for (int b = 0; b < kBranches; ++b) {
        // The valve law is linear in dpRef: size it so that it adds exactly (mx - dp) with this coolant.
        valve_[b]->setDpRef(1.0);
        const double unit = valve_[b]->dp(m, net.fluid());
        valveDp_[b] = unit > 0 ? (mx - dp[b]) / unit : 0.0;
        valve_[b]->setDpRef(valveDp_[b]);
    }
}

void Server::updateEdgeTemperatures(HydraulicNetwork& net, double supplyT) {
    net.edge(qdIn_).setTemperature(supplyT);
    net.edge(qdOut_).setTemperature(returnT_);
    for (int b = 0; b < kBranches; ++b) net.edge(branchEdge_[b]).setTemperature(chip(b).meanFluidT());
}

void Server::pullFlows(const HydraulicNetwork& net) {
    for (int b = 0; b < kBranches; ++b) {
        branchMdot_[b] = std::max(0.0, net.flow(branchEdge_[b]));
        chip(b).setMassFlow(branchMdot_[b]);
    }
}

void Server::setActivity(double u) {
    u_ = std::clamp(u, 0.0, 1.0);
    for (int b = 0; b < kBranches; ++b) chip(b).setActivity(u_);
    air_.setActivity(u_);
}

void Server::updateReturnT(double supplyT) {
    double q = 0.0;
    for (int b = 0; b < kBranches; ++b) q += chip(b).ledger().toCoolant;
    const double m = totalFlow();
    const double cp = fluid_->props(0.5 * (supplyT + returnT_)).cp;
    returnT_ = m > 0 ? supplyT + q / (m * cp) : supplyT;
}

void Server::step(double dt, double supplyT) {
    for (int b = 0; b < kBranches; ++b) {
        chip(b).setInletT(supplyT);
        chip(b).step(dt);
    }
    air_.step(dt);
    updateReturnT(supplyT);
}

void Server::solveSteady(double supplyT) {
    for (int b = 0; b < kBranches; ++b) {
        chip(b).setInletT(supplyT);
        chip(b).solveSteady();
    }
    air_.solveSteady();
    for (int k = 0; k < 3; ++k) updateReturnT(supplyT); // cp depends weakly on returnT
}

double Server::serverDp(const HydraulicNetwork& net) const {
    return net.edgeDp(qdIn_) + net.edgeDp(branchEdge_[0]) + net.edgeDp(qdOut_);
}

double Server::peakT() const {
    double m = -1e9;
    for (int b = 0; b < kBranches; ++b) m = std::max(m, chip(b).peakT());
    return m;
}

int Server::hottestBranch() const {
    int h = 0;
    for (int b = 1; b < kBranches; ++b)
        if (chip(b).peakT() > chip(h).peakT()) h = b;
    return h;
}

double Server::liquidPower() const {
    double s = 0.0;
    for (int b = 0; b < kBranches; ++b) s += chip(b).power();
    return s;
}

LiquidLedger Server::ledger() const {
    LiquidLedger l;
    for (int b = 0; b < kBranches; ++b) {
        const ChipLedger& c = chip(b).ledger();
        l.powerIn += c.powerIn;
        l.toCoolant += c.toCoolant;
        l.storageRate += c.storageRate;
    }
    return l;
}

// ---- ServerRig -----------------------------------------------------------------------------------
ServerRig::ServerRig(const ServerParams& p, std::shared_ptr<const Fluid> fluid, double totalFlowLpm, double supplyT)
    : fluid_(fluid), net_(std::make_unique<HydraulicNetwork>(fluid)), server_(p, fluid, 0),
      flowLpm_(totalFlowLpm), supplyT_(supplyT) {
    sup_ = net_->addNode("supply connector");
    ret_ = net_->addNode("return connector");
    net_->setReference(ret_, 150e3);
    server_.attach(*net_, sup_, ret_, 0.0);
    server_.autoBalance(*net_);
    for (int b = 0; b < Server::kBranches; ++b) server_.chip(b).setInletT(supplyT);
    solveHydraulics();
}

void ServerRig::setTotalFlowLpm(double lpm) { flowLpm_ = lpm; solveHydraulics(); }

double ServerRig::volumeFlow() const { return lpmToM3s(flowLpm_); }

void ServerRig::solveHydraulics() {
    const double rho = fluid_->props(supplyT_).rho;
    net_->setHydrostaticDensity(rho);
    net_->setInjection(sup_, rho * lpmToM3s(flowLpm_));
    server_.updateEdgeTemperatures(*net_, supplyT_);
    net_->solve();
    server_.pullFlows(*net_);
}

void ServerRig::step(double dt) {
    solveHydraulics();
    server_.step(dt, supplyT_);
}

void ServerRig::solveSteady() {
    double lastPeak = 1e9;
    for (int it = 0; it < 50; ++it) {
        solveHydraulics();
        server_.solveSteady(supplyT_);
        if (std::abs(server_.peakT() - lastPeak) < 1e-6) break;
        lastPeak = server_.peakT();
    }
}

double ServerRig::indicativePumpPower(double etaPump, double etaMotor) const {
    return totalDp() * volumeFlow() / (etaPump * etaMotor);
}

} // namespace mf
