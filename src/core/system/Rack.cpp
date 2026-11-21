#include "system/Rack.hpp"
#include "Units.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace mf {

const char* toString(ManifoldLayout l) {
    return l == ManifoldLayout::DirectReturn ? "direct-return" : "reverse-return";
}
const char* toString(PumpMode m) {
    switch (m) {
    case PumpMode::FixedSpeed: return "fixed-speed";
    case PumpMode::FixedFlow: return "fixed-flow";
    case PumpMode::PeakTemperature: return "peak-temperature";
    }
    return "?";
}

namespace {
constexpr double kPumpZ = 0.10; // m, pump and CDU sit in the bottom rack units
constexpr double kCduZ = 0.15;
} // namespace

Rack::Rack(const RackParams& rp, const ServerParams& sp, std::shared_ptr<const Fluid> baseFluid)
    : rp_(rp), sp_(sp), fluid_(std::make_shared<ScaledFluid>(std::move(baseFluid), 1.0)) {
    servers_.reserve(rp_.servers);
    for (int k = 0; k < rp_.servers; ++k) {
        servers_.emplace_back(sp_, fluid_, k);
        workloads_.push_back(std::make_shared<ConstantWorkload>(1.0));
    }
    returnSegT_.assign(rp_.servers, rp_.supplyT);
    returnT_ = rp_.supplyT;
    buildNetwork();
    for (auto& s : servers_) s.setActivity(workloads_[s.index()]->u(0.0));
    solveHydraulics();
}

double Rack::serverTapZ(int k) const {
    return (rp_.firstServerU - 1 + k * rp_.serverU + 0.5 * rp_.serverU) * kRackUnit;
}

void Rack::buildNetwork() {
    net_ = std::make_unique<HydraulicNetwork>(fluid_);
    const int n = rp_.servers;
    const double D = rp_.manifoldDiameter;
    const FluidProps ref = fluid_->props(rp_.supplyT);
    supplySeg_.assign(n, -1);
    returnSeg_.assign(n, -1);
    supplyNode_.assign(n, -1);
    returnNode_.assign(n, -1);

    pumpIn_ = net_->addNode("pump suction", kPumpZ);
    const int pumpOut = net_->addNode("pump discharge", kPumpZ);
    auto pump = std::make_unique<PumpEdge>("pump", rp_.pump);
    pump->setSpeed(rp_.pumpSpeed);
    pump_ = pump.get();
    pumpEdge_ = net_->addEdge(pumpIn_, pumpOut, std::move(pump));

    for (int k = 0; k < n; ++k) {
        supplyNode_[k] = net_->addNode("supply manifold tap " + std::to_string(k + 1), serverTapZ(k));
        returnNode_[k] = net_->addNode("return manifold tap " + std::to_string(k + 1), serverTapZ(k));
    }
    // Supply manifold rises from the pump; segment k ends at tap k.
    supplySeg_[0] = net_->addEdge(pumpOut, supplyNode_[0],
        std::make_unique<PipeEdge>("supply riser", D, rp_.connectLength + serverTapZ(0) - kPumpZ, rp_.connectK));
    for (int k = 1; k < n; ++k)
        supplySeg_[k] = net_->addEdge(supplyNode_[k - 1], supplyNode_[k],
            std::make_unique<PipeEdge>("supply manifold " + std::to_string(k), D, serverTapZ(k) - serverTapZ(k - 1), rp_.manifoldKPerTap));

    for (auto& s : servers_) s.attach(*net_, supplyNode_[s.index()], returnNode_[s.index()], serverTapZ(s.index()));

    const int cduIn = net_->addNode("CDU inlet", kCduZ);
    // Return manifold; returnSeg_[k] carries the flow leaving tap k towards the CDU.
    if (rp_.layout == ManifoldLayout::DirectReturn) {
        // Flow runs down the return manifold: first in, first out (shortest path for the lowest server).
        for (int k = n - 1; k >= 1; --k)
            returnSeg_[k] = net_->addEdge(returnNode_[k], returnNode_[k - 1],
                std::make_unique<PipeEdge>("return manifold " + std::to_string(k), D, serverTapZ(k) - serverTapZ(k - 1), rp_.manifoldKPerTap));
        returnSeg_[0] = net_->addEdge(returnNode_[0], cduIn,
            std::make_unique<PipeEdge>("return drop", D, rp_.connectLength + serverTapZ(0) - kCduZ, rp_.connectK));
    } else {
        // Reverse return: flow runs up and leaves at the top, so every server path has equal length.
        for (int k = 0; k < n - 1; ++k)
            returnSeg_[k] = net_->addEdge(returnNode_[k], returnNode_[k + 1],
                std::make_unique<PipeEdge>("return manifold " + std::to_string(k + 1), D, serverTapZ(k + 1) - serverTapZ(k), rp_.manifoldKPerTap));
        returnSeg_[n - 1] = net_->addEdge(returnNode_[n - 1], cduIn,
            std::make_unique<PipeEdge>("return riser", D, rp_.connectLength + serverTapZ(n - 1) - kCduZ, rp_.connectK));
    }
    const int filterOut = net_->addNode("filter outlet", kCduZ);
    filterEdge_ = net_->addEdge(cduIn, filterOut,
        std::make_unique<RefLoss>("filter", rp_.filterDpRefPa, lpmToM3s(rp_.filterQRefLpm), rp_.filterLaminarFraction, ref));
    hxEdge_ = net_->addEdge(filterOut, pumpIn_,
        std::make_unique<RefLoss>("CDU heat exchanger (placeholder)", rp_.hxDpRefPa, lpmToM3s(rp_.hxQRefLpm), rp_.hxLaminarFraction, ref));
    net_->setReference(pumpIn_, rp_.referencePressurePa);

    for (auto& s : servers_) s.autoBalance(*net_);
}

void Rack::setWorkload(int k, std::shared_ptr<Workload> w) { workloads_[k] = std::move(w); }
void Rack::setAllWorkloads(const std::shared_ptr<Workload>& w) {
    for (auto& x : workloads_) x = w;
}
void Rack::setViscosityScale(double s, bool rebalance) {
    fluid_->setScale(s);
    // Re-commission the balancing valves with the new coolant (as the server rig does), so a
    // viscosity comparison changes only the coolant, not the valve settings.
    if (rebalance)
        for (auto& sv : servers_) sv.autoBalance(*net_);
    solveHydraulics();
}
void Rack::setPumpMode(PumpMode m) { rp_.pumpMode = m; lastErr_ = 0.0; }
void Rack::setPumpSpeed(double s) {
    rp_.pumpSpeed = std::clamp(s, rp_.pump.minSpeed, rp_.pump.maxSpeed);
    pump_->setSpeed(rp_.pumpSpeed);
}
void Rack::setLayout(ManifoldLayout l) {
    if (l == rp_.layout) return;
    const double s = pump_->speed();
    rp_.layout = l;
    buildNetwork();
    pump_->setSpeed(s);
    solveHydraulics();
}
void Rack::setAutoBalance(bool on) {
    sp_.autoBalance = on;
    for (auto& s : servers_) { s.setBalanceEnabled(on); s.autoBalance(*net_); }
    solveHydraulics();
}

void Rack::updateReturnTemperatures() {
    const int n = serverCount();
    auto mix = [&](int a, int b) { // servers a..b inclusive
        double q = 0.0, m = 0.0;
        for (int k = a; k <= b; ++k) { q += servers_[k].ledger().toCoolant; m += servers_[k].totalFlow(); }
        const double cp = fluid_->props(rp_.supplyT + 4.0).cp;
        return m > 0 ? rp_.supplyT + q / (m * cp) : rp_.supplyT;
    };
    for (int k = 0; k < n; ++k)
        returnSegT_[k] = rp_.layout == ManifoldLayout::DirectReturn ? mix(k, n - 1) : mix(0, k);
    returnT_ = mix(0, n - 1);
}

void Rack::updateEdgeTemperatures() {
    for (int e : supplySeg_) net_->edge(e).setTemperature(rp_.supplyT);
    for (auto& s : servers_) s.updateEdgeTemperatures(*net_, rp_.supplyT);
    for (int k = 0; k < serverCount(); ++k) net_->edge(returnSeg_[k]).setTemperature(returnSegT_[k]);
    net_->edge(filterEdge_).setTemperature(returnT_);
    net_->edge(hxEdge_).setTemperature(0.5 * (returnT_ + rp_.supplyT));
    net_->edge(pumpEdge_).setTemperature(rp_.supplyT);
}

void Rack::solveHydraulics() {
    net_->setHydrostaticDensity(fluid_->props(rp_.supplyT).rho);
    updateEdgeTemperatures();
    if (rp_.pumpMode == PumpMode::FixedFlow) {
        const double target = rp_.targetFlowLpm;
        auto f = [&](double s) { pump_->setSpeed(s); net_->solve(); return totalFlowLpm() - target; };
        const double smin = rp_.pump.minSpeed, smax = rp_.pump.maxSpeed;
        double a = std::clamp(pump_->speed() * 0.95, smin, smax), b = std::clamp(pump_->speed() * 1.05, smin, smax);
        double fa = f(a), fb = f(b);
        while (fa > 0 && a > smin) { b = a; fb = fa; a = std::max(smin, a * 0.7); fa = f(a); }
        while (fb < 0 && b < smax) { a = b; fa = fb; b = std::min(smax, b * 1.4); fb = f(b); }
        if (fa > 0) { f(a); }
        else if (fb < 0) { f(b); }
        else {
            // Illinois (modified regula falsi).
            int side = 0;
            double s = b;
            for (int it = 0; it < 80; ++it) {
                s = (a * fb - b * fa) / (fb - fa);
                const double fs = f(s);
                if (std::abs(fs) < 1e-9 * target) break;
                if (fs * fb > 0) { b = s; fb = fs; if (side == -1) fa *= 0.5; side = -1; }
                else { a = s; fa = fs; if (side == +1) fb *= 0.5; side = +1; }
                if (std::abs(b - a) < 1e-12) break;
            }
            f(s);
        }
        rp_.pumpSpeed = pump_->speed();
    } else {
        pump_->setSpeed(rp_.pumpSpeed);
        net_->solve();
    }
    for (auto& s : servers_) s.pullFlows(*net_);
}

void Rack::pumpControlTick(double interval) {
    const double e = hotspot().t - rp_.targetPeakC;
    double s = rp_.pumpSpeed + rp_.piKp * (e - lastErr_) + rp_.piKi * e * interval;
    rp_.pumpSpeed = std::clamp(s, rp_.pump.minSpeed, rp_.pump.maxSpeed);
    lastErr_ = e;
}

void Rack::step(double dt) {
    if (!(dt > 0.0)) return; // a zero or negative step would divide by zero in the implicit solve
    t_ += dt;
    for (auto& s : servers_) s.setActivity(workloads_[s.index()]->u(t_));
    if (t_ >= nextHyd_ - 1e-12) {
        if (rp_.pumpMode == PumpMode::PeakTemperature) pumpControlTick(rp_.hydraulicInterval);
        solveHydraulics();
        nextHyd_ = t_ + rp_.hydraulicInterval;
    }
    for (auto& s : servers_) s.step(dt, rp_.supplyT);
    updateReturnTemperatures();
    pumpJ_ += pumpPower() * dt;
    fanJ_ += fanPower() * dt;
    itJ_ += itPower() * dt;
}

void Rack::steadyAtCurrentSpeed() {
    double lastFlow = -1.0, lastPeak = -1.0;
    for (int it = 0; it < 60; ++it) {
        solveHydraulics();
        for (auto& s : servers_) s.solveSteady(rp_.supplyT);
        updateReturnTemperatures();
        const double fl = totalMassFlow(), pk = hotspot().t;
        if (std::abs(fl - lastFlow) < 1e-10 && std::abs(pk - lastPeak) < 1e-6) break;
        lastFlow = fl;
        lastPeak = pk;
    }
}

void Rack::solveSteady() {
    for (auto& s : servers_) s.setActivity(workloads_[s.index()]->u(t_));
    if (rp_.pumpMode != PumpMode::PeakTemperature) {
        steadyAtCurrentSpeed();
        return;
    }
    // Find the pump speed whose steady hottest-chip temperature equals the target.
    auto peakAt = [&](double s) { rp_.pumpSpeed = s; steadyAtCurrentSpeed(); return hotspot().t; };
    double lo = rp_.pump.minSpeed, hi = rp_.pump.maxSpeed;
    if (peakAt(hi) > rp_.targetPeakC) return;  // cannot reach target even at max speed
    if (peakAt(lo) < rp_.targetPeakC) return;  // already below target at min speed
    for (int it = 0; it < 50 && hi - lo > 1e-6; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (peakAt(mid) > rp_.targetPeakC) lo = mid; else hi = mid;
    }
    peakAt(hi); // conservative side: peak <= target
}

double Rack::totalMassFlow() const { return net_->flow(pumpEdge_); }
double Rack::totalFlowLpm() const { return m3sToLpm(totalMassFlow() / fluid_->props(rp_.supplyT).rho); }
double Rack::pumpHead() const { return pump_->head(totalMassFlow() / fluid_->props(rp_.supplyT).rho); }
double Rack::pumpPower() const { return pump_->electricPower(totalMassFlow() / fluid_->props(rp_.supplyT).rho); }
double Rack::pumpEfficiency() const { return pump_->efficiency(totalMassFlow() / fluid_->props(rp_.supplyT).rho); }

double Rack::fanPower() const {
    double s = 0.0;
    for (const auto& sv : servers_) s += sv.air().fanPower();
    return s;
}
double Rack::airHeat() const {
    double s = 0.0;
    for (const auto& sv : servers_) s += sv.air().heatToAir();
    return s;
}
double Rack::itPower() const {
    double s = 0.0;
    for (const auto& sv : servers_) s += sv.totalPower();
    return s;
}

Hotspot Rack::hotspot() const {
    Hotspot h;
    h.t = -1e9;
    for (const auto& s : servers_)
        for (int b = 0; b < Server::kBranches; ++b) {
            const ChipGrid& c = s.chip(b);
            if (c.peakT() > h.t) h = Hotspot{s.index(), b, c.peakI(), c.peakJ(), c.peakT()};
        }
    return h;
}

LiquidLedger Rack::ledger() const {
    LiquidLedger l;
    for (const auto& s : servers_) l += s.ledger();
    return l;
}

double Rack::flowMaldistribution() const {
    double mn = 1e9, mx = -1e9, sum = 0.0;
    for (const auto& s : servers_) { mn = std::min(mn, s.totalFlow()); mx = std::max(mx, s.totalFlow()); sum += s.totalFlow(); }
    const double mean = sum / serverCount();
    return mean > 0 ? (mx - mn) / mean : 0.0;
}

} // namespace mf
