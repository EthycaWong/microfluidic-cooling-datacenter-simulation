#include <doctest/doctest.h>
#include "Units.hpp"
#include "config/Config.hpp"

#include <cmath>

using namespace mf;

namespace {
std::shared_ptr<WaterIAPWS> water() { static auto w = std::make_shared<WaterIAPWS>(); return w; }
}

TEST_CASE("Server auto-balance gives 0.5 L/min per branch at the isothermal baseline") {
    SimConfig c = defaultConfig();
    ServerRig rig(c.server, water(), 1.5, 35.0);
    for (int b = 0; b < Server::kBranches; ++b) rig.server().chip(b).setPower(0.0); // isothermal
    rig.solveSteady();
    const double rho = water()->props(35.0).rho;
    for (int b = 0; b < Server::kBranches; ++b) {
        CHECK(m3sToLpm(rig.server().branchFlow(b) / rho) == doctest::Approx(0.5).epsilon(0.01));
        CHECK(rig.server().branchDp(rig.network(), b) == doctest::Approx(rig.server().branchDp(rig.network(), 0)).epsilon(1e-9));
    }
}

TEST_CASE("Without balancing valves the low-resistance CPU cold plate takes extra flow") {
    SimConfig c = defaultConfig();
    c.server.autoBalance = false;
    ServerRig rig(c.server, water(), 1.5, 35.0);
    rig.server().setActivity(1.0);
    rig.solveSteady();
    CHECK(rig.server().branchFlow(2) > 1.3 * rig.server().branchFlow(0));
    CHECK(rig.server().totalFlow() == doctest::Approx(water()->props(35.0).rho * lpmToM3s(1.5)).epsilon(1e-9));
}

TEST_CASE("Rack at 12 L/min and full load: mixed return near 42.7 C, energy and mass balance") {
    SimConfig c = defaultConfig();
    Rack r(c.rack, c.server, water());
    r.solveSteady();
    CHECK(r.totalFlowLpm() == doctest::Approx(12.0).epsilon(1e-6));
    CHECK(r.returnT() == doctest::Approx(42.7).epsilon(0.2 / 42.7));
    CHECK(std::abs(r.ledger().relResidual()) < 1e-6);
    double sum = 0.0;
    for (int k = 0; k < r.serverCount(); ++k) sum += r.server(k).totalFlow();
    CHECK(sum == doctest::Approx(r.totalMassFlow()).epsilon(1e-9));
    // Mass conservation at every node.
    const HydraulicNetwork& net = r.network();
    std::vector<double> bal(net.nodeCount(), 0.0);
    for (int e = 0; e < net.edgeCount(); ++e) { bal[net.edgeFrom(e)] -= net.flow(e); bal[net.edgeTo(e)] += net.flow(e); }
    for (double b : bal) CHECK(std::abs(b) < 1e-9);
    // Pump operating point lies on its curve and equals the pressure rise across it (nodes 0 -> 1).
    CHECK(net.pressure(1) - net.pressure(0) == doctest::Approx(r.pumpHead()).epsilon(1e-6));
    const PumpCurve& pc = c.rack.pump;
    const double q = r.totalMassFlow() / water()->props(35.0).rho;
    const double s = r.pumpSpeed();
    CHECK(r.pumpHead() == doctest::Approx(s * s * pc.shutoffHeadPa - pc.shutoffHeadPa / (pc.maxFlowM3s * pc.maxFlowM3s) * q * q).epsilon(1e-6));
}

TEST_CASE("Reverse-return manifold shares flow more evenly than direct return") {
    SimConfig c = defaultConfig();
    RackParams direct = c.rack, reverse = c.rack;
    reverse.layout = ManifoldLayout::ReverseReturn;
    Rack a(direct, c.server, water()), b(reverse, c.server, water());
    a.solveSteady();
    b.solveSteady();
    CHECK(b.flowMaldistribution() < a.flowMaldistribution());
    // Direct return: the lowest server (shortest path) gets the most flow.
    CHECK(a.server(0).totalFlow() > a.server(7).totalFlow());
}

TEST_CASE("PI pump control holds the hottest chip at the target") {
    SimConfig c = defaultConfig();
    c.server.gpu.power.map.pattern = "hotspot";
    c.rack.pumpMode = PumpMode::PeakTemperature;
    c.rack.targetPeakC = 80.0;
    Rack r(c.rack, c.server, water());
    for (int n = 0; n < 3000; ++n) r.step(0.01);
    CHECK(r.hotspot().t == doctest::Approx(80.0).epsilon(0.1 / 80.0));
    CHECK(std::abs(r.ledger().relResidual()) < 1e-5);
}

TEST_CASE("Rack viscosity: same flow holds the same peak, pump power rises") {
    SimConfig c = defaultConfig();
    Rack a(c.rack, c.server, water()), b(c.rack, c.server, water());
    b.setViscosityScale(4.0);
    a.solveSteady();
    b.solveSteady();
    // Within 0.2 K: thicker coolant spreads flow slightly more evenly across servers (more laminar).
    CHECK(std::abs(b.hotspot().t - a.hotspot().t) < 0.2);
    CHECK(b.flowMaldistribution() < a.flowMaldistribution());
    CHECK(b.pumpHead() > 1.8 * a.pumpHead());
    CHECK(b.pumpPower() > 1.8 * a.pumpPower());
}

TEST_CASE("Config JSON overrides defaults with unit conversion") {
    SimConfig c = defaultConfig();
    applyConfigJson(c, R"({ "gpu": { "channels": { "width_um": 150, "count": 40 } },
                           "rack": { "layout": "reverse-return", "target_flow_lpm": 10, "pump": { "shutoff_kPa": 120 } },
                           "workload": { "type": "square", "period_s": 4 } })");
    CHECK(c.server.gpu.channels.width == doctest::Approx(150e-6));
    CHECK(c.server.gpu.channels.count == 40);
    CHECK(std::string(toString(c.rack.layout)) == "reverse-return");
    CHECK(c.rack.targetFlowLpm == 10.0);
    CHECK(c.rack.pump.shutoffHeadPa == doctest::Approx(120e3));
    CHECK(makeWorkload(c.workload)->u(1.0) == doctest::Approx(1.0));
    CHECK(makeWorkload(c.workload)->u(3.0) == doctest::Approx(0.2));
}
