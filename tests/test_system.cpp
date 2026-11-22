#include <doctest/doctest.h>
#include "Units.hpp"
#include "config/Config.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>

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

TEST_CASE("Job trace replay: idle before and after the job, phases, offset and time scale") {
    const std::string path = "job_trace_test.csv";
    {
        std::ofstream f(path);
        f << "t_s,u,phase,train_samples,infer_samples\n0.1,0.0,setup,0,0\n0.2,0.8,train,128,0\n0.3,0.6,train,256,0\n"
             "0.4,0.1,scan,256,32\n";
    }
    JobTraceWorkload w(path, 10.0, 2.0); // starts at t = 10 s, stretched x2
    CHECK(w.u(5.0) == 0.0);
    CHECK(w.at(5.0).phase == "waiting");
    CHECK(w.u(10.0 + 2.0 * 0.15) == doctest::Approx(0.8)); // trace t 0.15 -> window ending at 0.2
    CHECK(w.at(10.0 + 2.0 * 0.25).phase == "train");
    CHECK(w.at(10.0 + 2.0 * 0.35).inferSamples == 32);
    CHECK(w.durationSim() == doctest::Approx(0.8));
    CHECK(w.u(20.0) == 0.0);
    CHECK(w.at(20.0).phase == "done");
    std::remove(path.c_str());
}

TEST_CASE("Live job trace follows a file as it grows") {
    const std::string path = "job_trace_live_test.csv";
    std::remove(path.c_str());
    JobTraceWorkload w(path, 0.0, 1.0, true);
    CHECK(w.at(0.0).phase == "waiting"); // no file yet: servers wait on demand
    {
        std::ofstream f(path);
        f << "t_s,u,phase,train_samples,infer_samples\n0.1,0.7,train,128,0\n";
    }
    w.refresh();
    CHECK(w.u(123.0) == doctest::Approx(0.7)); // live: newest sample, whatever the simulation time
    {
        std::ofstream f(path, std::ios::app);
        f << "0.2,0.2,eval,128,0\n";
    }
    w.refresh();
    CHECK(w.at(0.0).phase == "eval");
    CHECK(w.size() == 2);
    std::remove(path.c_str());
}

TEST_CASE("Feedforward pump control limits overshoot on an on-demand load step") {
    SimConfig c = defaultConfig();
    c.server.gpu.power.map.pattern = "hotspot";
    c.rack.pumpMode = PumpMode::PeakTemperature;
    c.rack.targetPeakC = 80.0;
    Rack r(c.rack, c.server, water());
    r.setAllWorkloads(std::make_shared<StepWorkload>(0.0, 1.0, 5.0)); // idle, then a job arrives at t = 5 s
    r.solveSteady();                                                   // idle steady state (pump at its floor)
    CHECK(r.pumpSpeed() == doctest::Approx(c.rack.peakModeMinSpeed));
    double maxPeak = 0.0;
    for (int n = 0; n < 2000; ++n) { r.step(0.01); maxPeak = std::max(maxPeak, r.hotspot().t); } // 20 s
    CHECK(maxPeak < 86.0);                                               // feedback alone overshot past 100 C
    CHECK(r.hotspot().t == doctest::Approx(80.0).epsilon(1.0 / 80.0));
    CHECK(r.pumpSpeed() > 2.0 * c.rack.peakModeMinSpeed);               // the pump followed the load
}
