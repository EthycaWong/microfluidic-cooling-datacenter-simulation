#include <doctest/doctest.h>
#include "Units.hpp"
#include "hydraulics/HydraulicNetwork.hpp"

#include <cmath>

using namespace mf;

TEST_CASE("Shah & London rectangular duct limits") {
    CHECK(RectChannel::fanningFRe(0.0) == doctest::Approx(24.0));            // parallel plates
    CHECK(RectChannel::fanningFRe(1.0) == doctest::Approx(14.227).epsilon(2e-3)); // square duct
    CHECK(RectChannel::nu3WallFullyDeveloped(0.0) == doctest::Approx(8.235)); // plates, one side insulated
}

TEST_CASE("Fully developed channel pressure drop equals the analytic f*Re result and scales with viscosity") {
    ChannelGeometry g;
    g.kInlet = g.kOutlet = 0.0;
    g.hagenbach = false;
    RectChannel ch(g);
    FluidProps p{994.0, 7.2e-4, 4178.0, 0.62};
    const double q = lpmToM3s(0.5);
    const double u = q / g.count / g.area();
    const double dh = g.hydraulicDiameter();
    const double re = p.rho * u * dh / p.mu;
    const double analytic = 4.0 * RectChannel::fanningFRe(g.aspect()) / re * (g.length / dh) * 0.5 * p.rho * u * u;
    CHECK(ch.pressureDrop(q, p) == doctest::Approx(analytic).epsilon(1e-12));
    FluidProps p2 = p;
    p2.mu *= 2.0;
    CHECK(ch.pressureDrop(q, p2) / ch.pressureDrop(q, p) == doctest::Approx(2.0).epsilon(1e-12));
}

TEST_CASE("Baseline channel reproduces the hand calculation") {
    WaterIAPWS w;
    RectChannel ch(ChannelGeometry{});
    const ChannelState s = ch.state(lpmToM3s(0.5), w.props(35.0));
    CHECK(s.re == doctest::Approx(660).epsilon(0.02));
    CHECK(s.velocity == doctest::Approx(1.667).epsilon(1e-3));
    CHECK(s.dpFriction == doctest::Approx(9600).epsilon(0.02));
    CHECK(s.dpTotal > 12e3);
    CHECK(s.dpTotal < 15e3);
}

TEST_CASE("Two laminar branches split flow inversely to resistance, mass is conserved") {
    auto w = std::make_shared<WaterIAPWS>();
    const FluidProps ref = w->props(35.0);
    HydraulicNetwork net(w);
    const int a = net.addNode("a"), b = net.addNode("b");
    net.setReference(b, 0.0);
    const double mTot = 0.03;
    net.setInjection(a, mTot);
    // Purely laminar elements: dp = dpRef * Q/Qref, so resistance ratio 1 : 3
    net.addEdge(a, b, std::make_unique<RefLoss>("r1", 10e3, 1e-5, 1.0, ref));
    net.addEdge(a, b, std::make_unique<RefLoss>("r2", 30e3, 1e-5, 1.0, ref));
    REQUIRE(net.solve());
    CHECK(net.flow(0) + net.flow(1) == doctest::Approx(mTot).epsilon(1e-9));
    CHECK(net.flow(0) / net.flow(1) == doctest::Approx(3.0).epsilon(1e-5));
    CHECK(net.edgeDp(0) == doctest::Approx(net.edgeDp(1)).epsilon(1e-9));
}

namespace {
double loopFlow(double zMid) {
    auto w = std::make_shared<WaterIAPWS>();
    const FluidProps ref = w->props(35.0);
    HydraulicNetwork net(w);
    const int a = net.addNode("suction", 0.0), b = net.addNode("discharge", 0.0);
    const int c = net.addNode("high point", zMid), d = net.addNode("mid", 0.5 * zMid);
    net.setReference(a, 150e3);
    net.addEdge(a, b, std::make_unique<PumpEdge>("pump", PumpCurve{}));
    net.addEdge(b, c, std::make_unique<PipeEdge>("up", 0.02, 3.0, 2.0));
    net.addEdge(c, d, std::make_unique<PipeEdge>("down", 0.02, 3.0, 2.0));
    net.addEdge(d, a, std::make_unique<RefLoss>("load", 40e3, lpmToM3s(12), 0.5, ref));
    net.solve();
    return net.flow(0);
}
} // namespace

TEST_CASE("Elevation changes cancel around a closed loop") {
    const double flat = loopFlow(0.0);
    CHECK(flat > 0.0);
    CHECK(loopFlow(1.8) == doctest::Approx(flat).epsilon(1e-9));
    CHECK(loopFlow(5.0) == doctest::Approx(flat).epsilon(1e-9));
}

TEST_CASE("Pump curve and efficiency") {
    PumpEdge p("pump", PumpCurve{});
    CHECK(p.head(0.0) == doctest::Approx(100e3));
    CHECK(std::abs(p.head(lpmToM3s(20))) < 1.0);
    p.setSpeed(0.5);
    CHECK(p.head(0.0) == doctest::Approx(25e3)); // affinity law: H ~ s^2
    CHECK(p.efficiency(lpmToM3s(6)) == doctest::Approx(0.5)); // BEP flow scales with speed
}
