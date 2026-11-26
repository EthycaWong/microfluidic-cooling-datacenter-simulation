#include <doctest/doctest.h>
#include "cfd/DuctFlow.hpp"
#include "cfd/Presets.hpp"
#include "hydraulics/RectChannel.hpp"

#include <cmath>

using namespace mf;
using namespace mf::cfd;

TEST_CASE("Fully developed duct flow from Navier-Stokes matches the correlations used by the chip model") {
    const DuctResult sq = solveRectDuct(1.0, 1.0, 32);
    CHECK(sq.fRe == doctest::Approx(14.227).epsilon(0.005));         // Shah & London, square duct
    CHECK(sq.nuH1FourWalls == doctest::Approx(3.608).epsilon(0.01)); // H1, square duct
    const DuctResult ch = solveRectDuct(200e-6, 500e-6, 32);          // the proposal's microchannel
    CHECK(ch.fRe == doctest::Approx(RectChannel::fanningFRe(0.4)).epsilon(0.005));
    CHECK(ch.nuH1ThreeWalls == doctest::Approx(RectChannel::nu3WallFullyDeveloped(0.4)).epsilon(0.015));
}

TEST_CASE("Lid-driven cavity at Re = 100 reproduces Ghia, Ghia & Shin (1982)") {
    Case c = makeCase("lid-cavity", 100, 32);
    NavierStokes2D s(c.params);
    REQUIRE(s.runToSteady(1e-5, 1.0, 100000));
    const double ys[] = {0.1016, 0.2813, 0.5, 0.7344, 0.9531};
    const double ghia[] = {-0.06434, -0.15662, -0.20581, 0.00332, 0.68717};
    for (int k = 0; k < 5; ++k) {
        double u, v;
        s.velocityAt(0.5, ys[k], u, v);
        CHECK(std::abs(u - ghia[k]) < 0.02);
    }
    CHECK(s.maxDivergence() < 1e-5); // incompressible
}

TEST_CASE("Differentially heated cavity at Ra = 1e4 reproduces de Vahl Davis (1983)") {
    Case c = makeCase("heated-cavity", 1e4, 32);
    NavierStokes2D s(c.params);
    REQUIRE(s.runToSteady(1e-4, 1.0, 200000));
    const double nuHot = s.wallHeatFlux(0), nuCold = -s.wallHeatFlux(1); // L = k = dT = 1
    CHECK(nuHot == doctest::Approx(2.243).epsilon(0.03));
    CHECK(nuHot == doctest::Approx(nuCold).epsilon(0.01)); // steady energy balance: heat in = heat out
}

TEST_CASE("Rayleigh-Benard: conduction below the critical Rayleigh number, convection above it") {
    Case lo = makeCase("rayleigh-benard", 1000, 16), hi = makeCase("rayleigh-benard", 2e4, 16);
    NavierStokes2D a(lo.params), b(hi.params);
    // seed a small asymmetry so the unstable case can choose a roll direction
    a.params().t0 = 0.5; b.params().left.thermal = ThermalBC::FixedTemperature; b.params().left.T = 0.55;
    b.runUntil(1.0);
    b.params().left.thermal = ThermalBC::HeatFlux;
    a.runToSteady(1e-5, 1.0, 100000);
    b.runToSteady(1e-5, 1.0, 200000);
    CHECK(a.wallHeatFlux(2) == doctest::Approx(1.0).epsilon(0.02)); // Nu = 1: pure conduction
    CHECK(b.wallHeatFlux(2) > 1.5);                                // convection rolls carry more heat
}

TEST_CASE("Channel flow: fully developed pressure gradient, parabolic profile, mass conservation") {
    Case c = makeCase("poiseuille", 100, 12);
    NavierStokes2D s(c.params);
    REQUIRE(s.runToSteady(1e-5, c.timeScale, 200000));
    const int i1 = s.nx() / 2, i2 = 3 * s.nx() / 4;
    const double dpdx = (s.meanPressureColumn(i2) - s.meanPressureColumn(i1)) / ((i2 - i1) * s.dx());
    CHECK(dpdx == doctest::Approx(-12.0 * c.params.nu * 1.0 / 1.0).epsilon(0.03));
    double u, v;
    s.velocityAt(8.0, 0.5, u, v);
    CHECK(u == doctest::Approx(1.5).epsilon(0.03)); // u_max = 1.5 u_mean between plates
    double in = 0, out = 0;
    for (int j = 1; j <= s.ny(); ++j) { in += s.uFace(0, j); out += s.uFace(s.nx(), j); }
    CHECK(out == doctest::Approx(in).epsilon(1e-6));
}

TEST_CASE("Boundary-driven: walls at rest give no flow; reversing the lid mirrors the flow") {
    Case c = makeCase("lid-cavity", 100, 16);
    c.params.top.u = 0.0; // every wall u = v = 0, isothermal: the fluid must stay at rest
    NavierStokes2D rest(c.params);
    rest.runUntil(2.0);
    CHECK(rest.maxSpeed() < 1e-12);
    Case f = makeCase("lid-cavity", 100, 16), r = makeCase("lid-cavity", 100, 16);
    r.params.top.u = -1.0;
    NavierStokes2D a(f.params), b(r.params);
    a.runUntil(3.0);
    b.runUntil(3.0);
    for (int j = 1; j <= 16; j += 3) CHECK(b.uc(16 - 4 + 1, j) == doctest::Approx(-a.uc(4, j)).epsilon(1e-6));
}
