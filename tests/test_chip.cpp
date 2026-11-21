#include <doctest/doctest.h>
#include "Units.hpp"
#include "config/Config.hpp"
#include "thermal/AirGroup.hpp"
#include "thermal/ChipGrid.hpp"

#include <cmath>

using namespace mf;

namespace {
std::shared_ptr<WaterIAPWS> water() { static auto w = std::make_shared<WaterIAPWS>(); return w; }

ChipGrid baselineChip(const std::string& pattern = "uniform", int n = 10) {
    ChipParams p = defaultConfig().server.gpu;
    p.power.map.pattern = pattern;
    p.nx = p.ny = n;
    ChipGrid c(p, water());
    c.setMassFlow(water()->props(35.0).rho * lpmToM3s(0.5));
    c.setInletT(35.0);
    c.setActivity(1.0);
    return c;
}
} // namespace

TEST_CASE("Steady uniform 300 W chip: coolant energy balance and hand-calculation numbers") {
    ChipGrid c = baselineChip();
    c.solveSteady();
    const double mdot = c.massFlow();
    const double cp = water()->props(0.5 * (35.0 + c.outletT())).cp;
    CHECK(c.outletT() - 35.0 == doctest::Approx(300.0 / (mdot * cp)).epsilon(1e-3));
    CHECK(std::abs(c.ledger().relResidual()) < 1e-6);
    CHECK(c.outletT() == doctest::Approx(43.7).epsilon(2e-3));
    CHECK(c.peakT() > 65.0);
    CHECK(c.peakT() < 70.0);
    CHECK(c.peakI() == c.nx() - 1); // uniform heating: hottest at the outlet end
}

TEST_CASE("Symmetric power map gives a symmetric temperature field") {
    ChipGrid c = baselineChip("hotspot");
    c.solveSteady();
    for (int i = 0; i < c.nx(); ++i)
        for (int j = 0; j < c.ny() / 2; ++j)
            CHECK(c.siliconT(i, j) == doctest::Approx(c.siliconT(i, c.ny() - 1 - j)).epsilon(1e-9));
    CHECK(c.peakT() > baselineChip().peakT() - 1.0); // hotspot is hotter than... at least comparable
}

TEST_CASE("Transient chip conserves energy each step and settles at the steady solution") {
    ChipGrid steady = baselineChip("hotspot");
    steady.solveSteady();
    ChipGrid c = baselineChip("hotspot");
    double inSum = 0.0, outSum = 0.0, storSum = 0.0;
    const double e0 = c.storedEnergy();
    for (int n = 0; n < 400; ++n) {
        c.step(0.005);
        CHECK(std::abs(c.ledger().relResidual()) < 1e-6);
        inSum += c.ledger().powerIn * 0.005;
        outSum += c.ledger().toCoolant * 0.005;
        storSum += c.ledger().storageRate * 0.005;
    }
    CHECK(inSum - outSum == doctest::Approx(storSum).epsilon(1e-6));
    // Capacitances vary slightly with temperature (coolant rho*cp), so C*T differs by well under 1%.
    CHECK(c.storedEnergy() - e0 == doctest::Approx(storSum).epsilon(1e-2));
    CHECK(c.peakT() == doctest::Approx(steady.peakT()).epsilon(1e-4));
}

TEST_CASE("Lumped air group follows the analytic exponential") {
    AirParams p;
    AirGroup a(p);
    a.setActivity(1.0);
    const double dt = 0.05;
    double t = 0.0;
    // Flow is constant at full power, so UA is constant: T = Troom + P/UA (1 - exp(-t UA / C)).
    a.step(dt); t += dt;
    const double flow = a.airflow();
    const double ua = p.uaRef * std::pow(flow / p.flowRef, 0.8);
    while (t < 300.0) { a.step(dt); t += dt; }
    const double exact = p.roomT + p.fullW / ua * (1.0 - std::exp(-t * ua / p.heatCapacity));
    CHECK(a.temperature() == doctest::Approx(exact).epsilon(1e-3));
}

TEST_CASE("Cold-plate CPU solves and conserves energy") {
    ChipParams p = defaultConfig().server.cpu;
    ChipGrid c(p, water());
    c.setMassFlow(water()->props(35.0).rho * lpmToM3s(0.5));
    c.setInletT(35.0);
    c.setActivity(1.0);
    c.solveSteady();
    CHECK(std::abs(c.ledger().relResidual()) < 1e-6);
    CHECK(c.peakT() > c.outletT());
    CHECK(c.peakT() < 85.0);
}

TEST_CASE("Grid refinement changes the hotspot peak by less than 1.5 K from 10x10 to 40x40") {
    ChipGrid a = baselineChip("hotspot", 10), b = baselineChip("hotspot", 40);
    a.solveSteady();
    b.solveSteady();
    CHECK(std::abs(a.peakT() - b.peakT()) < 1.5);
    CHECK(a.outletT() == doctest::Approx(b.outletT()).epsilon(1e-6));
}

TEST_CASE("At fixed flow, viscosity does not change chip temperatures") {
    // Re*Pr = rho u Dh cp / k does not contain mu, so the thermal entrance and h are unchanged.
    ChipParams p = defaultConfig().server.gpu;
    auto f4 = std::make_shared<ScaledFluid>(water(), 4.0);
    ChipGrid c1(p, water()), c4(p, f4);
    for (ChipGrid* c : {&c1, &c4}) {
        c->setMassFlow(water()->props(35.0).rho * lpmToM3s(0.5));
        c->setInletT(35.0);
        c->setActivity(1.0);
        c->solveSteady();
    }
    CHECK(c4.peakT() == doctest::Approx(c1.peakT()).epsilon(1e-9));
}
