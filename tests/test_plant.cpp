#include <doctest/doctest.h>
#include "Units.hpp"
#include "config/Config.hpp"

#include <cmath>

using namespace mf;

namespace {
std::shared_ptr<WaterIAPWS> water() { static auto w = std::make_shared<WaterIAPWS>(); return w; }
SimConfig plantConfig(double outdoor, bool chiller = true) {
    SimConfig c = defaultConfig();
    c.plant.weather.meanC = outdoor;
    c.plant.chiller.enabled = chiller;
    return c;
}
} // namespace

TEST_CASE("Heat exchanger effectiveness limits") {
    CHECK(hx::effectivenessCounterflow(1.0, 1.0) == doctest::Approx(0.5));            // NTU/(1+NTU)
    CHECK(hx::effectivenessCounterflow(50.0, 0.5) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(hx::effectivenessCounterflow(2.0, 0.0) == doctest::Approx(1.0 - std::exp(-2.0)));
    CHECK(hx::effectivenessCrossflow(2.0, 1e-12, false) == doctest::Approx(1.0 - std::exp(-2.0)).epsilon(1e-6));
    // For the same NTU and Cr, counterflow beats crossflow.
    CHECK(hx::effectivenessCounterflow(1.5, 0.7) > hx::effectivenessCrossflow(1.5, 0.7, false));
}

TEST_CASE("Mild day: fans hold the supply setpoint, no chiller, heat balance closes") {
    SimConfig c = plantConfig(25.0);
    DataCenter dc(c.rack, c.server, c.plant, water());
    dc.solveSteady();
    const PlantState& p = dc.plant().state();
    CHECK(p.equipSupplyT == doctest::Approx(35.0).epsilon(1e-4));
    CHECK(p.fanSpeed < 1.0);
    CHECK(p.chillerHeat == 0.0);
    CHECK(p.dryCoolerHeat == doctest::Approx(dc.rack().ledger().toCoolant).epsilon(1e-6));
    CHECK(dc.rack().supplyT() == doctest::Approx(p.equipSupplyT).epsilon(1e-6));
    // Return - supply consistent with the rack heat and loop flow.
    const double cp = water()->props(39.0).cp;
    CHECK(dc.rack().returnT() - dc.rack().supplyT() ==
          doctest::Approx(dc.rack().ledger().toCoolant / (dc.rack().totalMassFlow() * cp)).epsilon(2e-3));
    CHECK(dc.evaporativeEquivalentLph() > 9.0); // ~6.4 kW / 2.43 MJ/kg
}

TEST_CASE("Hot day: the chiller holds the setpoint, without it the supply floats up") {
    SimConfig on = plantConfig(38.0, true), off = plantConfig(38.0, false);
    DataCenter a(on.rack, on.server, on.plant, water()), b(off.rack, off.server, off.plant, water());
    a.solveSteady();
    b.solveSteady();
    // At 38 C the economizer rule bypasses the dry cooler: its fans would cost more than they save.
    CHECK(a.plant().state().fanSpeed == 0.0);
    CHECK(a.plant().state().dryCoolerHeat == 0.0);
    CHECK(a.plant().state().chillerHeat > 0.0);
    CHECK(a.plant().state().equipSupplyT == doctest::Approx(35.0).epsilon(1e-4));
    CHECK(a.plant().state().dryCoolerHeat + a.plant().state().chillerHeat ==
          doctest::Approx(a.rack().ledger().toCoolant).epsilon(1e-6));
    CHECK(b.plant().state().equipSupplyT > 38.0);
    CHECK(b.rack().hotspot().t > a.rack().hotspot().t);
    CHECK(a.coolingPower() > b.coolingPower()); // chilling costs electricity
}

TEST_CASE("Warm day: fans at full speed plus a little trim chilling") {
    SimConfig c = plantConfig(30.0);
    DataCenter dc(c.rack, c.server, c.plant, water());
    dc.solveSteady();
    const PlantState& p = dc.plant().state();
    CHECK(p.fanSpeed == doctest::Approx(1.0));
    CHECK(p.dryCoolerHeat > 0.0);
    CHECK(p.chillerHeat > 0.0);
    CHECK(p.equipSupplyT == doctest::Approx(35.0).epsilon(1e-4));
}

TEST_CASE("Cold day: minimum fan speed already undershoots the setpoint") {
    SimConfig c = plantConfig(0.0);
    DataCenter dc(c.rack, c.server, c.plant, water());
    dc.solveSteady();
    CHECK(dc.plant().state().fanSpeed == doctest::Approx(c.plant.dryCooler.minFanSpeed));
    CHECK(dc.plant().state().equipSupplyT < 35.0);
}

TEST_CASE("Transient: after a load step the fan controller returns the supply to the setpoint") {
    SimConfig c = plantConfig(22.0);
    c.workload = WorkloadSpec{};
    DataCenter dc(c.rack, c.server, c.plant, water());
    dc.rack().setAllWorkloads(std::make_shared<StepWorkload>(0.3, 1.0, 20.0));
    dc.solveSteady(); // steady at u = 0.3 (t = 0)
    double maxSupply = 0.0;
    for (int n = 0; n < 30000; ++n) { // 600 s at 20 ms
        dc.step(0.02);
        maxSupply = std::max(maxSupply, dc.plant().supplyT());
    }
    CHECK(maxSupply > 35.05); // the step disturbs the supply...
    CHECK(dc.plant().supplyT() == doctest::Approx(35.0).epsilon(0.05 / 35.0)); // ...and the PI recovers it
    CHECK(std::abs(dc.plant().state().residual) < 0.01 * dc.rack().ledger().powerIn); // near steady again
}

TEST_CASE("A more viscous equipment coolant lowers the CDU conductance") {
    SimConfig c = plantConfig(25.0);
    DataCenter a(c.rack, c.server, c.plant, water()), b(c.rack, c.server, c.plant, water());
    b.rack().setViscosityScale(4.0);
    a.solveSteady();
    b.solveSteady();
    CHECK(b.plant().state().cduUA < 0.9 * a.plant().state().cduUA);
    CHECK(b.coolingPower() > a.coolingPower());
}

TEST_CASE("Disabling the plant restores the fixed 35 C supply boundary") {
    SimConfig c = plantConfig(40.0, false);
    DataCenter dc(c.rack, c.server, c.plant, water(), false);
    dc.solveSteady();
    CHECK(dc.rack().supplyT() == doctest::Approx(35.0));
    CHECK(dc.rack().returnT() == doctest::Approx(42.7).epsilon(0.2 / 42.7));
}
