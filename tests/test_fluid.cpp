#include <doctest/doctest.h>
#include "fluid/Fluid.hpp"

using namespace mf;

TEST_CASE("IAPWS-IF97 Region 1 matches the release's verification values") {
    // R7-97(2012) Table 5: T = 300 K, p = 3 MPa and p = 80 MPa
    CHECK(WaterIAPWS::if97SpecificVolume(300.0, 3e6) == doctest::Approx(0.100215168e-2).epsilon(1e-8));
    CHECK(WaterIAPWS::if97Cp(300.0, 3e6) == doctest::Approx(4173.01218).epsilon(1e-8));
    CHECK(WaterIAPWS::if97SpecificVolume(300.0, 80e6) == doctest::Approx(0.971180894e-3).epsilon(1e-8));
    CHECK(WaterIAPWS::if97Cp(300.0, 80e6) == doctest::Approx(4010.08987).epsilon(1e-8));
}

TEST_CASE("IAPWS 2008 viscosity matches the release's check values") {
    // R12-08 Table 4 (mu0 * mu1, no critical enhancement needed at these states)
    CHECK(WaterIAPWS::viscosity2008(298.15, 998.0) * 1e6 == doctest::Approx(889.735100).epsilon(1e-8));
    CHECK(WaterIAPWS::viscosity2008(298.15, 1200.0) * 1e6 == doctest::Approx(1437.649467).epsilon(1e-8));
}

TEST_CASE("IAPWS 2011 thermal conductivity matches the release's check values") {
    // R15-11 Table 4
    CHECK(WaterIAPWS::conductivity2011(298.15, 998.0) * 1e3 == doctest::Approx(607.712868).epsilon(1e-6));
    CHECK(WaterIAPWS::conductivity2011(298.15, 1200.0) * 1e3 == doctest::Approx(799.038144).epsilon(1e-6));
}

TEST_CASE("Tabulated water properties agree with direct evaluation") {
    WaterIAPWS w;
    for (double t : {5.0, 20.0, 35.05, 43.7, 60.0, 85.0}) {
        const FluidProps a = w.props(t), b = w.propsExact(t);
        CHECK(a.rho == doctest::Approx(b.rho).epsilon(1e-6));
        CHECK(a.mu == doctest::Approx(b.mu).epsilon(1e-4));
        CHECK(a.cp == doctest::Approx(b.cp).epsilon(1e-5));
        CHECK(a.k == doctest::Approx(b.k).epsilon(1e-5));
    }
}

TEST_CASE("ScaledFluid changes viscosity only") {
    auto w = std::make_shared<WaterIAPWS>();
    ScaledFluid s(w, 4.0);
    const FluidProps a = w->props(35.0), b = s.props(35.0);
    CHECK(b.mu == doctest::Approx(4.0 * a.mu));
    CHECK(b.rho == a.rho);
    CHECK(b.cp == a.cp);
    CHECK(b.k == a.k);
}

TEST_CASE("Melinder glycol correlation reduces to water at zero concentration") {
    WaterIAPWS w;
    GlycolMixture g(GlycolMixture::Kind::Propylene, 0.0), e(GlycolMixture::Kind::Ethylene, 0.0);
    for (double t : {20.0, 35.0, 45.0}) {
        for (const GlycolMixture* m : {&g, &e}) {
            const FluidProps a = m->propsExact(t), b = w.propsExact(t);
            CHECK(a.rho == doctest::Approx(b.rho).epsilon(2e-3));
            CHECK(a.mu == doctest::Approx(b.mu).epsilon(2.5e-2));
            CHECK(a.cp == doctest::Approx(b.cp).epsilon(5e-3));
            CHECK(a.k == doctest::Approx(b.k).epsilon(1e-2));
        }
    }
}

TEST_CASE("Glycol freezing points, concentration parsing and property trends") {
    auto pg25 = std::dynamic_pointer_cast<const GlycolMixture>(makeCoolant("PG25"));
    auto pg30m = std::dynamic_pointer_cast<const GlycolMixture>(makeCoolant("PG30m"));
    auto eg30m = std::dynamic_pointer_cast<const GlycolMixture>(makeCoolant("EG30m"));
    REQUIRE(pg25);
    CHECK(pg25->massFraction() == doctest::Approx(0.257).epsilon(0.01)); // 25 % by volume
    CHECK(pg30m->massFraction() == doctest::Approx(0.30));
    CHECK(pg30m->freezingPointC() == doctest::Approx(-12.8).epsilon(0.1)); // ~ -12.7 C in handbook tables
    CHECK(eg30m->freezingPointC() < pg30m->freezingPointC());             // EG depresses freezing more
    WaterIAPWS w;
    const FluidProps a = pg25->props(40.0), b = w.props(40.0);
    CHECK(a.mu > 2.0 * b.mu);
    CHECK(a.k < 0.8 * b.k);
    CHECK(a.cp < b.cp);
    CHECK(a.rho > b.rho);
    CHECK_THROWS(makeCoolant("brine"));
    CHECK(makeCoolant("water")->name() == "Water (IAPWS)");
}
