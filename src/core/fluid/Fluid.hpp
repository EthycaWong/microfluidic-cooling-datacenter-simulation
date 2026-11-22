#pragma once
#include <memory>
#include <string>
#include <vector>

namespace mf {

/// Thermophysical properties of a single-phase liquid at one state.
struct FluidProps {
    double rho; ///< density, kg/m^3
    double mu;  ///< dynamic viscosity, Pa s
    double cp;  ///< isobaric specific heat, J/(kg K)
    double k;   ///< thermal conductivity, W/(m K)
    double pr() const { return mu * cp / k; }
};

/// A coolant whose properties depend on temperature (pressure is fixed at the loop level).
class Fluid {
public:
    virtual ~Fluid() = default;
    virtual FluidProps props(double tC) const = 0;
    virtual std::string name() const = 0;
};

/// Pure water from IAPWS formulations, evaluated at a fixed loop pressure:
///  - density and cp: IAPWS-IF97 Region 1 (IAPWS R7-97(2012))
///  - viscosity:      IAPWS 2008 (R12-08), mu0 * mu1 terms, no critical enhancement
///  - conductivity:   IAPWS 2011 (R15-11), lambda0 * lambda1 terms, no critical enhancement
/// These pure-water formulations must not be used for glycol mixtures.
/// props() interpolates a table built once at construction (0.1 K spacing, 1-120 C).
class WaterIAPWS final : public Fluid {
public:
    explicit WaterIAPWS(double pressurePa = 0.3e6);
    FluidProps props(double tC) const override;
    std::string name() const override { return "Water (IAPWS)"; }

    /// Direct (non-tabulated) evaluation, exposed for verification against IAPWS check values.
    FluidProps propsExact(double tC) const;

    // Individual formulations (SI units, T in kelvin), public for unit tests.
    static double if97SpecificVolume(double T, double p);      ///< m^3/kg
    static double if97Cp(double T, double p);                  ///< J/(kg K)
    static double viscosity2008(double T, double rho);         ///< Pa s
    static double conductivity2011(double T, double rho);      ///< W/(m K)

private:
    double p_;
    static constexpr double kTmin = 1.0, kTmax = 120.0, kStep = 0.1;
    std::vector<FluidProps> table_;
};

/// Aqueous propylene- or ethylene-glycol coolant from the Melinder (2010) correlations (IIR
/// "Properties of Secondary Working Fluids for Indirect Systems"), using the polynomial
/// coefficients distributed by CoolProp. Mixture-specific data: pure-water IAPWS formulations are
/// never substituted. Valid for mass fraction 0-0.6 and temperatures above the freezing point.
class GlycolMixture final : public Fluid {
public:
    enum class Kind { Propylene, Ethylene };
    GlycolMixture(Kind kind, double massFraction);
    FluidProps props(double tC) const override;
    std::string name() const override;
    FluidProps propsExact(double tC) const;
    double freezingPointC() const;
    double massFraction() const { return x_; }
    Kind kind() const { return kind_; }
    /// Mass fraction from a volume fraction (how coolants are usually specified, e.g. PG25).
    static double massFromVolumeFraction(Kind kind, double volumeFraction);

private:
    Kind kind_;
    double x_;
    static constexpr double kTmin = -30.0, kTmax = 100.0, kStep = 0.1;
    std::vector<FluidProps> table_;
};

/// Plain IAPWS water at tC, cached. Component curves (quick disconnects, cold plates, HX, filter)
/// are specified for water; this is their reference state, whatever coolant is actually loaded.
FluidProps referenceWater(double tC);

/// "water", "PG25", "EG30" (volume %), or "PG25m" (mass %).
std::shared_ptr<const Fluid> makeCoolant(const std::string& spec);

/// Wraps another fluid and multiplies only its viscosity by a factor (1x, 2x, 4x ...).
/// Used for the "isolate viscosity" experiment: density, cp and k stay those of the base fluid.
class ScaledFluid final : public Fluid {
public:
    ScaledFluid(std::shared_ptr<const Fluid> base, double muScale)
        : base_(std::move(base)), scale_(muScale) {}
    FluidProps props(double tC) const override {
        FluidProps p = base_->props(tC);
        p.mu *= scale_;
        return p;
    }
    std::string name() const override { return base_->name() + " x" + std::to_string(scale_) + " mu"; }
    void setScale(double s) { scale_ = s; }
    void setBase(std::shared_ptr<const Fluid> b) { base_ = std::move(b); }
    const Fluid& base() const { return *base_; }
    double scale() const { return scale_; }

private:
    std::shared_ptr<const Fluid> base_;
    double scale_;
};

} // namespace mf
