#pragma once

namespace mf {

/// Memory, storage, NICs and VRMs of one server lumped as a single air-cooled node.
/// Heat leaves to room air (never to the coolant). Fans are controlled to hold an exhaust air
/// temperature rise; fan power follows a quadratic system curve (P_fan ~ V^3).
struct AirParams {
    double idleW = 120.0, fullW = 200.0;     ///< ASSUMED power range of the combined group
    double heatCapacity = 2000.0;            ///< J/K (boards, DIMMs, heatsinks)
    double uaRef = 10.0;                     ///< W/K component-to-air at the reference airflow
    double roomT = 27.0;                     ///< C, server inlet air
    double airDeltaT = 12.0;                 ///< K, fan-control target exhaust rise
    double fanDpRef = 150.0;                 ///< Pa at the reference airflow
    double flowRef = 0.0143;                 ///< m^3/s reference airflow (~30 CFM)
    double fanEta = 0.35;                    ///< fan + motor efficiency
    double minFlowFrac = 0.2;                ///< fans never drop below this fraction of flowRef
    double rhoAir = 1.16, cpAir = 1007.0;
};

class AirGroup {
public:
    explicit AirGroup(const AirParams& p = {}) : p_(p), t_(p.roomT) {}
    void setActivity(double u) { power_ = p_.idleW + u * (p_.fullW - p_.idleW); }
    void step(double dt);
    void solveSteady();

    double power() const { return power_; }
    double temperature() const { return t_; }
    double airflow() const { return flow_; }       ///< m^3/s
    double fanPower() const;                       ///< W
    double heatToAir() const { return toAir_; }     ///< W this step
    const AirParams& params() const { return p_; }

private:
    void updateFlow();
    AirParams p_;
    double power_ = 0.0, t_, flow_ = 0.0, toAir_ = 0.0;
};

} // namespace mf
