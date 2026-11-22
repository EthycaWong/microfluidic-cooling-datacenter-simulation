#pragma once
#include "fluid/Fluid.hpp"

#include <memory>
#include <string>

namespace mf {

/// Effectiveness of heat exchangers (epsilon-NTU method).
namespace hx {
/// Pure counterflow. cr = Cmin/Cmax in [0, 1].
double effectivenessCounterflow(double ntu, double cr);
/// Crossflow, Cmax mixed / Cmin unmixed (finned coil: air mixed, water unmixed when air is Cmax),
/// and the converse, chosen automatically from which stream is Cmin.
double effectivenessCrossflow(double ntu, double cr, bool cminIsMixed);
} // namespace hx

/// CDU liquid-to-liquid plate heat exchanger (equipment loop hot side, facility loop cold side).
/// Each side's conductance follows a Dittus-Boelter-type scaling from a design point:
///   hA = hA_ref (m/m_ref)^n (mu/mu_ref)^(1/3 - n) (k/k_ref)^(2/3) (cp/cp_ref)^(1/3)
/// so a more viscous coolant lowers hA (it lowers Re and raises Pr).
struct CduParams {
    double hAEquipRef = 2800.0;      ///< W/K, equipment side at design flow (ASSUMED sizing)
    double hAFacilityRef = 2800.0;   ///< W/K, facility side at design flow
    double equipRefLpm = 12.0, facilityRefLpm = 14.0;
    double exponent = 0.6;           ///< plate HX Nusselt ~ Re^0.6 (ASSUMED)
    double wallR = 0.0;              ///< K/W, plate wall (thin stainless, negligible by default)
    double facilityDpRefPa = 20e3;   ///< facility-side pressure drop at facilityRefLpm
};

/// Outdoor dry cooler: finned-tube coil with variable-speed fans, no evaporation.
struct DryCoolerParams {
    double uaRef = 2000.0;           ///< W/K at full fan speed and design water flow (ASSUMED sizing)
    double airShare = 0.75;          ///< fraction of total resistance on the air side at design
    double airExponent = 0.6;        ///< air-side hA ~ V^0.6
    double waterExponent = 0.8;      ///< water-side hA ~ m^0.8
    double maxAirflow = 2.0;         ///< m^3/s at full fan speed
    double maxFanPower = 450.0;      ///< W at full speed (fan power ~ speed^3)
    double minFanSpeed = 0.1;
    double coilDpRefPa = 30e3;       ///< water-side pressure drop at the facility design flow
    double rhoAir = 1.16, cpAir = 1007.0;
};

/// Optional trim chiller after the dry cooler, used only when fans at 100 % cannot hold the setpoint.
struct ChillerParams {
    bool enabled = true;
    double capacityW = 10e3;
    double carnotFraction = 0.5;     ///< COP = f * T_evap / (T_cond - T_evap)  (ASSUMED)
    double evapApproachK = 5.0;      ///< evaporator below the chilled-water target
    double condApproachK = 10.0;     ///< condenser above outdoor air
    double maxCop = 10.0;
};

struct FacilityParams {
    double flowLpm = 14.0;           ///< facility water flow (fixed-speed pump)
    double pipingDpRefPa = 10e3;     ///< piping + valves at flowLpm
    double pumpEta = 0.5, motorEta = 0.9;
};

/// Outdoor dry-bulb temperature: constant, or a daily sine wave.
struct WeatherParams {
    double meanC = 25.0;
    double amplitudeK = 0.0;         ///< 0 = constant
    double periodS = 86400.0;
    double phaseS = 0.0;             ///< time of the minimum (t = phase)
    double at(double t) const;
};

struct PlantParams {
    CduParams cdu;
    DryCoolerParams dryCooler;
    ChillerParams chiller;
    FacilityParams facility;
    WeatherParams weather;
    double supplySetpointC = 35.0;   ///< equipment-loop supply target (fans controlled to hold it)
    double equipSupplyVolumeL = 6.0; ///< coolant between CDU outlet and the servers (thermal mass)
    double equipReturnVolumeL = 8.0; ///< coolant between the servers and the CDU inlet
    double facilityVolumeL = 25.0;   ///< facility loop inventory (half before, half after the dry cooler)
    double fanKp = 0.08, fanKi = 0.0015; ///< PI on supply temperature error (tuned on a 30%->100% load step)
};

/// Everything the plant reports for one instant.
struct PlantState {
    double equipSupplyT = 35, equipReturnT = 35;   ///< at the rack (C)
    double cduHotInT = 35, cduHotOutT = 35;        ///< equipment side of the CDU
    double facilitySupplyT = 30, facilityReturnT = 30; ///< facility side into / out of the CDU
    double dryCoolerOutT = 30;
    double outdoorT = 25;
    double cduHeat = 0, dryCoolerHeat = 0, chillerHeat = 0; ///< W
    double cduUA = 0, cduEffectiveness = 0, dryCoolerUA = 0, dryCoolerEffectiveness = 0;
    double fanSpeed = 1, fanPower = 0, chillerPower = 0, chillerCop = 0, facilityPumpPower = 0;
    double residual = 0;  ///< W, heat into loop - heat out - storage (exact only at steady state)
    bool freezeRisk = false; ///< outdoor air within 2 K of the facility fluid's freezing point (dry cooler outdoors)
};

/// CDU + facility loop + dry cooler (+ trim chiller), coupled to the rack through the equipment
/// loop's mass flow, return temperature and supply temperature.
class CoolingPlant {
public:
    CoolingPlant(const PlantParams& p, std::shared_ptr<const Fluid> equipFluid,
                 std::shared_ptr<const Fluid> facilityFluid);

    /// Advance by dt with the rack's current loop flow (kg/s) and mixed return temperature.
    void step(double dt, double time, double equipMdot, double rackReturnT);
    /// Steady state for a given rack heat load (W) and loop flow. Fans are set to hold the
    /// setpoint (bisection); the chiller trims if fans at 100 % are not enough.
    void solveSteady(double time, double equipMdot, double rackHeatW);

    const PlantState& state() const { return s_; }
    double supplyT() const { return s_.equipSupplyT; }
    const PlantParams& params() const { return p_; }
    void setOutdoorMean(double c) { p_.weather.meanC = c; }
    void setSetpoint(double c) { p_.supplySetpointC = c; }
    void setChillerEnabled(bool on) { p_.chiller.enabled = on; }
    void setFanSpeed(double s);
    void setFacilityFluid(std::shared_ptr<const Fluid> f) { ff_ = std::move(f); }
    const Fluid& facilityFluid() const { return *ff_; }
    double facilityFreezingPointC() const;
    double facilityMdot() const;
    double coolingPower() const { return s_.fanPower + s_.chillerPower + s_.facilityPumpPower; }

    /// CDU conductance and effectiveness for the given flows and side temperatures.
    void cduPerformance(double mEquip, double tEquip, double tFac, double& ua, double& eps) const;
    void dryCoolerPerformance(double fanSpeed, double tWater, double tAir, double& ua, double& eps) const;

private:
    double chillerCop(double tChilled, double tOutdoor) const;
    double facilityPumpPower() const;

    PlantParams p_;
    std::shared_ptr<const Fluid> fe_, ff_;
    PlantState s_;
    double fanSpeed_ = 1.0, lastErr_ = 0.0;
    // lumped volumes: equipment supply / return, facility supply (after dry cooler) / return
    double tES_, tER_, tFS_, tFR_;
};

} // namespace mf
