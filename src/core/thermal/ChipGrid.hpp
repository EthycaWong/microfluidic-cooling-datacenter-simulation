#pragma once
#include "fluid/Fluid.hpp"
#include "hydraulics/RectChannel.hpp"
#include "thermal/PowerMap.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mf {

enum class CoolingMode { Microchannel, ColdPlate };

struct DieGeometry {
    double lx = 20e-3;        ///< heated footprint along the flow direction, m
    double ly = 20e-3;        ///< heated footprint across the flow, m
    double thickness = 0.8e-3;///< silicon thickness before etching, m
};

/// Conventional cold plate (copper base with internal fins) clamped on the die through a TIM.
/// Its convective conductance and pressure drop are lumped ASSUMPTIONS, typical of CPU cold plates.
struct ColdPlateParams {
    double timResistance = 8e-6;   ///< TIM area-specific resistance, m^2 K / W (0.08 cm^2 K/W)
    double baseThickness = 3e-3;   ///< copper base, m
    double kCopper = 390.0;        ///< W/(m K)
    double rhoCpCopper = 3.45e6;   ///< J/(m^3 K)
    double hARef = 25.0;           ///< total fin-side conductance at qRef, W/K
    double qRefLpm = 0.5;
    double hAExponent = 0.5;       ///< hA ~ Q^n
    double fluidVolume = 3e-6;     ///< coolant volume inside the plate, m^3
    double dpRefPa = 4000.0;       ///< pressure drop at qRef with reference water, Pa
    double laminarFraction = 0.6;  ///< share of dpRef that scales with viscosity
};

struct PowerSpec {
    double idleW = 60.0;
    double fullW = 300.0;
    PowerMapSpec map;
    /// P(u) = P_idle + u (P_full - P_idle); a linear ASSUMPTION needing calibration.
    double at(double u) const { return idleW + u * (fullW - idleW); }
};

struct ChipParams {
    std::string name = "GPU";
    CoolingMode mode = CoolingMode::Microchannel;
    DieGeometry die;
    ChannelGeometry channels;
    ColdPlateParams coldPlate;
    PowerSpec power;
    int nx = 10, ny = 10;
    double kSi300 = 148.0;         ///< silicon conductivity at 300 K, W/(m K)
    double kSiExponent = 1.3;      ///< k(T) = k300 (300/T)^n
    double rhoCpSi = 2329.0 * 705.0; ///< J/(m^3 K)
};

/// Per-step energy accounting for one chip.
struct ChipLedger {
    double powerIn = 0;      ///< W
    double toCoolant = 0;    ///< W, sum over columns of mdot cp (T_out - T_in)
    double storageRate = 0;  ///< W, dE/dt of all nodes
    double residual = 0;     ///< W, powerIn - toCoolant - storageRate
    double relResidual() const { return powerIn > 0 ? residual / powerIn : 0.0; }
};

/// Thermal resistance-capacitance network of one chip on an nx-by-ny grid (HotSpot-style).
///
/// Layers per cell: silicon (temperature of the heated active face), then for cold plates a copper
/// base node, then a coolant node. Coolant flows along +x; each grid column j carries the channels
/// (or plate passages) of its strip with first-order upwind advection, so coolant warms along the
/// flow. Integration is backward Euler (unconditionally stable) with a Gauss-Seidel solve.
/// Temperature-dependent coefficients (silicon k, coolant properties, h) are lagged one step.
class ChipGrid {
public:
    ChipGrid(ChipParams p, std::shared_ptr<const Fluid> fluid, double initialT = 35.0);

    void setFluid(std::shared_ptr<const Fluid> f) { fluid_ = std::move(f); }
    void setMassFlow(double mdot) { mdot_ = mdot; }   ///< kg/s through the whole chip
    void setInletT(double tC) { tin_ = tC; }
    void setActivity(double u) { power_ = p_.power.at(u); activity_ = u; }
    void setPower(double w) { power_ = w; }

    void step(double dt);
    /// Iterate to steady state for the current flow, inlet temperature and power.
    void solveSteady(double tolK = 1e-7);

    // ---- results ----
    const ChipParams& params() const { return p_; }
    int nx() const { return p_.nx; }
    int ny() const { return p_.ny; }
    double siliconT(int i, int j) const { return T_[idx(i, j, 0)]; }
    double copperT(int i, int j) const { return hasCopper() ? T_[idx(i, j, 1)] : siliconT(i, j); }
    double fluidT(int i, int j) const { return T_[idx(i, j, layers_ - 1)]; }
    double peakT() const { return peak_; }
    int peakI() const { return peakI_; }
    int peakJ() const { return peakJ_; }
    double meanSiliconT() const;
    double outletT() const;      ///< mixed outlet (columns carry equal flow)
    double meanFluidT() const;
    double inletT() const { return tin_; }
    double massFlow() const { return mdot_; }
    double power() const { return power_; }
    double activity() const { return activity_; }
    double storedEnergy() const;  ///< J relative to 0 C
    const ChipLedger& ledger() const { return ledger_; }
    double heatTransferCoeff(int i) const { return hRow_.empty() ? 0.0 : hRow_[i]; } ///< W/m^2K
    double cellPower(int i, int j) const { return power_ * weights_[j * p_.nx + i]; }
    /// Volumetric flow through the chip at the mean coolant temperature, m^3/s.
    double volumeFlow() const;
    /// Silicon footprint conductance to coolant summed over the chip (W/K), for diagnostics.
    double totalFluidConductance() const { return gTotal_; }
    int gsSweeps() const { return lastSweeps_; }

private:
    bool hasCopper() const { return p_.mode == CoolingMode::ColdPlate; }
    int idx(int i, int j, int l) const { return (j * p_.nx + i) * layers_ + l; }
    void build(double dt);
    int gaussSeidel(double omega, double tol, int maxSweeps);
    void updatePeak();

    struct Row { double diag = 0, rhs = 0; int n = 0; int col[6]; double val[6]; };
    void addG(int a, int b, double g);

    ChipParams p_;
    std::shared_ptr<const Fluid> fluid_;
    RectChannel channel_;
    int layers_;
    std::vector<double> T_, Told_, C_, weights_, hRow_;
    std::vector<Row> rows_;
    double mdot_ = 0.0, tin_ = 35.0, power_ = 0.0, activity_ = 0.0, cpChip_ = 4180.0;
    double peak_ = 35.0, gTotal_ = 0.0;
    int peakI_ = 0, peakJ_ = 0, lastSweeps_ = 0;
    ChipLedger ledger_;
};

} // namespace mf
