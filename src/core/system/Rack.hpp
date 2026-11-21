#pragma once
#include "hydraulics/HydraulicNetwork.hpp"
#include "system/Server.hpp"
#include "workload/Workload.hpp"

#include <memory>
#include <vector>

namespace mf {

enum class ManifoldLayout { DirectReturn, ReverseReturn };
enum class PumpMode { FixedSpeed, FixedFlow, PeakTemperature };

const char* toString(ManifoldLayout l);
const char* toString(PumpMode m);

struct RackParams {
    int servers = 8;
    int serverU = 4;
    int firstServerU = 7;            ///< U slot (1-based) of the lowest server; U1-U6 hold CDU + pump
    int totalU = 42;
    double manifoldDiameter = 20e-3; ///< supply/return manifold inner diameter, m
    double manifoldKPerTap = 0.5;    ///< tee run loss per manifold segment, velocity heads
    double connectLength = 0.5;      ///< pipe from CDU/pump to the manifold base, m (each side)
    double connectK = 1.5;           ///< elbows on that pipe
    // CDU placeholder (until the CDU heat exchanger milestone): fixed-curve resistances.
    double hxDpRefPa = 25e3, hxQRefLpm = 12.0, hxLaminarFraction = 0.5;
    double filterDpRefPa = 8e3, filterQRefLpm = 12.0, filterLaminarFraction = 0.3;
    double referencePressurePa = 150e3; ///< sealed expansion vessel at the pump suction
    PumpCurve pump;
    ManifoldLayout layout = ManifoldLayout::DirectReturn;
    PumpMode pumpMode = PumpMode::FixedFlow;
    double pumpSpeed = 1.0;          ///< used in FixedSpeed mode
    double targetFlowLpm = 12.0;     ///< FixedFlow mode
    double targetPeakC = 85.0;       ///< PeakTemperature mode (study assumption, not a hardware limit)
    double piKp = 0.02, piKi = 0.1;  ///< PI gains on pump speed per K and per K s
    double supplyT = 35.0;           ///< fixed supply boundary until the CDU + dry cooler milestone
    double hydraulicInterval = 0.05; ///< s between thermal-hydraulic coupling updates
};

/// Hottest silicon location in the rack.
struct Hotspot { int server = 0, branch = 0, i = 0, j = 0; double t = 0.0; };

/// 42U rack: pump -> supply manifold -> 8 servers in parallel -> return manifold -> CDU -> pump.
class Rack {
public:
    Rack(const RackParams& rp, const ServerParams& sp, std::shared_ptr<const Fluid> baseFluid);

    // ---- control ----
    void setWorkload(int server, std::shared_ptr<Workload> w);
    void setAllWorkloads(const std::shared_ptr<Workload>& w);
    /// Change coolant viscosity; by default the server balancing valves are re-sized for it.
    void setViscosityScale(double s, bool rebalance = true);
    double viscosityScale() const { return fluid_->scale(); }
    void setPumpMode(PumpMode m);
    void setPumpSpeed(double s);
    void setTargetFlowLpm(double lpm) { rp_.targetFlowLpm = lpm; }
    void setTargetPeak(double c) { rp_.targetPeakC = c; }
    void setSupplyT(double c) { rp_.supplyT = c; }
    void setLayout(ManifoldLayout l);
    void setAutoBalance(bool on);

    // ---- simulation ----
    void step(double dt);
    /// Steady state for the current workloads (activity at the current time) and pump mode.
    /// In PeakTemperature mode this finds the pump speed giving targetPeakC (bisection).
    void solveSteady();
    /// Solve hydraulics only (uses current temperatures).
    void solveHydraulics();

    // ---- results ----
    const RackParams& params() const { return rp_; }
    const ServerParams& serverParams() const { return sp_; }
    int serverCount() const { return static_cast<int>(servers_.size()); }
    Server& server(int k) { return servers_[k]; }
    const Server& server(int k) const { return servers_[k]; }
    const HydraulicNetwork& network() const { return *net_; }
    double time() const { return t_; }
    double serverTapZ(int k) const;
    double pumpSpeed() const { return pump_->speed(); }
    double totalMassFlow() const;           ///< kg/s through the pump
    double totalFlowLpm() const;            ///< at supply-temperature density
    double pumpHead() const;                ///< Pa (= loop pressure drop)
    double pumpPower() const;               ///< W electric
    double pumpEfficiency() const;
    double fanPower() const;                ///< W, all servers
    double airHeat() const;                 ///< W rejected to room air
    double itPower() const;                 ///< W, liquid + air groups
    double returnT() const { return returnT_; }
    double supplyT() const { return rp_.supplyT; }
    Hotspot hotspot() const;
    LiquidLedger ledger() const;
    double flowMaldistribution() const;     ///< (max - min) / mean of server flows
    double pumpEnergy() const { return pumpJ_; }   ///< J since start
    double fanEnergy() const { return fanJ_; }
    double itEnergy() const { return itJ_; }
    // manifold segments for display: supply segment k feeds tap k (k = 0 is the riser from the pump)
    double supplySegmentFlow(int k) const { return net_->flow(supplySeg_[k]); }
    double returnSegmentFlow(int k) const { return net_->flow(returnSeg_[k]); }
    double returnSegmentT(int k) const { return returnSegT_[k]; }
    int pressureNodeSupply(int k) const { return supplyNode_[k]; }
    int pressureNodeReturn(int k) const { return returnNode_[k]; }

private:
    void buildNetwork();
    void updateEdgeTemperatures();
    void updateReturnTemperatures();
    void pumpControlTick(double interval);
    void steadyAtCurrentSpeed();

    RackParams rp_;
    ServerParams sp_;
    std::shared_ptr<ScaledFluid> fluid_;
    std::unique_ptr<HydraulicNetwork> net_;
    std::vector<Server> servers_;
    std::vector<std::shared_ptr<Workload>> workloads_;
    PumpEdge* pump_ = nullptr;
    int pumpEdge_ = -1, filterEdge_ = -1, hxEdge_ = -1, pumpIn_ = -1;
    std::vector<int> supplySeg_, returnSeg_, supplyNode_, returnNode_;
    std::vector<double> returnSegT_;
    double t_ = 0.0, nextHyd_ = 0.0, returnT_ = 35.0, lastErr_ = 0.0;
    double pumpJ_ = 0.0, fanJ_ = 0.0, itJ_ = 0.0;
};

} // namespace mf
