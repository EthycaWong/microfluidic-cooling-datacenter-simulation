#pragma once
#include "hydraulics/HydraulicNetwork.hpp"
#include "thermal/AirGroup.hpp"
#include "thermal/ChipGrid.hpp"

#include <array>
#include <memory>
#include <string>

namespace mf {

struct ServerParams {
    ChipParams gpu;                 ///< used for both GPUs (microchannel mode)
    ChipParams cpu;                 ///< cold plate mode
    AirParams air;
    double hoseDiameter = 6e-3;     ///< branch hose inner diameter, m
    double hoseLength = 0.6;        ///< total branch hose length (supply + return side), m
    double hoseK = 3.0;             ///< branch fittings, velocity heads
    double qdDpRefPa = 3000.0;      ///< one quick disconnect at qdQRefLpm (ASSUMED, typical)
    double qdQRefLpm = 1.5;
    bool autoBalance = true;        ///< size balancing valves for equal branch flow at baseline
    double branchTargetLpm = 0.5;
    double balanceT = 35.0;         ///< isothermal temperature at which valves are sized
};

/// Energy accounting for the liquid-cooled part of a server (or rack).
struct LiquidLedger {
    double powerIn = 0, toCoolant = 0, storageRate = 0;
    double residual() const { return powerIn - toCoolant - storageRate; }
    double relResidual() const { return powerIn > 0 ? residual() / powerIn : 0.0; }
    LiquidLedger& operator+=(const LiquidLedger& o) {
        powerIn += o.powerIn; toCoolant += o.toCoolant; storageRate += o.storageRate;
        return *this;
    }
};

/// One 4U server: supply connection -> three parallel branches -> return connection.
///   branch 0: GPU 1 silicon microchannels
///   branch 1: GPU 2 silicon microchannels
///   branch 2: CPU conventional cold plate
/// Each branch is hose + component + balancing valve in series. Branch flows come from the
/// hydraulic network (equal pressure drop across parallel branches), never assumed equal.
class Server {
public:
    static constexpr int kBranches = 3;
    Server(const ServerParams& p, std::shared_ptr<const Fluid> fluid, int index);

    /// Create this server's nodes and edges between two manifold nodes.
    void attach(HydraulicNetwork& net, int supplyNode, int returnNode, double z);
    /// Size balancing valves so each branch carries branchTargetLpm at the same pressure drop
    /// (isothermal at balanceT). Requires attach() first.
    void autoBalance(HydraulicNetwork& net);
    void clearBalance();
    void setBalanceEnabled(bool on) { p_.autoBalance = on; }
    /// Push temperatures to edges (for mu(T), rho(T)) and pull branch flows into the chips.
    void updateEdgeTemperatures(HydraulicNetwork& net, double supplyT);
    void pullFlows(const HydraulicNetwork& net);

    void setActivity(double u);
    void step(double dt, double supplyT);
    void solveSteady(double supplyT);

    // ---- access / results ----
    int index() const { return index_; }
    ChipGrid& chip(int b) { return b < 2 ? gpu_[b] : cpu_; }
    const ChipGrid& chip(int b) const { return b < 2 ? gpu_[b] : cpu_; }
    const AirGroup& air() const { return air_; }
    double activity() const { return u_; }
    double branchFlow(int b) const { return branchMdot_[b]; }         ///< kg/s
    double branchDp(const HydraulicNetwork& net, int b) const { return net.edgeDp(branchEdge_[b]); }
    double totalFlow() const { return branchMdot_[0] + branchMdot_[1] + branchMdot_[2]; }
    double serverDp(const HydraulicNetwork& net) const;  ///< supply connector to return connector, Pa
    double returnT() const { return returnT_; }
    double peakT() const;
    int hottestBranch() const;
    double liquidPower() const;
    double totalPower() const { return liquidPower() + air_.power(); }
    LiquidLedger ledger() const;
    double valveDpRef(int b) const { return valveDp_[b]; }
    int qdSupplyEdge() const { return qdIn_; }
    int qdReturnEdge() const { return qdOut_; }
    int branchEdge(int b) const { return branchEdge_[b]; }
    std::string branchName(int b) const;

private:
    void updateReturnT(double supplyT);

    ServerParams p_;
    std::shared_ptr<const Fluid> fluid_;
    int index_;
    std::array<ChipGrid, 2> gpu_;
    ChipGrid cpu_;
    AirGroup air_;
    double u_ = 0.0;
    std::array<double, kBranches> branchMdot_{};
    std::array<double, kBranches> valveDp_{};
    std::array<int, kBranches> branchEdge_{-1, -1, -1};
    std::array<RefLoss*, kBranches> valve_{};
    int qdIn_ = -1, qdOut_ = -1, nodeS_ = -1, nodeR_ = -1;
    double returnT_ = 35.0;
};

/// Stand-alone server test rig: fixed supply temperature and fixed total flow (a temporary
/// boundary condition until the rack pump and CDU are attached).
class ServerRig {
public:
    ServerRig(const ServerParams& p, std::shared_ptr<const Fluid> fluid, double totalFlowLpm, double supplyT);
    void setTotalFlowLpm(double lpm);
    void step(double dt);
    void solveSteady();
    Server& server() { return server_; }
    const HydraulicNetwork& network() const { return *net_; }
    double supplyT() const { return supplyT_; }
    double totalDp() const { return server_.serverDp(*net_); }
    double volumeFlow() const;     ///< m^3/s at supply temperature
    /// Indicative pump electric power if one pump supplied only this server's pressure drop.
    double indicativePumpPower(double etaPump = 0.6, double etaMotor = 0.9) const;

private:
    void solveHydraulics();
    std::shared_ptr<const Fluid> fluid_;
    std::unique_ptr<HydraulicNetwork> net_;
    Server server_;
    int sup_ = -1, ret_ = -1;
    double flowLpm_, supplyT_;
};

} // namespace mf
