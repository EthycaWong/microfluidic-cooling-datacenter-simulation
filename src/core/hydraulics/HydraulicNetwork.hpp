#pragma once
#include "hydraulics/Edges.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mf {

/// Steady incompressible flow network solved for node pressures.
///
/// Unknowns are the static pressures at every node except the reference node (where a sealed
/// expansion vessel would fix the pressure). Each edge's mass flow follows from its pressure
/// difference by inverting its monotonic dp(mdot) law; Newton-Raphson drives the mass imbalance
/// at every node to zero. Mass flow (not volume flow) is conserved, so temperature-dependent
/// density cannot create or destroy coolant.
///
/// Hydrostatics: p_from - p_to = dp(mdot) + rho_h * g * (z_to - z_from), with one constant
/// density rho_h for the whole loop. Elevation changes therefore cancel around any closed loop
/// and do not change flows, while local static pressures still vary with height.
class HydraulicNetwork {
public:
    explicit HydraulicNetwork(std::shared_ptr<const Fluid> fluid);

    int addNode(const std::string& name, double z = 0.0);
    int addEdge(int from, int to, std::unique_ptr<Edge> e);
    void setReference(int node, double pressurePa);
    void setInjection(int node, double mdot); ///< external mass inflow, kg/s (negative = outflow)
    void setHydrostaticDensity(double rho) { rhoH_ = rho; }
    void setFluid(std::shared_ptr<const Fluid> f) { fluid_ = std::move(f); }

    /// Solve; returns true when the largest nodal mass imbalance is below tolerance.
    bool solve();

    int nodeCount() const { return static_cast<int>(nodes_.size()); }
    int edgeCount() const { return static_cast<int>(edges_.size()); }
    double pressure(int n) const { return nodes_[n].p; }
    double nodeZ(int n) const { return nodes_[n].z; }
    const std::string& nodeName(int n) const { return nodes_[n].name; }
    double flow(int e) const { return edges_[e].mdot; }        ///< kg/s, + means from -> to
    double edgeDp(int e) const;                                ///< loss across edge, Pa
    Edge& edge(int e) { return *edges_[e].e; }
    const Edge& edge(int e) const { return *edges_[e].e; }
    int edgeFrom(int e) const { return edges_[e].from; }
    int edgeTo(int e) const { return edges_[e].to; }
    double maxImbalance() const { return lastImbalance_; }     ///< kg/s
    int lastIterations() const { return lastIter_; }
    const Fluid& fluid() const { return *fluid_; }

private:
    struct Node { std::string name; double z; double p = 0.0; double inj = 0.0; };
    struct E { int from, to; std::unique_ptr<Edge> e; double mdot = 0.0; double g = 0.0; };

    double invert(E& e, double delta) const; ///< mdot with dp(mdot) = delta; also stores dmdot/ddelta
    void computeFlows(const std::vector<double>& p);
    double imbalance(std::vector<double>& r) const;

    std::shared_ptr<const Fluid> fluid_;
    std::vector<Node> nodes_;
    std::vector<E> edges_;
    int ref_ = -1;
    double rhoH_ = 994.0;
    double lastImbalance_ = 0.0;
    int lastIter_ = 0;
};

} // namespace mf
