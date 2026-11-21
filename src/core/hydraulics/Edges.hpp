#pragma once
#include "fluid/Fluid.hpp"
#include "hydraulics/RectChannel.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mf {

/// A two-port hydraulic element. dp(mdot) is the pressure LOSS from inlet to outlet at mass flow
/// mdot (kg/s); a pump returns a negative loss (a pressure rise). dp must increase monotonically
/// with mdot so the network problem has a unique solution.
class Edge {
public:
    explicit Edge(std::string name) : name_(std::move(name)) {}
    virtual ~Edge() = default;

    virtual double dp(double mdot, const Fluid& f) const = 0;
    /// d(dp)/d(mdot) by central difference (overridable).
    virtual double ddp(double mdot, const Fluid& f) const;
    /// Set the representative fluid temperature (C). Composite edges forward it.
    virtual void setTemperature(double tC) { tC_ = tC; }
    double temperature() const { return tC_; }
    const std::string& name() const { return name_; }

protected:
    std::string name_;
    double tC_ = 35.0;
};

/// Loss defined by one measured/assumed operating point and a laminar fraction:
///   dp = dpRef * [ f_lam * (Q/Qref) * (mu/muRef) + (1 - f_lam) * (Q|Q|/Qref^2) * (rho/rhoRef) ]
/// The laminar part scales with viscosity, the quadratic (inertial) part does not.
/// Used for quick disconnects, cold plates, valves, filters and the CDU placeholder.
class RefLoss : public Edge {
public:
    RefLoss(std::string name, double dpRefPa, double qRefM3s, double laminarFraction,
            const FluidProps& refProps);
    double dp(double mdot, const Fluid& f) const override;
    void setDpRef(double dpRefPa) { dpRef_ = dpRefPa; }
    double dpRef() const { return dpRef_; }
    double qRef() const { return qRef_; }
    double laminarFraction() const { return fLam_; }

private:
    double dpRef_, qRef_, fLam_, muRef_, rhoRef_;
};

/// Straight pipe with Churchill (1977) all-regime friction factor plus lumped fitting losses.
class PipeEdge : public Edge {
public:
    PipeEdge(std::string name, double diameter, double length, double kFittings = 0.0,
             double roughness = 1.5e-6);
    double dp(double mdot, const Fluid& f) const override;
    static double churchillDarcy(double re, double relRoughness);
    double diameter() const { return d_; }

private:
    double d_, l_, k_, eps_;
};

/// Silicon microchannel array (wraps RectChannel).
class ChannelEdge : public Edge {
public:
    ChannelEdge(std::string name, const ChannelGeometry& g) : Edge(std::move(name)), ch_(g) {}
    double dp(double mdot, const Fluid& f) const override;
    const RectChannel& channel() const { return ch_; }

private:
    RectChannel ch_;
};

/// Centrifugal pump with quadratic head curve and affinity-law speed scaling:
///   H(Q, s) = s^2 * H0 - b * Q|Q|            (Pa)
///   eta(Q, s) = etaMax * (1 - ((Q/s - Qbep)/Qbep)^2), clamped
/// Electric power = H * Q / (eta_pump * eta_motor).
struct PumpCurve {
    double shutoffHeadPa = 100e3;        ///< H0 at s = 1, Q = 0
    double maxFlowM3s = 20.0 / 60000.0;  ///< flow where H = 0 at s = 1 (defines b)
    double bepFlowM3s = 12.0 / 60000.0;  ///< best-efficiency flow at s = 1
    double etaMax = 0.50;                ///< peak hydraulic efficiency
    double etaMotor = 0.90;              ///< motor + drive efficiency
    double minSpeed = 0.05, maxSpeed = 2.0;
};

class PumpEdge : public Edge {
public:
    PumpEdge(std::string name, const PumpCurve& c) : Edge(std::move(name)), c_(c) {}
    double dp(double mdot, const Fluid& f) const override;
    void setSpeed(double s) { speed_ = s; }
    double speed() const { return speed_; }
    double head(double q) const;                  ///< Pa
    double efficiency(double q) const;            ///< hydraulic efficiency
    double electricPower(double q) const;         ///< W
    const PumpCurve& curve() const { return c_; }

private:
    PumpCurve c_;
    double speed_ = 1.0;
};

/// Several elements carrying the same mass flow in series; losses add.
class SeriesEdge : public Edge {
public:
    explicit SeriesEdge(std::string name) : Edge(std::move(name)) {}
    Edge* add(std::unique_ptr<Edge> e) {
        parts_.push_back(std::move(e));
        return parts_.back().get();
    }
    double dp(double mdot, const Fluid& f) const override;
    void setTemperature(double tC) override;
    const std::vector<std::unique_ptr<Edge>>& parts() const { return parts_; }

private:
    std::vector<std::unique_ptr<Edge>> parts_;
};

} // namespace mf
