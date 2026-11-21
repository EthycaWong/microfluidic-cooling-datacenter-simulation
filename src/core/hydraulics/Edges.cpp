#include "hydraulics/Edges.hpp"

#include <algorithm>
#include <cmath>

namespace mf {
namespace {
/// Q|Q| smoothed near zero so the derivative never vanishes (keeps Newton well conditioned).
double sq(double q, double scale) { return q * std::sqrt(q * q + 1e-8 * scale * scale); }
} // namespace

double Edge::ddp(double mdot, const Fluid& f) const {
    const double h = 1e-7 + 1e-6 * std::abs(mdot);
    return (dp(mdot + h, f) - dp(mdot - h, f)) / (2.0 * h);
}

// ---- RefLoss ---------------------------------------------------------------------------------
RefLoss::RefLoss(std::string name, double dpRefPa, double qRefM3s, double laminarFraction,
                 const FluidProps& ref)
    : Edge(std::move(name)), dpRef_(dpRefPa), qRef_(qRefM3s), fLam_(laminarFraction), muRef_(ref.mu),
      rhoRef_(ref.rho) {}

double RefLoss::dp(double mdot, const Fluid& f) const {
    const FluidProps p = f.props(tC_);
    const double q = mdot / p.rho;
    const double x = q / qRef_;
    const double lam = fLam_ * x * (p.mu / muRef_);
    const double tur = (1.0 - fLam_) * sq(x, 1.0) * (p.rho / rhoRef_);
    // A pure-quadratic element still gets a tiny linear term so dp'(0) > 0.
    return dpRef_ * (lam + tur + 1e-6 * x);
}

// ---- PipeEdge ----------------------------------------------------------------------------------
PipeEdge::PipeEdge(std::string name, double diameter, double length, double kFittings, double roughness)
    : Edge(std::move(name)), d_(diameter), l_(length), k_(kFittings), eps_(roughness) {}

double PipeEdge::churchillDarcy(double re, double rr) {
    if (re < 1.0) return 64.0 / std::max(re, 1e-12);
    const double a = std::pow(2.457 * std::log(1.0 / (std::pow(7.0 / re, 0.9) + 0.27 * rr)), 16.0);
    const double b = std::pow(37530.0 / re, 16.0);
    return 8.0 * std::pow(std::pow(8.0 / re, 12.0) + std::pow(a + b, -1.5), 1.0 / 12.0);
}

double PipeEdge::dp(double mdot, const Fluid& f) const {
    const FluidProps p = f.props(tC_);
    const double area = 0.25 * M_PI * d_ * d_;
    const double v = mdot / (p.rho * area);
    const double re = p.rho * std::abs(v) * d_ / p.mu;
    // Laminar limit written explicitly (dp = 32 mu L v / D^2) so it stays smooth through zero flow.
    double friction;
    if (re < 1.0)
        friction = 32.0 * p.mu * l_ * v / (d_ * d_);
    else
        friction = churchillDarcy(re, eps_ / d_) * (l_ / d_) * 0.5 * p.rho * v * std::abs(v);
    const double minor = k_ * 0.5 * p.rho * sq(v, 1.0);
    return friction + minor;
}

// ---- ChannelEdge -------------------------------------------------------------------------------
double ChannelEdge::dp(double mdot, const Fluid& f) const {
    const FluidProps p = f.props(tC_);
    return ch_.pressureDrop(mdot / p.rho, p);
}

// ---- PumpEdge ------------------------------------------------------------------------------------
double PumpEdge::head(double q) const {
    const double b = c_.shutoffHeadPa / (c_.maxFlowM3s * c_.maxFlowM3s);
    return speed_ * speed_ * c_.shutoffHeadPa - b * sq(q, c_.maxFlowM3s);
}

double PumpEdge::dp(double mdot, const Fluid& f) const {
    const double rho = f.props(tC_).rho;
    return -head(mdot / rho);
}

double PumpEdge::efficiency(double q) const {
    const double s = std::max(speed_, 1e-6);
    const double r = (q / s - c_.bepFlowM3s) / c_.bepFlowM3s;
    return std::clamp(c_.etaMax * (1.0 - r * r), 0.05, c_.etaMax);
}

double PumpEdge::electricPower(double q) const {
    const double hyd = std::max(0.0, head(q) * q);
    return hyd / (efficiency(q) * c_.etaMotor);
}

// ---- SeriesEdge ----------------------------------------------------------------------------------
double SeriesEdge::dp(double mdot, const Fluid& f) const {
    double s = 0.0;
    for (const auto& e : parts_) s += e->dp(mdot, f);
    return s;
}

void SeriesEdge::setTemperature(double tC) {
    tC_ = tC;
    for (auto& e : parts_) e->setTemperature(tC);
}

} // namespace mf
