#include "hydraulics/RectChannel.hpp"

#include <algorithm>
#include <cmath>

namespace mf {

double RectChannel::fanningFRe(double a) {
    // Shah & London (1978), rectangular ducts, 0 <= a <= 1 (short side / long side).
    return 24.0 * (1.0 - 1.3553 * a + 1.9467 * a * a - 1.7012 * std::pow(a, 3) + 0.9564 * std::pow(a, 4) -
                   0.2537 * std::pow(a, 5));
}

double RectChannel::hagenbachK(double a) {
    // Shah & London (1978) curve fit for K(inf) in rectangular ducts.
    return 0.6796 + 1.2197 * a + 3.3089 * a * a - 9.5921 * std::pow(a, 3) + 8.9089 * std::pow(a, 4) -
           2.9959 * std::pow(a, 5);
}

double RectChannel::nu3WallFullyDeveloped(double a) {
    // Shah & London (1978): H1 boundary, three walls heated, one wall of length "width" adiabatic.
    a = std::clamp(a, 0.0, 1.0);
    return 8.235 * (1.0 - 1.883 * a + 3.767 * a * a - 5.814 * std::pow(a, 3) + 5.361 * std::pow(a, 4) -
                    2.0 * std::pow(a, 5));
}

double RectChannel::nuCircularDevelopingH(double xs) {
    // Shah & London (1978) local Nusselt number, thermally developing, hydrodynamically developed,
    // circular duct, constant heat flux. x* = x / (D Re Pr).
    xs = std::max(xs, 1e-6);
    if (xs <= 0.0015) return 1.302 * std::cbrt(1.0 / xs) - 0.5;
    return 4.364 + 8.68 * std::pow(1e3 * xs, -0.506) * std::exp(-41.0 * xs);
}

ChannelState RectChannel::state(double q, const FluidProps& p) const {
    ChannelState s{};
    const double dh = g_.hydraulicDiameter();
    const double qc = q / g_.count;
    s.velocity = qc / g_.area();
    const double u = std::abs(s.velocity);
    s.re = p.rho * u * dh / p.mu;
    const double dynHead = 0.5 * p.rho * s.velocity * u; // signed velocity head
    // Fully developed: dp = (4 fRe / Re) (L/Dh) rho u^2 / 2 = 2 fRe mu L u / Dh^2  (linear in mu and u)
    s.dpFriction = 2.0 * fanningFRe(g_.aspect()) * p.mu * g_.length * s.velocity / (dh * dh);
    s.dpDeveloping = g_.hagenbach ? hagenbachK(g_.aspect()) * dynHead : 0.0;
    s.dpMinor = (g_.kInlet + g_.kOutlet) * dynHead;
    s.dpTotal = s.dpFriction + s.dpDeveloping + s.dpMinor;
    s.lhPlus = s.re > 0 ? g_.length / (dh * s.re) : INFINITY;
    s.entranceLengthHydro = 0.05 * s.re * dh;
    s.entranceLengthThermal = 0.05 * s.re * p.pr() * dh;
    return s;
}

double RectChannel::localH(double x, double qTotal, const FluidProps& p) const {
    const double dh = g_.hydraulicDiameter();
    const double u = std::abs(qTotal) / g_.count / g_.area();
    const double nuFd = nu3WallFullyDeveloped(g_.width / g_.depth);
    // Re * Pr = rho u Dh cp / k: independent of viscosity, so viscosity does not move the thermal
    // entrance region at fixed flow.
    const double rePr = p.rho * u * dh * p.cp / p.k;
    if (rePr <= 0.0) return nuFd * p.k / dh;
    const double xs = std::max(x, 1e-9) / (dh * rePr);
    const double enhancement = std::max(1.0, nuCircularDevelopingH(xs) / 4.364);
    return nuFd * enhancement * p.k / dh;
}

} // namespace mf
