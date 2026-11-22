#include "fluid/Fluid.hpp"
#include "Units.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace mf {
namespace {

// Melinder (2010) aqueous glycol correlations as distributed by CoolProp (MIT licence),
// files dev/incompressible_liquids/json/MPG.json and MEG.json.
// f(T, x) = sum_i sum_j C[i][j] (T - Tbase)^i (x - xbase)^j   (i: temperature power, j: mass-fraction power)
// Viscosity is exp() of the same polynomial (Pa s). Units: K, kg/m3, J/(kg K), W/(m K).
struct Coeffs {
    double tBase, xBase;
    double rho[4][6], cp[4][6], k[4][6], lnMu[4][6];
    double tFreeze[6]; // K, polynomial in (x - xbase) only
};

const Coeffs kMPG = {
    305.8583, 0.307031,
    {{1018.0, 76.04, -24.98, -155.0, -113.1, 234.2},
     {-0.5406, -0.945, 0.27, 2.829, -2.221, 0.0},
     {-0.002666, 0.005541, -0.004018, -0.007175, 0.0, 0.0},
     {1.347e-05, -1.343e-05, 3.376e-05, 0.0, 0.0, 0.0}},
    {{3882.0, -1304.0, -1598.0, 353.9, 5000.0, -4959.0},
     {2.699, 5.07, 0.9534, 31.02, -71.35, 0.0},
     {-0.001659, -0.004752, 0.1167, -0.295, 0.0, 0.0},
     {-1.032e-05, 0.0001522, -0.000487, 0.0, 0.0, 0.0}},
    {{0.4513, -0.4795, 0.2076, -0.09083, -0.05952, 0.2104},
     {0.0007955, -0.001678, 0.001563, -0.002518, -0.003605, 0.0},
     {3.482e-08, 8.941e-06, -4.615e-05, 6.543e-05, 0.0, 0.0},
     {-5.966e-09, 1.493e-08, 9.897e-08, 0.0, 0.0, 0.0}},
    {{-6.224055, 3.328, 0.5453, -3.9, -1.587, 35.64},
     {-0.03045, -0.03984, -0.00086, 0.1054, 0.04475, 0.0},
     {0.0002525, 0.0004332, -0.0001593, -0.001589, 0.0, 0.0},
     {-1.399e-06, -1.86e-06, -4.465e-07, 0.0, 0.0, 0.0}},
    {259.9, -66.31, -109.4, -228.3, -340.9, 146.5}};

const Coeffs kMEG = {
    304.878, 0.308462,
    {{1034.0, 131.1, 0.749, -106.2, -96.23, 489.1},
     {-0.4781, -0.6876, 0.7855, 1.229, -7.211, 0.0},
     {-0.002692, 0.004805, -0.003995, -0.01153, 0.0, 0.0},
     {4.725e-06, 1.69e-06, 4.982e-05, 0.0, 0.0, 0.0}},
    {{3737.0, -1799.0, -993.3, 2610.0, 1537.0, -16180.0},
     {2.93, 10.46, 3.516, -1.189, -42.72, 0.0},
     {-0.004675, -0.04147, 0.05109, -0.1643, 0.0, 0.0},
     {-1.389e-05, 1.847e-05, -0.0007138, 0.0, 0.0, 0.0}},
    {{0.472, -0.4286, 0.1747, 0.03017, -0.1322, 0.2678},
     {0.0008903, -0.001473, 0.0006814, -0.002412, 0.002555, 0.0},
     {-1.058e-06, 1.059e-05, -3.612e-05, 4.004e-05, 0.0, 0.0},
     {-2.789e-09, -1.142e-08, 2.365e-08, 0.0, 0.0, 0.0}},
    {{-6.437255, 2.471, 0.03328, 1.659, 3.089, -18.65},
     {-0.0255, -0.01171, 0.01086, 0.003157, 0.01831, 0.0},
     {0.0001782, 0.0001052, 0.0001051, 0.0004063, 0.0, 0.0},
     {-7.669e-07, -1.634e-06, -6.475e-06, 0.0, 0.0, 0.0}},
    {257.9, -80.8, -133.4, -72.93, 100.6, 11.4}};

double poly(const double c[4][6], double ti, double xi) {
    double s = 0.0, tp = 1.0;
    for (int i = 0; i < 4; ++i) {
        double row = 0.0, xp = 1.0;
        for (int j = 0; j < 6; ++j) { row += c[i][j] * xp; xp *= xi; }
        s += row * tp;
        tp *= ti;
    }
    return s;
}

const Coeffs& coeffs(GlycolMixture::Kind k) { return k == GlycolMixture::Kind::Propylene ? kMPG : kMEG; }

} // namespace

GlycolMixture::GlycolMixture(Kind kind, double massFraction) : kind_(kind), x_(massFraction) {
    if (x_ < 0.0 || x_ > 0.6) throw std::runtime_error("glycol mass fraction must be within 0 - 0.6 (Melinder range)");
    const int n = static_cast<int>(std::lround((kTmax - kTmin) / kStep)) + 1;
    table_.reserve(n);
    for (int i = 0; i < n; ++i) table_.push_back(propsExact(kTmin + i * kStep));
}

double GlycolMixture::massFromVolumeFraction(Kind kind, double phi) {
    // Pure-component densities at 20 C; volume contraction on mixing is neglected (~1 % effect).
    const double rhoG = kind == Kind::Propylene ? 1036.0 : 1113.0, rhoW = 998.2;
    return phi * rhoG / (phi * rhoG + (1.0 - phi) * rhoW);
}

FluidProps GlycolMixture::propsExact(double tC) const {
    const Coeffs& c = coeffs(kind_);
    const double ti = toKelvin(tC) - c.tBase, xi = x_ - c.xBase;
    return FluidProps{poly(c.rho, ti, xi), std::exp(poly(c.lnMu, ti, xi)), poly(c.cp, ti, xi), poly(c.k, ti, xi)};
}

double GlycolMixture::freezingPointC() const {
    const Coeffs& c = coeffs(kind_);
    double s = 0.0, xp = 1.0;
    for (int j = 0; j < 6; ++j) { s += c.tFreeze[j] * xp; xp *= (x_ - c.xBase); }
    return s - kZeroCelsius;
}

FluidProps GlycolMixture::props(double tC) const {
    const double x = (std::clamp(tC, kTmin, kTmax) - kTmin) / kStep;
    const int i = std::min(static_cast<int>(x), static_cast<int>(table_.size()) - 2);
    const double f = x - i;
    const FluidProps& a = table_[i];
    const FluidProps& b = table_[i + 1];
    return FluidProps{a.rho + f * (b.rho - a.rho), a.mu + f * (b.mu - a.mu), a.cp + f * (b.cp - a.cp),
                      a.k + f * (b.k - a.k)};
}

std::string GlycolMixture::name() const {
    char b[64];
    std::snprintf(b, sizeof b, "%s glycol %.1f%% mass (Melinder)", kind_ == Kind::Propylene ? "Propylene" : "Ethylene",
                  100.0 * x_);
    return b;
}

std::shared_ptr<const Fluid> makeCoolant(const std::string& spec) {
    // "water" | "PG25" | "PG40" | "EG30" ... number = volume percent (as coolants are specified);
    // "PG25m" / "EG30m" = mass percent.
    if (spec.empty() || spec == "water" || spec == "Water") return std::make_shared<WaterIAPWS>();
    if (spec.size() >= 3 && (spec[0] == 'P' || spec[0] == 'E') && spec[1] == 'G') {
        const auto kind = spec[0] == 'P' ? GlycolMixture::Kind::Propylene : GlycolMixture::Kind::Ethylene;
        const bool mass = spec.back() == 'm';
        const double pct = std::stod(spec.substr(2, spec.size() - 2 - (mass ? 1 : 0)));
        const double w = mass ? pct / 100.0 : GlycolMixture::massFromVolumeFraction(kind, pct / 100.0);
        return std::make_shared<GlycolMixture>(kind, w);
    }
    throw std::runtime_error("unknown coolant '" + spec + "' (use water, PG25, EG30, PG25m ...)");
}

} // namespace mf
