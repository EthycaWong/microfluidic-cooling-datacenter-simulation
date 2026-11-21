#include "fluid/Fluid.hpp"
#include "Units.hpp"

#include <algorithm>
#include <cmath>

namespace mf {
namespace {

// ---- IAPWS-IF97 Region 1 (R7-97(2012), Table 2) ------------------------------------------------
constexpr int kI1[34] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                         2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32};
constexpr int kJ1[34] = {-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1,
                         3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41};
constexpr double kN1[34] = {
    0.14632971213167,     -0.84548187169114,    -0.37563603672040e1,  0.33855169168385e1,
    -0.95791963387872,    0.15772038513228,     -0.16616417199501e-1, 0.81214629983568e-3,
    0.28319080123804e-3,  -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1,
    -0.21841717175414e-1, -0.52838357969930e-4, -0.47184321073267e-3, -0.30001780793026e-3,
    0.47661393906987e-4,  -0.44141845330846e-5, -0.72694996297594e-15, -0.31679644845054e-4,
    -0.28270797985312e-5, -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6,
    -0.14341729937924e-12, -0.40516996860117e-6, -0.12734301741641e-8, -0.17427420029190e-9,
    -0.68762131295531e-18, 0.14478307828521e-19, 0.26335781662795e-22, -0.11947622640071e-22,
    0.18228094581404e-23, -0.93537087292458e-25};
constexpr double kR = 461.526;    // J/(kg K), specific gas constant of water (IF97)
constexpr double kPstar1 = 16.53e6;
constexpr double kTstar1 = 1386.0;

// ---- IAPWS 2008 viscosity (R12-08) ---------------------------------------------------------------
constexpr double kTc = 647.096, kRhoc = 322.0;
constexpr double kH0[4] = {1.67752, 2.20462, 0.6366564, -0.241605};
// H1[i][j], i = 0..5 (powers of 1/Tbar - 1), j = 0..6 (powers of rhobar - 1)
constexpr double kH1[6][7] = {
    {5.20094e-1, 2.22531e-1, -2.81378e-1, 1.61913e-1, -3.25372e-2, 0.0, 0.0},
    {8.50895e-2, 9.99115e-1, -9.06851e-1, 2.57399e-1, 0.0, 0.0, 0.0},
    {-1.08374, 1.88797, -7.72479e-1, 0.0, 0.0, 0.0, 0.0},
    {-2.89555e-1, 1.26613, -4.89837e-1, 0.0, 6.98452e-2, 0.0, -4.35673e-3},
    {0.0, 0.0, -2.57040e-1, 0.0, 0.0, 8.72102e-3, 0.0},
    {0.0, 1.20573e-1, 0.0, 0.0, 0.0, 0.0, -5.93264e-4}};

// ---- IAPWS 2011 thermal conductivity (R15-11) ----------------------------------------------------
constexpr double kL0[5] = {2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4};
constexpr double kL1[5][6] = {
    {1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258},
    {2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245},
    {2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816},
    {-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0.0, 0.0},
    {-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842}};

} // namespace

double WaterIAPWS::if97SpecificVolume(double T, double p) {
    const double pi = p / kPstar1, tau = kTstar1 / T;
    double gPi = 0.0;
    for (int i = 0; i < 34; ++i)
        gPi += -kN1[i] * kI1[i] * std::pow(7.1 - pi, kI1[i] - 1) * std::pow(tau - 1.222, kJ1[i]);
    return kR * T / p * pi * gPi;
}

double WaterIAPWS::if97Cp(double T, double p) {
    const double pi = p / kPstar1, tau = kTstar1 / T;
    double gTauTau = 0.0;
    for (int i = 0; i < 34; ++i)
        gTauTau += kN1[i] * std::pow(7.1 - pi, kI1[i]) * kJ1[i] * (kJ1[i] - 1) *
                   std::pow(tau - 1.222, kJ1[i] - 2);
    return -kR * tau * tau * gTauTau;
}

double WaterIAPWS::viscosity2008(double T, double rho) {
    const double Tb = T / kTc, rb = rho / kRhoc;
    double s0 = 0.0;
    for (int i = 0; i < 4; ++i) s0 += kH0[i] / std::pow(Tb, i);
    const double mu0 = 100.0 * std::sqrt(Tb) / s0;
    double s1 = 0.0;
    for (int i = 0; i < 6; ++i) {
        const double a = std::pow(1.0 / Tb - 1.0, i);
        double inner = 0.0;
        for (int j = 0; j < 7; ++j) inner += kH1[i][j] * std::pow(rb - 1.0, j);
        s1 += a * inner;
    }
    const double mu1 = std::exp(rb * s1);
    return mu0 * mu1 * 1e-6; // reference viscosity 1 uPa s
}

double WaterIAPWS::conductivity2011(double T, double rho) {
    const double Tb = T / kTc, rb = rho / kRhoc;
    double s0 = 0.0;
    for (int k = 0; k < 5; ++k) s0 += kL0[k] / std::pow(Tb, k);
    const double l0 = std::sqrt(Tb) / s0;
    double s1 = 0.0;
    for (int i = 0; i < 5; ++i) {
        const double a = std::pow(1.0 / Tb - 1.0, i);
        double inner = 0.0;
        for (int j = 0; j < 6; ++j) inner += kL1[i][j] * std::pow(rb - 1.0, j);
        s1 += a * inner;
    }
    const double l1 = std::exp(rb * s1);
    return l0 * l1 * 1e-3; // reference conductivity 1 mW/(m K)
}

WaterIAPWS::WaterIAPWS(double pressurePa) : p_(pressurePa) {
    const int n = static_cast<int>(std::lround((kTmax - kTmin) / kStep)) + 1;
    table_.reserve(n);
    for (int i = 0; i < n; ++i) table_.push_back(propsExact(kTmin + i * kStep));
}

FluidProps WaterIAPWS::propsExact(double tC) const {
    const double T = toKelvin(tC);
    const double rho = 1.0 / if97SpecificVolume(T, p_);
    return FluidProps{rho, viscosity2008(T, rho), if97Cp(T, p_), conductivity2011(T, rho)};
}

FluidProps WaterIAPWS::props(double tC) const {
    const double x = (std::clamp(tC, kTmin, kTmax) - kTmin) / kStep;
    const int i = std::min(static_cast<int>(x), static_cast<int>(table_.size()) - 2);
    const double f = x - i;
    const FluidProps& a = table_[i];
    const FluidProps& b = table_[i + 1];
    return FluidProps{a.rho + f * (b.rho - a.rho), a.mu + f * (b.mu - a.mu), a.cp + f * (b.cp - a.cp),
                      a.k + f * (b.k - a.k)};
}

} // namespace mf
