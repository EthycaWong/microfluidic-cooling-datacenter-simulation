#include "cfd/DuctFlow.hpp"

#include <algorithm>
#include <cmath>

namespace mf::cfd {
namespace {

/// Solve -lap(phi) = src on an nx x ny cell grid (spacing h), phi = 0 on Dirichlet walls,
/// zero gradient on the top wall when topNeumann. Jacobi-preconditioned CG.
std::vector<double> poisson(int nx, int ny, double h, const std::vector<double>& src, bool topNeumann) {
    const size_t N = static_cast<size_t>(nx) * ny;
    auto apply = [&](const std::vector<double>& x, std::vector<double>& out) {
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const double c = x[j * nx + i];
                const double w = i > 0 ? x[j * nx + i - 1] : -c;          // Dirichlet wall: ghost = -c
                const double e = i < nx - 1 ? x[j * nx + i + 1] : -c;
                const double s = j > 0 ? x[(j - 1) * nx + i] : -c;
                const double n = j < ny - 1 ? x[(j + 1) * nx + i] : (topNeumann ? c : -c);
                out[j * nx + i] = -(e + w + n + s - 4 * c) / (h * h);
            }
    };
    std::vector<double> diag(N);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            double d = 4;
            if (i == 0) d += 1;
            if (i == nx - 1) d += 1;
            if (j == 0) d += 1;
            if (j == ny - 1) d += topNeumann ? -1 : 1;
            diag[j * nx + i] = d / (h * h);
        }
    std::vector<double> x(N, 0.0), r = src, z(N), d(N), q(N);
    for (size_t k = 0; k < N; ++k) { z[k] = r[k] / diag[k]; d[k] = z[k]; }
    auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0;
        for (size_t k = 0; k < N; ++k) s += a[k] * b[k];
        return s;
    };
    const double bn = std::sqrt(dot(src, src));
    double rz = dot(r, z);
    for (int it = 0; it < 20 * static_cast<int>(N); ++it) {
        if (std::sqrt(dot(r, r)) < 1e-12 * bn) break;
        apply(d, q);
        const double a = rz / dot(d, q);
        for (size_t k = 0; k < N; ++k) { x[k] += a * d[k]; r[k] -= a * q[k]; z[k] = r[k] / diag[k]; }
        const double rz2 = dot(r, z);
        const double b = rz2 / rz;
        rz = rz2;
        for (size_t k = 0; k < N; ++k) d[k] = z[k] + b * d[k];
    }
    return x;
}

} // namespace

DuctResult solveRectDuct(double width, double height, int nShort) {
    DuctResult res;
    const double shortSide = std::min(width, height);
    const double h = shortSide / nShort;
    const int nx = std::max(2, static_cast<int>(std::lround(width / h)));
    const int ny = std::max(2, static_cast<int>(std::lround(height / h)));
    const double hx = width / nx; // equal to h up to rounding
    res.nx = nx;
    res.ny = ny;
    const size_t N = static_cast<size_t>(nx) * ny;
    const double area = width * height, perim = 2 * (width + height), dh = 4 * area / perim;

    // Momentum: -lap(w) = G with G = -dp/dz / mu = 1.
    const std::vector<double> w = poisson(nx, ny, hx, std::vector<double>(N, 1.0), false);
    double wm = 0, wmax = 0;
    for (double x : w) { wm += x; wmax = std::max(wmax, x); }
    wm /= N;
    // Wall shear balance: tau_mean * P = (-dp/dz) A  ->  Fanning f Re = Dh^2 G / (2 w_mean) with G = 1.
    res.fRe = dh * dh / (2.0 * wm);
    res.maxOverMean = wmax / wm;

    // Energy, H1: lap(theta) = S w, theta = 0 on heated walls. Solve -lap(theta) = -w (S = 1).
    auto nusselt = [&](bool topInsulated, std::vector<double>* field) {
        std::vector<double> src(N);
        for (size_t k = 0; k < N; ++k) src[k] = -w[k];
        const std::vector<double> th = poisson(nx, ny, hx, src, topInsulated);
        double num = 0, den = 0;
        for (size_t k = 0; k < N; ++k) { num += w[k] * th[k]; den += w[k]; }
        const double thetaB = num / den; // bulk (Tb - Tw), negative for heating
        const double heatedPerimeter = topInsulated ? perim - width : perim;
        // q'' P_h = rho cp w_m A dTb/dz  ->  with S = (dTb/dz)/alpha = 1: q''/k = w_m A / P_h
        const double qOverK = wm * area / heatedPerimeter;
        if (field) {
            field->resize(N);
            for (size_t k = 0; k < N; ++k) (*field)[k] = th[k] / thetaB;
        }
        return qOverK * dh / (-thetaB);
    };
    res.nuH1FourWalls = nusselt(false, nullptr);
    res.nuH1ThreeWalls = nusselt(true, &res.theta);
    res.w.resize(N);
    for (size_t k = 0; k < N; ++k) res.w[k] = w[k] / wm;
    return res;
}

} // namespace mf::cfd
