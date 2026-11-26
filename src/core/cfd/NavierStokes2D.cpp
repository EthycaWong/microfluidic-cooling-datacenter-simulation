#include "cfd/NavierStokes2D.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace mf::cfd {

NavierStokes2D::NavierStokes2D(const Params& p) : p_(p) {
    if (p_.nx < 2 || p_.ny < 2) throw std::runtime_error("NavierStokes2D: grid must be at least 2x2");
    reset();
}

void NavierStokes2D::reset() {
    dx_ = p_.lx / p_.nx;
    dy_ = p_.ly / p_.ny;
    const size_t nc = static_cast<size_t>(p_.nx + 2) * (p_.ny + 2);
    u_.assign(static_cast<size_t>(p_.nx + 1) * (p_.ny + 2), 0.0);
    v_.assign(static_cast<size_t>(p_.nx + 2) * (p_.ny + 1), 0.0);
    pr_.assign(nc, 0.0);
    t_.assign(nc, p_.t0);
    f_ = u_;
    g_ = v_;
    time_ = 0.0;
    dt_ = 0.0;
    steps_ = 0;
    applyVelocityBC();
    applyTemperatureBC();
}

void NavierStokes2D::setSide(int which, const Side& s) {
    Side* sides[4] = {&p_.left, &p_.right, &p_.bottom, &p_.top};
    *sides[std::clamp(which, 0, 3)] = s;
}

bool NavierStokes2D::anyOutflow() const {
    return p_.left.velocity == VelocityBC::Outflow || p_.right.velocity == VelocityBC::Outflow ||
           p_.bottom.velocity == VelocityBC::Outflow || p_.top.velocity == VelocityBC::Outflow;
}

double NavierStokes2D::inflowProfile(const Side& s, double frac) const {
    return s.parabolic ? 6.0 * frac * (1.0 - frac) : 1.0; // mean of the parabola = 1
}

// ---- boundary conditions ------------------------------------------------------------------------
void NavierStokes2D::applyVelocityBC() {
    const int nx = p_.nx, ny = p_.ny;
    const Side &L = p_.left, &R = p_.right, &B = p_.bottom, &T = p_.top;
    for (int j = 0; j <= ny + 1; ++j) {
        const double fr = std::clamp((j - 0.5) / ny, 0.0, 1.0);
        if (j >= 1 && j <= ny) {
            // Outflow normal velocities are left alone: the pressure projection sets them (mass balance).
            if (L.velocity == VelocityBC::Wall) U(0, j) = L.u * inflowProfile(L, fr);
            if (R.velocity == VelocityBC::Wall) U(nx, j) = R.u * inflowProfile(R, fr);
        }
    }
    for (int j = 0; j <= ny; ++j) {
        if (L.velocity == VelocityBC::Wall) V(0, j) = 2.0 * L.v - V(1, j);  // tangential: ghost mirror
        else V(0, j) = V(1, j);
        if (R.velocity == VelocityBC::Wall) V(nx + 1, j) = 2.0 * R.v - V(nx, j);
        else V(nx + 1, j) = V(nx, j);
    }
    for (int i = 0; i <= nx + 1; ++i) {
        const double fr = std::clamp((i - 0.5) / nx, 0.0, 1.0);
        if (i >= 1 && i <= nx) {
            if (B.velocity == VelocityBC::Wall) V(i, 0) = B.v * inflowProfile(B, fr);
            if (T.velocity == VelocityBC::Wall) V(i, ny) = T.v * inflowProfile(T, fr);
        }
    }
    for (int i = 0; i <= nx; ++i) {
        if (B.velocity == VelocityBC::Wall) U(i, 0) = 2.0 * B.u - U(i, 1);
        else U(i, 0) = U(i, 1);
        if (T.velocity == VelocityBC::Wall) U(i, ny + 1) = 2.0 * T.u - U(i, ny);
        else U(i, ny + 1) = U(i, ny);
    }
}

void NavierStokes2D::applyTemperatureBC() {
    const int nx = p_.nx, ny = p_.ny;
    auto ghost = [&](const Side& s, double inside, double h) {
        if (s.velocity == VelocityBC::Outflow) return inside;               // zero gradient
        if (s.thermal == ThermalBC::FixedTemperature) return 2.0 * s.T - inside;
        return inside + s.q * h / p_.k;                                      // flux into the fluid
    };
    for (int j = 1; j <= ny; ++j) {
        Tm(0, j) = ghost(p_.left, Tm(1, j), dx_);
        Tm(nx + 1, j) = ghost(p_.right, Tm(nx, j), dx_);
    }
    for (int i = 1; i <= nx; ++i) {
        Tm(i, 0) = ghost(p_.bottom, Tm(i, 1), dy_);
        Tm(i, ny + 1) = ghost(p_.top, Tm(i, ny), dy_);
    }
}

// ---- time step -----------------------------------------------------------------------------------
double NavierStokes2D::computeDt() const {
    double umax = 1e-30, vmax = 1e-30;
    for (double x : u_) umax = std::max(umax, std::abs(x));
    for (double x : v_) vmax = std::max(vmax, std::abs(x));
    const double inv = 1.0 / (dx_ * dx_) + 1.0 / (dy_ * dy_);
    double dt = std::min({0.5 / (p_.nu * inv), 0.5 / (p_.alpha * inv), dx_ / umax, dy_ / vmax});
    return std::min(p_.tau * dt, p_.maxDt);
}

double NavierStokes2D::step() {
    const int nx = p_.nx, ny = p_.ny;
    applyVelocityBC();
    applyTemperatureBC();
    const double dt = computeDt();
    double umax = 0, vmax = 0;
    for (double x : u_) umax = std::max(umax, std::abs(x));
    for (double x : v_) vmax = std::max(vmax, std::abs(x));
    const double gamma = std::clamp(std::max(umax * dt / dx_, vmax * dt / dy_), p_.gammaMin, 1.0);
    const double idx = 1.0 / dx_, idy = 1.0 / dy_, idx2 = idx * idx, idy2 = idy * idy;

    // ---- energy equation (explicit, donor-cell / central blend for convection) ----
    std::vector<double> tn = t_;
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j) {
            const double lap = (Tm(i + 1, j) - 2 * Tm(i, j) + Tm(i - 1, j)) * idx2 +
                               (Tm(i, j + 1) - 2 * Tm(i, j) + Tm(i, j - 1)) * idy2;
            const double duT = idx * (U(i, j) * 0.5 * (Tm(i, j) + Tm(i + 1, j)) - U(i - 1, j) * 0.5 * (Tm(i - 1, j) + Tm(i, j))) +
                               gamma * idx * (std::abs(U(i, j)) * 0.5 * (Tm(i, j) - Tm(i + 1, j)) -
                                              std::abs(U(i - 1, j)) * 0.5 * (Tm(i - 1, j) - Tm(i, j)));
            const double dvT = idy * (V(i, j) * 0.5 * (Tm(i, j) + Tm(i, j + 1)) - V(i, j - 1) * 0.5 * (Tm(i, j - 1) + Tm(i, j))) +
                               gamma * idy * (std::abs(V(i, j)) * 0.5 * (Tm(i, j) - Tm(i, j + 1)) -
                                              std::abs(V(i, j - 1)) * 0.5 * (Tm(i, j - 1) - Tm(i, j)));
            tn[i * (ny + 2) + j] = Tm(i, j) + dt * (p_.alpha * lap - duT - dvT);
        }
    t_.swap(tn);
    applyTemperatureBC();

    // ---- momentum predictor F, G (without pressure) ----
    for (int i = 1; i <= nx - 1; ++i)
        for (int j = 1; j <= ny; ++j) {
            const double lap = (U(i + 1, j) - 2 * U(i, j) + U(i - 1, j)) * idx2 + (U(i, j + 1) - 2 * U(i, j) + U(i, j - 1)) * idy2;
            const double a = U(i, j) + U(i + 1, j), b = U(i - 1, j) + U(i, j);
            const double du2dx = idx * (0.25 * a * a - 0.25 * b * b) +
                                 gamma * idx * (0.25 * std::abs(a) * (U(i, j) - U(i + 1, j)) - 0.25 * std::abs(b) * (U(i - 1, j) - U(i, j)));
            const double c = V(i, j) + V(i + 1, j), d = V(i, j - 1) + V(i + 1, j - 1);
            const double duvdy = idy * (0.25 * c * (U(i, j) + U(i, j + 1)) - 0.25 * d * (U(i, j - 1) + U(i, j))) +
                                 gamma * idy * (0.25 * std::abs(c) * (U(i, j) - U(i, j + 1)) - 0.25 * std::abs(d) * (U(i, j - 1) - U(i, j)));
            const double tf = 0.5 * (Tm(i, j) + Tm(i + 1, j));
            f_[i * (ny + 2) + j] = U(i, j) + dt * (p_.nu * lap - du2dx - duvdy - p_.beta * (tf - p_.tRef) * p_.gx);
        }
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny - 1; ++j) {
            const double lap = (V(i + 1, j) - 2 * V(i, j) + V(i - 1, j)) * idx2 + (V(i, j + 1) - 2 * V(i, j) + V(i, j - 1)) * idy2;
            const double a = V(i, j) + V(i, j + 1), b = V(i, j - 1) + V(i, j);
            const double dv2dy = idy * (0.25 * a * a - 0.25 * b * b) +
                                 gamma * idy * (0.25 * std::abs(a) * (V(i, j) - V(i, j + 1)) - 0.25 * std::abs(b) * (V(i, j - 1) - V(i, j)));
            const double c = U(i, j) + U(i, j + 1), d = U(i - 1, j) + U(i - 1, j + 1);
            const double duvdx = idx * (0.25 * c * (V(i, j) + V(i + 1, j)) - 0.25 * d * (V(i - 1, j) + V(i, j))) +
                                 gamma * idx * (0.25 * std::abs(c) * (V(i, j) - V(i + 1, j)) - 0.25 * std::abs(d) * (V(i - 1, j) - V(i, j)));
            const double tf = 0.5 * (Tm(i, j) + Tm(i, j + 1));
            g_[i * (ny + 1) + j] = V(i, j) + dt * (p_.nu * lap - dv2dy - duvdx - p_.beta * (tf - p_.tRef) * p_.gy);
        }
    // Boundary faces: prescribed (walls) or extrapolated predictor (outflow, corrected by pressure below).
    for (int j = 1; j <= ny; ++j) {
        f_[0 * (ny + 2) + j] = p_.left.velocity == VelocityBC::Outflow ? f_[1 * (ny + 2) + j] : U(0, j);
        f_[nx * (ny + 2) + j] = p_.right.velocity == VelocityBC::Outflow ? f_[(nx - 1) * (ny + 2) + j] : U(nx, j);
    }
    for (int i = 1; i <= nx; ++i) {
        g_[i * (ny + 1) + 0] = p_.bottom.velocity == VelocityBC::Outflow ? g_[i * (ny + 1) + 1] : V(i, 0);
        g_[i * (ny + 1) + ny] = p_.top.velocity == VelocityBC::Outflow ? g_[i * (ny + 1) + ny - 1] : V(i, ny);
    }

    // ---- pressure Poisson: lap(p) = div(F, G) / dt ----
    std::vector<double> rhs(pr_.size(), 0.0);
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j)
            rhs[i * (ny + 2) + j] = ((f_[i * (ny + 2) + j] - f_[(i - 1) * (ny + 2) + j]) * idx +
                                     (g_[i * (ny + 1) + j] - g_[i * (ny + 1) + j - 1]) * idy) / dt;
    solvePressure(rhs);
    applyPressureBC(pr_);

    // ---- projection: u = F - dt grad(p) ----
    const int i0 = p_.left.velocity == VelocityBC::Outflow ? 0 : 1;
    const int i1 = p_.right.velocity == VelocityBC::Outflow ? nx : nx - 1;
    for (int i = i0; i <= i1; ++i)
        for (int j = 1; j <= ny; ++j) U(i, j) = f_[i * (ny + 2) + j] - dt * idx * (P(i + 1, j) - P(i, j));
    const int j0 = p_.bottom.velocity == VelocityBC::Outflow ? 0 : 1;
    const int j1 = p_.top.velocity == VelocityBC::Outflow ? ny : ny - 1;
    for (int i = 1; i <= nx; ++i)
        for (int j = j0; j <= j1; ++j) V(i, j) = g_[i * (ny + 1) + j] - dt * idy * (P(i, j + 1) - P(i, j));
    applyVelocityBC();

    time_ += dt;
    dt_ = dt;
    ++steps_;
    return dt;
}

void NavierStokes2D::runUntil(double t, int maxSteps) {
    for (int n = 0; n < maxSteps && time_ < t; ++n) step();
}

bool NavierStokes2D::runToSteady(double tol, double tScale, int maxSteps) {
    std::vector<double> u0, v0, t0;
    for (int n = 0; n < maxSteps; ++n) {
        u0 = u_; v0 = v_; t0 = t_;
        const double dt = step();
        double du = 0, us = 1e-30, dv = 0, vs = 1e-30, dT = 0, ts = 1e-30;
        for (size_t k = 0; k < u_.size(); ++k) { du = std::max(du, std::abs(u_[k] - u0[k])); us = std::max(us, std::abs(u_[k])); }
        for (size_t k = 0; k < v_.size(); ++k) { dv = std::max(dv, std::abs(v_[k] - v0[k])); vs = std::max(vs, std::abs(v_[k])); }
        double tmin = 1e300, tmax = -1e300;
        for (size_t k = 0; k < t_.size(); ++k) { dT = std::max(dT, std::abs(t_[k] - t0[k])); tmin = std::min(tmin, t_[k]); tmax = std::max(tmax, t_[k]); }
        ts = std::max(tmax - tmin, 1e-12);
        // floor at a tiny fraction of the viscous velocity scale nu/L, so a fluid at rest counts as steady
        const double speed = std::max({us, vs, 1e-3 * p_.nu / std::min(p_.lx, p_.ly)});
        const double rate = std::max({du / speed, dv / speed, dT / ts}) / dt * tScale;
        if (n > 20 && rate < tol) return true;
    }
    return false;
}

// ---- pressure solve: (-lap) p = -rhs, Jacobi-preconditioned conjugate gradients -------------------
void NavierStokes2D::applyPressureBC(std::vector<double>& p) const {
    const int nx = p_.nx, ny = p_.ny;
    auto at = [&](int i, int j) -> double& { return p[i * (ny + 2) + j]; };
    for (int j = 1; j <= ny; ++j) {
        at(0, j) = p_.left.velocity == VelocityBC::Outflow ? -at(1, j) : at(1, j);
        at(nx + 1, j) = p_.right.velocity == VelocityBC::Outflow ? -at(nx, j) : at(nx, j);
    }
    for (int i = 1; i <= nx; ++i) {
        at(i, 0) = p_.bottom.velocity == VelocityBC::Outflow ? -at(i, 1) : at(i, 1);
        at(i, ny + 1) = p_.top.velocity == VelocityBC::Outflow ? -at(i, ny) : at(i, ny);
    }
}

void NavierStokes2D::laplacian(const std::vector<double>& x, std::vector<double>& out) const {
    const int nx = p_.nx, ny = p_.ny;
    const double idx2 = 1.0 / (dx_ * dx_), idy2 = 1.0 / (dy_ * dy_);
    const double sl = p_.left.velocity == VelocityBC::Outflow ? -1.0 : 1.0;
    const double sr = p_.right.velocity == VelocityBC::Outflow ? -1.0 : 1.0;
    const double sb = p_.bottom.velocity == VelocityBC::Outflow ? -1.0 : 1.0;
    const double st = p_.top.velocity == VelocityBC::Outflow ? -1.0 : 1.0;
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j) {
            const double c = x[i * (ny + 2) + j];
            const double w = i == 1 ? sl * c : x[(i - 1) * (ny + 2) + j];
            const double e = i == nx ? sr * c : x[(i + 1) * (ny + 2) + j];
            const double s = j == 1 ? sb * c : x[i * (ny + 2) + j - 1];
            const double n = j == ny ? st * c : x[i * (ny + 2) + j + 1];
            out[i * (ny + 2) + j] = -((e - 2 * c + w) * idx2 + (n - 2 * c + s) * idy2);
        }
}

void NavierStokes2D::solvePressure(const std::vector<double>& rhs) {
    const int nx = p_.nx, ny = p_.ny;
    const size_t N = pr_.size();
    const bool neumann = !anyOutflow();
    const double idx2 = 1.0 / (dx_ * dx_), idy2 = 1.0 / (dy_ * dy_);
    std::vector<double> b(N, 0.0), diag(N, 1.0);
    double mean = 0.0;
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j) { b[i * (ny + 2) + j] = -rhs[i * (ny + 2) + j]; mean += b[i * (ny + 2) + j]; }
    mean /= nx * ny;
    auto sideCoef = [](bool outflow) { return outflow ? 3.0 : 1.0; }; // ghost = -c (Dirichlet) or +c (Neumann)
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j) {
            if (neumann) b[i * (ny + 2) + j] -= mean; // compatibility condition for the pure-Neumann problem
            double dxx = (i == 1 ? sideCoef(p_.left.velocity == VelocityBC::Outflow) : 2.0) - 2.0 +
                         (i == nx ? sideCoef(p_.right.velocity == VelocityBC::Outflow) : 2.0);
            double dyy = (j == 1 ? sideCoef(p_.bottom.velocity == VelocityBC::Outflow) : 2.0) - 2.0 +
                         (j == ny ? sideCoef(p_.top.velocity == VelocityBC::Outflow) : 2.0);
            diag[i * (ny + 2) + j] = std::max(1e-30, dxx * idx2 + dyy * idy2);
        }
    // interior mask helpers
    auto dotI = [&](const std::vector<double>& a, const std::vector<double>& c) {
        double s = 0.0;
        for (int i = 1; i <= nx; ++i)
            for (int j = 1; j <= ny; ++j) s += a[i * (ny + 2) + j] * c[i * (ny + 2) + j];
        return s;
    };
    std::vector<double> x = pr_, r(N, 0.0), z(N, 0.0), d(N, 0.0), q(N, 0.0);
    laplacian(x, q);
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j) {
            const size_t k = i * (ny + 2) + j;
            r[k] = b[k] - q[k];
            z[k] = r[k] / diag[k];
            d[k] = z[k];
        }
    const double bnorm = std::sqrt(std::max(dotI(b, b), 1e-300));
    double rz = dotI(r, z);
    int it = 0;
    const int maxIt = 4 * nx * ny + 100;
    for (; it < maxIt; ++it) {
        if (std::sqrt(dotI(r, r)) <= p_.pressureTol * bnorm) break;
        laplacian(d, q);
        const double dq = dotI(d, q);
        if (std::abs(dq) < 1e-300) break;
        const double a = rz / dq;
        for (int i = 1; i <= nx; ++i)
            for (int j = 1; j <= ny; ++j) {
                const size_t k = i * (ny + 2) + j;
                x[k] += a * d[k];
                r[k] -= a * q[k];
                z[k] = r[k] / diag[k];
            }
        const double rzNew = dotI(r, z);
        const double beta = rzNew / rz;
        rz = rzNew;
        for (int i = 1; i <= nx; ++i)
            for (int j = 1; j <= ny; ++j) {
                const size_t k = i * (ny + 2) + j;
                d[k] = z[k] + beta * d[k];
            }
    }
    if (neumann) { // pin the free constant: zero mean pressure
        double m = 0.0;
        for (int i = 1; i <= nx; ++i)
            for (int j = 1; j <= ny; ++j) m += x[i * (ny + 2) + j];
        m /= nx * ny;
        for (int i = 1; i <= nx; ++i)
            for (int j = 1; j <= ny; ++j) x[i * (ny + 2) + j] -= m;
    }
    pr_.swap(x);
    cgIter_ = it;
}

// ---- diagnostics ---------------------------------------------------------------------------------
double NavierStokes2D::maxSpeed() const {
    double m = 0.0;
    for (int i = 1; i <= p_.nx; ++i)
        for (int j = 1; j <= p_.ny; ++j) m = std::max(m, std::hypot(uc(i, j), vc(i, j)));
    return m;
}

double NavierStokes2D::maxDivergence() const {
    double m = 0.0;
    for (int i = 1; i <= p_.nx; ++i)
        for (int j = 1; j <= p_.ny; ++j)
            m = std::max(m, std::abs((U(i, j) - U(i - 1, j)) / dx_ + (V(i, j) - V(i, j - 1)) / dy_));
    return m;
}

double NavierStokes2D::kineticEnergy() const {
    double e = 0.0;
    for (int i = 1; i <= p_.nx; ++i)
        for (int j = 1; j <= p_.ny; ++j) e += 0.5 * (uc(i, j) * uc(i, j) + vc(i, j) * vc(i, j)) * dx_ * dy_;
    return e;
}

double NavierStokes2D::wallHeatFlux(int which) const {
    const Side* sides[4] = {&p_.left, &p_.right, &p_.bottom, &p_.top};
    const Side& s = *sides[std::clamp(which, 0, 3)];
    if (s.thermal == ThermalBC::HeatFlux || s.velocity == VelocityBC::Outflow) return s.thermal == ThermalBC::HeatFlux ? s.q : 0.0;
    const int nx = p_.nx, ny = p_.ny;
    double sum = 0.0;
    int n = 0;
    // second-order one-sided derivative at the wall: dT/dn ~ (-8 Tw + 9 T1 - T2) / (3 h)
    auto flux = [&](double tw, double t1, double t2, double h) { return -p_.k * (-8.0 * tw + 9.0 * t1 - t2) / (3.0 * h); };
    if (which == 0) for (int j = 1; j <= ny; ++j, ++n) sum += flux(s.T, Tm(1, j), Tm(2, j), dx_);
    if (which == 1) for (int j = 1; j <= ny; ++j, ++n) sum += flux(s.T, Tm(nx, j), Tm(nx - 1, j), dx_);
    if (which == 2) for (int i = 1; i <= nx; ++i, ++n) sum += flux(s.T, Tm(i, 1), Tm(i, 2), dy_);
    if (which == 3) for (int i = 1; i <= nx; ++i, ++n) sum += flux(s.T, Tm(i, ny), Tm(i, ny - 1), dy_);
    return n ? sum / n : 0.0;
}

double NavierStokes2D::bulkTemperatureColumn(int i) const {
    double num = 0.0, den = 0.0;
    for (int j = 1; j <= p_.ny; ++j) { num += std::abs(uc(i, j)) * tc(i, j); den += std::abs(uc(i, j)); }
    return den > 0 ? num / den : tc(i, p_.ny / 2);
}

double NavierStokes2D::meanPressureColumn(int i) const {
    double s = 0.0;
    for (int j = 1; j <= p_.ny; ++j) s += P(i, j);
    return s / p_.ny;
}

void NavierStokes2D::velocityAt(double x, double y, double& u, double& v) const {
    const int nx = p_.nx, ny = p_.ny;
    // u lives at (i dx, (j - 0.5) dy), i = 0..nx, j = 0..ny+1
    double fi = std::clamp(x / dx_, 0.0, static_cast<double>(nx) - 1e-9);
    double fj = std::clamp(y / dy_ + 0.5, 0.0, static_cast<double>(ny + 1) - 1e-9);
    int i = static_cast<int>(fi), j = static_cast<int>(fj);
    double a = fi - i, b = fj - j;
    u = (1 - a) * (1 - b) * U(i, j) + a * (1 - b) * U(i + 1, j) + (1 - a) * b * U(i, j + 1) + a * b * U(i + 1, j + 1);
    // v lives at ((i - 0.5) dx, j dy), i = 0..nx+1, j = 0..ny
    fi = std::clamp(x / dx_ + 0.5, 0.0, static_cast<double>(nx + 1) - 1e-9);
    fj = std::clamp(y / dy_, 0.0, static_cast<double>(ny) - 1e-9);
    i = static_cast<int>(fi); j = static_cast<int>(fj);
    a = fi - i; b = fj - j;
    v = (1 - a) * (1 - b) * V(i, j) + a * (1 - b) * V(i + 1, j) + (1 - a) * b * V(i, j + 1) + a * b * V(i + 1, j + 1);
}

double NavierStokes2D::temperatureAt(double x, double y) const {
    const double fi = std::clamp(x / dx_ + 0.5, 0.0, static_cast<double>(p_.nx + 1) - 1e-9);
    const double fj = std::clamp(y / dy_ + 0.5, 0.0, static_cast<double>(p_.ny + 1) - 1e-9);
    const int i = static_cast<int>(fi), j = static_cast<int>(fj);
    const double a = fi - i, b = fj - j;
    return (1 - a) * (1 - b) * Tm(i, j) + a * (1 - b) * Tm(i + 1, j) + (1 - a) * b * Tm(i, j + 1) + a * b * Tm(i + 1, j + 1);
}

std::vector<double> NavierStokes2D::streamFunction() const {
    const int nx = p_.nx, ny = p_.ny;
    std::vector<double> psi(static_cast<size_t>(nx + 1) * (ny + 1), 0.0);
    auto at = [&](int i, int j) -> double& { return psi[i * (ny + 1) + j]; };
    for (int i = 1; i <= nx; ++i) at(i, 0) = at(i - 1, 0) - V(i, 0) * dx_;   // d(psi)/dx = -v
    for (int i = 0; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j) at(i, j) = at(i, j - 1) + U(i, j) * dy_; // d(psi)/dy = u
    return psi;
}

} // namespace mf::cfd
