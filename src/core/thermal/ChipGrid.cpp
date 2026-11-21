#include "thermal/ChipGrid.hpp"
#include "Units.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mf {

ChipGrid::ChipGrid(ChipParams p, std::shared_ptr<const Fluid> fluid, double initialT)
    : p_(std::move(p)), fluid_(std::move(fluid)), channel_(p_.channels), tin_(initialT) {
    if (p_.nx < 1 || p_.ny < 1) throw std::runtime_error("ChipGrid: grid must be at least 1x1");
    if (p_.mode == CoolingMode::Microchannel && p_.channels.depth >= p_.die.thickness)
        throw std::runtime_error("ChipGrid: channel depth must be less than die thickness");
    layers_ = hasCopper() ? 3 : 2;
    const size_t n = static_cast<size_t>(p_.nx) * p_.ny * layers_;
    T_.assign(n, initialT);
    Told_ = T_;
    C_.assign(n, 0.0);
    rows_.resize(n);
    weights_ = buildPowerWeights(p_.power.map, p_.nx, p_.ny);
    hRow_.assign(p_.nx, 0.0);
    power_ = p_.power.at(0.0);
    peak_ = initialT;
    build(1.0); // fills the capacitances so storedEnergy() is meaningful before the first step
}

void ChipGrid::addG(int a, int b, double g) {
    Row& ra = rows_[a];
    Row& rb = rows_[b];
    ra.diag += g;
    rb.diag += g;
    ra.col[ra.n] = b; ra.val[ra.n++] = -g;
    rb.col[rb.n] = a; rb.val[rb.n++] = -g;
}

void ChipGrid::build(double dt) {
    const int nx = p_.nx, ny = p_.ny;
    const double dx = p_.die.lx / nx, dy = p_.die.ly / ny, A = dx * dy;
    for (auto& r : rows_) r = Row{};

    auto kSi = [&](double tC) { return p_.kSi300 * std::pow(300.0 / toKelvin(tC), p_.kSiExponent); };
    auto harm = [](double a, double b) { return 2.0 * a * b / (a + b); };

    // Coolant properties: one cp per chip (at mean coolant T) keeps the energy ledger exact.
    const FluidProps pm = fluid_->props(meanFluidT());
    cpChip_ = pm.cp;
    const double qTotal = mdot_ / pm.rho;
    const double mcol = mdot_ / ny;
    const double adv = mcol * cpChip_;
    gTotal_ = 0.0;

    const ChannelGeometry& g = p_.channels;
    const bool micro = !hasCopper();
    // Silicon conduction thicknesses.
    const double tFloor = micro ? p_.die.thickness - g.depth : p_.die.thickness;
    const double pitch = p_.die.ly / g.count;
    const double ribFrac = micro ? g.rib / pitch : 0.0;
    const double tAlongX = micro ? tFloor + g.depth * ribFrac : tFloor; // ribs conduct along the flow
    const double tAcrossY = tFloor;                                      // channels cut the ribs across
    const double nchCol = static_cast<double>(g.count) / ny;

    // Cold plate lumped convection, shared uniformly by cell area.
    const ColdPlateParams& cp = p_.coldPlate;
    const double qRef = lpmToM3s(cp.qRefLpm);
    const double hAcell = cp.hARef * std::pow(std::max(std::abs(qTotal), 1e-3 * qRef) / qRef, cp.hAExponent) / (nx * ny);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int s = idx(i, j, 0), f = idx(i, j, layers_ - 1);
            const double ks = kSi(T_[s]);
            const FluidProps pf = fluid_->props(T_[f]);

            // --- capacitances ---
            C_[s] = p_.rhoCpSi * A * tAlongX;
            if (micro)
                C_[f] = pf.rho * cpChip_ * nchCol * g.width * g.depth * dx;
            else
                C_[f] = pf.rho * cpChip_ * cp.fluidVolume / (nx * ny);

            // --- lateral silicon conduction (to +x and +y neighbours) ---
            if (i + 1 < nx) {
                const int s2 = idx(i + 1, j, 0);
                addG(s, s2, harm(ks, kSi(T_[s2])) * dy * tAlongX / dx);
            }
            if (j + 1 < ny) {
                const int s2 = idx(i, j + 1, 0);
                addG(s, s2, harm(ks, kSi(T_[s2])) * dx * tAcrossY / dy);
            }

            // --- vertical path to the coolant ---
            if (micro) {
                const double x = (i + 0.5) * dx;
                const double h = channel_.localH(x, qTotal, pf);
                if (j == 0) hRow_[i] = h;
                const double m = std::sqrt(2.0 * h / (ks * g.rib));
                const double mH = m * g.depth;
                const double eta = mH > 1e-9 ? std::tanh(mH) / mH : 1.0;
                const double gConv = h * (g.width + 2.0 * g.depth * eta) * dx * nchCol;
                const double gFloor = ks * A / tFloor;
                const double gsf = 1.0 / (1.0 / gConv + 1.0 / gFloor);
                addG(s, f, gsf);
                gTotal_ += gsf;
            } else {
                const int c = idx(i, j, 1);
                C_[c] = cp.rhoCpCopper * A * cp.baseThickness;
                const double gsc = 1.0 / (tFloor / (ks * A) + cp.timResistance / A +
                                          0.5 * cp.baseThickness / (cp.kCopper * A));
                addG(s, c, gsc);
                const double gcf = 1.0 / (0.5 * cp.baseThickness / (cp.kCopper * A) + 1.0 / hAcell);
                addG(c, f, gcf);
                gTotal_ += 1.0 / (1.0 / gsc + 1.0 / gcf);
                if (i + 1 < nx) addG(c, idx(i + 1, j, 1), cp.kCopper * dy * cp.baseThickness / dx);
                if (j + 1 < ny) addG(c, idx(i, j + 1, 1), cp.kCopper * dx * cp.baseThickness / dy);
                if (j == 0) hRow_[i] = hAcell / A;
            }

            // --- heat input on the silicon node ---
            rows_[s].rhs += power_ * weights_[j * nx + i];

            // --- upwind advection into the coolant node ---
            Row& rf = rows_[f];
            rf.diag += adv;
            if (i == 0) {
                rf.rhs += adv * tin_;
            } else {
                rf.col[rf.n] = idx(i - 1, j, layers_ - 1);
                rf.val[rf.n++] = -adv;
            }
        }
    }
    // --- storage terms ---
    for (size_t n = 0; n < rows_.size(); ++n) {
        const double cdt = C_[n] / dt;
        rows_[n].diag += cdt;
        rows_[n].rhs += cdt * Told_[n];
    }
}

int ChipGrid::gaussSeidel(double omega, double tol, int maxSweeps) {
    int sweep = 0;
    for (; sweep < maxSweeps; ++sweep) {
        double maxd = 0.0;
        for (size_t n = 0; n < rows_.size(); ++n) {
            const Row& r = rows_[n];
            double s = r.rhs;
            for (int k = 0; k < r.n; ++k) s -= r.val[k] * T_[r.col[k]];
            const double tNew = T_[n] + omega * (s / r.diag - T_[n]);
            maxd = std::max(maxd, std::abs(tNew - T_[n]));
            T_[n] = tNew;
        }
        if (maxd < tol) { ++sweep; break; }
    }
    return sweep;
}

void ChipGrid::step(double dt) {
    if (!(dt > 0.0)) return;
    Told_ = T_;
    build(dt);
    lastSweeps_ = gaussSeidel(1.0, 1e-9, 2000);

    // Energy ledger. Summing the discrete equations, all internal conductances cancel, leaving
    // P = sum(C dT/dt) + sum_j mcol cp (T_out,j - T_in): the residual measures solver error only.
    double stor = 0.0;
    for (size_t n = 0; n < T_.size(); ++n) stor += C_[n] * (T_[n] - Told_[n]) / dt;
    double toCool = 0.0;
    const double mcol = mdot_ / p_.ny;
    for (int j = 0; j < p_.ny; ++j) toCool += mcol * cpChip_ * (fluidT(p_.nx - 1, j) - tin_);
    ledger_.powerIn = power_;
    ledger_.toCoolant = toCool;
    ledger_.storageRate = stor;
    ledger_.residual = power_ - toCool - stor;
    updatePeak();
}

void ChipGrid::solveSteady(double tolK) {
    // Pseudo-transient continuation with very large steps: each step is an (almost) steady solve;
    // repeating it converges the lagged temperature-dependent coefficients.
    for (int outer = 0; outer < 500; ++outer) {
        Told_ = T_;
        build(1e4);
        lastSweeps_ = gaussSeidel(1.3, 1e-11, 50000);
        double maxd = 0.0;
        for (size_t n = 0; n < T_.size(); ++n) maxd = std::max(maxd, std::abs(T_[n] - Told_[n]));
        if (maxd < tolK) break;
    }
    step(1e4); // final consistent step so the ledger reflects the converged state
}

void ChipGrid::updatePeak() {
    peak_ = -1e9;
    for (int j = 0; j < p_.ny; ++j)
        for (int i = 0; i < p_.nx; ++i)
            if (siliconT(i, j) > peak_) { peak_ = siliconT(i, j); peakI_ = i; peakJ_ = j; }
}

double ChipGrid::meanSiliconT() const {
    double s = 0.0;
    for (int j = 0; j < p_.ny; ++j)
        for (int i = 0; i < p_.nx; ++i) s += siliconT(i, j);
    return s / (p_.nx * p_.ny);
}

double ChipGrid::outletT() const {
    double s = 0.0;
    for (int j = 0; j < p_.ny; ++j) s += fluidT(p_.nx - 1, j);
    return s / p_.ny;
}

double ChipGrid::meanFluidT() const {
    double s = 0.0;
    for (int j = 0; j < p_.ny; ++j)
        for (int i = 0; i < p_.nx; ++i) s += fluidT(i, j);
    return s / (p_.nx * p_.ny);
}

double ChipGrid::storedEnergy() const {
    double e = 0.0;
    for (size_t n = 0; n < T_.size(); ++n) e += C_[n] * T_[n];
    return e;
}

double ChipGrid::volumeFlow() const { return mdot_ / fluid_->props(meanFluidT()).rho; }

} // namespace mf
