#include "hydraulics/HydraulicNetwork.hpp"
#include "Units.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mf {
namespace {

/// Dense Gaussian elimination with partial pivoting; A is n x n row-major, solved in place into b.
bool solveDense(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r)
            if (std::abs(A[r * n + c]) > std::abs(A[piv * n + c])) piv = r;
        if (std::abs(A[piv * n + c]) < 1e-300) return false;
        if (piv != c) {
            for (int k = 0; k < n; ++k) std::swap(A[c * n + k], A[piv * n + k]);
            std::swap(b[c], b[piv]);
        }
        for (int r = c + 1; r < n; ++r) {
            const double f = A[r * n + c] / A[c * n + c];
            if (f == 0.0) continue;
            for (int k = c; k < n; ++k) A[r * n + k] -= f * A[c * n + k];
            b[r] -= f * b[c];
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        double s = b[r];
        for (int k = r + 1; k < n; ++k) s -= A[r * n + k] * b[k];
        b[r] = s / A[r * n + r];
    }
    return true;
}

} // namespace

HydraulicNetwork::HydraulicNetwork(std::shared_ptr<const Fluid> fluid) : fluid_(std::move(fluid)) {}

int HydraulicNetwork::addNode(const std::string& name, double z) {
    nodes_.push_back(Node{name, z});
    return static_cast<int>(nodes_.size()) - 1;
}

int HydraulicNetwork::addEdge(int from, int to, std::unique_ptr<Edge> e) {
    edges_.push_back(E{from, to, std::move(e)});
    return static_cast<int>(edges_.size()) - 1;
}

void HydraulicNetwork::setReference(int node, double p) {
    ref_ = node;
    nodes_[node].p = p;
}

void HydraulicNetwork::setInjection(int node, double mdot) { nodes_[node].inj = mdot; }

double HydraulicNetwork::edgeDp(int e) const { return edges_[e].e->dp(edges_[e].mdot, *fluid_); }

double HydraulicNetwork::invert(E& e, double delta) const {
    const Fluid& f = *fluid_;
    auto fn = [&](double m) { return e.e->dp(m, f) - delta; };
    const double tol = 1e-10 * (std::abs(delta) + 1.0);

    // Bracket the root around the previous solution (dp is monotonic increasing).
    double m = e.mdot;
    double fm = fn(m);
    if (std::abs(fm) <= tol) {
        e.g = 1.0 / std::max(e.e->ddp(m, f), 1e-12);
        return m;
    }
    double step = std::max(std::abs(m) * 0.5, 1e-3);
    double lo, hi;
    if (fm < 0) {
        lo = m;
        hi = m + step;
        while (fn(hi) < 0) { lo = hi; step *= 2; hi += step; if (step > 1e6) break; }
    } else {
        hi = m;
        lo = m - step;
        while (fn(lo) > 0) { hi = lo; step *= 2; lo -= step; if (step > 1e6) break; }
    }
    // Safeguarded Newton.
    m = 0.5 * (lo + hi);
    for (int it = 0; it < 100; ++it) {
        const double fv = fn(m);
        if (std::abs(fv) <= tol) break;
        if (fv < 0) lo = m; else hi = m;
        const double d = e.e->ddp(m, f);
        double mn = (d > 0) ? m - fv / d : 0.5 * (lo + hi);
        if (!(mn > lo && mn < hi)) mn = 0.5 * (lo + hi);
        m = mn;
        if (hi - lo < 1e-15) break;
    }
    e.g = 1.0 / std::max(e.e->ddp(m, f), 1e-12);
    return m;
}

void HydraulicNetwork::computeFlows(const std::vector<double>& p) {
    for (auto& e : edges_) {
        const double dz = nodes_[e.to].z - nodes_[e.from].z;
        const double delta = p[e.from] - p[e.to] - rhoH_ * kGravity * dz;
        e.mdot = invert(e, delta);
    }
}

double HydraulicNetwork::imbalance(std::vector<double>& r) const {
    std::fill(r.begin(), r.end(), 0.0);
    for (int n = 0; n < nodeCount(); ++n) r[n] = -nodes_[n].inj;
    for (const auto& e : edges_) {
        r[e.from] += e.mdot;
        r[e.to] -= e.mdot;
    }
    double mx = 0.0;
    for (int n = 0; n < nodeCount(); ++n)
        if (n != ref_) mx = std::max(mx, std::abs(r[n]));
    return mx;
}

bool HydraulicNetwork::solve() {
    if (ref_ < 0) throw std::runtime_error("HydraulicNetwork: no reference node");
    const int nn = nodeCount();
    // Map non-reference nodes to unknown indices.
    std::vector<int> idx(nn, -1);
    int nu = 0;
    for (int n = 0; n < nn; ++n)
        if (n != ref_) idx[n] = nu++;

    std::vector<double> p(nn), r(nn), A, b;
    for (int n = 0; n < nn; ++n) p[n] = nodes_[n].p;
    computeFlows(p);
    double res = imbalance(r);
    const double tol = 1e-11;
    int it = 0;
    for (; it < 100 && res > tol; ++it) {
        A.assign(static_cast<size_t>(nu) * nu, 0.0);
        b.assign(nu, 0.0);
        for (const auto& e : edges_) {
            const int a = idx[e.from], c = idx[e.to];
            if (a >= 0) { A[a * nu + a] += e.g; if (c >= 0) A[a * nu + c] -= e.g; }
            if (c >= 0) { A[c * nu + c] += e.g; if (a >= 0) A[c * nu + a] -= e.g; }
        }
        for (int n = 0; n < nn; ++n)
            if (idx[n] >= 0) b[idx[n]] = -r[n];
        if (!solveDense(A, b, nu)) break;

        // Backtracking line search on the max imbalance.
        const std::vector<double> p0 = p;
        std::vector<double> m0(edges_.size());
        for (size_t k = 0; k < edges_.size(); ++k) m0[k] = edges_[k].mdot;
        double lambda = 1.0;
        for (int ls = 0; ls < 30; ++ls) {
            for (int n = 0; n < nn; ++n)
                if (idx[n] >= 0) p[n] = p0[n] + lambda * b[idx[n]];
            computeFlows(p);
            const double rn = imbalance(r);
            if (rn < res || ls == 29) { res = rn; break; }
            lambda *= 0.5;
            for (size_t k = 0; k < edges_.size(); ++k) edges_[k].mdot = m0[k];
        }
    }
    for (int n = 0; n < nn; ++n) nodes_[n].p = p[n];
    lastImbalance_ = res;
    lastIter_ = it;
    return res <= tol * 10;
}

} // namespace mf
