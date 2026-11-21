#include "thermal/PowerMap.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace mf {
namespace {
double overlap1d(double a0, double a1, double b0, double b1) {
    return std::max(0.0, std::min(a1, b1) - std::max(a0, b0));
}
} // namespace

std::vector<double> buildPowerWeights(const PowerMapSpec& s, int nx, int ny) {
    std::vector<double> w(static_cast<size_t>(nx) * ny, 1.0 / (nx * ny));
    if (s.pattern == "uniform") return w;

    if (s.pattern == "hotspot") {
        const double x0 = s.centerX - 0.5 * s.sizeX, x1 = s.centerX + 0.5 * s.sizeX;
        const double y0 = s.centerY - 0.5 * s.sizeY, y1 = s.centerY + 0.5 * s.sizeY;
        const double blockArea = overlap1d(x0, x1, 0, 1) * overlap1d(y0, y1, 0, 1);
        if (blockArea <= 0.0 || blockArea >= 1.0) return w;
        const double cellArea = 1.0 / (nx * ny);
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                const double ov = overlap1d(double(i) / nx, double(i + 1) / nx, x0, x1) *
                                  overlap1d(double(j) / ny, double(j + 1) / ny, y0, y1);
                w[j * nx + i] = s.coreShare * ov / blockArea + (1.0 - s.coreShare) * (cellArea - ov) / (1.0 - blockArea);
            }
        return w;
    }

    if (s.pattern == "custom") {
        if (s.customNx <= 0 || s.customNy <= 0 ||
            s.custom.size() != static_cast<size_t>(s.customNx) * s.customNy)
            throw std::runtime_error("custom power map needs customNx*customNy weights");
        // Area-weighted resampling of the custom map onto the simulation grid.
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                double acc = 0.0;
                for (int cj = 0; cj < s.customNy; ++cj)
                    for (int ci = 0; ci < s.customNx; ++ci) {
                        const double ov =
                            overlap1d(double(i) / nx, double(i + 1) / nx, double(ci) / s.customNx, double(ci + 1) / s.customNx) *
                            overlap1d(double(j) / ny, double(j + 1) / ny, double(cj) / s.customNy, double(cj + 1) / s.customNy);
                        acc += ov * s.custom[cj * s.customNx + ci] * s.customNx * s.customNy;
                    }
                w[j * nx + i] = acc;
            }
        const double sum = std::accumulate(w.begin(), w.end(), 0.0);
        for (double& v : w) v /= sum;
        return w;
    }
    throw std::runtime_error("unknown power map pattern: " + s.pattern);
}

} // namespace mf
