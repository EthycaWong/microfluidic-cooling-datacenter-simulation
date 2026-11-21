#pragma once
#include <string>
#include <vector>

namespace mf {

/// How a chip's total power is spread over its grid. These spatial patterns are explicit study
/// ASSUMPTIONS, not measurements; a measured power map can replace them via "custom" weights.
struct PowerMapSpec {
    std::string pattern = "uniform"; ///< "uniform", "hotspot" or "custom"
    // hotspot: a rectangular core block (normalised 0..1 chip coordinates, x = flow direction)
    double centerX = 0.5, centerY = 0.5;
    double sizeX = 0.5, sizeY = 0.5;
    double coreShare = 0.5;          ///< fraction of total power dissipated inside the block
    std::vector<double> custom;      ///< nx*ny weights (row-major, j*nx + i), any scale
    int customNx = 0, customNy = 0;
};

/// Cell weights summing to 1 for an nx-by-ny grid (index j*nx + i, i along flow).
/// Hotspot weights use exact cell/block overlap areas, so they are consistent across resolutions.
std::vector<double> buildPowerWeights(const PowerMapSpec& spec, int nx, int ny);

} // namespace mf
