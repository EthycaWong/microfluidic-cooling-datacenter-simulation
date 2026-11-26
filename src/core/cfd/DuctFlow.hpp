#pragma once
#include <vector>

namespace mf::cfd {

/// Fully developed laminar flow and heat transfer in a rectangular duct, from first principles.
///
/// Momentum (Navier-Stokes reduced to fully developed flow, w = w(x, y), u = v = 0):
///     mu (d2w/dx2 + d2w/dy2) = dp/dz,  w = 0 on the walls (no-slip)
/// Energy, H1 condition (axially uniform heat input, peripherally uniform wall temperature):
///     alpha lap(theta) = w dTb/dz,  theta = 0 on heated walls, d(theta)/dn = 0 on an insulated wall
/// Both Poisson problems are solved on an n x m cell-centred grid with conjugate gradients.
/// Results are the Fanning friction factor x Reynolds number and the Nusselt numbers (based on
/// the hydraulic diameter), which the chip model otherwise takes from Shah & London's tables.
struct DuctResult {
    double fRe = 0;           ///< Fanning f * Re
    double nuH1FourWalls = 0; ///< all four walls heated
    double nuH1ThreeWalls = 0;///< top wall (length = width) insulated, as under a bonded cover
    double maxOverMean = 0;   ///< w_max / w_mean
    int nx = 0, ny = 0;
    std::vector<double> w;    ///< normalised velocity w / w_mean, row-major (j * nx + i)
    std::vector<double> theta;///< normalised 3-wall temperature (T - Tw)/(Tb - Tw)
};

/// width = horizontal side (with the insulated top wall), height = vertical side; cells on the
/// shorter side = nShort (the longer side gets proportionally more, square cells).
DuctResult solveRectDuct(double width, double height, int nShort = 64);

} // namespace mf::cfd
