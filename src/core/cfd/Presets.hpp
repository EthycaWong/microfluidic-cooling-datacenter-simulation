#pragma once
#include "cfd/NavierStokes2D.hpp"

#include <string>
#include <vector>

namespace mf::cfd {

/// A ready-made boundary-driven case. Every wall's (u, v), temperature and flux remain editable.
struct Case {
    std::string name, description;
    Params params;
    double velocityScale = 1.0;  ///< for colouring, steadiness checks and Re
    double timeScale = 1.0;      ///< reference time for runToSteady
    double tCold = 0.0, tHot = 1.0;
    std::string reynoldsNote;    ///< e.g. "Re = 100" or "Ra = 1e5, Pr = 0.71"
};

/// Names: "lid-cavity", "heated-cavity", "rayleigh-benard", "mixed-convection",
///        "microchannel", "poiseuille". `value` is Re (lid, mixed, poiseuille), Ra (heated, RB) or
///        the mean inlet velocity in m/s (microchannel); <= 0 picks the default. n = cells on the short side.
Case makeCase(const std::string& name, double value = -1.0, int n = -1);
std::vector<std::string> caseNames();

} // namespace mf::cfd
