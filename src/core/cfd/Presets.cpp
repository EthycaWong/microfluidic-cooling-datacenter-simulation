#include "cfd/Presets.hpp"
#include "fluid/Fluid.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace mf::cfd {
namespace {
std::string fmt(const char* f, double a) { char b[96]; std::snprintf(b, sizeof b, f, a); return b; }
Side fixedT(double t) { Side s; s.thermal = ThermalBC::FixedTemperature; s.T = t; return s; }
} // namespace

std::vector<std::string> caseNames() {
    return {"lid-cavity", "heated-cavity", "rayleigh-benard", "mixed-convection", "microchannel", "poiseuille"};
}

Case makeCase(const std::string& name, double value, int n) {
    Case c;
    c.name = name;
    Params& p = c.params;
    if (name == "lid-cavity") {
        // Square box, all walls no-slip (u = v = 0) except the lid moving at u = 1. Benchmark: Ghia et al. (1982).
        const double re = value > 0 ? value : 100.0;
        p.nx = p.ny = n > 0 ? n : 64;
        p.nu = 1.0 / re;
        p.alpha = p.nu;
        p.top.u = 1.0;
        c.reynoldsNote = "Re = " + fmt("%.0f", re);
        c.description = "Lid-driven cavity: walls u = v = 0, lid u = 1 (Ghia et al. 1982 benchmark)";
    } else if (name == "heated-cavity") {
        // Hot left wall, cold right wall, insulated top/bottom, all no-slip; buoyancy drives the flow.
        // Benchmark: de Vahl Davis (1983). Non-dimensional: L = 1, alpha = 1, nu = Pr, g beta dT = Ra Pr.
        const double ra = value > 0 ? value : 1e5;
        const double pr = 0.71;
        p.nx = p.ny = n > 0 ? n : 64;
        p.alpha = 1.0;
        p.nu = pr;
        p.beta = 1.0;
        p.gy = -ra * pr;
        p.tRef = 0.5;
        p.t0 = 0.5;
        p.left = fixedT(1.0);
        p.right = fixedT(0.0);
        c.velocityScale = std::sqrt(ra);
        c.timeScale = 1.0;
        c.reynoldsNote = "Ra = " + fmt("%.0e", ra) + ", Pr = 0.71";
        c.description = "Differentially heated cavity: hot left (T = 1), cold right (T = 0), u = v = 0 on all walls";
    } else if (name == "rayleigh-benard") {
        // Hot floor, cold lid, insulated sides: convection rolls above the critical Rayleigh number (~1708).
        const double ra = value > 0 ? value : 5e4;
        const double pr = 0.71;
        p.lx = 2.0;
        p.ly = 1.0;
        p.ny = n > 0 ? n : 40;
        p.nx = 2 * p.ny;
        p.alpha = 1.0;
        p.nu = pr;
        p.beta = 1.0;
        p.gy = -ra * pr;
        p.tRef = 0.5;
        p.t0 = 0.5;
        p.bottom = fixedT(1.0);
        p.top = fixedT(0.0);
        c.velocityScale = std::sqrt(ra);
        c.reynoldsNote = "Ra = " + fmt("%.0e", ra) + ", Pr = 0.71";
        c.description = "Rayleigh-Benard: hot floor, cold lid, insulated sides, all walls u = v = 0";
    } else if (name == "mixed-convection") {
        // Moving lid over a hot floor and under it a cold lid: shear and buoyancy together.
        const double re = value > 0 ? value : 100.0;
        const double ri = 1.0; // Richardson number Gr/Re^2
        p.nx = p.ny = n > 0 ? n : 64;
        p.nu = 1.0 / re;
        p.alpha = p.nu / 0.71;
        p.beta = 1.0;
        p.gy = -ri;            // g beta dT L / U^2 = Ri with U = L = dT = 1
        p.tRef = 0.5;
        p.t0 = 0.5;
        p.top = fixedT(0.0);
        p.top.u = 1.0;
        p.bottom = fixedT(1.0);
        c.reynoldsNote = "Re = " + fmt("%.0f", re) + ", Ri = 1, Pr = 0.71";
        c.description = "Mixed convection: cold moving lid (u = 1) over a hot floor, side walls u = v = 0";
    } else if (name == "poiseuille") {
        const double re = value > 0 ? value : 100.0;
        p.lx = 10.0;
        p.ly = 1.0;
        p.ny = n > 0 ? n : 20;
        p.nx = 10 * p.ny;
        p.nu = 2.0 / re;       // Re on hydraulic diameter 2h with mean velocity 1
        p.alpha = p.nu / 0.71;
        p.left = fixedT(0.0);
        p.left.u = 1.0;
        p.right.velocity = VelocityBC::Outflow;
        p.bottom = fixedT(1.0);
        p.top = fixedT(1.0);
        c.timeScale = 10.0;
        c.reynoldsNote = "Re = " + fmt("%.0f", re);
        c.description = "Channel between hot plates: uniform cold inflow, outflow right, walls u = v = 0";
    } else if (name == "microchannel") {
        // Plan-view slice of one silicon microchannel: the 200 um gap between two ribs, 20 mm long.
        // Ribs are silicon at a fixed hot temperature; water enters cold at 35 C. Real water properties.
        const double um = value > 0 ? value : 1.667; // mean velocity at 0.5 L/min per GPU
        const WaterIAPWS water;
        const FluidProps w = water.props(40.0);
        p.lx = 20e-3;
        p.ly = 200e-6;
        p.ny = n > 0 ? n : 16;
        p.nx = 25 * p.ny;      // cells 3.1x longer than tall (flow is nearly parallel)
        p.nu = w.mu / w.rho;
        p.alpha = w.k / (w.rho * w.cp);
        p.k = w.k;
        p.t0 = 35.0;
        p.gammaMin = 0.5;
        p.left = fixedT(35.0);
        p.left.u = um;
        p.right.velocity = VelocityBC::Outflow;
        p.bottom = fixedT(65.0);
        p.top = fixedT(65.0);
        c.velocityScale = um;
        c.timeScale = p.lx / um;
        c.tCold = 35.0;
        c.tHot = 65.0;
        c.reynoldsNote = "Re_Dh = " + fmt("%.0f", um * 2 * p.ly / p.nu) + ", water at 40 C";
        c.description = "Microchannel slice: 200 um gap between hot silicon ribs (65 C), water in at 35 C";
    } else {
        throw std::runtime_error("unknown CFD case " + name);
    }
    if (c.name != "microchannel") {
        double tmin = 1e9, tmax = -1e9;
        for (const Side* s : {&p.left, &p.right, &p.bottom, &p.top})
            if (s->thermal == ThermalBC::FixedTemperature) { tmin = std::min(tmin, s->T); tmax = std::max(tmax, s->T); }
        if (tmax > tmin) { c.tCold = tmin; c.tHot = tmax; }
    }
    return c;
}

} // namespace mf::cfd
