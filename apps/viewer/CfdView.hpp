#pragma once
#include "Scenes.hpp"
#include "cfd/Presets.hpp"

#include <memory>
#include <vector>

namespace viewer {

/// State of the interactive 2-D Navier-Stokes scene.
struct CfdApp {
    std::unique_ptr<mf::cfd::NavierStokes2D> solver;
    mf::cfd::Case kase;
    int caseIndex = 1;            ///< index into mf::cfd::caseNames()
    float value = -1.0f;          ///< Re / Ra / inlet velocity (slider)
    int wall = 3;                 ///< selected wall: 0 left, 1 right, 2 bottom, 3 top
    int colorMode = 0;            ///< 0 temperature, 1 speed, 2 pressure, 3 vorticity
    bool arrows = true, tracers = true, paused = false;
    float budgetMs = 14.0f;       ///< solver time per frame
    struct Particle { double x, y; float age; };
    std::vector<Particle> particles;
    double simPerFrame = 0.0;
};

void cfdLoadCase(CfdApp& c, int index, double value);
void cfdUpdate(CfdApp& c);
void drawCfdField(CfdApp& c, const ViewerState& s, Rectangle area);
/// Controls; returns the y position reached.
float drawCfdPanel(CfdApp& c, ViewerState& s, float x, float y, float w);
void drawCfdReadouts(CfdApp& c, const ViewerState& s, Rectangle area);

} // namespace viewer
