#pragma once
#include "Scenes.hpp"
#include "workload/Workload.hpp"

#include <memory>
#include <vector>

namespace viewer {

struct App {
    mf::DataCenter* dc = nullptr;
    mf::Rack* rack = nullptr;
    std::vector<std::shared_ptr<mf::ConstantWorkload>> manual;
    ViewerState s;
    bool requestSteady = false;
    bool requestFrame = false;
};

/// Install the workloads for the selected profile (manual sliders or square waves).
void applyProfile(App& app);
/// Left control panel; applies any changed control to the rack.
void drawControlPanel(App& app, Rectangle area);
/// Right readout panel (system numbers + details of the selected scene).
void drawReadouts(App& app, Rectangle area);

} // namespace viewer
