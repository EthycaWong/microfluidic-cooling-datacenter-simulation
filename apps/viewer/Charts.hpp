#pragma once
#include "Scenes.hpp"

namespace viewer {

/// Append one history sample of every quantity the dashboard plots (call every ~0.1 s of sim time).
void recordSample(ViewerState& s, const mf::DataCenter& dc);

/// Full-area dashboard: six live panels (temperatures, IT power, cooling electricity by component,
/// flow and pressure, per-server state, selected chip heat map) plus an energy and efficiency summary.
void drawDashboard(ViewerState& s, const mf::DataCenter& dc, Rectangle area);

} // namespace viewer
