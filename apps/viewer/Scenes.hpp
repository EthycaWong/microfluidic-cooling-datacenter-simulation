#pragma once
#include "raylib.h"
#include "system/DataCenter.hpp"

#include <deque>
#include <string>
#include <vector>

namespace viewer {

enum class SceneId { Rack = 0, Server = 1, Chip = 2 };

/// A text label anchored at a 3D point, drawn after the 3D pass.
struct Label {
    Vector3 pos;
    std::string text;
    Color color;
    bool strong = false;
};

struct HistorySample { double t, peak, pumpW, flowLpm, returnT; };

/// Everything the UI can change, plus scene/camera state.
struct ViewerState {
    SceneId scene = SceneId::Rack;
    int server = 7;            ///< selected server (0-based, 7 = top)
    int branch = 0;            ///< selected chip in the server (0 GPU1, 1 GPU2, 2 CPU)
    bool paused = false;
    float simSpeed = 1.0f;     ///< simulated seconds per real second
    bool plant = true;         ///< CDU + dry cooler model on (off = fixed supply temperature)
    float outdoor = 25.0f;
    float setpoint = 35.0f;
    bool chiller = true;
    int profile = 0;           ///< 0 manual, 1 square wave, 2 staggered square
    float load[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    float loadAll = 1.0f;
    float viscosity = 1.0f;
    float supplyT = 35.0f;
    int pumpMode = 1;          ///< 0 speed, 1 flow, 2 peak temperature
    float pumpSpeed = 0.94f;
    float targetFlow = 12.0f;
    float targetPeak = 85.0f;
    int layout = 0;            ///< 0 direct, 1 reverse return
    bool balance = true;
    bool labels = true;
    // camera orbit
    float yaw = 0.75f, pitch = 0.32f, distance = 3.6f;
    Vector3 target{0, 0.95f, 0};
    float anim = 0.0f;         ///< animation clock (s), advances with simulated time
    std::vector<Label> labelList;
    std::deque<HistorySample> history;
    Font font{};
};

// Temperature colour scales shared by every scene and the legend.
constexpr double kSiliconLo = 35.0, kSiliconHi = 90.0;
constexpr double kCoolantLo = 28.0, kCoolantHi = 48.0;
Color heatColor(double t, double lo, double hi);

void drawRack3D(const mf::Rack& r, ViewerState& s);
/// Outdoor dry cooler, trim chiller and facility piping next to the rack (plant enabled).
void drawPlant3D(const mf::DataCenter& dc, ViewerState& s);
void drawServer3D(const mf::Rack& r, ViewerState& s);
void drawChip3D(const mf::Rack& r, ViewerState& s);
/// Bounding box of server k in the rack scene, for mouse picking.
BoundingBox serverBox(const mf::Rack& r, int k);
/// Default camera framing for each scene.
void frameScene(ViewerState& s);

void drawLabels(const ViewerState& s, Camera3D cam);
void drawLegend(const ViewerState& s, Rectangle area);
void drawHistory(const ViewerState& s, Rectangle area);
void drawText(const ViewerState& s, const std::string& t, float x, float y, float size, Color c);

} // namespace viewer
