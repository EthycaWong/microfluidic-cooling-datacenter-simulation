#include "Scenes.hpp"
#include "Units.hpp"

#include "raymath.h"
#include "rlgl.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace viewer {
using mf::Rack;
using mf::Server;
using mf::ChipGrid;
using mf::Hotspot;

namespace {

constexpr float kU = static_cast<float>(mf::kRackUnit);
constexpr float kPlinth = 0.10f;           // rack floor height, m
constexpr float kRackW = 0.60f, kRackD = 1.00f;
constexpr float kSrvW = 0.44f, kSrvD = 0.78f;

const Color kFrame{70, 76, 86, 255};
const Color kChassis{46, 50, 58, 255};
const Color kBoard{24, 62, 48, 255};
const Color kCopper{196, 120, 70, 255};
const Color kInk{226, 232, 240, 255};
const Color kDim{150, 160, 175, 255};

std::string fmt(const char* f, double a) { char b[64]; std::snprintf(b, sizeof b, f, a); return b; }
std::string fmt2(const char* f, double a, double c) { char b[96]; std::snprintf(b, sizeof b, f, a, c); return b; }

float uBottom(int u) { return kPlinth + (u - 1) * kU; } // bottom of 1-based rack unit u

Color lerpColor(Color a, Color b, float t) {
    return Color{static_cast<unsigned char>(a.r + (b.r - a.r) * t), static_cast<unsigned char>(a.g + (b.g - a.g) * t),
                 static_cast<unsigned char>(a.b + (b.b - a.b) * t), 255};
}

void pipe(Vector3 a, Vector3 b, float r, Color c) { DrawCylinderEx(a, b, r, r, 10, c); }

/// Polyline hose.
void hose(const std::vector<Vector3>& pts, float r, Color c) {
    for (size_t i = 0; i + 1 < pts.size(); ++i) {
        pipe(pts[i], pts[i + 1], r, c);
        DrawSphere(pts[i + 1], r, c);
    }
}

/// Small moving dots along a segment to show flow direction and rate.
void flowDots(Vector3 a, Vector3 b, double flowFrac, float anim, Color c, int n = 3, float r = 0.008f) {
    if (flowFrac <= 1e-4) return;
    const float len = Vector3Distance(a, b);
    const float speed = 0.35f * static_cast<float>(flowFrac); // visual m/s
    for (int i = 0; i < n; ++i) {
        float f = std::fmod(anim * speed / std::max(len, 0.01f) + static_cast<float>(i) / n, 1.0f);
        DrawSphere(Vector3Lerp(a, b, f), r, c);
    }
}

void label(ViewerState& s, Vector3 p, const std::string& t, Color c = kInk, bool strong = false) {
    if (s.labels || strong) s.labelList.push_back(Label{p, t, c, strong});
}

const char* branchName(int b) { return b == 0 ? "GPU 1" : b == 1 ? "GPU 2" : "CPU"; }

} // namespace

// ------------------------------------------------------------------------------------------------
Color heatColor(double t, double lo, double hi) {
    // Turbo-like perceptual ramp: dark blue -> blue -> cyan -> green -> yellow -> orange -> red.
    static const Color stops[] = {{35, 23, 110, 255},  {58, 110, 230, 255}, {30, 190, 220, 255}, {70, 220, 120, 255},
                                  {225, 220, 60, 255}, {250, 140, 40, 255}, {200, 30, 30, 255}};
    const int n = sizeof(stops) / sizeof(stops[0]);
    double x = std::clamp((t - lo) / (hi - lo), 0.0, 1.0) * (n - 1);
    const int i = std::min(static_cast<int>(x), n - 2);
    return lerpColor(stops[i], stops[i + 1], static_cast<float>(x - i));
}

void drawText(const ViewerState& s, const std::string& t, float x, float y, float size, Color c) {
    DrawTextEx(s.font, t.c_str(), Vector2{x, y}, size, 0.5f, c);
}

BoundingBox serverBox(const Rack& r, int k) {
    const int u0 = r.params().firstServerU + k * r.params().serverU;
    const float y0 = uBottom(u0) + 0.002f, y1 = uBottom(u0 + r.params().serverU) - 0.002f;
    return BoundingBox{{-kSrvW / 2, y0, -kSrvD / 2}, {kSrvW / 2, y1, kSrvD / 2}};
}

void frameScene(ViewerState& s) {
    switch (s.scene) {
    case SceneId::Rack: s.target = {0.95f, 0.85f, -0.9f}; s.distance = 6.5f; s.yaw = 0.42f; s.pitch = 0.30f; break;
    case SceneId::Server: s.target = {0.04f, 0.02f, 0.02f}; s.distance = 1.45f; s.yaw = 0.55f; s.pitch = 0.80f; break;
    case SceneId::Chip: s.target = {0.05f, 0.18f, 0}; s.distance = 3.6f; s.yaw = 0.55f; s.pitch = 0.42f; break;
    case SceneId::Charts: break;
    }
}

// ================================================================================================
// RACK
// ================================================================================================
void drawRack3D(const Rack& r, ViewerState& s) {
    const auto& rp = r.params();
    const int n = r.serverCount();
    const float top = uBottom(rp.totalU + 1);

    // Floor and frame.
    DrawPlane({0, 0, 0}, {6, 6}, Color{22, 25, 31, 255});
    DrawGrid(24, 0.25f);
    DrawCube({0, kPlinth / 2, 0}, kRackW, kPlinth, kRackD, kFrame);
    DrawCube({0, top + 0.02f, 0}, kRackW, 0.04f, kRackD, kFrame);
    for (float x : {-kRackW / 2 + 0.02f, kRackW / 2 - 0.02f})
        for (float z : {-kRackD / 2 + 0.02f, kRackD / 2 - 0.02f})
            DrawCube({x, (kPlinth + top) / 2, z}, 0.04f, top - kPlinth, 0.04f, kFrame);
    // Mounting rails with U ticks (front).
    for (int u = 1; u <= rp.totalU + 1; ++u) {
        const float y = uBottom(u);
        DrawLine3D({-0.25f, y, kRackD / 2 - 0.03f}, {-0.24f, y, kRackD / 2 - 0.03f}, kDim);
    }

    // ---- CDU + pump in U1..U(first-1) ----
    const float cduTop = uBottom(rp.firstServerU) - 0.003f;
    DrawCube({0, (kPlinth + cduTop) / 2, 0}, kSrvW, cduTop - kPlinth - 0.004f, kSrvD, Color{40, 58, 72, 255});
    DrawCubeWires({0, (kPlinth + cduTop) / 2, 0}, kSrvW, cduTop - kPlinth - 0.004f, kSrvD, Color{90, 120, 140, 255});
    // Pump (front right), impeller spins with speed.
    const Vector3 pc{0.12f, kPlinth + 0.09f, kSrvD / 2 + 0.005f};
    DrawCylinderEx(pc, Vector3Add(pc, {0, 0, 0.02f}), 0.07f, 0.07f, 24, Color{60, 80, 96, 255});
    const float ang = s.anim * static_cast<float>(r.pumpSpeed()) * 18.0f;
    for (int b = 0; b < 5; ++b) {
        const float a = ang + b * 2.0f * PI / 5;
        DrawLine3D(Vector3Add(pc, {0, 0, 0.022f}), Vector3Add(pc, {0.06f * cosf(a), 0.06f * sinf(a), 0.022f}), kInk);
    }
    label(s, Vector3Add(pc, {0.12f, 0.02f, 0.05f}), fmt2("Pump %.0f%% speed, %.1f W", 100 * r.pumpSpeed(), r.pumpPower()));
    // Heat exchanger plates + filter + expansion vessel (front left).
    for (int i = 0; i < 9; ++i)
        DrawCube({-0.16f + i * 0.012f, kPlinth + 0.12f, kSrvD / 2 + 0.01f}, 0.006f, 0.16f, 0.02f,
                 lerpColor(heatColor(r.returnT(), kCoolantLo, kCoolantHi), heatColor(r.supplyT(), kCoolantLo, kCoolantHi), i / 8.0f));
    DrawCylinder({-0.02f, kPlinth + 0.05f, kSrvD / 2 + 0.03f}, 0.022f, 0.022f, 0.12f, 12, Color{120, 130, 140, 255});
    label(s, {-0.30f, kPlinth + 0.04f, kSrvD / 2 + 0.05f}, "CDU (heat exchanger + filter)");

    // ---- servers ----
    const Hotspot hs = r.hotspot();
    for (int k = 0; k < n; ++k) {
        const Server& sv = r.server(k);
        BoundingBox b = serverBox(r, k);
        const Vector3 c = Vector3Scale(Vector3Add(b.min, b.max), 0.5f);
        const Vector3 sz = Vector3Subtract(b.max, b.min);
        const bool sel = (k == s.server);
        DrawCube(c, sz.x, sz.y, sz.z, kChassis);
        DrawCubeWires(c, sz.x, sz.y, sz.z, sel ? Color{255, 255, 255, 255} : Color{80, 86, 96, 255});
        // Front bezel glows with the server's hottest chip, rear panel with its return temperature.
        const Color glow = heatColor(sv.peakT(), kSiliconLo, kSiliconHi);
        DrawCube({c.x, c.y, b.max.z + 0.004f}, sz.x - 0.01f, sz.y - 0.012f, 0.006f, glow);
        for (int v = 0; v < 6; ++v) // vent slots
            DrawCube({c.x - 0.15f + v * 0.06f, c.y, b.max.z + 0.008f}, 0.035f, sz.y - 0.05f, 0.002f, ColorAlpha(BLACK, 0.35f));
        // Rear: dark panel with a thin exhaust stripe in the server's return-coolant colour.
        DrawCube({c.x, c.y, b.min.z - 0.003f}, sz.x - 0.01f, sz.y - 0.012f, 0.004f, Color{36, 40, 47, 255});
        DrawCube({c.x, b.max.y - 0.012f, b.min.z - 0.006f}, sz.x - 0.04f, 0.008f, 0.003f, heatColor(sv.returnT(), kCoolantLo, kCoolantHi));
        std::string t = "S" + std::to_string(k + 1) + "  " + fmt("%.1f \xC2\xB0""C", sv.peakT()) +
                        fmt("  %.0f%% load", 100 * sv.activity());
        if (hs.server == k) t += "  (hottest)";
        label(s, {b.max.x + 0.03f, c.y, b.max.z}, t, sel ? Color{255, 255, 255, 255} : kInk, hs.server == k);
    }

    // ---- network switch + spare ----
    const int netU = rp.firstServerU + n * rp.serverU;
    for (int u = netU; u < netU + 2 && u <= rp.totalU; ++u) {
        const float y = uBottom(u) + kU / 2;
        DrawCube({0, y, 0.1f}, kSrvW, kU - 0.004f, 0.55f, Color{34, 38, 44, 255});
        for (int p = 0; p < 16; ++p) {
            const bool on = std::fmod(s.anim * (1.3f + 0.17f * p) + p * 0.37f, 1.0f) < 0.55f;
            DrawCube({-0.19f + p * 0.022f, y, 0.379f}, 0.012f, 0.012f, 0.004f, on ? Color{60, 230, 120, 255} : Color{30, 60, 40, 255});
        }
    }
    label(s, {kSrvW / 2 + 0.03f, uBottom(netU) + kU, 0.38f}, "Network switches (2U)", kDim);
    for (int u = netU + 2; u <= rp.totalU; ++u)
        DrawCube({0, uBottom(u) + kU / 2, kRackD / 2 - 0.05f}, kSrvW, kU - 0.003f, 0.01f, Color{28, 30, 34, 255});
    if (netU + 2 <= rp.totalU) label(s, {kSrvW / 2 + 0.03f, uBottom(netU + 2) + kU, 0.45f}, "Spare (blank panels)", kDim);

    // ---- rear manifolds ----
    // Manifolds run up the rear corners, beside (not behind) the servers, so they stay visible.
    const float zM = -kRackD / 2 + 0.06f;
    const float xS = 0.255f, xR = -0.255f;
    const float rM = 0.016f;
    const float yPump = kPlinth + 0.20f;
    const double fRef = r.totalMassFlow() > 0 ? r.totalMassFlow() : 1.0;
    const double flowNom = 0.2; // kg/s reference for dot speed
    const Color supC = heatColor(r.supplyT(), kCoolantLo, kCoolantHi);
    // Supply: riser from CDU then segments up to each tap.
    float yPrev = yPump;
    for (int k = 0; k < n; ++k) {
        const float y = static_cast<float>(r.serverTapZ(k)) + kPlinth;
        const Vector3 a{xS, yPrev, zM}, b{xS, y, zM};
        pipe(a, b, rM, supC);
        flowDots(a, b, r.supplySegmentFlow(k) / flowNom, s.anim, Color{230, 245, 255, 255});
        yPrev = y;
    }
    hose({{xS, yPump, zM}, {xS, yPump, -kSrvD / 2 + 0.02f}, {0.12f, yPump, -kSrvD / 2 + 0.02f}}, 0.012f, supC);
    // Return manifold.
    for (int k = 0; k < n; ++k) {
        const float y = static_cast<float>(r.serverTapZ(k)) + kPlinth;
        const Color rc = heatColor(r.returnSegmentT(k), kCoolantLo, kCoolantHi);
        Vector3 a{xR, y, zM}, b;
        if (rp.layout == mf::ManifoldLayout::DirectReturn)
            b = (k == 0) ? Vector3{xR, yPump, zM} : Vector3{xR, static_cast<float>(r.serverTapZ(k - 1)) + kPlinth, zM};
        else
            b = (k == n - 1) ? Vector3{xR, y + 0.10f, zM} : Vector3{xR, static_cast<float>(r.serverTapZ(k + 1)) + kPlinth, zM};
        pipe(a, b, rM, rc);
        flowDots(a, b, r.returnSegmentFlow(k) / flowNom, s.anim, Color{255, 235, 220, 255});
    }
    if (rp.layout == mf::ManifoldLayout::ReverseReturn) {
        const float yTop = static_cast<float>(r.serverTapZ(n - 1)) + kPlinth + 0.10f;
        const Color rc = heatColor(r.returnT(), kCoolantLo, kCoolantHi);
        hose({{xR, yTop, zM}, {xR - 0.07f, yTop, zM}, {xR - 0.07f, yPump, zM}, {xR, yPump, zM}}, 0.012f, rc);
        flowDots({xR - 0.07f, yTop, zM}, {xR - 0.07f, yPump, zM}, fRef / flowNom, s.anim, Color{255, 235, 220, 255});
    }
    hose({{xR, yPump, zM}, {xR, yPump, -kSrvD / 2 + 0.02f}, {-0.12f, yPump, -kSrvD / 2 + 0.02f}}, 0.012f,
         heatColor(r.returnT(), kCoolantLo, kCoolantHi));
    // Server hoses (QD pairs).
    for (int k = 0; k < n; ++k) {
        const Server& sv = r.server(k);
        const float y = static_cast<float>(r.serverTapZ(k)) + kPlinth;
        pipe({xS, y, zM}, {0.06f, y, -kSrvD / 2}, 0.007f, supC);
        pipe({xR, y, zM}, {-0.06f, y, -kSrvD / 2}, 0.007f, heatColor(sv.returnT(), kCoolantLo, kCoolantHi));
        const double rho = r.network().fluid().props(r.supplyT()).rho;
        label(s, {xS + 0.04f, y, zM - 0.02f}, fmt2("%.3f L/min  %.1f kPa", mf::m3sToLpm(sv.totalFlow() / rho),
                                                    sv.serverDp(r.network()) / 1e3), kDim);
    }
    label(s, {xS + 0.03f, yPrev + 0.08f, zM}, fmt("Supply manifold %.1f \xC2\xB0""C", r.supplyT()), supC, true);
    label(s, {xR - 0.03f, yPrev + 0.16f, zM}, fmt("Return manifold %.2f \xC2\xB0""C", r.returnT()),
          heatColor(r.returnT(), kCoolantLo, kCoolantHi), true);
}

// ================================================================================================
// PLANT (outdoor dry cooler, trim chiller, facility loop)
// ================================================================================================
void drawPlant3D(const mf::DataCenter& dc, ViewerState& s) {
    const mf::PlantState& p = dc.plant().state();
    // Rear building wall (translucent) behind the rack; outdoor pad beyond it.
    const float wallZ = -0.85f;
    DrawCube({1.5f, 0.01f, -1.75f}, 2.6f, 0.02f, 1.5f, Color{38, 46, 38, 255});
    label(s, {-0.6f, 2.25f, wallZ}, "building wall  |  outdoors behind", kDim);

    // Dry cooler: frame, coil coloured by facility water, fans on top spinning with speed.
    const Vector3 dcC{1.45f, 0.45f, -1.75f};
    const Color coilIn = heatColor(p.facilityReturnT, kCoolantLo, kCoolantHi);
    const Color coilOut = heatColor(p.dryCoolerOutT, kCoolantLo, kCoolantHi);
    DrawCube(dcC, 1.3f, 0.7f, 0.9f, Color{58, 64, 72, 255});
    DrawCubeWires(dcC, 1.3f, 0.7f, 0.9f, Color{110, 118, 130, 255});
    for (int i = 0; i < 14; ++i) { // coil fins, warm (inlet) to cool (outlet)
        const float x = dcC.x - 0.6f + i * 0.092f;
        const Color c = lerpColor(coilIn, coilOut, i / 13.0f);
        DrawCube({x, dcC.y, dcC.z + 0.46f}, 0.05f, 0.6f, 0.01f, c);
        DrawCube({x, dcC.y, dcC.z - 0.46f}, 0.05f, 0.6f, 0.01f, c);
    }
    for (int f = 0; f < 2; ++f) {
        const Vector3 c{dcC.x - 0.32f + f * 0.64f, dcC.y + 0.36f, dcC.z};
        DrawCylinder(c, 0.27f, 0.27f, 0.03f, 28, Color{30, 34, 40, 255});
        const float a0 = s.anim * static_cast<float>(p.fanSpeed) * 14.0f + f;
        for (int b = 0; b < 6; ++b) {
            const float a = a0 + b * 2 * PI / 6;
            DrawLine3D({c.x, c.y + 0.035f, c.z}, {c.x + 0.25f * cosf(a), c.y + 0.035f, c.z + 0.25f * sinf(a)}, kInk);
        }
        for (int k = 0; k < 4; ++k) { // warm air leaving upwards
            const float h = std::fmod(s.anim * (0.2f + 0.5f * static_cast<float>(p.fanSpeed)) + k * 0.25f, 1.0f);
            DrawSphere({c.x + 0.1f * cosf(k * 1.7f), c.y + 0.08f + 0.6f * h, c.z + 0.1f * sinf(k * 1.7f)}, 0.012f,
                       ColorAlpha(Color{255, 190, 120, 255}, 1.0f - h));
        }
    }
    label(s, {dcC.x - 0.4f, dcC.y + 1.05f, dcC.z}, fmt("Outdoor dry cooler   air %.1f \xC2\xB0""C", p.outdoorT), kInk, true);
    label(s, {dcC.x - 0.4f, dcC.y + 0.92f, dcC.z}, fmt2("fans %.0f%%  %.0f W", 100 * p.fanSpeed, p.fanPower) +
              fmt("   rejects %.0f W, evaporates no water", p.dryCoolerHeat), kDim);

    // Trim chiller beside it.
    const Vector3 chC{2.35f, 0.2f, -1.45f};
    const bool on = p.chillerHeat > 1.0;
    DrawCube(chC, 0.5f, 0.36f, 0.6f, on ? Color{40, 70, 120, 255} : Color{50, 54, 62, 255});
    DrawCubeWires(chC, 0.5f, 0.36f, 0.6f, on ? Color{120, 180, 255, 255} : Color{90, 96, 108, 255});
    label(s, {chC.x - 0.2f, chC.y + 0.3f, chC.z}, on ? fmt2("Trim chiller ON  %.0f W  COP %.1f", p.chillerPower, p.chillerCop)
                                                    : std::string(dc.plant().params().chiller.enabled ? "Trim chiller (standby)" : "Trim chiller (off)"),
          on ? Color{140, 200, 255, 255} : kDim);

    // Facility pipes: CDU (rack bottom, rear) -> through the wall -> dry cooler, and back.
    const float yS = 0.14f, yR = 0.24f, x0 = 0.18f, zOut = -0.42f;
    const Color fsC = heatColor(p.facilitySupplyT, kCoolantLo, kCoolantHi), frC = heatColor(p.facilityReturnT, kCoolantLo, kCoolantHi);
    const std::vector<Vector3> ret = {{x0, yR, zOut}, {0.55f, yR, zOut}, {0.55f, yR, dcC.z - 0.15f}, {dcC.x - 0.66f, yR, dcC.z - 0.15f}};
    const std::vector<Vector3> sup = {{dcC.x - 0.66f, yS, dcC.z + 0.15f}, {0.42f, yS, dcC.z + 0.15f}, {0.42f, yS, zOut}, {x0, yS, zOut}};
    hose(ret, 0.022f, frC);
    hose(sup, 0.022f, fsC);
    for (size_t i = 0; i + 1 < ret.size(); ++i) flowDots(ret[i], ret[i + 1], 1.0, s.anim, Color{255, 235, 220, 255}, 2, 0.01f);
    for (size_t i = 0; i + 1 < sup.size(); ++i) flowDots(sup[i], sup[i + 1], 1.0, s.anim, Color{230, 245, 255, 255}, 2, 0.01f);
    label(s, {0.6f, yR + 0.1f, -1.2f}, fmt("facility return %.1f \xC2\xB0""C", p.facilityReturnT), frC);
    label(s, {0.45f, yS - 0.06f, -1.3f}, fmt("facility supply %.1f \xC2\xB0""C", p.facilitySupplyT), fsC);
    label(s, {-0.30f, kPlinth - 0.03f, kSrvD / 2 + 0.05f}, fmt2("UA %.0f W/K, effectiveness %.2f", p.cduUA, p.cduEffectiveness), kDim);
    // Translucent wall last, so the equipment behind it still shows through.
    DrawCube({0.9f, 1.1f, wallZ}, 3.2f, 2.2f, 0.04f, ColorAlpha(Color{120, 130, 150, 255}, 0.12f));
    DrawCubeWires({0.9f, 1.1f, wallZ}, 3.2f, 2.2f, 0.04f, ColorAlpha(Color{150, 160, 180, 255}, 0.45f));
}

// ================================================================================================
// SERVER
// ================================================================================================
namespace {

/// Draw a chip's silicon temperature map as a grid of tiles on a horizontal plane.
void drawTileMap(const ChipGrid& c, Vector3 center, float sx, float sz, float h, int layer) {
    const int nx = c.nx(), ny = c.ny();
    const float dx = sx / nx, dz = sz / ny;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double t = layer == 0 ? c.siliconT(i, j) : layer == 1 ? c.copperT(i, j) : c.fluidT(i, j);
            const Color col = layer == 2 ? heatColor(t, kCoolantLo, kCoolantHi) : heatColor(t, kSiliconLo, kSiliconHi);
            const Vector3 p{center.x - sx / 2 + (i + 0.5f) * dx, center.y, center.z - sz / 2 + (j + 0.5f) * dz};
            DrawCube(p, dx * 0.96f, h, dz * 0.96f, col);
        }
}

} // namespace

void drawServer3D(const Rack& r, ViewerState& s) {
    const Server& sv = r.server(s.server);
    const float H = 4 * kU;
    const double rho = r.network().fluid().props(r.supplyT()).rho;
    DrawPlane({0, -0.002f, 0}, {3, 3}, Color{22, 25, 31, 255});
    DrawGrid(30, 0.05f);
    // Chassis: floor, low side walls (cut away), front bezel.
    DrawCube({0, 0.003f, 0}, kSrvW, 0.006f, kSrvD, kChassis);
    for (float x : {-kSrvW / 2, kSrvW / 2}) DrawCube({x, H / 2, 0}, 0.004f, H, kSrvD, ColorAlpha(kChassis, 0.55f));
    DrawCubeWires({0, H / 2, 0}, kSrvW, H, kSrvD, Color{90, 96, 108, 255});
    DrawCube({0, H / 2, kSrvD / 2}, kSrvW, H, 0.006f, ColorAlpha(heatColor(sv.peakT(), kSiliconLo, kSiliconHi), 0.35f));

    // Motherboard.
    DrawCube({0, 0.009f, 0}, kSrvW - 0.02f, 0.004f, kSrvD - 0.05f, kBoard);

    // Fans (front), spin with airflow.
    const float fanSpeed = static_cast<float>(sv.air().airflow() / sv.air().params().flowRef);
    for (int f = 0; f < 4; ++f) {
        const Vector3 c{-0.15f + f * 0.1f, H / 2, kSrvD / 2 - 0.05f};
        DrawCylinderEx(Vector3Add(c, {0, 0, -0.012f}), Vector3Add(c, {0, 0, 0.012f}), 0.04f, 0.04f, 20, Color{40, 44, 50, 255});
        const float a0 = s.anim * fanSpeed * 25.0f + f;
        for (int b = 0; b < 5; ++b) {
            const float a = a0 + b * 2 * PI / 5;
            DrawLine3D(Vector3Add(c, {0, 0, 0.014f}), Vector3Add(c, {0.036f * cosf(a), 0.036f * sinf(a), 0.014f}), kDim);
        }
    }
    label(s, {0.22f, H + 0.02f, kSrvD / 2 - 0.05f}, fmt2("Fans %.1f W, airflow %.1f L/s", sv.air().fanPower(), sv.air().airflow() * 1000));

    // Air-cooled group: DIMM banks each side of the CPU.
    const Color dimm = heatColor(sv.air().temperature(), kSiliconLo, kSiliconHi);
    for (int side : {-1, 1})
        for (int d = 0; d < 8; ++d)
            DrawCube({side * (0.09f + d * 0.011f), 0.04f, 0.16f}, 0.004f, 0.06f, 0.13f, dimm);
    label(s, {-0.16f, 0.09f, 0.16f}, fmt2("Memory/storage/other %.0f W  %.1f \xC2\xB0""C (air)", sv.air().power(), sv.air().temperature()));

    // Rear connectors.
    const Vector3 qdS{-0.05f, 0.05f, -kSrvD / 2}, qdR{0.05f, 0.05f, -kSrvD / 2};
    const Color supC = heatColor(r.supplyT(), kCoolantLo, kCoolantHi);
    DrawCube(qdS, 0.03f, 0.03f, 0.03f, supC);
    DrawCube(qdR, 0.03f, 0.03f, 0.03f, heatColor(sv.returnT(), kCoolantLo, kCoolantHi));
    label(s, Vector3Add(qdS, {-0.04f, 0.04f, -0.02f}), fmt("Supply in %.1f \xC2\xB0""C", r.supplyT()), supC, true);
    label(s, Vector3Add(qdR, {0.04f, 0.04f, -0.02f}), fmt2("Return out %.2f \xC2\xB0""C  (server dp %.1f kPa)", sv.returnT(),
                                                         sv.serverDp(r.network()) / 1e3),
          heatColor(sv.returnT(), kCoolantLo, kCoolantHi), true);

    // Components: GPU1, GPU2 (rear), CPU (front centre). Dies are drawn enlarged for visibility.
    struct Slot { Vector3 c; float disp; };
    const Slot slots[3] = {{{-0.11f, 0.02f, -0.12f}, 0.10f}, {{0.11f, 0.02f, -0.12f}, 0.10f}, {{0.0f, 0.02f, 0.16f}, 0.08f}};
    for (int b = 0; b < Server::kBranches; ++b) {
        const ChipGrid& c = sv.chip(b);
        const Slot& sl = slots[b];
        const bool sel = (b == s.branch);
        // module board / package
        DrawCube({sl.c.x, 0.014f, sl.c.z}, sl.disp + 0.05f, 0.006f, sl.disp + 0.05f, b < 2 ? Color{30, 34, 40, 255} : Color{50, 54, 60, 255});
        const float dieY = 0.022f;
        if (b == 2) { // cold plate on CPU
            DrawCube({sl.c.x, dieY + 0.012f, sl.c.z}, sl.disp + 0.012f, 0.022f, sl.disp + 0.012f, kCopper);
            drawTileMap(c, {sl.c.x, dieY + 0.024f, sl.c.z}, sl.disp, sl.disp, 0.003f, 1);
        } else {
            drawTileMap(c, {sl.c.x, dieY, sl.c.z}, sl.disp, sl.disp, 0.006f, 0);
            // microchannel inlet/outlet headers along the flow direction (x)
            DrawCube({sl.c.x - sl.disp / 2 - 0.008f, dieY, sl.c.z}, 0.012f, 0.012f, sl.disp, supC);
            DrawCube({sl.c.x + sl.disp / 2 + 0.008f, dieY, sl.c.z}, 0.012f, 0.012f, sl.disp, heatColor(c.outletT(), kCoolantLo, kCoolantHi));
        }
        if (sel) DrawCubeWires({sl.c.x, dieY + 0.01f, sl.c.z}, sl.disp + 0.03f, 0.04f, sl.disp + 0.03f, WHITE);
        // hoses: supply connector -> inlet side, outlet side -> return connector
        const float hy = 0.06f;
        const Vector3 in{sl.c.x - sl.disp / 2 - 0.01f, dieY + 0.01f, sl.c.z};
        const Vector3 out{sl.c.x + sl.disp / 2 + 0.01f, dieY + 0.01f, sl.c.z};
        const float lane = (b - 1) * 0.012f;
        hose({qdS, {qdS.x + lane, hy, qdS.z + 0.04f}, {in.x, hy, in.z + lane}, in}, 0.005f, supC);
        hose({out, {out.x, hy + 0.012f, out.z + lane}, {qdR.x + lane, hy + 0.012f, qdR.z + 0.04f}, qdR}, 0.005f,
             heatColor(c.outletT(), kCoolantLo, kCoolantHi));
        flowDots({in.x, hy, in.z + lane}, in, sv.branchFlow(b) / 0.0083, s.anim, WHITE, 2, 0.004f);
        std::string t = std::string(branchName(b)) + fmt("  %.0f W", c.power()) +
                        fmt("  %.3f L/min", mf::m3sToLpm(sv.branchFlow(b) / rho)) + fmt("  peak %.1f \xC2\xB0""C", c.peakT());
        label(s, {sl.c.x, 0.11f, sl.c.z}, t, sel ? WHITE : kInk, sel);
        label(s, {sl.c.x, 0.085f, sl.c.z}, b < 2 ? "silicon microchannels" : "conventional cold plate", kDim);
    }
}

// ================================================================================================
// CHIP
// ================================================================================================
void drawChip3D(const Rack& r, ViewerState& s) {
    const Server& sv = r.server(s.server);
    const ChipGrid& c = sv.chip(s.branch);
    const auto& p = c.params();
    const float scale = 1.0f / static_cast<float>(std::max(p.die.lx, p.die.ly)); // die -> ~1 unit
    const float sx = static_cast<float>(p.die.lx) * scale, sz = static_cast<float>(p.die.ly) * scale;
    const int nx = c.nx(), ny = c.ny();
    DrawPlane({0, -0.06f, 0}, {4, 4}, Color{22, 25, 31, 255});
    DrawGrid(20, 0.1f);

    // Layer 1: silicon active face (heat source), coloured by temperature.
    const float ySi = 0.0f;
    drawTileMap(c, {0, ySi, 0}, sx, sz, 0.05f, 0);
    label(s, {-sx / 2 - 0.05f, ySi, sz / 2 + 0.05f}, "Silicon (heated active face)", kInk, true);

    // Hotspot marker.
    const float dx = sx / nx, dz = sz / ny;
    const Vector3 hp{-sx / 2 + (c.peakI() + 0.5f) * dx, ySi + 0.03f, -sz / 2 + (c.peakJ() + 0.5f) * dz};
    DrawLine3D(hp, Vector3Add(hp, {0, 0.25f, 0}), WHITE);
    DrawSphere(Vector3Add(hp, {0, 0.25f, 0}), 0.015f, WHITE);
    label(s, Vector3Add(hp, {0.02f, 0.28f, 0}), fmt("Peak %.2f \xC2\xB0""C", c.peakT()) +
              "  cell (" + std::to_string(c.peakI()) + "," + std::to_string(c.peakJ()) + ")", WHITE, true);

    const float gap = 0.40f;
    if (p.mode == mf::CoolingMode::Microchannel) {
        // Layer 2 (exploded upward): the etched microchannels, each coloured along its length by
        // the coolant temperature of its grid column; silicon ribs between them.
        const auto& g = p.channels;
        const float yCh = ySi + gap;
        const float chW = static_cast<float>(g.width) * scale, chH = static_cast<float>(g.depth) * scale * 1.6f;
        const float pitch = sz / g.count;
        for (int k = 0; k < g.count; ++k) {
            const float z = -sz / 2 + (k + 0.5f) * pitch;
            const int j = std::min(ny - 1, static_cast<int>(k * ny / g.count));
            for (int i = 0; i < nx; ++i) {
                const Color col = ColorAlpha(heatColor(c.fluidT(i, j), kCoolantLo, kCoolantHi), 0.9f);
                DrawCube({-sx / 2 + (i + 0.5f) * dx, yCh, z}, dx * 1.001f, chH, chW, col);
            }
            if (k + 1 < g.count)
                DrawCube({0, yCh, z + pitch / 2}, sx, chH, pitch - chW, ColorAlpha(Color{150, 156, 168, 255}, 0.18f));
        }
        DrawCubeWires({0, yCh, 0}, sx, chH, sz, Color{120, 128, 140, 255});
        label(s, {-sx / 2 - 0.05f, yCh, sz / 2 + 0.05f},
              std::to_string(g.count) + " microchannels " + fmt2("%.0f x %.0f um", g.width * 1e6, g.depth * 1e6), kInk, true);
        // Cover lid (wireframe) above.
        DrawCubeWires({0, yCh + 0.10f, 0}, sx, 0.02f, sz, Color{110, 120, 135, 255});
        label(s, {sx / 2, yCh + 0.12f, -sz / 2}, "bonded cover (adiabatic)", kDim);
        // Headers and flow arrows.
        const Color cin = heatColor(c.inletT(), kCoolantLo, kCoolantHi), cout = heatColor(c.outletT(), kCoolantLo, kCoolantHi);
        DrawCube({-sx / 2 - 0.06f, yCh, 0}, 0.08f, chH * 1.6f, sz, ColorAlpha(cin, 0.85f));
        DrawCube({sx / 2 + 0.06f, yCh, 0}, 0.08f, chH * 1.6f, sz, ColorAlpha(cout, 0.85f));
        for (float z : {-0.3f, 0.0f, 0.3f}) {
            DrawCylinderEx({-sx / 2 - 0.35f, yCh, z * sz}, {-sx / 2 - 0.14f, yCh, z * sz}, 0.012f, 0.012f, 8, cin);
            DrawCylinderEx({-sx / 2 - 0.14f, yCh, z * sz}, {-sx / 2 - 0.1f, yCh, z * sz}, 0.03f, 0.0f, 8, cin);
            DrawCylinderEx({sx / 2 + 0.1f, yCh, z * sz}, {sx / 2 + 0.31f, yCh, z * sz}, 0.012f, 0.012f, 8, cout);
            DrawCylinderEx({sx / 2 + 0.31f, yCh, z * sz}, {sx / 2 + 0.35f, yCh, z * sz}, 0.03f, 0.0f, 8, cout);
        }
        label(s, {-sx / 2 - 0.3f, yCh + 0.06f, 0}, fmt("Inlet %.1f \xC2\xB0""C", c.inletT()), cin, true);
        label(s, {sx / 2 + 0.2f, yCh + 0.06f, 0}, fmt("Outlet %.2f \xC2\xB0""C", c.outletT()), cout, true);
        // Dotted lines linking layers at the corners (exploded-view guides).
        for (float x : {-sx / 2, sx / 2})
            for (float z : {-sz / 2, sz / 2}) DrawLine3D({x, ySi + 0.025f, z}, {x, yCh - chH / 2, z}, ColorAlpha(kDim, 0.6f));
    } else {
        // Cold plate: TIM + copper base layer, then the coolant passages.
        const float yCu = ySi + gap, yF = yCu + gap;
        drawTileMap(c, {0, yCu, 0}, sx, sz, 0.06f, 1);
        drawTileMap(c, {0, yF, 0}, sx, sz, 0.03f, 2);
        label(s, {-sx / 2 - 0.05f, yCu, sz / 2 + 0.05f}, "TIM + copper base", kInk, true);
        label(s, {-sx / 2 - 0.05f, yF, sz / 2 + 0.05f}, "cold-plate coolant passages", kInk, true);
        const Color cin = heatColor(c.inletT(), kCoolantLo, kCoolantHi), cout = heatColor(c.outletT(), kCoolantLo, kCoolantHi);
        DrawCylinderEx({-sx / 2 - 0.35f, yF, 0}, {-sx / 2 - 0.06f, yF, 0}, 0.02f, 0.02f, 10, cin);
        DrawCylinderEx({sx / 2 + 0.06f, yF, 0}, {sx / 2 + 0.35f, yF, 0}, 0.02f, 0.02f, 10, cout);
        label(s, {-sx / 2 - 0.3f, yF + 0.06f, 0}, fmt("Inlet %.1f \xC2\xB0""C", c.inletT()), cin, true);
        label(s, {sx / 2 + 0.2f, yF + 0.06f, 0}, fmt("Outlet %.2f \xC2\xB0""C", c.outletT()), cout, true);
        for (float x : {-sx / 2, sx / 2})
            for (float z : {-sz / 2, sz / 2}) DrawLine3D({x, ySi + 0.025f, z}, {x, yF, z}, ColorAlpha(kDim, 0.6f));
    }
    label(s, {0, ySi - 0.04f, sz / 2 + 0.12f}, "flow direction  \xE2\x86\x92", kDim, true);
}

// ================================================================================================
// 2D overlays
// ================================================================================================
void drawLabels(const ViewerState& s, Camera3D cam) {
    const Vector3 fwd = Vector3Normalize(Vector3Subtract(cam.target, cam.position));
    for (const Label& l : s.labelList) {
        if (Vector3DotProduct(Vector3Subtract(l.pos, cam.position), fwd) <= 0) continue;
        const Vector2 p = GetWorldToScreen(l.pos, cam);
        const float size = l.strong ? 17.0f : 15.0f;
        const Vector2 m = MeasureTextEx(s.font, l.text.c_str(), size, 0.5f);
        DrawRectangleRounded({p.x - 4, p.y - 2, m.x + 8, m.y + 4}, 0.3f, 4, Color{12, 14, 18, 170});
        drawText(s, l.text, p.x, p.y, size, l.color);
    }
}

void drawLegend(const ViewerState& s, Rectangle a) {
    auto bar = [&](float y, double lo, double hi, const char* name) {
        drawText(s, name, a.x, y, 15, kInk);
        const float w = a.width, h = 12, y0 = y + 20;
        for (int i = 0; i < 64; ++i) {
            const double t = lo + (hi - lo) * i / 63.0;
            DrawRectangleRec({a.x + w * i / 64.0f, y0, w / 64.0f + 1, h}, heatColor(t, lo, hi));
        }
        drawText(s, fmt("%.0f", lo), a.x, y0 + h + 2, 14, kDim);
        drawText(s, fmt("%.0f \xC2\xB0""C", hi), a.x + w - 40, y0 + h + 2, 14, kDim);
    };
    bar(a.y, kSiliconLo, kSiliconHi, "Silicon / components");
    bar(a.y + 50, kCoolantLo, kCoolantHi, "Coolant");
}

void drawHistory(const ViewerState& s, Rectangle a) {
    DrawRectangleRounded(a, 0.04f, 6, Color{14, 16, 20, 200});
    drawText(s, "Hottest chip (\xC2\xB0""C, red), pump power (W, blue), server 1 activity u (green): last 60 s", a.x + 10, a.y + 6, 15, kInk);
    if (s.history.size() < 2) return;
    const Rectangle p{a.x + 44, a.y + 30, a.width - 92, a.height - 44};
    DrawRectangleLinesEx(p, 1, Color{60, 66, 76, 255});
    const double t1 = s.history.back().t, t0 = t1 - 60.0;
    double pmax = 10.0;
    for (const auto& h : s.history) pmax = std::max(pmax, h.pumpW * 1.15);
    const double tLo = 30, tHi = 100;
    for (double tl : {40.0, 60.0, 85.0}) {
        const float y = p.y + p.height * static_cast<float>(1 - (tl - tLo) / (tHi - tLo));
        DrawLine(static_cast<int>(p.x), static_cast<int>(y), static_cast<int>(p.x + p.width), static_cast<int>(y),
                 tl == 85.0 ? Color{200, 60, 60, 120} : Color{50, 56, 66, 255});
        drawText(s, fmt("%.0f", tl), a.x + 10, y - 8, 13, kDim);
    }
    drawText(s, fmt("%.0f W", pmax), p.x + p.width + 6, p.y - 2, 13, Color{110, 170, 255, 255});
    Vector2 prevT{}, prevP{}, prevU{};
    bool first = true;
    for (const auto& h : s.history) {
        if (h.t < t0) continue;
        const float x = p.x + p.width * static_cast<float>((h.t - t0) / 60.0);
        const Vector2 vt{x, p.y + p.height * static_cast<float>(1 - std::clamp((h.peak - tLo) / (tHi - tLo), 0.0, 1.0))};
        const Vector2 vp{x, p.y + p.height * static_cast<float>(1 - std::clamp(h.pumpW / pmax, 0.0, 1.0))};
        const Vector2 vu{x, p.y + p.height * static_cast<float>(1 - std::clamp(h.u, 0.0, 1.0))};
        if (!first) {
            DrawLineEx(prevU, vu, 1.5f, Color{110, 220, 140, 200});
            DrawLineEx(prevT, vt, 2.0f, Color{255, 90, 80, 255});
            DrawLineEx(prevP, vp, 2.0f, Color{110, 170, 255, 255});
        }
        prevT = vt; prevP = vp; prevU = vu; first = false;
    }
}

} // namespace viewer
