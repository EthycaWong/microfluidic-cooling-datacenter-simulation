// mfsim_viewer: interactive 3D view of the simulated rack, server and chip (raylib + raygui).
//
//   mfsim_viewer [--config configs/rack_baseline.json] [--scene rack|server|chip] [--server 1-8] [--chip 0-2]
//                [--screenshot out.png --frames 90]   (render, save a screenshot, exit)
#include "Scenes.hpp"
#include "UiPanel.hpp"
#include "Units.hpp"
#include "config/Config.hpp"

#include "raylib.h"
#include "raymath.h"
#include "raygui.h"

extern "C" void ViewerLoadDarkStyle(void);

#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

using namespace viewer;

namespace {

Font loadUiFont() {
    // Latin-1 plus a few symbols (degree, micro, arrows, superscript two).
    std::vector<int> cps;
    for (int c = 32; c < 256; ++c) cps.push_back(c);
    for (int c : {0x2192, 0x2190, 0x0394, 0x03BC, 0x2013, 0x00B2}) cps.push_back(c);
    const char* candidates[] = {"/System/Library/Fonts/Supplemental/Arial.ttf",
                                "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
                                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "C:/Windows/Fonts/arial.ttf"};
    for (const char* p : candidates) {
        if (!FileExists(p)) continue;
        Font f = LoadFontEx(p, 34, cps.data(), static_cast<int>(cps.size()));
        if (f.glyphCount > 0) {
            SetTextureFilter(f.texture, TEXTURE_FILTER_BILINEAR);
            return f;
        }
    }
    return GetFontDefault();
}

Camera3D cameraFor(const ViewerState& s) {
    Camera3D c{};
    c.target = s.target;
    c.position = Vector3Add(s.target, {s.distance * cosf(s.pitch) * sinf(s.yaw), s.distance * sinf(s.pitch),
                                        s.distance * cosf(s.pitch) * cosf(s.yaw)});
    c.up = {0, 1, 0};
    c.fovy = 40.0f;
    c.projection = CAMERA_PERSPECTIVE;
    return c;
}

std::string arg(int argc, char** argv, const char* key, const char* def) {
    for (int i = 1; i + 1 < argc; ++i)
        if (std::strcmp(argv[i], key) == 0) return argv[i + 1];
    return def;
}

} // namespace

int main(int argc, char** argv) {
    const std::string cfgPath = arg(argc, argv, "--config", "");
    const std::string shot = arg(argc, argv, "--screenshot", "");
    const int shotFrames = std::stoi(arg(argc, argv, "--frames", "90"));
    const std::string sceneArg = arg(argc, argv, "--scene", "rack");

    mf::SimConfig cfg = mf::defaultConfig();
    if (!cfgPath.empty()) cfg = mf::loadConfig(cfgPath);
    else if (FileExists("configs/rack_baseline.json")) cfg = mf::loadConfig("configs/rack_baseline.json");

    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE | FLAG_WINDOW_HIGHDPI);
    InitWindow(1600, 960, "Microfluidic cooling - data center simulation");
    SetTargetFPS(60);
    ViewerLoadDarkStyle();
    Font font = loadUiFont();
    GuiSetFont(font);
    GuiSetStyle(DEFAULT, TEXT_SIZE, 15);

    const bool plantOn = arg(argc, argv, "--plant", "on") != "off";
    mf::DataCenter dc(cfg.rack, cfg.server, cfg.plant, std::make_shared<mf::WaterIAPWS>(), plantOn);
    mf::Rack& rack = dc.rack();
    rack.setViscosityScale(cfg.viscosityScale);
    App app;
    app.dc = &dc;
    app.rack = &rack;
    app.s.plant = plantOn;
    app.s.outdoor = static_cast<float>(cfg.plant.weather.meanC);
    app.s.setpoint = static_cast<float>(cfg.plant.supplySetpointC);
    app.s.chiller = cfg.plant.chiller.enabled;
    const std::string outdoorArg = arg(argc, argv, "--outdoor", "");
    if (!outdoorArg.empty()) { app.s.outdoor = std::stof(outdoorArg); dc.plant().setOutdoorMean(app.s.outdoor); }
    app.s.font = font;
    for (int k = 0; k < rack.serverCount(); ++k) app.manual.push_back(std::make_shared<mf::ConstantWorkload>(1.0));
    app.s.viscosity = static_cast<float>(cfg.viscosityScale);
    app.s.supplyT = static_cast<float>(cfg.rack.supplyT);
    app.s.targetFlow = static_cast<float>(cfg.rack.targetFlowLpm);
    app.s.targetPeak = static_cast<float>(cfg.rack.targetPeakC);
    app.s.layout = cfg.rack.layout == mf::ManifoldLayout::ReverseReturn ? 1 : 0;
    app.s.pumpMode = static_cast<int>(cfg.rack.pumpMode);
    app.s.server = std::clamp(std::stoi(arg(argc, argv, "--server", "8")) - 1, 0, rack.serverCount() - 1);
    app.s.branch = std::clamp(std::stoi(arg(argc, argv, "--chip", "0")), 0, 2);
    app.s.scene = sceneArg == "chip" ? SceneId::Chip : sceneArg == "server" ? SceneId::Server : SceneId::Rack;
    const std::string loadArg = arg(argc, argv, "--loads", "");
    if (!loadArg.empty()) { // e.g. --loads 1,1,1,1,0.2,0.2,0.2,0.2
        int k = 0;
        for (const char* p = loadArg.c_str(); *p && k < 8; ++k) {
            app.s.load[k] = std::strtof(p, const_cast<char**>(&p));
            if (*p == ',') ++p;
        }
    }
    for (int k = 0; k < 8; ++k) app.manual[k]->set(app.s.load[k]);
    applyProfile(app);
    frameScene(app.s);
    dc.solveSteady(); // start from a converged state
    app.s.pumpSpeed = static_cast<float>(rack.pumpSpeed());

    const float panelW = 360.0f;
    double lastSample = -1.0, lastClick = -1.0;
    int frame = 0;
    while (!WindowShouldClose()) {
        const float sw = static_cast<float>(GetScreenWidth()), sh = static_cast<float>(GetScreenHeight());
        const Vector2 mouse = GetMousePosition();
        const Rectangle readoutArea{sw - 380, 12, 368, 750};
        const bool overUi = mouse.x < panelW || CheckCollisionPointRec(mouse, readoutArea);

        // ---- camera input ----
        if (!overUi) {
            if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
                const Vector2 d = GetMouseDelta();
                app.s.yaw -= d.x * 0.006f;
                app.s.pitch = std::clamp(app.s.pitch + d.y * 0.006f, -0.2f, 1.45f);
            }
            if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
                const Vector2 d = GetMouseDelta();
                const float k = app.s.distance * 0.0015f;
                const Vector3 right{cosf(app.s.yaw), 0, -sinf(app.s.yaw)};
                app.s.target = Vector3Add(app.s.target, Vector3Scale(right, -d.x * k));
                app.s.target.y += d.y * k;
            }
            const float wheel = GetMouseWheelMove();
            if (wheel != 0) app.s.distance = std::clamp(app.s.distance * (1.0f - 0.1f * wheel), 0.25f, 14.0f);
            if (app.s.scene == SceneId::Rack && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
                const Ray ray = GetScreenToWorldRay(mouse, cameraFor(app.s));
                for (int k = 0; k < rack.serverCount(); ++k)
                    if (GetRayCollisionBox(ray, serverBox(rack, k)).hit) {
                        if (k == app.s.server && GetTime() - lastClick < 0.35) {
                            app.s.scene = SceneId::Server;
                            frameScene(app.s);
                        }
                        app.s.server = k;
                        lastClick = GetTime();
                        break;
                    }
            }
        }
        if (IsKeyPressed(KEY_SPACE)) app.s.paused = !app.s.paused;
        if (IsKeyPressed(KEY_ONE)) { app.s.scene = SceneId::Rack; frameScene(app.s); }
        if (IsKeyPressed(KEY_TWO)) { app.s.scene = SceneId::Server; frameScene(app.s); }
        if (IsKeyPressed(KEY_THREE)) { app.s.scene = SceneId::Chip; frameScene(app.s); }
        if (app.requestFrame) { frameScene(app.s); app.requestFrame = false; }
        if (app.requestSteady) { dc.solveSteady(); app.requestSteady = false; }

        // ---- physics: fixed implicit sub-steps decoupled from the frame rate ----
        if (!app.s.paused && GetFrameTime() > 0.0f) {
            const double simDt = std::min<double>(GetFrameTime(), 0.05) * app.s.simSpeed;
            const int n = std::clamp(static_cast<int>(std::ceil(simDt / 0.02)), 1, 40);
            for (int i = 0; i < n; ++i) dc.step(simDt / n);
            app.s.anim += static_cast<float>(simDt);
        }
        if (app.s.pumpMode == 2) app.s.pumpSpeed = static_cast<float>(rack.pumpSpeed());
        if (rack.time() - lastSample >= 0.1) {
            app.s.history.push_back({rack.time(), rack.hotspot().t, rack.pumpPower(), rack.totalFlowLpm(), rack.returnT()});
            while (!app.s.history.empty() && app.s.history.front().t < rack.time() - 62.0) app.s.history.pop_front();
            lastSample = rack.time();
        }

        // ---- draw ----
        BeginDrawing();
        ClearBackground(Color{18, 20, 26, 255});
        const Camera3D cam = cameraFor(app.s);
        app.s.labelList.clear();
        BeginMode3D(cam);
        switch (app.s.scene) {
        case SceneId::Rack:
            drawRack3D(rack, app.s);
            if (dc.plantEnabled()) drawPlant3D(dc, app.s);
            break;
        case SceneId::Server: drawServer3D(rack, app.s); break;
        case SceneId::Chip: drawChip3D(rack, app.s); break;
        }
        EndMode3D();
        drawLabels(app.s, cam);
        const char* titles[] = {"42U rack + CDU, closed loop to an outdoor dry cooler (no evaporation)",
                                "Server: supply \xE2\x86\x92 GPU 1 | GPU 2 | CPU in parallel \xE2\x86\x92 return  (dies enlarged)",
                                "Chip: exploded view, silicon heat map and microchannel coolant"};
        drawText(app.s, titles[static_cast<int>(app.s.scene)], panelW + 20, 16, 20, Color{230, 236, 245, 255});
        drawControlPanel(app, {0, 0, panelW, sh});
        drawReadouts(app, readoutArea);
        drawHistory(app.s, {panelW + 20, sh - 160, sw - panelW - 420, 146});
        drawLegend(app.s, {sw - 360, sh - 130, 330, 110});
        EndDrawing();

        if (!shot.empty() && ++frame >= shotFrames) {
            Image img = LoadImageFromScreen();
            ExportImage(img, shot.c_str());
            UnloadImage(img);
            break;
        }
    }
    if (font.texture.id != GetFontDefault().texture.id) UnloadFont(font);
    CloseWindow();
    return 0;
}
