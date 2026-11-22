#include "UiPanel.hpp"
#include "Units.hpp"
#include "hydraulics/RectChannel.hpp"

#include "raygui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace viewer {
namespace {

const Color kInk{226, 232, 240, 255};
const Color kDim{150, 160, 175, 255};
const Color kHead{120, 190, 255, 255};

std::string f1(const char* f, double a) { char b[96]; std::snprintf(b, sizeof b, f, a); return b; }

void heading(const ViewerState& s, const char* t, float x, float& y, float w) {
    y += 6;
    drawText(s, t, x, y, 16, kHead);
    DrawLine(static_cast<int>(x), static_cast<int>(y + 20), static_cast<int>(x + w), static_cast<int>(y + 20), Color{60, 70, 84, 255});
    y += 26;
}

/// Labelled slider; returns true when the value changed.
bool slider(const ViewerState& s, const char* name, const std::string& value, float x, float& y, float w, float* v,
            float lo, float hi) {
    drawText(s, name, x, y + 3, 15, kInk);
    const float before = *v;
    GuiSlider({x + 112, y, w - 112 - 62, 20}, nullptr, nullptr, v, lo, hi);
    drawText(s, value, x + w - 56, y + 3, 15, kDim);
    y += 27;
    return *v != before;
}

void row(const ViewerState& s, const char* k, const std::string& v, float x, float& y, float w, Color c = kInk) {
    drawText(s, k, x, y, 15, kDim);
    const Vector2 m = MeasureTextEx(s.font, v.c_str(), 15, 0.5f);
    drawText(s, v, x + w - m.x, y, 15, c);
    y += 20;
}

const char* chipName(int b) { return b == 0 ? "GPU 1" : b == 1 ? "GPU 2" : "CPU"; }

} // namespace

void applyProfile(App& app) {
    mf::Rack& r = *app.rack;
    for (int k = 0; k < r.serverCount(); ++k) {
        switch (app.s.profile) {
        case 0: r.setWorkload(k, app.manual[k]); break;
        case 1: r.setWorkload(k, std::make_shared<mf::SquareWorkload>(0.2, 1.0, 20.0, 0.5, 0.0)); break;
        default: r.setWorkload(k, std::make_shared<mf::SquareWorkload>(0.2, 1.0, 20.0, 0.5, 2.5 * k)); break;
        }
    }
}

void drawControlPanel(App& app, Rectangle a) {
    ViewerState& s = app.s;
    mf::Rack& r = *app.rack;
    DrawRectangleRec(a, Color{16, 18, 23, 235});
    DrawLine(static_cast<int>(a.x + a.width), 0, static_cast<int>(a.x + a.width), static_cast<int>(a.height), Color{50, 58, 70, 255});
    const float x = a.x + 16, w = a.width - 32;
    float y = a.y + 14;
    drawText(s, "Microfluidic cooling", x, y, 22, kInk);
    y += 26;
    drawText(s, "chip \xE2\x86\x92 server \xE2\x86\x92 42U rack, live physics", x, y, 14, kDim);
    y += 24;

    // ---- view ----
    int scene = static_cast<int>(s.scene);
    GuiToggleGroup({x, y, (w - 8) / 3, 28}, "Rack;Server;Chip", &scene);
    if (scene != static_cast<int>(s.scene)) { s.scene = static_cast<SceneId>(scene); frameScene(s); }
    y += 36;
    drawText(s, "Server", x, y + 4, 15, kDim);
    int srv = s.server;
    GuiToggleGroup({x + 60, y, (w - 60 - 14) / 8, 24}, "1;2;3;4;5;6;7;8", &srv);
    s.server = srv;
    y += 30;
    drawText(s, "Chip", x, y + 4, 15, kDim);
    int br = s.branch;
    GuiToggleGroup({x + 60, y, (w - 60 - 4) / 3, 24}, "GPU 1;GPU 2;CPU", &br);
    s.branch = br;
    y += 30;

    // ---- simulation ----
    heading(s, "Simulation", x, y, w);
    if (GuiButton({x, y, 90, 26}, s.paused ? "Resume" : "Pause")) s.paused = !s.paused;
    if (GuiButton({x + 98, y, 120, 26}, "Steady state")) app.requestSteady = true;
    if (GuiButton({x + 226, y, w - 226, 26}, "Reset view")) app.requestFrame = true;
    y += 34;
    slider(s, "Speed", f1("%.1fx", s.simSpeed), x, y, w, &s.simSpeed, 0.1f, 60.0f);

    // ---- workload ----
    heading(s, "Workload (activity u)", x, y, w);
    int prof = s.profile;
    GuiComboBox({x, y, w, 26}, "Manual sliders;Square wave, in sync;Square wave, staggered", &prof);
    if (prof != s.profile) { s.profile = prof; applyProfile(app); }
    y += 34;
    if (s.profile == 0) {
        if (slider(s, "All servers", f1("%.0f%%", 100 * s.loadAll), x, y, w, &s.loadAll, 0.0f, 1.0f))
            for (int k = 0; k < 8; ++k) s.load[k] = s.loadAll;
        for (int k = 7; k >= 0; --k) {
            char name[24];
            std::snprintf(name, sizeof name, "  Server %d", k + 1);
            const float before = s.load[k];
            drawText(s, name, x, y + 1, 14, kDim);
            GuiSlider({x + 112, y, w - 112 - 62, 16}, nullptr, nullptr, &s.load[k], 0.0f, 1.0f);
            drawText(s, f1("%.0f%%", 100 * s.load[k]), x + w - 56, y + 1, 14, kDim);
            if (s.load[k] != before) app.manual[k]->set(s.load[k]);
            y += 19;
        }
        for (int k = 0; k < 8; ++k) app.manual[k]->set(s.load[k]);
    } else {
        drawText(s, "u alternates 20% / 100% every 10 s", x, y, 14, kDim);
        y += 22;
    }

    // ---- coolant ----
    heading(s, "Coolant", x, y, w);
    if (slider(s, "Viscosity", f1("%.2fx", s.viscosity), x, y, w, &s.viscosity, 0.5f, 5.0f)) r.setViscosityScale(s.viscosity);
    drawText(s, "viscosity only; rho, cp, k stay water's", x, y, 13, kDim);
    y += 20;

    // ---- cooling plant ----
    heading(s, "Cooling plant (CDU + dry cooler)", x, y, w);
    int plantMode = s.plant ? 1 : 0;
    GuiToggleGroup({x, y, (w - 4) / 2, 26}, "Fixed supply;CDU + dry cooler", &plantMode);
    if ((plantMode == 1) != s.plant) { s.plant = plantMode == 1; app.dc->setPlantEnabled(s.plant); }
    y += 34;
    if (s.plant) {
        if (slider(s, "Supply setpoint", f1("%.1f \xC2\xB0""C", s.setpoint), x, y, w, &s.setpoint, 20.0f, 45.0f))
            app.dc->plant().setSetpoint(s.setpoint);
        if (slider(s, "Outdoor air", f1("%.0f \xC2\xB0""C", s.outdoor), x, y, w, &s.outdoor, -10.0f, 45.0f))
            app.dc->plant().setOutdoorMean(s.outdoor);
        bool ch = s.chiller;
        GuiToggle({x, y, w, 24}, ch ? "Trim chiller allowed" : "Trim chiller off (dry cooler only)", &ch);
        if (ch != s.chiller) { s.chiller = ch; app.dc->plant().setChillerEnabled(ch); }
        y += 30;
    } else {
        if (slider(s, "Supply temp", f1("%.1f \xC2\xB0""C", s.supplyT), x, y, w, &s.supplyT, 20.0f, 45.0f))
            app.dc->setFixedSupply(s.supplyT);
    }

    // ---- pump ----
    heading(s, "Pump control", x, y, w);
    int mode = s.pumpMode;
    GuiToggleGroup({x, y, (w - 8) / 3, 26}, "Fixed speed;Fixed flow;Hold peak T", &mode);
    if (mode != s.pumpMode) {
        s.pumpMode = mode;
        if (mode == 0) { s.pumpSpeed = static_cast<float>(r.pumpSpeed()); r.setPumpMode(mf::PumpMode::FixedSpeed); r.setPumpSpeed(s.pumpSpeed); }
        if (mode == 1) r.setPumpMode(mf::PumpMode::FixedFlow);
        if (mode == 2) r.setPumpMode(mf::PumpMode::PeakTemperature);
    }
    y += 34;
    if (s.pumpMode == 0) {
        if (slider(s, "Pump speed", f1("%.0f%%", 100 * s.pumpSpeed), x, y, w, &s.pumpSpeed, 0.1f, 2.0f)) r.setPumpSpeed(s.pumpSpeed);
    } else if (s.pumpMode == 1) {
        if (slider(s, "Rack flow", f1("%.1f L/m", s.targetFlow), x, y, w, &s.targetFlow, 2.0f, 24.0f)) r.setTargetFlowLpm(s.targetFlow);
    } else {
        if (slider(s, "Peak limit", f1("%.0f \xC2\xB0""C", s.targetPeak), x, y, w, &s.targetPeak, 60.0f, 95.0f)) r.setTargetPeak(s.targetPeak);
    }

    // ---- manifold ----
    heading(s, "Rack plumbing", x, y, w);
    int lay = s.layout;
    GuiToggleGroup({x, y, (w - 4) / 2, 26}, "Direct return;Reverse return", &lay);
    if (lay != s.layout) {
        s.layout = lay;
        r.setLayout(lay == 0 ? mf::ManifoldLayout::DirectReturn : mf::ManifoldLayout::ReverseReturn);
    }
    y += 34;
    bool bal = s.balance;
    GuiToggle({x, y, (w - 8) / 2, 26}, bal ? "Balancing valves: ON" : "Balancing valves: OFF", &bal);
    if (bal != s.balance) { s.balance = bal; r.setAutoBalance(bal); }
    GuiToggle({x + (w + 8) / 2, y, (w - 8) / 2, 26}, s.labels ? "Labels: ON" : "Labels: OFF", &s.labels);
    y += 34;
    drawText(s, "Drag: orbit   Wheel: zoom   Right-drag: pan", x, a.y + a.height - 44, 13, kDim);
    drawText(s, "Click a server, double-click to open it", x, a.y + a.height - 26, 13, kDim);
}

void drawReadouts(App& app, Rectangle a) {
    const ViewerState& s = app.s;
    const mf::Rack& r = *app.rack;
    DrawRectangleRounded(a, 0.03f, 6, Color{14, 16, 20, 215});
    const float x = a.x + 14, w = a.width - 28;
    float y = a.y + 10;
    const mf::Hotspot h = r.hotspot();
    const mf::LiquidLedger l = r.ledger();
    drawText(s, f1("t = %.1f s", r.time()) + (s.paused ? "  (paused)" : ""), x, y, 17, kInk);
    y += 26;
    row(s, "IT power", f1("%.0f W", r.itPower()), x, y, w);
    row(s, "  liquid (GPUs + CPUs)", f1("%.0f W", l.powerIn), x, y, w);
    row(s, "  air group", f1("%.0f W", r.itPower() - l.powerIn), x, y, w);
    row(s, "Rack flow", f1("%.2f L/min", r.totalFlowLpm()), x, y, w);
    row(s, "Pump speed / head", f1("%.0f%%", 100 * r.pumpSpeed()) + f1("  /  %.1f kPa", r.pumpHead() / 1e3), x, y, w);
    row(s, "Pump electric", f1("%.1f W", r.pumpPower()), x, y, w, Color{110, 170, 255, 255});
    row(s, "Fans (8 servers)", f1("%.1f W", r.fanPower()), x, y, w);
    row(s, "Supply / return", f1("%.1f", r.supplyT()) + f1("  /  %.2f \xC2\xB0""C", r.returnT()), x, y, w);
    row(s, "Hottest chip", f1("%.2f \xC2\xB0""C", h.t), x, y, w, heatColor(h.t, kSiliconLo, kSiliconHi));
    row(s, "  where", "S" + std::to_string(h.server + 1) + " " + chipName(h.branch) + "  cell (" +
                          std::to_string(h.i) + "," + std::to_string(h.j) + ")", x, y, w);
    row(s, "Server flow spread", f1("%.2f %%", 100 * r.flowMaldistribution()), x, y, w);
    row(s, "Energy residual", f1("%.1e W", l.residual()), x, y, w);
    row(s, "Pump / fan energy", f1("%.2f", r.pumpEnergy() / 1e3) + f1(" / %.2f kJ", r.fanEnergy() / 1e3), x, y, w);
    row(s, "Viscosity / layout", f1("%.2fx", r.viscosityScale()) + "  /  " + mf::toString(r.params().layout), x, y, w);
    if (app.dc->plantEnabled()) {
        const mf::PlantState& p = app.dc->plant().state();
        row(s, "Outdoor / facility sup-ret", f1("%.1f", p.outdoorT) + f1("  /  %.1f", p.facilitySupplyT) +
                                                 f1("-%.1f \xC2\xB0""C", p.facilityReturnT), x, y, w);
        row(s, "Dry-cooler fans", f1("%.0f%%", 100 * p.fanSpeed) + f1("  %.1f W", p.fanPower), x, y, w);
        row(s, "Trim chiller", f1("%.1f W", p.chillerPower) + (p.chillerHeat > 1 ? f1("  (COP %.1f)", p.chillerCop) : ""), x, y, w);
        row(s, "Facility pump", f1("%.1f W", p.facilityPumpPower), x, y, w);
    } else {
        row(s, "Plant", "off: supply held fixed", x, y, w, kDim);
    }
    row(s, "Cooling electricity", f1("%.1f W", app.dc->coolingPower()), x, y, w, Color{110, 170, 255, 255});
    row(s, "Partial PUE (modelled)", f1("%.4f", app.dc->partialPue()), x, y, w);
    row(s, "Water evaporated", f1("0 L/h  (tower: %.1f L/h)", app.dc->evaporativeEquivalentLph()), x, y, w);

    // Selected chip details.
    const mf::Server& sv = r.server(s.server);
    const mf::ChipGrid& c = sv.chip(s.branch);
    const double rho = r.network().fluid().props(r.supplyT()).rho;
    y += 8;
    drawText(s, "Selected: server " + std::to_string(s.server + 1) + ", " + chipName(s.branch), x, y, 16, kHead);
    y += 24;
    row(s, "Power", f1("%.0f W", c.power()) + f1("  (u = %.0f%%)", 100 * c.activity()), x, y, w);
    row(s, "Flow", f1("%.4f L/min", mf::m3sToLpm(sv.branchFlow(s.branch) / rho)), x, y, w);
    row(s, "Branch dp", f1("%.2f kPa", sv.branchDp(r.network(), s.branch) / 1e3) +
                            f1("  (valve %.1f)", sv.valveDpRef(s.branch) / 1e3), x, y, w);
    row(s, "Inlet / outlet", f1("%.1f", c.inletT()) + f1("  /  %.2f \xC2\xB0""C", c.outletT()), x, y, w);
    row(s, "Peak / mean silicon", f1("%.2f", c.peakT()) + f1("  /  %.2f \xC2\xB0""C", c.meanSiliconT()), x, y, w);
    if (c.params().mode == mf::CoolingMode::Microchannel) {
        mf::RectChannel ch(c.params().channels);
        const mf::FluidProps p = r.network().fluid().props(c.meanFluidT());
        const mf::ChannelState st = ch.state(c.volumeFlow(), p);
        row(s, "Channel velocity / Re", f1("%.2f m/s", st.velocity) + f1("  /  %.0f", st.re), x, y, w);
        row(s, "h inlet / outlet", f1("%.0f", c.heatTransferCoeff(0)) + f1(" / %.0f W/m\xC2\xB2K", c.heatTransferCoeff(c.nx() - 1)), x, y, w);
        row(s, "Channel dp (f/dev/hdr)", f1("%.1f", st.dpFriction / 1e3) + f1("/%.1f", st.dpDeveloping / 1e3) +
                                             f1("/%.1f kPa", st.dpMinor / 1e3), x, y, w);
    } else {
        row(s, "Cold plate hA", f1("%.1f W/K", c.heatTransferCoeff(0) * c.params().die.lx * c.params().die.ly), x, y, w);
    }
}

} // namespace viewer
