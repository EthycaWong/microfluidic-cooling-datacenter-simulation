#include "CfdView.hpp"

#include "raygui.h"
#include "raymath.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>

namespace viewer {
using namespace mf::cfd;

namespace {
const Color kInk{226, 232, 240, 255};
const Color kDim{150, 160, 175, 255};
const Color kHead{120, 190, 255, 255};
const char* kCaseLabels = "Lid-driven cavity;Heated cavity;Rayleigh-Benard;Mixed convection;Microchannel slice;Channel flow";
std::mt19937 rng(42);

std::string f1(const char* f, double a) { char b[96]; std::snprintf(b, sizeof b, f, a); return b; }

Side& sideRef(Params& p, int k) { Side* s[4] = {&p.left, &p.right, &p.bottom, &p.top}; return *s[k]; }

bool isRaCase(const std::string& n) { return n == "heated-cavity" || n == "rayleigh-benard"; }

void spawn(CfdApp& c, CfdApp::Particle& p, bool atInlet) {
    const Params& q = c.solver->params();
    std::uniform_real_distribution<double> ux(0.0, 1.0);
    if (atInlet && q.left.u > 0 && q.left.velocity == VelocityBC::Wall) p.x = 0.01 * q.lx * ux(rng);
    else p.x = q.lx * ux(rng);
    p.y = q.ly * (0.02 + 0.96 * ux(rng));
    p.age = 0.0f;
}
} // namespace

void cfdLoadCase(CfdApp& c, int index, double value) {
    const auto names = caseNames();
    c.caseIndex = std::clamp(index, 0, static_cast<int>(names.size()) - 1);
    const int n = names[c.caseIndex] == "microchannel" ? 10 : names[c.caseIndex] == "rayleigh-benard" ? 32 : 48;
    c.kase = makeCase(names[c.caseIndex], value, n);
    c.solver = std::make_unique<NavierStokes2D>(c.kase.params);
    if (c.value <= 0 || value <= 0) {
        const std::string& nm = names[c.caseIndex];
        c.value = static_cast<float>(nm == "microchannel" ? c.kase.velocityScale
                                     : isRaCase(nm) ? std::log10(-c.kase.params.gy / c.kase.params.nu)
                                                    : (nm == "poiseuille" ? 2.0 / c.kase.params.nu : 1.0 / c.kase.params.nu));
    }
    c.particles.assign(500, {});
    for (auto& p : c.particles) spawn(c, p, false);
    if (c.kase.name == "microchannel") c.colorMode = 0;
}

void cfdUpdate(CfdApp& c) {
    if (!c.solver || c.paused) { c.simPerFrame = 0.0; return; }
    const auto t0 = std::chrono::steady_clock::now();
    double sim = 0.0;
    do {
        sim += c.solver->step();
    } while (std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() < c.budgetMs);
    c.simPerFrame = sim;
    // Tracers ride the computed velocity field (midpoint rule); show the flow pattern.
    const Params& q = c.solver->params();
    for (auto& p : c.particles) {
        double u1, v1, u2, v2;
        c.solver->velocityAt(p.x, p.y, u1, v1);
        c.solver->velocityAt(p.x + 0.5 * sim * u1, p.y + 0.5 * sim * v1, u2, v2);
        p.x += sim * u2;
        p.y += sim * v2;
        p.age += 1.0f;
        if (p.x < 0 || p.x > q.lx || p.y < 0 || p.y > q.ly || p.age > 600) spawn(c, p, true);
    }
}

void drawCfdField(CfdApp& c, const ViewerState& s, Rectangle a) {
    if (!c.solver) return;
    NavierStokes2D& sv = *c.solver;
    const Params& q = sv.params();
    drawText(s, "2-D Navier-Stokes + energy (Boussinesq)", a.x, a.y, 18, kInk);
    drawText(s, c.kase.description, a.x, a.y + 22, 14, kInk);
    drawText(s, c.kase.reynoldsNote + f1("   grid %.0f", sv.nx()) + f1(" x %.0f", sv.ny()) + "   staggered MAC grid, projection method",
             a.x, a.y + 40, 13, kDim);
    // Fit the domain, stretching thin channels vertically for visibility.
    const float top = a.y + 64, availW = a.width - 40, availH = a.height - 150;
    double stretch = 1.0;
    const double aspect = q.lx / q.ly;
    if (aspect > 4.0) stretch = aspect / 4.0;
    float scale = static_cast<float>(std::min(availW / q.lx, availH / (q.ly * stretch)));
    const float W = static_cast<float>(q.lx * scale), H = static_cast<float>(q.ly * stretch * scale);
    const float ox = a.x + 20 + (availW - W) / 2, oy = top + (availH - H) / 2;
    auto toScreen = [&](double x, double y) { return Vector2{ox + static_cast<float>(x * scale), oy + H - static_cast<float>(y * stretch * scale)}; };

    // colour range
    double lo = c.kase.tCold, hi = c.kase.tHot;
    const int nx = sv.nx(), ny = sv.ny();
    std::vector<double> val(static_cast<size_t>(nx) * ny);
    for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j) {
            double v = 0;
            switch (c.colorMode) {
            case 0: v = sv.tc(i, j); break;
            case 1: v = std::hypot(sv.uc(i, j), sv.vc(i, j)); break;
            case 2: v = sv.pc(i, j); break;
            default: {
                const int ip = std::min(i + 1, nx), im = std::max(i - 1, 1), jp = std::min(j + 1, ny), jm = std::max(j - 1, 1);
                v = (sv.vc(ip, j) - sv.vc(im, j)) / ((ip - im) * sv.dx()) - (sv.uc(i, jp) - sv.uc(i, jm)) / ((jp - jm) * sv.dy());
            }
            }
            val[(j - 1) * nx + (i - 1)] = v;
        }
    if (c.colorMode != 0) {
        lo = 1e300; hi = -1e300;
        for (double v : val) { lo = std::min(lo, v); hi = std::max(hi, v); }
        if (c.colorMode == 3) { const double m = std::max(std::abs(lo), std::abs(hi)); lo = -m; hi = m; }
        if (hi - lo < 1e-12) hi = lo + 1;
    }
    const float cw = W / nx, ch = H / ny;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            DrawRectangle(static_cast<int>(ox + i * cw), static_cast<int>(oy + H - (j + 1) * ch), static_cast<int>(cw + 1.5f),
                          static_cast<int>(ch + 1.5f), heatColor(val[j * nx + i], lo, hi));

    // velocity arrows
    if (c.arrows) {
        const int stepI = std::max(1, nx / 24), stepJ = std::max(1, ny / (aspect > 4 ? 8 : 24));
        const double vmax = std::max(sv.maxSpeed(), 1e-12);
        const float len = std::min(cw * stepI, ch * stepJ) * 0.9f;
        for (int i = stepI / 2 + 1; i <= nx; i += stepI)
            for (int j = stepJ / 2 + 1; j <= ny; j += stepJ) {
                const double u = sv.uc(i, j), v = sv.vc(i, j);
                const double m = std::hypot(u, v);
                if (m < 0.02 * vmax) continue;
                const Vector2 p0 = toScreen((i - 0.5) * sv.dx(), (j - 0.5) * sv.dy());
                const float sx = static_cast<float>(u / vmax), sy = static_cast<float>(-v * stretch / vmax);
                const Vector2 p1{p0.x + sx * len, p0.y + sy * len};
                DrawLineEx(p0, p1, 1.4f, Color{255, 255, 255, 200});
                const Vector2 d = Vector2Normalize(Vector2Subtract(p1, p0));
                const Vector2 nrm{-d.y, d.x};
                DrawTriangle(p1, Vector2Add(Vector2Subtract(p1, Vector2Scale(d, 5)), Vector2Scale(nrm, 2.5f)),
                             Vector2Subtract(Vector2Subtract(p1, Vector2Scale(d, 5)), Vector2Scale(nrm, 2.5f)), Color{255, 255, 255, 220});
            }
    }
    if (c.tracers)
        for (const auto& p : c.particles) DrawCircleV(toScreen(p.x, p.y), 1.6f, Color{15, 15, 20, 200});

    // walls: colour by thermal condition, arrow for moving walls, dashed for outflow
    auto wallColor = [&](const Side& sd) {
        if (sd.velocity == VelocityBC::Outflow) return Color{200, 200, 210, 255};
        if (sd.thermal == ThermalBC::FixedTemperature)
            return sd.T >= 0.5 * (c.kase.tCold + c.kase.tHot) ? Color{240, 70, 60, 255} : Color{70, 140, 255, 255};
        return Color{150, 150, 160, 255};
    };
    const Vector2 c00 = toScreen(0, 0), c10 = toScreen(q.lx, 0), c01 = toScreen(0, q.ly), c11 = toScreen(q.lx, q.ly);
    const Vector2 seg[4][2] = {{c00, c01}, {c10, c11}, {c00, c10}, {c01, c11}};
    for (int k = 0; k < 4; ++k) {
        const Side& sd = sideRef(const_cast<Params&>(q), k);
        const float th = k == c.wall ? 7.0f : 5.0f;
        if (sd.velocity == VelocityBC::Outflow) {
            for (int d = 0; d < 20; d += 2)
                DrawLineEx(Vector2Lerp(seg[k][0], seg[k][1], d / 20.0f), Vector2Lerp(seg[k][0], seg[k][1], (d + 1) / 20.0f), th, wallColor(sd));
        } else {
            DrawLineEx(seg[k][0], seg[k][1], th, wallColor(sd));
        }
        if (k == c.wall) DrawLineEx(seg[k][0], seg[k][1], 1.5f, WHITE);
        if (sd.velocity == VelocityBC::Wall && (std::abs(sd.u) > 1e-12 || std::abs(sd.v) > 1e-12)) {
            const Vector2 mid = Vector2Lerp(seg[k][0], seg[k][1], 0.5f);
            const double m = std::hypot(sd.u, sd.v);
            const Vector2 dir{static_cast<float>(sd.u / m), static_cast<float>(-sd.v / m)};
            const Vector2 off = k == 0 ? Vector2{-18, 0} : k == 1 ? Vector2{18, 0} : k == 2 ? Vector2{0, 18} : Vector2{0, -18};
            const Vector2 a0 = Vector2Add(mid, off), a1 = Vector2Add(a0, Vector2Scale(dir, 40));
            DrawLineEx(a0, a1, 3, Color{255, 220, 90, 255});
            DrawCircleV(a1, 5, Color{255, 220, 90, 255});
            drawText(s, f1("u=%.3g", sd.u) + f1(" v=%.3g", sd.v), a0.x + 6, a0.y + (k == 2 ? 4 : -18), 13, Color{255, 220, 90, 255});
        }
    }
    if (stretch > 1.0) drawText(s, f1("vertical scale stretched x%.0f", stretch), ox, oy + H + 8, 13, kDim);
    // colour bar
    const char* names[] = {"temperature", "speed (m/s)", "kinematic pressure p/rho", "vorticity (1/s)"};
    const float bx = a.x + 20, by = a.y + a.height - 70, bw = a.width - 60;
    for (int k = 0; k < 64; ++k) DrawRectangleRec({bx + bw * k / 64.0f, by, bw / 64.0f + 1, 12}, heatColor(lo + (hi - lo) * k / 63.0, lo, hi));
    drawText(s, names[c.colorMode], bx, by - 18, 14, kInk);
    drawText(s, f1("%.4g", lo), bx, by + 14, 13, kDim);
    drawText(s, f1("%.4g", hi), bx + bw - 60, by + 14, 13, kDim);
    drawText(s, "walls: red = fixed hot, blue = fixed cold, grey = insulated, dashed = outflow, yellow arrow = moving wall / inflow",
             bx, by + 32, 13, kDim);
}

float drawCfdPanel(CfdApp& c, ViewerState& s, float x, float y, float w) {
    auto heading = [&](const char* t) {
        y += 6;
        drawText(s, t, x, y, 16, kHead);
        DrawLine(static_cast<int>(x), static_cast<int>(y + 20), static_cast<int>(x + w), static_cast<int>(y + 20), Color{60, 70, 84, 255});
        y += 26;
    };
    auto slider = [&](const char* name, const std::string& value, float* v, float lo, float hi) {
        drawText(s, name, x, y + 3, 15, kInk);
        const float before = *v;
        GuiSlider({x + 100, y, w - 100 - 70, 20}, nullptr, nullptr, v, lo, hi);
        drawText(s, value, x + w - 64, y + 3, 14, kDim);
        y += 27;
        return *v != before;
    };
    Params& q = c.solver->params();
    heading("Case");
    int idx = c.caseIndex;
    GuiComboBox({x, y, w, 26}, kCaseLabels, &idx);
    if (idx != c.caseIndex) { c.value = -1; cfdLoadCase(c, idx, -1); }
    y += 34;
    const std::string nm = c.kase.name;
    bool changed = false;
    if (nm == "microchannel") changed = slider("Inlet speed", f1("%.2f m/s", c.value), &c.value, 0.1f, 3.0f);
    else if (isRaCase(nm)) changed = slider("log10 Ra", f1("1e%.1f", c.value), &c.value, 3.0f, 6.0f);
    else changed = slider("Re", f1("%.0f", c.value), &c.value, 10.0f, nm == "poiseuille" ? 400.0f : 1000.0f);
    if (GuiButton({x, y, (w - 8) / 2, 26}, "Apply + restart")) {
        const double v = isRaCase(nm) ? std::pow(10.0, c.value) : c.value;
        const int keep = c.caseIndex;
        cfdLoadCase(c, keep, v);
    }
    (void)changed;
    if (GuiButton({x + (w + 8) / 2, y, (w - 8) / 2, 26}, "Reset field")) c.solver->reset();
    y += 34;

    heading("Boundary conditions (live)");
    GuiToggleGroup({x, y, (w - 12) / 4, 24}, "Left;Right;Bottom;Top", &c.wall);
    y += 32;
    Side sd = sideRef(q, c.wall);
    const float vs = static_cast<float>(c.kase.velocityScale);
    float u = static_cast<float>(sd.u), v = static_cast<float>(sd.v);
    bool sc = false;
    sc |= slider("Wall u", f1("%.3g", u), &u, -2 * vs, 2 * vs);
    sc |= slider("Wall v", f1("%.3g", v), &v, -2 * vs, 2 * vs);
    if (GuiButton({x, y, (w - 8) / 2, 22}, "u = v = 0 (no-slip)")) { u = 0; v = 0; sc = true; }
    bool outflow = sd.velocity == VelocityBC::Outflow;
    GuiToggle({x + (w + 8) / 2, y, (w - 8) / 2, 22}, outflow ? "Outflow: ON" : "Outflow: off", &outflow);
    y += 30;
    int th = sd.thermal == ThermalBC::HeatFlux ? 2 : (sd.T >= 0.5 * (c.kase.tCold + c.kase.tHot) ? 0 : 1);
    const int thBefore = th;
    GuiToggleGroup({x, y, (w - 8) / 3, 24}, "Hot (fixed);Cold (fixed);Insulated", &th);
    y += 30;
    float tv = static_cast<float>(sd.T);
    bool tChanged = false;
    if (sd.thermal == ThermalBC::FixedTemperature)
        tChanged = slider("Wall T", f1("%.2f", tv), &tv, static_cast<float>(c.kase.tCold - 0.5 * (c.kase.tHot - c.kase.tCold)),
                          static_cast<float>(c.kase.tHot + 0.5 * (c.kase.tHot - c.kase.tCold)));
    if (sc || tChanged || th != thBefore || outflow != (sd.velocity == VelocityBC::Outflow)) {
        sd.u = u;
        sd.v = v;
        sd.velocity = outflow ? VelocityBC::Outflow : VelocityBC::Wall;
        if (th != thBefore) {
            if (th == 2) { sd.thermal = ThermalBC::HeatFlux; sd.q = 0.0; }
            else { sd.thermal = ThermalBC::FixedTemperature; tv = static_cast<float>(th == 0 ? c.kase.tHot : c.kase.tCold); }
        }
        if (sd.thermal == ThermalBC::FixedTemperature) sd.T = tv;
        c.solver->setSide(c.wall, sd);
    }
    drawText(s, "Velocity BCs apply to both components; normal", x, y, 13, kDim);
    y += 16;
    drawText(s, "values add inflow. Changes act on the running flow.", x, y, 13, kDim);
    y += 22;

    heading("Display");
    GuiToggleGroup({x, y, (w - 12) / 4, 24}, "Temp;Speed;Press.;Vort.", &c.colorMode);
    y += 30;
    GuiToggle({x, y, (w - 8) / 3, 22}, c.arrows ? "Arrows ON" : "Arrows off", &c.arrows);
    GuiToggle({x + (w - 8) / 3 + 4, y, (w - 8) / 3, 22}, c.tracers ? "Tracers ON" : "Tracers off", &c.tracers);
    GuiToggle({x + 2 * (w - 8) / 3 + 8, y, (w - 8) / 3, 22}, c.paused ? "Paused" : "Running", &c.paused);
    y += 30;
    slider("CPU / frame", f1("%.0f ms", c.budgetMs), &c.budgetMs, 2.0f, 30.0f);
    return y;
}

void drawCfdReadouts(CfdApp& c, const ViewerState& s, Rectangle a) {
    DrawRectangleRounded(a, 0.03f, 6, Color{14, 16, 20, 215});
    const float x = a.x + 14, w = a.width - 28;
    float y = a.y + 10;
    auto row = [&](const std::string& k, const std::string& v, Color col = kInk) {
        drawText(s, k, x, y, 15, kDim);
        const float tw = MeasureTextEx(s.font, v.c_str(), 15, 0.5f).x;
        drawText(s, v, x + w - tw, y, 15, col);
        y += 20;
    };
    NavierStokes2D& sv = *c.solver;
    const Params& q = sv.params();
    drawText(s, "2-D Navier-Stokes solver", x, y, 17, kInk);
    y += 26;
    row("case", c.kase.reynoldsNote);
    row("simulated time", f1("%.4g s", sv.time()));
    row("steps", f1("%.0f", static_cast<double>(sv.steps())));
    row("time step", f1("%.3g s", sv.lastDt()));
    row("sim time / frame", f1("%.3g s", c.simPerFrame));
    row("max speed", f1("%.4g m/s", sv.maxSpeed()));
    row("max |div u| (mass)", f1("%.1e 1/s", sv.maxDivergence()), sv.maxDivergence() < 1e-3 ? Color{120, 230, 150, 255} : Color{255, 140, 120, 255});
    row("kinetic energy", f1("%.4g", sv.kineticEnergy()));
    row("pressure CG iterations", f1("%.0f", sv.lastPressureIterations()));
    y += 8;
    drawText(s, "Heat into the fluid through each wall", x, y, 16, kHead);
    y += 24;
    const char* nms[4] = {"left", "right", "bottom", "top"};
    const double dT = std::max(c.kase.tHot - c.kase.tCold, 1e-12);
    for (int k = 0; k < 4; ++k) {
        const Side& sd = k == 0 ? q.left : k == 1 ? q.right : k == 2 ? q.bottom : q.top;
        const double L = (k < 2) ? q.ly : q.lx;
        const double qf = sv.wallHeatFlux(k);
        std::string v = f1("%.4g W/m\xC2\xB2", qf);
        if (sd.thermal == ThermalBC::FixedTemperature && c.kase.name != "microchannel") v += f1("  (Nu %.3f)", qf * (k < 2 ? q.lx : q.ly) / (q.k * dT));
        (void)L;
        row(nms[k], v);
    }
    y += 8;
    drawText(s, "Equations solved", x, y, 16, kHead);
    y += 24;
    const char* eq[] = {"du/dt + (u.grad)u = -grad p/rho + nu lap u", "                  - beta (T - Tref) g",
                        "div u = 0", "dT/dt + u.grad T = alpha lap T"};
    for (const char* e : eq) { drawText(s, e, x, y, 14, kInk); y += 18; }
    y += 6;
    drawText(s, f1("nu %.3g m\xC2\xB2/s", q.nu) + f1("   alpha %.3g", q.alpha) + f1("   beta g %.3g", q.beta * std::hypot(q.gx, q.gy)),
             x, y, 13, kDim);
}

} // namespace viewer
