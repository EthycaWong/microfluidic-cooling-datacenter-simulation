#include "Charts.hpp"
#include "Units.hpp"

#include "raygui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace viewer {
namespace {

const Color kInk{226, 232, 240, 255};
const Color kDim{150, 160, 175, 255};
const Color kGrid{44, 50, 60, 255};
const Color kPanel{14, 16, 20, 225};

std::string fmt(const char* f, double a) { char b[96]; std::snprintf(b, sizeof b, f, a); return b; }

struct Series {
    const char* name;
    Color color;
    double HistorySample::*field;
    double scale = 1.0;
};

/// Visible samples, decimated to at most ~maxPts for drawing.
std::vector<const HistorySample*> visible(const ViewerState& s, double t0, size_t maxPts = 700) {
    std::vector<const HistorySample*> v;
    for (const auto& h : s.history)
        if (h.t >= t0) v.push_back(&h);
    if (v.size() <= maxPts) return v;
    std::vector<const HistorySample*> d;
    const double stride = static_cast<double>(v.size()) / maxPts;
    for (size_t i = 0; i < maxPts; ++i) d.push_back(v[static_cast<size_t>(i * stride)]);
    d.push_back(v.back());
    return d;
}

double niceStep(double range) {
    const double raw = range / 4.0;
    const double mag = std::pow(10.0, std::floor(std::log10(std::max(raw, 1e-9))));
    for (double m : {1.0, 2.0, 2.5, 5.0, 10.0})
        if (raw <= m * mag) return m * mag;
    return 10 * mag;
}

struct Frame { Rectangle plot; double t0, t1, lo, hi; };

Frame panelFrame(const ViewerState& s, Rectangle r, const char* title, double t0, double t1, double lo, double hi,
                 const char* unit) {
    DrawRectangleRounded(r, 0.03f, 6, kPanel);
    drawText(s, title, r.x + 10, r.y + 6, 15, kInk);
    Frame f{{r.x + 46, r.y + 30, r.width - 58, r.height - 58}, t0, t1, lo, hi};
    if (f.hi - f.lo < 1e-9) { f.hi += 1; f.lo -= 1; }
    const double step = niceStep(f.hi - f.lo);
    f.lo = std::floor(f.lo / step) * step;
    f.hi = std::ceil(f.hi / step) * step;
    for (double v = f.lo; v <= f.hi + 1e-9; v += step) {
        const float y = f.plot.y + f.plot.height * static_cast<float>(1 - (v - f.lo) / (f.hi - f.lo));
        DrawLine(static_cast<int>(f.plot.x), static_cast<int>(y), static_cast<int>(f.plot.x + f.plot.width), static_cast<int>(y), kGrid);
        drawText(s, std::abs(v) >= 100 || step >= 1 ? fmt("%.0f", v) : fmt("%.1f", v), r.x + 6, y - 7, 12, kDim);
    }
    drawText(s, unit, r.x + 6, r.y + r.height - 22, 12, kDim);
    const double span = t1 - t0;
    const std::string win = span >= 120 ? fmt("last %.0f min", span / 60) : fmt("last %.0f s", span);
    drawText(s, win, r.x + r.width - 12 - MeasureTextEx(s.font, win.c_str(), 12, 0.5f).x, r.y + 8, 12, kDim);
    DrawRectangleLinesEx(f.plot, 1, kGrid);
    return f;
}

Vector2 toScreen(const Frame& f, double t, double v) {
    return {f.plot.x + f.plot.width * static_cast<float>((t - f.t0) / (f.t1 - f.t0)),
            f.plot.y + f.plot.height * static_cast<float>(1 - std::clamp((v - f.lo) / (f.hi - f.lo), -0.02, 1.02))};
}

void legend(const ViewerState& s, Rectangle r, std::vector<std::pair<std::string, Color>> items) {
    float size = 12;
    auto width = [&] {
        float w = 0;
        for (auto& it : items) w += MeasureTextEx(s.font, it.first.c_str(), size, 0.5f).x + 22;
        return w;
    };
    if (width() > r.width - 50) // drop the live values first
        for (auto& it : items) { const auto p = it.first.rfind(' '); if (p != std::string::npos) it.first.resize(p); }
    while (width() > r.width - 50 && size > 9) size -= 1;
    float x = r.x + r.width - 8;
    for (auto it = items.rbegin(); it != items.rend(); ++it) {
        const float w = MeasureTextEx(s.font, it->first.c_str(), size, 0.5f).x;
        x -= w + 22;
        DrawRectangle(static_cast<int>(x), static_cast<int>(r.y + 11), 10, 4, it->second);
        drawText(s, it->first, x + 13, r.y + 5, size, kDim);
    }
}

void lineChart(const ViewerState& s, Rectangle r, const char* title, const std::vector<Series>& ser,
               const std::vector<const HistorySample*>& pts, double t0, double t1, const char* unit,
               double refLine = NAN, double fixLo = NAN, double fixHi = NAN) {
    double lo = 1e300, hi = -1e300;
    for (const auto* p : pts)
        for (const auto& sr : ser) {
            const double v = p->*sr.field * sr.scale;
            if (std::isnan(v)) continue;
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
    if (lo > hi) { lo = 0; hi = 1; }
    if (pts.empty()) { lo = 0; hi = 1; }
    if (!std::isnan(refLine)) { lo = std::min(lo, refLine); hi = std::max(hi, refLine); }
    if (!std::isnan(fixLo)) lo = std::min(lo, fixLo);
    if (!std::isnan(fixHi)) hi = std::max(hi, fixHi);
    const double pad = 0.08 * (hi - lo + 1e-9);
    const double lo2 = (!std::isnan(fixLo) && lo >= fixLo) ? fixLo : lo - pad;
    Frame f = panelFrame(s, r, title, t0, t1, lo2, hi + pad, unit);
    if (!std::isnan(refLine)) {
        const Vector2 a = toScreen(f, t0, refLine), b = toScreen(f, t1, refLine);
        DrawLineEx(a, b, 1.0f, Color{220, 70, 70, 160});
    }
    for (const auto& sr : ser) {
        Vector2 prev{};
        bool first = true;
        for (const auto* p : pts) {
            const double val = p->*sr.field * sr.scale;
            if (std::isnan(val)) { first = true; continue; }
            const Vector2 v = toScreen(f, p->t, val);
            if (!first) DrawLineEx(prev, v, 1.8f, sr.color);
            prev = v;
            first = false;
        }
    }
    std::vector<std::pair<std::string, Color>> items;
    for (const auto& sr : ser) {
        std::string label = sr.name;
        if (!pts.empty() && !std::isnan(pts.back()->*sr.field)) label += fmt(" %.1f", pts.back()->*sr.field * sr.scale);
        items.push_back({label, sr.color});
    }
    legend(s, {r.x, r.y + r.height - 22, r.width, 18}, items);
}

void stackedChart(const ViewerState& s, Rectangle r, const char* title, const std::vector<Series>& ser,
                  const std::vector<const HistorySample*>& pts, double t0, double t1, const char* unit) {
    double hi = 1.0;
    for (const auto* p : pts) {
        double sum = 0;
        for (const auto& sr : ser) sum += std::max(0.0, p->*sr.field * sr.scale);
        hi = std::max(hi, sum);
    }
    Frame f = panelFrame(s, r, title, t0, t1, 0.0, hi * 1.1, unit);
    for (size_t k = 1; k < pts.size(); ++k) {
        const auto* a = pts[k - 1];
        const auto* b = pts[k];
        double baseA = 0, baseB = 0;
        for (const auto& sr : ser) {
            const double va = std::max(0.0, a->*sr.field * sr.scale), vb = std::max(0.0, b->*sr.field * sr.scale);
            const Vector2 p0 = toScreen(f, a->t, baseA), p1 = toScreen(f, b->t, baseB);
            const Vector2 q0 = toScreen(f, a->t, baseA + va), q1 = toScreen(f, b->t, baseB + vb);
            DrawTriangle(p0, q1, q0, sr.color);   // two triangles per segment (counter-clockwise)
            DrawTriangle(p0, p1, q1, sr.color);
            baseA += va;
            baseB += vb;
        }
    }
    std::vector<std::pair<std::string, Color>> items;
    for (const auto& sr : ser) {
        std::string label = sr.name;
        if (!pts.empty()) label += fmt(" %.0f", pts.back()->*sr.field * sr.scale);
        items.push_back({label, sr.color});
    }
    legend(s, {r.x, r.y + r.height - 22, r.width, 18}, items);
}

void serverPanel(const ViewerState& s, Rectangle r, const mf::Rack& rack) {
    DrawRectangleRounded(r, 0.03f, 6, kPanel);
    drawText(s, "Per server now: hottest chip (bar), flow and activity", r.x + 10, r.y + 6, 15, kInk);
    const int n = rack.serverCount();
    const float x0 = r.x + 46, w = r.width - 58, y0 = r.y + 32, h = r.height - 76;
    const double lo = 30, hi = 95;
    for (double v : {40.0, 60.0, 85.0}) {
        const float y = y0 + h * static_cast<float>(1 - (v - lo) / (hi - lo));
        DrawLine(static_cast<int>(x0), static_cast<int>(y), static_cast<int>(x0 + w), static_cast<int>(y),
                 v == 85.0 ? Color{220, 70, 70, 160} : kGrid);
        drawText(s, fmt("%.0f", v), r.x + 10, y - 7, 12, kDim);
    }
    const double rho = rack.network().fluid().props(rack.supplyT()).rho;
    const float bw = w / n;
    for (int k = 0; k < n; ++k) {
        const mf::Server& sv = rack.server(k);
        const double t = sv.peakT();
        const float bh = h * static_cast<float>(std::clamp((t - lo) / (hi - lo), 0.0, 1.0));
        const float bx = x0 + k * bw + bw * 0.15f;
        DrawRectangle(static_cast<int>(bx), static_cast<int>(y0 + h - bh), static_cast<int>(bw * 0.7f), static_cast<int>(bh),
                      heatColor(t, kSiliconLo, kSiliconHi));
        drawText(s, fmt("%.0f", t), bx + 2, y0 + h - bh - 15, 12, kInk);
        // activity strip under the bar
        DrawRectangle(static_cast<int>(bx), static_cast<int>(y0 + h + 4), static_cast<int>(bw * 0.7f * sv.activity()), 5,
                      Color{110, 220, 140, 255});
        drawText(s, "S" + std::to_string(k + 1), bx, y0 + h + 12, 12, kDim);
        drawText(s, fmt("%.2f", mf::m3sToLpm(sv.totalFlow() / rho)), bx - 2, y0 + h + 26, 11, Color{120, 180, 255, 255});
    }
    drawText(s, "\xC2\xB0""C", r.x + 10, y0 - 4, 12, kDim);
    drawText(s, "L/min", r.x + 6, y0 + h + 26, 11, Color{120, 180, 255, 255});
}

void chipPanel(const ViewerState& s, Rectangle r, const mf::Rack& rack) {
    DrawRectangleRounded(r, 0.03f, 6, kPanel);
    const mf::ChipGrid& c = rack.server(s.server).chip(s.branch);
    const char* names[] = {"GPU 1", "GPU 2", "CPU"};
    drawText(s, "Server " + std::to_string(s.server + 1) + " " + names[s.branch] + ": silicon map + coolant along flow",
             r.x + 10, r.y + 6, 15, kInk);
    const int nx = c.nx(), ny = c.ny();
    const float side = std::min(r.height - 70, r.width * 0.5f);
    const float cx = r.x + 14, cy = r.y + 32, cell = side / std::max(nx, ny);
    double tmin = 1e9, tmax = -1e9;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) { tmin = std::min(tmin, c.siliconT(i, j)); tmax = std::max(tmax, c.siliconT(i, j)); }
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double t = c.siliconT(i, j);
            const Rectangle cr{cx + i * cell, cy + (ny - 1 - j) * cell, cell - 1, cell - 1};
            DrawRectangleRec(cr, heatColor(t, kSiliconLo, kSiliconHi));
            if (cell > 26) drawText(s, fmt("%.0f", t), cr.x + 3, cr.y + cell / 2 - 7, 11, Color{10, 10, 14, 230});
        }
    DrawRectangleLinesEx({cx + c.peakI() * cell - 1, cy + (ny - 1 - c.peakJ()) * cell - 1, cell + 1, cell + 1}, 2, WHITE);
    drawText(s, "flow \xE2\x86\x92", cx, cy + side + 4, 12, kDim);
    drawText(s, fmt("min %.1f", tmin) + fmt("  max %.1f \xC2\xB0""C", tmax), cx + 60, cy + side + 4, 12, kDim);

    // coolant temperature and heat-transfer coefficient along the channel (column-averaged)
    const Rectangle pr{cx + side + 44, cy + 4, r.x + r.width - (cx + side + 44) - 12, side - 30};
    DrawRectangleLinesEx(pr, 1, kGrid);
    double flo = 1e9, fhi = -1e9;
    std::vector<double> fl(nx), si(nx);
    for (int i = 0; i < nx; ++i) {
        double a = 0, b = 0;
        for (int j = 0; j < ny; ++j) { a += c.fluidT(i, j); b += c.siliconT(i, j); }
        fl[i] = a / ny; si[i] = b / ny;
        flo = std::min({flo, fl[i], c.inletT()}); fhi = std::max({fhi, si[i]});
    }
    auto P = [&](int i, double v) {
        return Vector2{pr.x + pr.width * (i + 0.5f) / nx,
                       pr.y + pr.height * static_cast<float>(1 - (v - flo) / std::max(fhi - flo, 1e-6))};
    };
    for (int i = 1; i < nx; ++i) {
        DrawLineEx(P(i - 1, fl[i - 1]), P(i, fl[i]), 2, Color{90, 170, 255, 255});
        DrawLineEx(P(i - 1, si[i - 1]), P(i, si[i]), 2, Color{255, 140, 80, 255});
    }
    drawText(s, fmt("%.0f", fhi), pr.x - 26, pr.y - 6, 11, kDim);
    drawText(s, fmt("%.0f", flo), pr.x - 26, pr.y + pr.height - 8, 11, kDim);
    drawText(s, "silicon (column mean)", pr.x + 4, pr.y + 4, 11, Color{255, 140, 80, 255});
    drawText(s, "coolant", pr.x + 4, pr.y + 18, 11, Color{90, 170, 255, 255});
    drawText(s, "inlet \xE2\x86\x92 outlet", pr.x + pr.width - 80, pr.y + pr.height + 4, 11, kDim);
    float ty = pr.y + pr.height + 22;
    drawText(s, fmt("power %.0f W", c.power()) + fmt("   in %.1f", c.inletT()) + fmt(" / out %.2f \xC2\xB0""C", c.outletT()),
             pr.x, ty, 12, kInk);
    ty += 16;
    drawText(s, fmt("h inlet %.0f", c.heatTransferCoeff(0)) + fmt(" / outlet %.0f W/m\xC2\xB2K", c.heatTransferCoeff(nx - 1)),
             pr.x, ty, 12, kDim);
}

void summaryPanel(const ViewerState& s, Rectangle r, const mf::DataCenter& dc, const std::vector<const HistorySample*>& pts) {
    DrawRectangleRounded(r, 0.03f, 6, kPanel);
    drawText(s, "Energy and efficiency since start", r.x + 10, r.y + 6, 15, kInk);
    const mf::Rack& rack = dc.rack();
    const double itJ = rack.itEnergy(), fanJ = rack.fanEnergy(), coolJ = dc.coolingEnergy();
    const double hours = std::max(rack.time(), 1e-9) / 3600.0;
    float y = r.y + 30, col = r.x + 12;
    const float colW = r.width * 0.31f;
    int n = 0;
    auto row = [&](const std::string& k, const std::string& v, Color c = kInk) {
        if (n == 6) { y = r.y + 30; col += colW + 10; }
        drawText(s, k, col, y, 13, kDim);
        const float w = MeasureTextEx(s.font, v.c_str(), 13, 0.5f).x;
        drawText(s, v, col + colW - w, y, 13, c);
        y += 18;
        ++n;
    };
    row("simulated time", fmt("%.1f min", rack.time() / 60));
    row("IT energy", fmt("%.3f kWh", itJ / 3.6e6));
    row("server fans (IT)", fmt("%.1f Wh", fanJ / 3600));
    row("cooling electricity", fmt("%.1f Wh", coolJ / 3600), Color{110, 170, 255, 255});
    row("  of which rack pump", fmt("%.1f Wh", rack.pumpEnergy() / 3600));
    row("  dry-cooler fans", fmt("%.1f Wh", dc.plantFanEnergy() / 3600));
    row("  chiller", fmt("%.1f Wh", dc.chillerEnergy() / 3600));
    const double pue = (itJ + fanJ) > 0 ? (itJ + fanJ + coolJ) / (itJ + fanJ) : 1.0;
    row("partial PUE (cumulative)", fmt("%.4f", pue));
    row("average IT power", fmt("%.2f kW", itJ / 3.6e6 / std::max(hours, 1e-9)));
    row("water evaporated", "0 L");
    row("evaporative tower would use", fmt("%.1f L", dc.evaporativeEquivalentLph() * hours));

    // mini pPUE trend on the right
    const Rectangle pr{r.x + r.width * 0.66f, r.y + 34, r.width * 0.31f, r.height - 70};
    DrawRectangleLinesEx(pr, 1, kGrid);
    drawText(s, "instant pPUE", pr.x, pr.y - 16 + 0, 11, kDim);
    double lo = 1e9, hi = -1e9;
    for (const auto* p : pts) { lo = std::min(lo, p->pue); hi = std::max(hi, p->pue); }
    if (pts.size() > 1 && hi > lo) {
        Vector2 prev{};
        for (size_t k = 0; k < pts.size(); ++k) {
            const Vector2 v{pr.x + pr.width * static_cast<float>(k) / (pts.size() - 1),
                            pr.y + pr.height * static_cast<float>(1 - (pts[k]->pue - lo) / (hi - lo))};
            if (k) DrawLineEx(prev, v, 1.5f, Color{200, 160, 255, 255});
            prev = v;
        }
        drawText(s, fmt("%.3f", hi), pr.x + 2, pr.y + 2, 11, kDim);
        drawText(s, fmt("%.3f", lo), pr.x + 2, pr.y + pr.height - 14, 11, kDim);
    }
}

} // namespace

void recordSample(ViewerState& s, const mf::DataCenter& dc) {
    const mf::Rack& r = dc.rack();
    const mf::PlantState& p = dc.plant().state();
    const bool plant = dc.plantEnabled();
    HistorySample h{};
    h.t = r.time();
    h.peak = r.hotspot().t;
    h.pumpW = r.pumpPower();
    h.flowLpm = r.totalFlowLpm();
    h.returnT = r.returnT();
    h.u = r.server(0).activity();
    h.supplyT = r.supplyT();
    h.facSupplyT = plant ? p.facilitySupplyT : NAN;
    h.facReturnT = plant ? p.facilityReturnT : NAN;
    h.outdoorT = plant ? p.outdoorT : NAN;
    h.itKW = r.itPower() / 1000.0;
    h.liquidKW = r.ledger().powerIn / 1000.0;
    h.serverFanW = r.fanPower();
    h.facPumpW = plant ? p.facilityPumpPower : 0.0;
    h.dcFanW = plant ? p.fanPower : 0.0;
    h.chillerW = plant ? p.chillerPower : 0.0;
    h.coolingW = dc.coolingPower();
    h.headKPa = r.pumpHead() / 1000.0;
    h.pumpPct = 100.0 * r.pumpSpeed();
    h.pue = dc.partialPue();
    double meanU = 0;
    for (int k = 0; k < r.serverCount(); ++k) meanU += r.server(k).activity();
    h.meanU = meanU / r.serverCount();
    s.history.push_back(h);
    while (!s.history.empty() && s.history.front().t < h.t - 3700.0) s.history.pop_front();
}

void drawDashboard(ViewerState& s, const mf::DataCenter& dc, Rectangle a) {
    const double spans[] = {60.0, 600.0, 3600.0};
    const double t1 = s.history.empty() ? 0.0 : s.history.back().t;
    const double t0 = t1 - spans[s.chartWindow];
    drawText(s, "Live charts", a.x, a.y, 20, kInk);
    GuiToggleGroup({a.x + a.width - 3 * 82, a.y - 2, 78, 24}, "1 min;10 min;60 min", &s.chartWindow);
    const auto pts = visible(s, t0);
    const bool plant = dc.plantEnabled();

    const float gx = 12, gy = 10, top = a.y + 30;
    const float cw = (a.width - gx) / 2, ch = (a.height - 30 - 3 * gy) / 4;
    auto cell = [&](int col, int row) { return Rectangle{a.x + col * (cw + gx), top + row * (ch + gy), cw, ch}; };

    std::vector<Series> temps = {{"chip max", Color{255, 90, 80, 255}, &HistorySample::peak},
                                 {"return", Color{255, 170, 60, 255}, &HistorySample::returnT},
                                 {"supply", Color{90, 170, 255, 255}, &HistorySample::supplyT}};
    if (plant) {
        temps.push_back({"fac. sup", Color{80, 220, 220, 255}, &HistorySample::facSupplyT});
        temps.push_back({"outdoor", Color{170, 170, 170, 255}, &HistorySample::outdoorT});
    }
    lineChart(s, cell(0, 0), "Temperatures", temps, pts, t0, t1, "\xC2\xB0""C", 85.0);

    lineChart(s, cell(1, 0), "IT power and activity",
              {{"IT", Color{230, 230, 240, 255}, &HistorySample::itKW},
               {"to liquid", Color{255, 140, 80, 255}, &HistorySample::liquidKW},
               {"activity x10", Color{110, 220, 140, 255}, &HistorySample::meanU, 10.0}},
              pts, t0, t1, "kW", NAN, 0.0);

    stackedChart(s, cell(0, 1), "Cooling electricity by component (stacked)",
                 {{"rack pump", Color{90, 150, 255, 255}, &HistorySample::pumpW},
                  {"fac. pump", Color{80, 210, 210, 255}, &HistorySample::facPumpW},
                  {"dry-cooler fans", Color{120, 220, 120, 255}, &HistorySample::dcFanW},
                  {"chiller", Color{200, 120, 255, 255}, &HistorySample::chillerW}},
                 pts, t0, t1, "W");

    lineChart(s, cell(1, 1), "Rack hydraulics",
              {{"flow L/min", Color{90, 170, 255, 255}, &HistorySample::flowLpm},
               {"head kPa/10", Color{255, 200, 90, 255}, &HistorySample::headKPa, 0.1},
               {"pump %/10", Color{200, 160, 255, 255}, &HistorySample::pumpPct, 0.1}},
              pts, t0, t1, "", NAN, 0.0);

    serverPanel(s, cell(0, 2), dc.rack());
    chipPanel(s, cell(1, 2), dc.rack());
    Rectangle bottom = cell(0, 3);
    bottom.width = a.width;
    summaryPanel(s, bottom, dc, pts);
}

} // namespace viewer
