// mfsim: headless experiments for the microfluidic-cooled chip / server / rack model.
//
//   mfsim info            [--config F]                      hand-calculation diagnostics for one chip
//   mfsim steady          [--config F] [--level L]          steady state summary
//   mfsim run             [--config F] [--level L] [--duration s] [--dt s] [--every s] [--out f.csv]
//   mfsim resolution      [--config F]                      grid-resolution sensitivity (chip)
//   mfsim sweep-viscosity [--config F] [--level L] [--scales 1,2,4] [--out f.csv]
//   mfsim uneven-load     [--config F]                      half the rack loaded, direct vs reverse return
//   mfsim climate-sweep   [--config F] [--temps ...]        plant: outdoor temperature vs cooling energy
//
// Common: --viscosity X (viscosity multiplier), --level chip|server|rack (overrides the config).
#include "Units.hpp"
#include "config/Config.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace mf;

namespace {

struct Args {
    std::string cmd;
    std::map<std::string, std::string> kv;
    std::string get(const std::string& k, const std::string& def = "") const {
        auto it = kv.find(k);
        return it == kv.end() ? def : it->second;
    }
    double num(const std::string& k, double def) const {
        auto it = kv.find(k);
        return it == kv.end() ? def : std::stod(it->second);
    }
    bool has(const std::string& k) const { return kv.count(k) > 0; }
};

Args parse(int argc, char** argv) {
    Args a;
    if (argc > 1) a.cmd = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string k = argv[i];
        if (k.rfind("--", 0) != 0) continue;
        k = k.substr(2);
        if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) a.kv[k] = argv[++i];
        else a.kv[k] = "1";
    }
    return a;
}

SimConfig configFrom(const Args& a) {
    SimConfig c = a.has("config") ? loadConfig(a.get("config")) : defaultConfig();
    if (a.has("level")) c.level = a.get("level");
    if (a.has("viscosity")) c.viscosityScale = a.num("viscosity", 1.0);
    if (a.has("outdoor")) c.plant.weather.meanC = a.num("outdoor", 25.0);
    if (a.has("chiller")) c.plant.chiller.enabled = a.get("chiller") != "off";
    return c;
}

std::vector<double> parseList(const std::string& s) {
    std::vector<double> v;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) v.push_back(std::stod(item));
    return v;
}

std::shared_ptr<ScaledFluid> makeFluid(double scale) {
    return std::make_shared<ScaledFluid>(std::make_shared<WaterIAPWS>(), scale);
}

void applyWorkloads(Rack& r, const SimConfig& c) {
    r.setAllWorkloads(makeWorkload(c.workload));
    for (size_t k = 0; k < c.serverWorkloads.size() && static_cast<int>(k) < r.serverCount(); ++k)
        r.setWorkload(static_cast<int>(k), makeWorkload(c.serverWorkloads[k]));
}

std::ofstream openCsv(const std::string& path) {
    std::filesystem::path p(path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    return f;
}

const char* chipName(int b) { return b == 0 ? "GPU1" : b == 1 ? "GPU2" : "CPU"; }

// ---------------------------------------------------------------------------------------------------
int cmdInfo(const Args& a) {
    SimConfig c = configFrom(a);
    auto fluid = makeFluid(c.viscosityScale);
    const ChipParams& g = c.server.gpu;
    const FluidProps p = fluid->props(c.chipInletC);
    RectChannel ch(g.channels);
    const ChannelState s = ch.state(lpmToM3s(c.chipFlowLpm), p);
    std::printf("== One microchannel GPU: geometry and hand-check quantities ==\n");
    std::printf("Channels: %d x %.0f um wide x %.0f um deep x %.1f mm long, ribs %.0f um\n", g.channels.count,
                g.channels.width * 1e6, g.channels.depth * 1e6, g.channels.length * 1e3, g.channels.rib * 1e6);
    std::printf("Width check: %d x (channel + rib) = %.2f mm vs die %.2f mm\n", g.channels.count,
                g.channels.count * (g.channels.width + g.channels.rib) * 1e3, g.die.ly * 1e3);
    std::printf("Floor under channels: %.0f um\n", (g.die.thickness - g.channels.depth) * 1e6);
    std::printf("Coolant at %.1f C: rho %.2f kg/m3, mu %.4e Pa s, cp %.1f J/kgK, k %.4f W/mK, Pr %.3f\n",
                c.chipInletC, p.rho, p.mu, p.cp, p.k, p.pr());
    std::printf("Flow %.3f L/min: velocity %.3f m/s, Dh %.1f um, Re %.1f (laminar < ~2300)\n", c.chipFlowLpm,
                s.velocity, g.channels.hydraulicDiameter() * 1e6, s.re);
    std::printf("Entrance lengths: hydrodynamic %.2f mm, thermal %.2f mm (channel %.1f mm)\n",
                s.entranceLengthHydro * 1e3, s.entranceLengthThermal * 1e3, g.channels.length * 1e3);
    std::printf("Pressure drop: friction %.0f Pa + developing %.0f Pa + headers %.0f Pa = %.0f Pa\n", s.dpFriction,
                s.dpDeveloping, s.dpMinor, s.dpTotal);
    std::printf("Hydraulic power through the channels: %.3f W\n", s.dpTotal * lpmToM3s(c.chipFlowLpm));
    const double mdot = p.rho * lpmToM3s(c.chipFlowLpm);
    std::printf("Coolant rise at full power %.0f W: %.2f K\n", g.power.fullW, g.power.fullW / (mdot * p.cp));
    std::printf("h: fully developed (3 walls) %.0f W/m2K\n",
                RectChannel::nu3WallFullyDeveloped(g.channels.width / g.channels.depth) * p.k / g.channels.hydraulicDiameter());
    std::printf("fRe(Fanning) %.3f, K(inf) %.3f, Nu_fd(3 walls) %.3f\n", RectChannel::fanningFRe(g.channels.aspect()),
                RectChannel::hagenbachK(g.channels.aspect()), RectChannel::nu3WallFullyDeveloped(g.channels.width / g.channels.depth));
    return 0;
}

// ---------------------------------------------------------------------------------------------------
void printRackSummary(const Rack& r) {
    const Hotspot h = r.hotspot();
    const LiquidLedger l = r.ledger();
    std::printf("Pump: %s, speed %.3f, flow %.3f L/min, head %.2f kPa, eta %.2f, electric %.2f W\n",
                toString(r.params().pumpMode), r.pumpSpeed(), r.totalFlowLpm(), r.pumpHead() / 1e3,
                r.pumpEfficiency(), r.pumpPower());
    std::printf("Layout %s, viscosity x%.2f, supply %.2f C, mixed return %.3f C\n", toString(r.params().layout),
                r.viscosityScale(), r.supplyT(), r.returnT());
    std::printf("IT %.1f W (liquid %.1f W, air group %.1f W, of which %.1f W is reaching room air now)\n", r.itPower(),
                l.powerIn, r.itPower() - l.powerIn, r.airHeat());
    std::printf("Fans %.2f W, liquid energy residual %.2e W\n", r.fanPower(), l.residual());
    std::printf("Hottest silicon: %.2f C in server %d %s cell (%d,%d); server flow spread %.2f%%\n", h.t, h.server + 1,
                chipName(h.branch), h.i, h.j, 100.0 * r.flowMaldistribution());
    std::printf("  srv  U-slot  activity  flow L/min  GPU1/GPU2/CPU L/min   dp kPa  return C  peak C\n");
    const double rho = r.network().fluid().props(r.supplyT()).rho;
    for (int k = r.serverCount() - 1; k >= 0; --k) {
        const Server& s = r.server(k);
        const int u0 = r.params().firstServerU + k * r.params().serverU;
        std::printf("  %3d  U%02d-%02d  %8.2f  %10.4f  %5.3f/%5.3f/%5.3f   %6.2f  %8.3f  %6.2f\n", k + 1, u0,
                    u0 + r.params().serverU - 1, s.activity(), m3sToLpm(s.totalFlow() / rho),
                    m3sToLpm(s.branchFlow(0) / rho), m3sToLpm(s.branchFlow(1) / rho), m3sToLpm(s.branchFlow(2) / rho),
                    s.serverDp(r.network()) / 1e3, s.returnT(), s.peakT());
    }
}

std::unique_ptr<DataCenter> makeDataCenter(const SimConfig& c, double viscosity) {
    auto dc = std::make_unique<DataCenter>(c.rack, c.server, c.plant, std::make_shared<WaterIAPWS>(), true);
    dc->rack().setViscosityScale(viscosity);
    applyWorkloads(dc->rack(), c);
    return dc;
}

void printPlantSummary(const DataCenter& dc) {
    const PlantState& p = dc.plant().state();
    std::printf("Plant: outdoor %.1f C, supply setpoint %.1f C -> equipment supply %.2f C\n", p.outdoorT,
                dc.plant().params().supplySetpointC, p.equipSupplyT);
    std::printf("  CDU: UA %.0f W/K, effectiveness %.3f, equipment %.2f -> %.2f C, facility %.2f -> %.2f C, %.0f W\n",
                p.cduUA, p.cduEffectiveness, p.cduHotInT, p.cduHotOutT, p.facilitySupplyT, p.facilityReturnT, p.cduHeat);
    std::printf("  Dry cooler: fans %.0f%% = %.1f W, UA %.0f W/K, effectiveness %.3f, rejects %.0f W, outlet %.2f C\n",
                100 * p.fanSpeed, p.fanPower, p.dryCoolerUA, p.dryCoolerEffectiveness, p.dryCoolerHeat, p.dryCoolerOutT);
    std::printf("  Trim chiller: %s%.0f W removed, COP %.2f, %.1f W electric\n",
                dc.plant().params().chiller.enabled ? "" : "(disabled) ", p.chillerHeat, p.chillerCop, p.chillerPower);
    std::printf("  Electricity: equipment pump %.1f W + facility pump %.1f W + dry-cooler fans %.1f W + chiller %.1f W = %.1f W"
                " (server fans %.1f W counted as IT)\n",
                dc.rack().pumpPower(), p.facilityPumpPower, p.fanPower, p.chillerPower, dc.coolingPower(), dc.serverFanPower());
    std::printf("  Partial PUE of the modelled boundary %.4f; heat to outdoor air %.0f W; water evaporated 0 L/h "
                "(an evaporative tower would use ~%.1f L/h); plant residual %.1e W\n",
                dc.partialPue(), dc.heatToOutdoor(), dc.evaporativeEquivalentLph(), p.residual);
}

int cmdSteady(const Args& a) {
    SimConfig c = configFrom(a);
    auto fluid = makeFluid(c.viscosityScale);
    if (c.level == "chip") {
        ChipGrid chip(c.server.gpu, fluid);
        chip.setMassFlow(fluid->props(c.chipInletC).rho * lpmToM3s(c.chipFlowLpm));
        chip.setInletT(c.chipInletC);
        chip.setActivity(makeWorkload(c.workload)->u(0.0));
        chip.solveSteady();
        std::printf("Chip steady: power %.1f W, flow %.3f L/min, inlet %.2f C, outlet %.3f C\n", chip.power(),
                    c.chipFlowLpm, c.chipInletC, chip.outletT());
        std::printf("Peak silicon %.2f C at cell (%d,%d) (i along flow), mean silicon %.2f C, residual %.2e W\n",
                    chip.peakT(), chip.peakI(), chip.peakJ(), chip.meanSiliconT(), chip.ledger().residual);
        std::printf("Silicon temperature map (rows = across flow j, columns = along flow i, inlet on the left):\n");
        for (int j = chip.ny() - 1; j >= 0; --j) {
            for (int i = 0; i < chip.nx(); ++i) std::printf("%6.1f", chip.siliconT(i, j));
            std::printf("\n");
        }
        return 0;
    }
    if (c.level == "server") {
        ServerRig rig(c.server, fluid, c.serverFlowLpm, c.server.balanceT);
        rig.server().setActivity(makeWorkload(c.workload)->u(0.0));
        rig.solveSteady();
        Server& s = rig.server();
        const double rho = fluid->props(rig.supplyT()).rho;
        std::printf("Server steady: total %.3f L/min, server dp %.2f kPa, return %.3f C, indicative pump %.3f W\n",
                    c.serverFlowLpm, rig.totalDp() / 1e3, s.returnT(), rig.indicativePumpPower());
        for (int b = 0; b < Server::kBranches; ++b)
            std::printf("  %-4s flow %.4f L/min, branch dp %.2f kPa (valve %.2f kPa), power %.0f W, outlet %.3f C, peak %.2f C\n",
                        chipName(b), m3sToLpm(s.branchFlow(b) / rho), s.branchDp(rig.network(), b) / 1e3,
                        s.valveDpRef(b) / 1e3, s.chip(b).power(), s.chip(b).outletT(), s.chip(b).peakT());
        std::printf("  Air group %.0f W -> room air, component %.1f C, fans %.2f W\n", s.air().power(),
                    s.air().temperature(), s.air().fanPower());
        std::printf("  Liquid energy residual %.2e W\n", s.ledger().residual());
        return 0;
    }
    if (c.level == "plant") {
        auto dc = makeDataCenter(c, c.viscosityScale);
        dc->solveSteady();
        printRackSummary(dc->rack());
        printPlantSummary(*dc);
        return 0;
    }
    Rack rack(c.rack, c.server, std::make_shared<WaterIAPWS>());
    rack.setViscosityScale(c.viscosityScale);
    applyWorkloads(rack, c);
    rack.solveSteady();
    printRackSummary(rack);
    return 0;
}

// ---------------------------------------------------------------------------------------------------
int cmdRun(const Args& a) {
    SimConfig c = configFrom(a);
    const double duration = a.num("duration", 30.0), dt = a.num("dt", 0.005), every = a.num("every", 0.1);
    const std::string out = a.get("out", "results/run_" + c.level + ".csv");
    auto f = openCsv(out);
    const int steps = static_cast<int>(std::lround(duration / dt));
    const int stride = std::max(1, static_cast<int>(std::lround(every / dt)));
    const auto t0 = std::chrono::steady_clock::now();
    auto fluid = makeFluid(c.viscosityScale);
    auto wl = makeWorkload(c.workload);

    if (c.level == "chip") {
        ChipGrid chip(c.server.gpu, fluid);
        chip.setMassFlow(fluid->props(c.chipInletC).rho * lpmToM3s(c.chipFlowLpm));
        chip.setInletT(c.chipInletC);
        f << "t_s,activity,power_W,flow_lpm,inlet_C,outlet_C,peak_C,peak_i,peak_j,mean_si_C,to_coolant_W,storage_W,residual_W\n";
        for (int n = 1; n <= steps; ++n) {
            const double t = n * dt;
            chip.setActivity(wl->u(t));
            chip.step(dt);
            if (n % stride == 0) {
                const auto& l = chip.ledger();
                f << t << ',' << chip.activity() << ',' << chip.power() << ',' << c.chipFlowLpm << ',' << c.chipInletC << ','
                  << chip.outletT() << ',' << chip.peakT() << ',' << chip.peakI() << ',' << chip.peakJ() << ','
                  << chip.meanSiliconT() << ',' << l.toCoolant << ',' << l.storageRate << ',' << l.residual << '\n';
            }
        }
        std::printf("chip run: %d steps, final peak %.2f C, outlet %.3f C\n", steps, chip.peakT(), chip.outletT());
    } else if (c.level == "server") {
        ServerRig rig(c.server, fluid, c.serverFlowLpm, c.server.balanceT);
        Server& s = rig.server();
        const double rho = fluid->props(rig.supplyT()).rho;
        f << "t_s,activity,liquid_W,air_W,gpu1_lpm,gpu2_lpm,cpu_lpm,server_dp_kPa,return_C,gpu1_peak_C,gpu2_peak_C,"
             "cpu_peak_C,air_component_C,fan_W,indicative_pump_W,residual_W\n";
        for (int n = 1; n <= steps; ++n) {
            const double t = n * dt;
            s.setActivity(wl->u(t));
            rig.step(dt);
            if (n % stride == 0)
                f << t << ',' << s.activity() << ',' << s.liquidPower() << ',' << s.air().power() << ','
                  << m3sToLpm(s.branchFlow(0) / rho) << ',' << m3sToLpm(s.branchFlow(1) / rho) << ','
                  << m3sToLpm(s.branchFlow(2) / rho) << ',' << rig.totalDp() / 1e3 << ',' << s.returnT() << ','
                  << s.chip(0).peakT() << ',' << s.chip(1).peakT() << ',' << s.chip(2).peakT() << ','
                  << s.air().temperature() << ',' << s.air().fanPower() << ',' << rig.indicativePumpPower() << ','
                  << s.ledger().residual() << '\n';
        }
        std::printf("server run: %d steps, final peak %.2f C, return %.3f C\n", steps, s.peakT(), s.returnT());
    } else if (c.level == "plant") {
        auto dc = makeDataCenter(c, c.viscosityScale);
        if (a.has("start-steady")) dc->solveSteady();
        Rack& r = dc->rack();
        f << "t_s,outdoor_C,it_W,liquid_W,flow_lpm,pump_speed,supply_C,return_C,facility_supply_C,facility_return_C,"
             "cdu_W,cdu_eps,drycooler_W,fan_speed,drycooler_fan_W,chiller_W_heat,chiller_W,facility_pump_W,equip_pump_W,"
             "server_fans_W,cooling_W,partial_pue,peak_C,plant_residual_W,liquid_residual_W\n";
        for (int n = 1; n <= steps; ++n) {
            dc->step(dt);
            if (n % stride == 0) {
                const PlantState& p = dc->plant().state();
                f << r.time() << ',' << p.outdoorT << ',' << r.itPower() << ',' << r.ledger().powerIn << ','
                  << r.totalFlowLpm() << ',' << r.pumpSpeed() << ',' << p.equipSupplyT << ',' << r.returnT() << ','
                  << p.facilitySupplyT << ',' << p.facilityReturnT << ',' << p.cduHeat << ',' << p.cduEffectiveness << ','
                  << p.dryCoolerHeat << ',' << p.fanSpeed << ',' << p.fanPower << ',' << p.chillerHeat << ','
                  << p.chillerPower << ',' << p.facilityPumpPower << ',' << r.pumpPower() << ',' << r.fanPower() << ','
                  << dc->coolingPower() << ',' << dc->partialPue() << ',' << r.hotspot().t << ',' << p.residual << ','
                  << r.ledger().residual() << '\n';
            }
        }
        std::printf("plant run: %d steps of %.3f s\n", steps, dt);
        printRackSummary(r);
        printPlantSummary(*dc);
        std::printf("Energy over run: IT %.1f kJ, cooling %.2f kJ (dry-cooler fans %.2f kJ, chiller %.2f kJ)\n",
                    r.itEnergy() / 1e3, dc->coolingEnergy() / 1e3, dc->plantFanEnergy() / 1e3, dc->chillerEnergy() / 1e3);
    } else {
        Rack r(c.rack, c.server, std::make_shared<WaterIAPWS>());
        r.setViscosityScale(c.viscosityScale);
        applyWorkloads(r, c);
        f << "t_s,it_W,liquid_W,air_W,flow_lpm,pump_speed,pump_head_kPa,pump_W,fan_W,supply_C,return_C,peak_C,"
             "peak_server,peak_chip,peak_i,peak_j,flow_spread,residual_W,pump_J,fan_J,it_J";
        for (int k = 0; k < r.serverCount(); ++k) f << ",s" << k + 1 << "_u,s" << k + 1 << "_lpm,s" << k + 1 << "_peak_C";
        f << '\n';
        const double rho = r.network().fluid().props(r.supplyT()).rho;
        for (int n = 1; n <= steps; ++n) {
            r.step(dt);
            if (n % stride == 0) {
                const Hotspot h = r.hotspot();
                const LiquidLedger l = r.ledger();
                f << r.time() << ',' << r.itPower() << ',' << l.powerIn << ',' << r.airHeat() << ',' << r.totalFlowLpm()
                  << ',' << r.pumpSpeed() << ',' << r.pumpHead() / 1e3 << ',' << r.pumpPower() << ',' << r.fanPower()
                  << ',' << r.supplyT() << ',' << r.returnT() << ',' << h.t << ',' << h.server + 1 << ','
                  << chipName(h.branch) << ',' << h.i << ',' << h.j << ',' << r.flowMaldistribution() << ','
                  << l.residual() << ',' << r.pumpEnergy() << ',' << r.fanEnergy() << ',' << r.itEnergy();
                for (int k = 0; k < r.serverCount(); ++k)
                    f << ',' << r.server(k).activity() << ',' << m3sToLpm(r.server(k).totalFlow() / rho) << ','
                      << r.server(k).peakT();
                f << '\n';
            }
        }
        std::printf("rack run: %d steps of %.3f s\n", steps, dt);
        printRackSummary(r);
        std::printf("Energy over run: IT %.1f kJ, pump %.3f kJ, fans %.3f kJ\n", r.itEnergy() / 1e3,
                    r.pumpEnergy() / 1e3, r.fanEnergy() / 1e3);
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("wrote %s (%.2f s wall time)\n", out.c_str(), secs);
    return 0;
}

// ---------------------------------------------------------------------------------------------------
int cmdResolution(const Args& a) {
    SimConfig c = configFrom(a);
    auto fluid = makeFluid(c.viscosityScale);
    std::printf("Grid-resolution sensitivity, one GPU at full power, %.2f L/min, inlet %.1f C, power map '%s'\n",
                c.chipFlowLpm, c.chipInletC, c.server.gpu.power.map.pattern.c_str());
    std::printf("  grid     peak C   hotspot (x,y) mm   outlet C   mean Si C   change vs previous\n");
    double prev = NAN;
    auto f = openCsv(a.get("out", "results/resolution.csv"));
    f << "n,peak_C,peak_x_mm,peak_y_mm,outlet_C,mean_si_C\n";
    for (int n : {5, 10, 20, 40}) {
        ChipParams p = c.server.gpu;
        p.nx = p.ny = n;
        ChipGrid chip(p, fluid);
        chip.setMassFlow(fluid->props(c.chipInletC).rho * lpmToM3s(c.chipFlowLpm));
        chip.setInletT(c.chipInletC);
        chip.setActivity(1.0);
        chip.solveSteady();
        const double x = (chip.peakI() + 0.5) * p.die.lx / n * 1e3, y = (chip.peakJ() + 0.5) * p.die.ly / n * 1e3;
        std::printf("  %2dx%-2d  %7.2f   (%5.2f,%5.2f)     %7.3f   %8.2f    %s\n", n, n, chip.peakT(), x, y,
                    chip.outletT(), chip.meanSiliconT(),
                    std::isnan(prev) ? "-" : (std::to_string(chip.peakT() - prev).substr(0, 6) + " K").c_str());
        f << n << ',' << chip.peakT() << ',' << x << ',' << y << ',' << chip.outletT() << ',' << chip.meanSiliconT() << '\n';
        prev = chip.peakT();
    }
    return 0;
}

// ---------------------------------------------------------------------------------------------------
/// Bisection for x in [lo, hi] such that g(x) = target, g monotonic. Returns NAN if not bracketed.
double bisect(const std::function<double(double)>& g, double lo, double hi, double target, bool decreasing) {
    double glo = g(lo), ghi = g(hi);
    auto above = [&](double v) { return decreasing ? v > target : v < target; };
    if (above(ghi) || !above(glo)) return NAN;
    for (int it = 0; it < 60 && (hi - lo) > 1e-9 * hi; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (above(g(mid))) lo = mid; else hi = mid;
    }
    return hi;
}

int cmdSweepViscosity(const Args& a) {
    SimConfig c = configFrom(a);
    const std::vector<double> scales = parseList(a.get("scales", "1,2,4"));
    const double limit = c.rack.targetPeakC;
    const std::string out = a.get("out", "results/viscosity_" + c.level + ".csv");
    auto f = openCsv(out);
    std::printf("Viscosity sensitivity (%s level). Only viscosity changes; rho, cp, k stay water's.\n", c.level.c_str());
    std::printf("Power map '%s', peak limit %.1f C (study assumption).\n\n", c.server.gpu.power.map.pattern.c_str(), limit);

    if (c.level == "chip") {
        f << "mu_scale,case,flow_lpm,dp_kPa,hydraulic_W,peak_C,outlet_C\n";
        double dp1 = NAN;
        std::printf(" mu   case               flow L/min   dp kPa   hydraulic W   peak C   outlet C\n");
        for (double s : scales) {
            auto fluid = makeFluid(s);
            ChipParams p = c.server.gpu;
            ChannelEdge edge("chip", p.channels);
            auto steadyAt = [&](double lpm, ChipGrid& chip) {
                chip.setMassFlow(fluid->props(c.chipInletC).rho * lpmToM3s(lpm));
                chip.setInletT(c.chipInletC);
                chip.setActivity(1.0);
                chip.solveSteady();
                edge.setTemperature(chip.meanFluidT());
                return edge.dp(fluid->props(c.chipInletC).rho * lpmToM3s(lpm), *fluid);
            };
            auto report = [&](const char* name, double lpm) {
                ChipGrid chip(p, fluid);
                const double dp = steadyAt(lpm, chip);
                std::printf(" %3.1f  %-18s %10.4f %8.2f %12.4f %8.2f %9.3f\n", s, name, lpm, dp / 1e3, dp * lpmToM3s(lpm),
                            chip.peakT(), chip.outletT());
                f << s << ',' << name << ',' << lpm << ',' << dp / 1e3 << ',' << dp * lpmToM3s(lpm) << ',' << chip.peakT()
                  << ',' << chip.outletT() << '\n';
                return dp;
            };
            const double dpA = report("fixed-flow", c.chipFlowLpm);
            if (std::isnan(dp1)) dp1 = dpA;
            // (b) same pressure difference as the 1x case: flow drops.
            const double qB = bisect([&](double lpm) { ChipGrid ch(p, fluid); return steadyAt(lpm, ch); }, 1e-3, 5.0, dp1, false);
            if (!std::isnan(qB)) report("fixed-dp", qB);
            // (c) minimum flow that keeps the peak at the limit.
            const double qC = bisect([&](double lpm) { ChipGrid ch(p, fluid); steadyAt(lpm, ch); return ch.peakT(); },
                                     0.01, 5.0, limit, true);
            if (!std::isnan(qC)) report("min-flow-for-limit", qC);
        }
    } else if (c.level == "server") {
        f << "mu_scale,case,flow_lpm,server_dp_kPa,indicative_pump_W,peak_C,gpu1_lpm,gpu2_lpm,cpu_lpm\n";
        std::printf(" mu   case               flow L/min   dp kPa   pump W*   peak C   GPU1/GPU2/CPU L/min\n");
        double dp1 = NAN;
        for (double s : scales) {
            auto fluid = makeFluid(s);
            auto at = [&](double lpm) {
                auto rig = std::make_unique<ServerRig>(c.server, fluid, lpm, c.server.balanceT);
                rig->server().setActivity(1.0);
                rig->solveSteady();
                return rig;
            };
            auto report = [&](const char* name, double lpm) {
                auto rig = at(lpm);
                Server& sv = rig->server();
                const double rho = fluid->props(rig->supplyT()).rho;
                std::printf(" %3.1f  %-18s %10.4f %8.2f %8.3f %8.2f   %5.3f/%5.3f/%5.3f\n", s, name, lpm,
                            rig->totalDp() / 1e3, rig->indicativePumpPower(), sv.peakT(), m3sToLpm(sv.branchFlow(0) / rho),
                            m3sToLpm(sv.branchFlow(1) / rho), m3sToLpm(sv.branchFlow(2) / rho));
                f << s << ',' << name << ',' << lpm << ',' << rig->totalDp() / 1e3 << ',' << rig->indicativePumpPower() << ','
                  << sv.peakT() << ',' << m3sToLpm(sv.branchFlow(0) / rho) << ',' << m3sToLpm(sv.branchFlow(1) / rho) << ','
                  << m3sToLpm(sv.branchFlow(2) / rho) << '\n';
                return rig->totalDp();
            };
            const double dpA = report("fixed-flow", c.serverFlowLpm);
            if (std::isnan(dp1)) dp1 = dpA;
            const double qB = bisect([&](double lpm) { return at(lpm)->totalDp(); }, 0.01, 10.0, dp1, false);
            if (!std::isnan(qB)) report("fixed-dp", qB);
            const double qC = bisect([&](double lpm) { return at(lpm)->server().peakT(); }, 0.05, 10.0, limit, true);
            if (!std::isnan(qC)) report("min-flow-for-limit", qC);
        }
        std::printf("* indicative pump power = server dp x flow / (0.6 x 0.9); the real pump arrives at rack level.\n");
    } else if (c.level == "plant") {
        f << "mu_scale,case,pump_speed,flow_lpm,pump_head_kPa,equip_pump_W,cdu_UA,supply_C,drycooler_fan_W,chiller_W,"
             "facility_pump_W,cooling_W,partial_pue,peak_C\n";
        std::printf("Outdoor %.1f C, supply setpoint %.1f C, chiller %s.\n", c.plant.weather.meanC,
                    c.plant.supplySetpointC, c.plant.chiller.enabled ? "enabled" : "disabled");
        std::printf(" mu   case                 flow L/min  head kPa  pump W  CDU UA  supply C  fans W  chiller W  cooling W  pPUE    peak C\n");
        double speed1 = NAN;
        for (double s : scales) {
            auto run = [&](const char* name, PumpMode mode, double speed) {
                SimConfig cc = c;
                cc.rack.pumpMode = mode;
                cc.rack.targetPeakC = limit;
                if (mode == PumpMode::FixedSpeed) cc.rack.pumpSpeed = speed;
                auto dc = makeDataCenter(cc, s);
                dc->solveSteady();
                const Rack& r = dc->rack();
                const PlantState& p = dc->plant().state();
                std::printf(" %3.1f  %-20s %10.3f %9.2f %7.2f %7.0f %9.2f %7.1f %10.1f %10.1f  %.4f %7.2f\n", s, name,
                            r.totalFlowLpm(), r.pumpHead() / 1e3, r.pumpPower(), p.cduUA, p.equipSupplyT, p.fanPower,
                            p.chillerPower, dc->coolingPower(), dc->partialPue(), r.hotspot().t);
                f << s << ',' << name << ',' << r.pumpSpeed() << ',' << r.totalFlowLpm() << ',' << r.pumpHead() / 1e3 << ','
                  << r.pumpPower() << ',' << p.cduUA << ',' << p.equipSupplyT << ',' << p.fanPower << ',' << p.chillerPower
                  << ',' << p.facilityPumpPower << ',' << dc->coolingPower() << ',' << dc->partialPue() << ','
                  << r.hotspot().t << '\n';
                return r.pumpSpeed();
            };
            const double sp = run("fixed-flow", PumpMode::FixedFlow, 1.0);
            if (std::isnan(speed1)) speed1 = sp;
            run("fixed-pump-speed", PumpMode::FixedSpeed, speed1);
            run("hold-peak-at-limit", PumpMode::PeakTemperature, 1.0);
        }
    } else {
        f << "mu_scale,case,pump_speed,flow_lpm,pump_head_kPa,pump_W,pump_eta,peak_C,return_C,flow_spread\n";
        std::printf(" mu   case                 speed   flow L/min  head kPa   pump W   eta   peak C  return C  spread\n");
        double speed1 = NAN;
        for (double s : scales) {
            auto run = [&](const char* name, PumpMode mode, double speed) {
                Rack r(c.rack, c.server, std::make_shared<WaterIAPWS>());
                r.setViscosityScale(s);
                applyWorkloads(r, c);
                r.setPumpMode(mode);
                if (mode == PumpMode::FixedSpeed) r.setPumpSpeed(speed);
                r.setTargetPeak(limit);
                r.solveSteady();
                const Hotspot h = r.hotspot();
                std::printf(" %3.1f  %-20s %6.3f %10.3f %9.2f %8.2f  %4.2f %7.2f %8.3f  %5.2f%%\n", s, name, r.pumpSpeed(),
                            r.totalFlowLpm(), r.pumpHead() / 1e3, r.pumpPower(), r.pumpEfficiency(), h.t, r.returnT(),
                            100 * r.flowMaldistribution());
                f << s << ',' << name << ',' << r.pumpSpeed() << ',' << r.totalFlowLpm() << ',' << r.pumpHead() / 1e3 << ','
                  << r.pumpPower() << ',' << r.pumpEfficiency() << ',' << h.t << ',' << r.returnT() << ','
                  << r.flowMaldistribution() << '\n';
                return r.pumpSpeed();
            };
            const double sp = run("fixed-flow", PumpMode::FixedFlow, 1.0);
            if (std::isnan(speed1)) speed1 = sp;
            run("fixed-pump-speed", PumpMode::FixedSpeed, speed1);
            run("hold-peak-at-limit", PumpMode::PeakTemperature, 1.0);
        }
        std::printf("fixed-pump-speed uses the speed that gives %.1f L/min at %.1fx viscosity.\n", c.rack.targetFlowLpm,
                    scales.front());
    }
    std::printf("\nwrote %s\n", out.c_str());
    return 0;
}

// ---------------------------------------------------------------------------------------------------
int cmdUnevenLoad(const Args& a) {
    SimConfig c = configFrom(a);
    auto f = openCsv(a.get("out", "results/uneven_load.csv"));
    f << "layout,pattern,server,activity,flow_lpm,peak_C,return_C\n";
    struct Pattern { const char* name; std::vector<double> u; };
    const std::vector<Pattern> patterns = {
        {"bottom-half-busy", {1, 1, 1, 1, 0, 0, 0, 0}},
        {"top-half-busy", {0, 0, 0, 0, 1, 1, 1, 1}},
        {"all-busy", {1, 1, 1, 1, 1, 1, 1, 1}},
    };
    for (ManifoldLayout layout : {ManifoldLayout::DirectReturn, ManifoldLayout::ReverseReturn}) {
        for (const auto& pat : patterns) {
            RackParams rp = c.rack;
            rp.layout = layout;
            Rack r(rp, c.server, std::make_shared<WaterIAPWS>());
            r.setViscosityScale(c.viscosityScale);
            for (int k = 0; k < r.serverCount(); ++k)
                r.setWorkload(k, std::make_shared<ConstantWorkload>(pat.u[k % pat.u.size()]));
            r.solveSteady();
            const Hotspot h = r.hotspot();
            std::printf("\n== %s, %s ==  hottest: server %d %s %.2f C, flow spread %.2f%%, pump %.2f W\n",
                        toString(layout), pat.name, h.server + 1, chipName(h.branch), h.t,
                        100 * r.flowMaldistribution(), r.pumpPower());
            const double rho = r.network().fluid().props(r.supplyT()).rho;
            for (int k = r.serverCount() - 1; k >= 0; --k) {
                const Server& s = r.server(k);
                std::printf("  server %d (top=8): u=%.0f  flow %.4f L/min  peak %.2f C  return %.2f C\n", k + 1,
                            s.activity(), m3sToLpm(s.totalFlow() / rho), s.peakT(), s.returnT());
                f << toString(layout) << ',' << pat.name << ',' << k + 1 << ',' << s.activity() << ','
                  << m3sToLpm(s.totalFlow() / rho) << ',' << s.peakT() << ',' << s.returnT() << '\n';
            }
        }
    }
    return 0;
}

int cmdClimateSweep(const Args& a) {
    SimConfig c = configFrom(a);
    const std::vector<double> temps = parseList(a.get("temps", "0,10,20,25,30,35,40,45"));
    auto f = openCsv(a.get("out", "results/climate.csv"));
    f << "chiller,outdoor_C,supply_C,fan_speed,drycooler_fan_W,chiller_W,chiller_cop,facility_pump_W,equip_pump_W,"
         "cooling_W,partial_pue,peak_C,evap_equiv_Lph\n";
    std::printf("Outdoor-temperature sweep: rack at full load, supply setpoint %.1f C, pump %s.\n",
                c.plant.supplySetpointC, toString(c.rack.pumpMode));
    for (bool chill : {true, false}) {
        std::printf("\n== trim chiller %s ==\n", chill ? "enabled" : "disabled");
        std::printf(" outdoor C  supply C  fans   fan W  chiller W  COP   cooling W   pPUE    peak C   evap-equiv L/h\n");
        for (double t : temps) {
            SimConfig cc = c;
            cc.plant.weather.meanC = t;
            cc.plant.weather.amplitudeK = 0.0;
            cc.plant.chiller.enabled = chill;
            auto dc = makeDataCenter(cc, c.viscosityScale);
            dc->solveSteady();
            const PlantState& p = dc->plant().state();
            std::printf(" %8.1f  %8.2f  %4.0f%% %7.1f %10.1f %5.2f %10.1f  %.4f %8.2f %10.2f\n", t, p.equipSupplyT,
                        100 * p.fanSpeed, p.fanPower, p.chillerPower, p.chillerCop, dc->coolingPower(), dc->partialPue(),
                        dc->rack().hotspot().t, dc->evaporativeEquivalentLph());
            f << (chill ? "on" : "off") << ',' << t << ',' << p.equipSupplyT << ',' << p.fanSpeed << ',' << p.fanPower << ','
              << p.chillerPower << ',' << p.chillerCop << ',' << p.facilityPumpPower << ',' << dc->rack().pumpPower() << ','
              << dc->coolingPower() << ',' << dc->partialPue() << ',' << dc->rack().hotspot().t << ','
              << dc->evaporativeEquivalentLph() << '\n';
        }
    }
    return 0;
}

int usage() {
    std::printf(
        "mfsim - microfluidic cooling data-center simulation (headless)\n\n"
        "  mfsim info            [--config F]\n"
        "  mfsim steady          [--config F] [--level chip|server|rack]\n"
        "  mfsim run             [--config F] [--level L] [--duration 30] [--dt 0.005] [--every 0.1] [--out f.csv]\n"
        "  mfsim resolution      [--config F]\n"
        "  mfsim sweep-viscosity [--config F] [--level L] [--scales 1,2,4] [--out f.csv]\n"
        "  mfsim uneven-load     [--config F]\n"
        "  mfsim climate-sweep   [--config F] [--temps 0,10,20,30,40]\n\n"
        "Levels: chip | server | rack (fixed 35 C supply) | plant (rack + CDU + dry cooler + trim chiller)\n"
        "Common options: --viscosity X, --outdoor C, --chiller on|off; run: --start-steady\n");
    return 1;
}

} // namespace

int main(int argc, char** argv) {
    const Args a = parse(argc, argv);
    try {
        if (a.cmd == "info") return cmdInfo(a);
        if (a.cmd == "steady") return cmdSteady(a);
        if (a.cmd == "run") return cmdRun(a);
        if (a.cmd == "resolution") return cmdResolution(a);
        if (a.cmd == "sweep-viscosity") return cmdSweepViscosity(a);
        if (a.cmd == "uneven-load") return cmdUnevenLoad(a);
        if (a.cmd == "climate-sweep") return cmdClimateSweep(a);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
    return usage();
}
