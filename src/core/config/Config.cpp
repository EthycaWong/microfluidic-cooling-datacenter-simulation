#include "config/Config.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace mf {
using nlohmann::json;

namespace {

/// Read j[key] into out (multiplied by scale) when present.
template <class T>
void get(const json& j, const char* key, T& out, double scale = 1.0) {
    if (!j.contains(key)) return;
    if constexpr (std::is_same_v<T, double>) out = j.at(key).get<double>() * scale;
    else out = j.at(key).get<T>();
}

void readPowerMap(const json& j, PowerMapSpec& m) {
    get(j, "pattern", m.pattern);
    get(j, "center_x", m.centerX);
    get(j, "center_y", m.centerY);
    get(j, "size_x", m.sizeX);
    get(j, "size_y", m.sizeY);
    get(j, "core_share", m.coreShare);
    get(j, "custom_nx", m.customNx);
    get(j, "custom_ny", m.customNy);
    if (j.contains("custom")) m.custom = j.at("custom").get<std::vector<double>>();
}

void readChip(const json& j, ChipParams& c) {
    if (j.contains("die")) {
        const json& d = j["die"];
        get(d, "lx_mm", c.die.lx, 1e-3);
        get(d, "ly_mm", c.die.ly, 1e-3);
        get(d, "thickness_mm", c.die.thickness, 1e-3);
    }
    if (j.contains("channels")) {
        const json& d = j["channels"];
        get(d, "count", c.channels.count);
        get(d, "width_um", c.channels.width, 1e-6);
        get(d, "depth_um", c.channels.depth, 1e-6);
        get(d, "length_mm", c.channels.length, 1e-3);
        get(d, "rib_um", c.channels.rib, 1e-6);
        get(d, "k_inlet", c.channels.kInlet);
        get(d, "k_outlet", c.channels.kOutlet);
        get(d, "hagenbach", c.channels.hagenbach);
    }
    if (j.contains("grid")) {
        get(j["grid"], "nx", c.nx);
        get(j["grid"], "ny", c.ny);
    }
    if (j.contains("power")) {
        const json& p = j["power"];
        get(p, "idle_W", c.power.idleW);
        get(p, "full_W", c.power.fullW);
        if (p.contains("map")) readPowerMap(p["map"], c.power.map);
    }
    if (j.contains("silicon")) {
        get(j["silicon"], "k300_W_per_mK", c.kSi300);
        get(j["silicon"], "k_exponent", c.kSiExponent);
    }
    if (j.contains("cold_plate")) {
        const json& p = j["cold_plate"];
        ColdPlateParams& cp = c.coldPlate;
        get(p, "tim_m2K_per_W", cp.timResistance);
        get(p, "base_mm", cp.baseThickness, 1e-3);
        get(p, "hA_ref_W_per_K", cp.hARef);
        get(p, "q_ref_lpm", cp.qRefLpm);
        get(p, "hA_exponent", cp.hAExponent);
        get(p, "fluid_volume_ml", cp.fluidVolume, 1e-6);
        get(p, "dp_ref_kPa", cp.dpRefPa, 1e3);
        get(p, "laminar_fraction", cp.laminarFraction);
    }
    if (j.contains("cooling")) {
        const std::string m = j["cooling"].get<std::string>();
        c.mode = (m == "cold-plate") ? CoolingMode::ColdPlate : CoolingMode::Microchannel;
    }
}

void readWorkload(const json& j, WorkloadSpec& w) {
    get(j, "type", w.type);
    get(j, "u", w.u);
    get(j, "before", w.before);
    get(j, "after", w.after);
    get(j, "t_step_s", w.tStep);
    get(j, "low", w.low);
    get(j, "high", w.high);
    get(j, "period_s", w.period);
    get(j, "duty", w.duty);
    get(j, "phase_s", w.phase);
    get(j, "file", w.file);
}

} // namespace

std::shared_ptr<Workload> makeWorkload(const WorkloadSpec& w) {
    if (w.type == "constant") return std::make_shared<ConstantWorkload>(w.u);
    if (w.type == "step") return std::make_shared<StepWorkload>(w.before, w.after, w.tStep);
    if (w.type == "square") return std::make_shared<SquareWorkload>(w.low, w.high, w.period, w.duty, w.phase);
    if (w.type == "trace") return std::make_shared<TraceWorkload>(w.file);
    throw std::runtime_error("unknown workload type: " + w.type);
}

SimConfig defaultConfig() {
    SimConfig c;
    // GPU: silicon microchannels, the proposal's candidate geometry.
    c.server.gpu.name = "GPU";
    c.server.gpu.mode = CoolingMode::Microchannel;
    c.server.gpu.power.idleW = 60.0;
    c.server.gpu.power.fullW = 300.0;
    // CPU: conventional cold plate. Die size and idle power are ASSUMPTIONS.
    c.server.cpu.name = "CPU";
    c.server.cpu.mode = CoolingMode::ColdPlate;
    c.server.cpu.die = DieGeometry{25e-3, 25e-3, 0.8e-3};
    c.server.cpu.power.idleW = 50.0;
    c.server.cpu.power.fullW = 200.0;
    return c;
}

void applyConfigJson(SimConfig& c, const std::string& text) {
    const json j = json::parse(text, nullptr, true, /*ignore_comments=*/true);
    get(j, "level", c.level);
    if (j.contains("fluid")) {
        get(j["fluid"], "viscosity_scale", c.viscosityScale);
        get(j["fluid"], "coolant", c.coolant);
        get(j["fluid"], "facility_coolant", c.facilityCoolant);
    }
    if (j.contains("gpu")) readChip(j["gpu"], c.server.gpu);
    if (j.contains("cpu")) readChip(j["cpu"], c.server.cpu);
    if (j.contains("air")) {
        const json& a = j["air"];
        AirParams& p = c.server.air;
        get(a, "idle_W", p.idleW);
        get(a, "full_W", p.fullW);
        get(a, "room_C", p.roomT);
        get(a, "air_dT_K", p.airDeltaT);
        get(a, "fan_dp_ref_Pa", p.fanDpRef);
        get(a, "flow_ref_m3s", p.flowRef);
        get(a, "fan_eta", p.fanEta);
        get(a, "ua_ref_W_per_K", p.uaRef);
        get(a, "heat_capacity_J_per_K", p.heatCapacity);
    }
    if (j.contains("server")) {
        const json& s = j["server"];
        ServerParams& p = c.server;
        get(s, "hose_d_mm", p.hoseDiameter, 1e-3);
        get(s, "hose_length_m", p.hoseLength);
        get(s, "hose_K", p.hoseK);
        get(s, "qd_dp_ref_kPa", p.qdDpRefPa, 1e3);
        get(s, "qd_q_ref_lpm", p.qdQRefLpm);
        get(s, "auto_balance", p.autoBalance);
        get(s, "branch_target_lpm", p.branchTargetLpm);
        get(s, "rig_total_flow_lpm", c.serverFlowLpm);
    }
    if (j.contains("rack")) {
        const json& r = j["rack"];
        RackParams& p = c.rack;
        get(r, "servers", p.servers);
        get(r, "server_U", p.serverU);
        get(r, "first_server_U", p.firstServerU);
        get(r, "manifold_d_mm", p.manifoldDiameter, 1e-3);
        get(r, "manifold_K_per_tap", p.manifoldKPerTap);
        get(r, "connect_length_m", p.connectLength);
        get(r, "connect_K", p.connectK);
        get(r, "hx_dp_ref_kPa", p.hxDpRefPa, 1e3);
        get(r, "hx_q_ref_lpm", p.hxQRefLpm);
        get(r, "hx_laminar_fraction", p.hxLaminarFraction);
        get(r, "filter_dp_ref_kPa", p.filterDpRefPa, 1e3);
        get(r, "filter_q_ref_lpm", p.filterQRefLpm);
        get(r, "filter_laminar_fraction", p.filterLaminarFraction);
        get(r, "reference_pressure_kPa", p.referencePressurePa, 1e3);
        get(r, "pump_speed", p.pumpSpeed);
        get(r, "target_flow_lpm", p.targetFlowLpm);
        get(r, "target_peak_C", p.targetPeakC);
        get(r, "pi_kp", p.piKp);
        get(r, "pi_ki", p.piKi);
        get(r, "supply_C", p.supplyT);
        get(r, "hydraulic_interval_s", p.hydraulicInterval);
        if (r.contains("layout"))
            p.layout = r["layout"].get<std::string>() == "reverse-return" ? ManifoldLayout::ReverseReturn
                                                                          : ManifoldLayout::DirectReturn;
        if (r.contains("pump_mode")) {
            const std::string m = r["pump_mode"].get<std::string>();
            p.pumpMode = m == "fixed-speed" ? PumpMode::FixedSpeed
                       : m == "peak-temperature" ? PumpMode::PeakTemperature : PumpMode::FixedFlow;
        }
        if (r.contains("pump")) {
            const json& q = r["pump"];
            get(q, "shutoff_kPa", p.pump.shutoffHeadPa, 1e3);
            get(q, "max_flow_lpm", p.pump.maxFlowM3s, 1.0 / 60000.0);
            get(q, "bep_flow_lpm", p.pump.bepFlowM3s, 1.0 / 60000.0);
            get(q, "eta_max", p.pump.etaMax);
            get(q, "eta_motor", p.pump.etaMotor);
            get(q, "min_speed", p.pump.minSpeed);
            get(q, "max_speed", p.pump.maxSpeed);
        }
    }
    if (j.contains("plant")) {
        const json& q = j["plant"];
        PlantParams& p = c.plant;
        get(q, "supply_setpoint_C", p.supplySetpointC);
        get(q, "equip_supply_volume_L", p.equipSupplyVolumeL);
        get(q, "equip_return_volume_L", p.equipReturnVolumeL);
        get(q, "facility_volume_L", p.facilityVolumeL);
        get(q, "fan_kp", p.fanKp);
        get(q, "fan_ki", p.fanKi);
        if (q.contains("cdu")) {
            const json& d = q["cdu"];
            get(d, "hA_equip_W_per_K", p.cdu.hAEquipRef);
            get(d, "hA_facility_W_per_K", p.cdu.hAFacilityRef);
            get(d, "equip_ref_lpm", p.cdu.equipRefLpm);
            get(d, "facility_ref_lpm", p.cdu.facilityRefLpm);
            get(d, "exponent", p.cdu.exponent);
            get(d, "wall_R_K_per_W", p.cdu.wallR);
            get(d, "facility_dp_ref_kPa", p.cdu.facilityDpRefPa, 1e3);
        }
        if (q.contains("dry_cooler")) {
            const json& d = q["dry_cooler"];
            get(d, "UA_W_per_K", p.dryCooler.uaRef);
            get(d, "air_share", p.dryCooler.airShare);
            get(d, "air_exponent", p.dryCooler.airExponent);
            get(d, "water_exponent", p.dryCooler.waterExponent);
            get(d, "max_airflow_m3s", p.dryCooler.maxAirflow);
            get(d, "max_fan_W", p.dryCooler.maxFanPower);
            get(d, "min_fan_speed", p.dryCooler.minFanSpeed);
            get(d, "coil_dp_ref_kPa", p.dryCooler.coilDpRefPa, 1e3);
        }
        if (q.contains("chiller")) {
            const json& d = q["chiller"];
            get(d, "enabled", p.chiller.enabled);
            get(d, "capacity_W", p.chiller.capacityW);
            get(d, "carnot_fraction", p.chiller.carnotFraction);
            get(d, "evap_approach_K", p.chiller.evapApproachK);
            get(d, "cond_approach_K", p.chiller.condApproachK);
            get(d, "max_cop", p.chiller.maxCop);
        }
        if (q.contains("facility")) {
            const json& d = q["facility"];
            get(d, "flow_lpm", p.facility.flowLpm);
            get(d, "piping_dp_ref_kPa", p.facility.pipingDpRefPa, 1e3);
            get(d, "pump_eta", p.facility.pumpEta);
            get(d, "motor_eta", p.facility.motorEta);
        }
        if (q.contains("weather")) {
            const json& d = q["weather"];
            get(d, "outdoor_mean_C", p.weather.meanC);
            get(d, "daily_amplitude_K", p.weather.amplitudeK);
            get(d, "period_s", p.weather.periodS);
            get(d, "phase_s", p.weather.phaseS);
        }
    }
    if (j.contains("chip_rig")) {
        get(j["chip_rig"], "flow_lpm", c.chipFlowLpm);
        get(j["chip_rig"], "inlet_C", c.chipInletC);
    }
    if (j.contains("workload")) readWorkload(j["workload"], c.workload);
    if (j.contains("server_workloads")) {
        c.serverWorkloads.clear();
        for (const auto& w : j["server_workloads"]) {
            WorkloadSpec s = c.workload;
            readWorkload(w, s);
            c.serverWorkloads.push_back(s);
        }
    }
}

SimConfig loadConfig(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open config " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    SimConfig c = defaultConfig();
    applyConfigJson(c, ss.str());
    return c;
}

} // namespace mf
