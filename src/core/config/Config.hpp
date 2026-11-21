#pragma once
#include "system/Rack.hpp"
#include "system/Server.hpp"
#include "workload/Workload.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mf {

/// Description of a workload, as read from JSON.
struct WorkloadSpec {
    std::string type = "constant"; ///< constant | step | square | trace
    double u = 1.0;
    double before = 0.0, after = 1.0, tStep = 1.0;
    double low = 0.2, high = 1.0, period = 10.0, duty = 0.5, phase = 0.0;
    std::string file;
};
std::shared_ptr<Workload> makeWorkload(const WorkloadSpec& w);

/// Everything needed to build a chip, server or rack simulation. Defaults in code are the
/// proposal's baseline values; JSON files override any subset (keys carry their units).
struct SimConfig {
    std::string level = "rack";       ///< chip | server | rack
    double viscosityScale = 1.0;
    ServerParams server;
    RackParams rack;
    double chipFlowLpm = 0.5;         ///< chip rig flow
    double chipInletC = 35.0;
    double serverFlowLpm = 1.5;       ///< server rig total flow
    WorkloadSpec workload;            ///< applied to every server unless overridden
    std::vector<WorkloadSpec> serverWorkloads;
};

SimConfig defaultConfig();
SimConfig loadConfig(const std::string& path);           ///< defaults + overrides from file
void applyConfigJson(SimConfig& c, const std::string& jsonText);

} // namespace mf
