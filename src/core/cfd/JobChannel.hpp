#pragma once
#include "cfd/NavierStokes2D.hpp"
#include "fluid/Fluid.hpp"
#include "thermal/ChipGrid.hpp"
#include "workload/Workload.hpp"

#include <memory>
#include <string>

namespace mf::cfd {

struct JobChannelParams {
    std::string coolant = "water";   ///< water | PG25 | EG25 ... (real properties)
    double inletC = 35.0;
    double velocity = 1.667;         ///< mean channel velocity, m/s (1.667 = 0.5 L/min per GPU)
    int cellsAcross = 12;            ///< Navier-Stokes cells across the 200 um gap
    double dt = 0.05;                ///< energy time step while following the job, s
    ChipParams gpu;                  ///< geometry, power range and power map of the GPU
};

/// One instant of the job-driven run.
struct JobChannelSample {
    double t = 0, u = 0, powerW = 0;
    std::string phase;
    double qMean = 0;                ///< mean wall heat flux, W/m^2
    double wallMeanC = 0, wallMaxC = 0, outletC = 0;
    double heatInW = 0, heatOutW = 0;  ///< one channel, 3-D equivalent (x channel depth)
    double chipPeakC = 0, chipOutletC = 0; ///< reduced chip model, same trace
    int sweeps = 0;
};

/// The CNN workload driving the Navier-Stokes microchannel slice.
///
/// The recorded device activity u(t) sets the GPU power P = P_idle + u (P_full - P_idle). One
/// channel's share, P / N_channels, enters the fluid through the two silicon rib walls of the
/// slice, shaped along the flow like the central row of the GPU's power map (the channel that
/// crosses the hotspot). The velocity field comes from the full Navier-Stokes equations (solved
/// once: forced convection, buoyancy negligible); the energy equation is then marched implicitly
/// through the job. The reduced chip model (ChipGrid) is stepped with the same trace for comparison.
class JobChannel {
public:
    JobChannel(const JobChannelParams& p, std::shared_ptr<const JobTraceWorkload> job);

    /// Solve the velocity field (full Navier-Stokes to steady state); returns steps used.
    long solveFlow();
    /// Advance both models by one energy step; returns the sample.
    JobChannelSample step();
    double time() const { return t_; }
    double duration() const { return job_->durationSim(); }
    NavierStokes2D& solver() { return *ns_; }
    const ChipGrid& chip() const { return *chip_; }
    const Fluid& fluid() const { return *fluid_; }
    double reynolds() const;
    double channelHeatPerWatt() const { return 1.0 / p_.gpu.channels.count; }

private:
    void setWallFlux(double powerW);

    JobChannelParams p_;
    std::shared_ptr<const JobTraceWorkload> job_;
    std::shared_ptr<const Fluid> fluid_;
    FluidProps props_;
    std::unique_ptr<NavierStokes2D> ns_;
    std::unique_ptr<ChipGrid> chip_;
    std::vector<double> shape_;
    double t_ = 0.0;
};

} // namespace mf::cfd
