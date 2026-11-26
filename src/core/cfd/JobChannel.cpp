#include "cfd/JobChannel.hpp"
#include "Units.hpp"
#include "thermal/PowerMap.hpp"

#include <algorithm>
#include <cmath>

namespace mf::cfd {

JobChannel::JobChannel(const JobChannelParams& p, std::shared_ptr<const JobTraceWorkload> job)
    : p_(p), job_(std::move(job)), fluid_(makeCoolant(p.coolant)) {
    props_ = fluid_->props(p_.inletC + 5.0);
    const ChannelGeometry& g = p_.gpu.channels;

    // Hotspot shape along the flow for the channel through the middle of the die.
    const int gx = p_.gpu.nx, gy = p_.gpu.ny;
    const std::vector<double> w = buildPowerWeights(p_.gpu.power.map, gx, gy);
    const int jc = gy / 2;
    shape_.resize(gx);
    double mean = 0.0;
    for (int i = 0; i < gx; ++i) { shape_[i] = w[jc * gx + i] * gx * gy; mean += shape_[i]; }
    // shape_ is relative to the die-average flux, so the central (hotspot) row may average above 1.
    (void)mean;

    Params q;
    q.lx = g.length;
    q.ly = g.width;
    q.ny = std::max(4, p_.cellsAcross);
    q.nx = std::max(10, static_cast<int>(std::lround(q.ny * g.length / g.width / 3.0)));
    q.nu = props_.mu / props_.rho;
    q.alpha = props_.k / (props_.rho * props_.cp);
    q.k = props_.k;
    q.t0 = p_.inletC;
    q.gammaMin = 0.5;
    q.left.thermal = ThermalBC::FixedTemperature;
    q.left.T = p_.inletC;
    q.left.u = p_.velocity;
    q.right.velocity = VelocityBC::Outflow;
    for (Side* s : {&q.bottom, &q.top}) {
        s->thermal = ThermalBC::HeatFlux;
        s->qShape = shape_;
    }
    ns_ = std::make_unique<NavierStokes2D>(q);
    setWallFlux(p_.gpu.power.at(0.0));

    // Reduced model with the same coolant, flow and power map.
    chip_ = std::make_unique<ChipGrid>(p_.gpu, fluid_, p_.inletC);
    chip_->setMassFlow(props_.rho * p_.velocity * g.area() * g.count);
    chip_->setInletT(p_.inletC);
    chip_->setActivity(0.0);
}

double JobChannel::reynolds() const { return p_.velocity * 2.0 * p_.gpu.channels.width * props_.rho / props_.mu; }

void JobChannel::setWallFlux(double powerW) {
    // One channel's heat enters through its two rib walls (height = channel depth, length = channel length).
    const ChannelGeometry& g = p_.gpu.channels;
    const double qAvg = powerW / g.count / (2.0 * g.length * g.depth);
    Params& q = ns_->params();
    q.bottom.q = qAvg;
    q.top.q = qAvg;
}

long JobChannel::solveFlow() {
    ns_->runToSteady(1e-4, p_.gpu.channels.length / p_.velocity, 400000);
    return ns_->steps();
}

JobChannelSample JobChannel::step() {
    const double dt = p_.dt;
    t_ += dt;
    const TracePoint tp = job_->at(t_);
    const double power = p_.gpu.power.at(tp.u);
    setWallFlux(power);
    JobChannelSample s;
    s.sweeps = ns_->stepEnergyImplicit(dt, 1e-6, 4000);
    chip_->setActivity(tp.u);
    chip_->step(dt);

    const ChannelGeometry& g = p_.gpu.channels;
    s.t = t_;
    s.u = tp.u;
    s.phase = tp.phase;
    s.powerW = power;
    s.qMean = ns_->params().bottom.q;
    s.wallMeanC = 0.5 * (ns_->wallTemperatureMean(2) + ns_->wallTemperatureMean(3));
    s.wallMaxC = std::max(ns_->wallTemperatureMax(2), ns_->wallTemperatureMax(3));
    s.outletC = ns_->bulkTemperatureColumn(ns_->nx());
    s.heatInW = ns_->wallHeatInput() * g.depth;
    s.heatOutW = ns_->advectedHeat(props_.rho * props_.cp) * g.depth;
    s.chipPeakC = chip_->peakT();
    s.chipOutletC = chip_->outletT();
    return s;
}

} // namespace mf::cfd
