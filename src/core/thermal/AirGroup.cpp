#include "thermal/AirGroup.hpp"

#include <algorithm>
#include <cmath>

namespace mf {

void AirGroup::updateFlow() {
    const double needed = power_ / (p_.rhoAir * p_.cpAir * p_.airDeltaT);
    flow_ = std::max(needed, p_.minFlowFrac * p_.flowRef);
}

void AirGroup::step(double dt) {
    updateFlow();
    const double ua = p_.uaRef * std::pow(flow_ / p_.flowRef, 0.8);
    const double cdt = p_.heatCapacity / dt;
    t_ = (cdt * t_ + power_ + ua * p_.roomT) / (cdt + ua); // backward Euler
    toAir_ = ua * (t_ - p_.roomT);
}

void AirGroup::solveSteady() {
    updateFlow();
    const double ua = p_.uaRef * std::pow(flow_ / p_.flowRef, 0.8);
    t_ = p_.roomT + power_ / ua;
    toAir_ = power_;
}

double AirGroup::fanPower() const {
    const double r = flow_ / p_.flowRef;
    return p_.fanDpRef * r * r * flow_ / p_.fanEta;
}

} // namespace mf
