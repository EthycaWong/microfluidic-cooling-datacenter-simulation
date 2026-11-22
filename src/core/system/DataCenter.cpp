#include "system/DataCenter.hpp"

#include <cmath>

namespace mf {

namespace {
constexpr double kLatentHeat = 2.43e6; // J/kg, water evaporating near 30 C
}

DataCenter::DataCenter(const RackParams& rp, const ServerParams& sp, const PlantParams& pp,
                       std::shared_ptr<const Fluid> baseFluid, bool plantEnabled)
    : rack_(rp, sp, baseFluid), plant_(pp, rack_.fluid(), std::make_shared<WaterIAPWS>()), plantOn_(plantEnabled),
      fixedSupply_(rp.supplyT) {}

void DataCenter::setPlantEnabled(bool on) {
    plantOn_ = on;
    if (!on) rack_.setSupplyT(fixedSupply_);
}

void DataCenter::step(double dt) {
    if (!(dt > 0.0)) return;
    rack_.step(dt);
    if (plantOn_) {
        plant_.step(dt, rack_.time(), rack_.totalMassFlow(), rack_.returnT());
        rack_.setSupplyT(plant_.supplyT());
        fanJ_ += plant_.state().fanPower * dt;
        chJ_ += plant_.state().chillerPower * dt;
    }
    coolJ_ += coolingPower() * dt;
}

void DataCenter::solveSteady() {
    if (!plantOn_) {
        rack_.setSupplyT(fixedSupply_);
        rack_.solveSteady();
        return;
    }
    for (int it = 0; it < 40; ++it) {
        rack_.solveSteady();
        plant_.solveSteady(rack_.time(), rack_.totalMassFlow(), rack_.ledger().toCoolant);
        const double t = plant_.supplyT();
        const bool done = std::abs(t - rack_.supplyT()) < 1e-6;
        rack_.setSupplyT(t);
        if (done) break;
    }
}

double DataCenter::coolingPower() const {
    double p = rack_.pumpPower();
    if (plantOn_) p += plant_.coolingPower();
    return p;
}

double DataCenter::partialPue() const {
    const double it = rack_.itPower() + rack_.fanPower();
    return it > 0 ? (it + coolingPower()) / it : 1.0;
}

double DataCenter::heatToOutdoor() const {
    if (!plantOn_) return 0.0;
    const PlantState& s = plant_.state();
    return s.dryCoolerHeat + s.chillerHeat + s.chillerPower;
}

double DataCenter::evaporativeEquivalentLph() const { return heatToOutdoor() / kLatentHeat * 3600.0; }

} // namespace mf
