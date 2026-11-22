#pragma once
#include "system/CoolingPlant.hpp"
#include "system/Rack.hpp"

#include <memory>

namespace mf {

/// Rack + cooling plant. With the plant enabled, the rack's supply temperature comes from the CDU
/// (which depends on the dry cooler, the outdoor air and the trim chiller); with it disabled, the
/// supply is held at RackParams::supplyT as in the rack-only milestone.
class DataCenter {
public:
    DataCenter(const RackParams& rp, const ServerParams& sp, const PlantParams& pp,
               std::shared_ptr<const Fluid> baseFluid, bool plantEnabled = true,
               std::shared_ptr<const Fluid> facilityFluid = nullptr);

    Rack& rack() { return rack_; }
    const Rack& rack() const { return rack_; }
    CoolingPlant& plant() { return plant_; }
    const CoolingPlant& plant() const { return plant_; }
    bool plantEnabled() const { return plantOn_; }
    void setPlantEnabled(bool on);
    /// Supply temperature used while the plant is disabled.
    void setFixedSupply(double c) { fixedSupply_ = c; if (!plantOn_) rack_.setSupplyT(c); }

    void step(double dt);
    void solveSteady();

    // ---- energy (W, and J since start) ----
    double serverFanPower() const { return rack_.fanPower(); }
    /// Cooling electricity outside the IT equipment: equipment pump + facility pump + dry-cooler
    /// fans + chiller. Zero plant terms when the plant is disabled.
    double coolingPower() const;
    /// Partial PUE of the modelled boundary: (IT incl. server fans + cooling) / (IT incl. server fans).
    double partialPue() const;
    /// Heat rejected to outdoor air (dry cooler + chiller condenser), W.
    double heatToOutdoor() const;
    /// Water an evaporative tower would have evaporated to reject the same heat, L/h. The closed loop
    /// with a dry cooler evaporates none of it.
    double evaporativeEquivalentLph() const;
    double coolingEnergy() const { return coolJ_; }
    double plantFanEnergy() const { return fanJ_; }
    double chillerEnergy() const { return chJ_; }

private:
    Rack rack_;
    CoolingPlant plant_;
    bool plantOn_;
    double fixedSupply_;
    double coolJ_ = 0, fanJ_ = 0, chJ_ = 0;
};

} // namespace mf
