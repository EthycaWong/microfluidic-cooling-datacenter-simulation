#pragma once
#include "fluid/Fluid.hpp"

namespace mf {

/// Geometry of an array of identical parallel rectangular microchannels etched into silicon.
struct ChannelGeometry {
    double width = 200e-6;   ///< channel width (across the chip), m
    double depth = 500e-6;   ///< etched channel depth, m
    double length = 20e-3;   ///< channel length (flow direction), m
    double rib = 200e-6;     ///< silicon rib (fin) thickness between channels, m
    int count = 50;          ///< number of parallel channels
    double kInlet = 0.5;     ///< header-to-channel contraction loss, velocity heads (channel velocity)
    double kOutlet = 1.0;    ///< channel-to-header expansion loss, velocity heads
    bool hagenbach = true;   ///< include the developing-flow (Hagenbach) excess pressure drop K(inf)

    double area() const { return width * depth; }                               ///< one channel, m^2
    double hydraulicDiameter() const { return 2.0 * width * depth / (width + depth); }
    double aspect() const { return width < depth ? width / depth : depth / width; } ///< short/long, 0..1
};

/// Per-state results of the channel model, useful for diagnostics and the UI.
struct ChannelState {
    double velocity;   ///< mean channel velocity, m/s
    double re;         ///< Reynolds number
    double dpFriction; ///< fully developed friction part, Pa (scales with mu * Q)
    double dpDeveloping; ///< Hagenbach excess, Pa
    double dpMinor;    ///< header contraction + expansion, Pa
    double dpTotal;    ///< Pa
    double lhPlus;     ///< L / (Dh Re), dimensionless hydrodynamic length
    double entranceLengthHydro;   ///< ~0.05 Re Dh, m
    double entranceLengthThermal; ///< ~0.05 Re Pr Dh, m
};

/// Laminar single-phase pressure drop and heat transfer for rectangular microchannels.
///
/// Friction:  Shah & London fully developed f*Re(alpha) for rectangular ducts, plus the Hagenbach
///            excess pressure drop K(inf)(alpha) for hydrodynamically developing flow, plus header
///            contraction/expansion minor losses.
/// Heat:      fully developed Nusselt number for three heated walls (top wall adiabatic, bonded
///            cover), from Shah & London; thermal-entrance enhancement from the Shah & London local
///            Nusselt solution for thermally developing flow (H boundary condition) applied as a
///            ratio to its fully developed value. This is an approximation for rectangular ducts.
class RectChannel {
public:
    explicit RectChannel(const ChannelGeometry& g) : g_(g) {}
    const ChannelGeometry& geometry() const { return g_; }

    /// Fanning friction factor times Reynolds number, fully developed laminar flow.
    static double fanningFRe(double aspect);
    /// Hagenbach factor K(inf), velocity heads, for fully developed velocity profile at the exit.
    static double hagenbachK(double aspect);
    /// Fully developed Nu for three heated walls, one short wall adiabatic. alpha = width/depth.
    static double nu3WallFullyDeveloped(double alphaWidthOverDepth);
    /// Local Nu for thermally developing laminar flow in a circular duct, H boundary condition.
    static double nuCircularDevelopingH(double xStar);

    /// Pressure drop across the whole channel array at total volumetric flow q (m^3/s).
    ChannelState state(double q, const FluidProps& p) const;
    double pressureDrop(double q, const FluidProps& p) const { return state(q, p).dpTotal; }

    /// Local heat transfer coefficient (W/m^2K) at distance x from the channel inlet.
    double localH(double x, double qTotal, const FluidProps& p) const;

private:
    ChannelGeometry g_;
};

} // namespace mf
