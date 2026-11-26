#pragma once
#include <string>
#include <vector>

namespace mf::cfd {

/// Velocity boundary condition on one side of the rectangle.
///  Wall:    both velocity components prescribed (u, v); default 0, 0 = no-slip, stationary wall.
///           A tangential value makes a moving wall (lid-driven flow); a normal value makes an inflow.
///  Outflow: zero normal gradient for velocity and temperature, pressure fixed at 0.
enum class VelocityBC { Wall, Outflow };
/// Temperature boundary condition: fixed temperature (hot / cold wall) or fixed heat flux (0 = insulated).
enum class ThermalBC { FixedTemperature, HeatFlux };

struct Side {
    VelocityBC velocity = VelocityBC::Wall;
    double u = 0.0, v = 0.0;      ///< prescribed velocity components on the boundary, m/s
    bool parabolic = false;       ///< shape a normal (inflow) velocity as a parabola with this mean
    ThermalBC thermal = ThermalBC::HeatFlux;
    double T = 0.0;               ///< fixed temperature (FixedTemperature), C or K
    double q = 0.0;               ///< heat flux INTO the fluid, W/m^2 (HeatFlux; 0 = adiabatic)
};

/// Physical and numerical set-up. SI units; non-dimensional cases simply use L = 1 etc.
struct Params {
    int nx = 64, ny = 64;
    double lx = 1.0, ly = 1.0;      ///< domain size, m
    double nu = 0.01;               ///< kinematic viscosity, m^2/s
    double alpha = 0.01;            ///< thermal diffusivity, m^2/s
    double k = 1.0;                 ///< thermal conductivity (only for heat-flux walls and Nu), W/(m K)
    double gx = 0.0, gy = 0.0;      ///< gravity, m/s^2 (e.g. gy = -9.81)
    double beta = 0.0;              ///< thermal expansion coefficient (Boussinesq), 1/K
    double tRef = 0.0;              ///< Boussinesq reference temperature
    double t0 = 0.0;                ///< initial temperature
    double tau = 0.5;               ///< time-step safety factor on the CFL / diffusion limits
    double gammaMin = 0.0;          ///< minimum donor-cell (upwind) blending, 0..1
    double maxDt = 1e30;
    double pressureTol = 1e-7;      ///< relative residual for the pressure Poisson solve
    Side left, right, bottom, top;
};

/// Unsteady incompressible Navier-Stokes + energy equations in 2-D (Boussinesq buoyancy):
///   du/dt + (u.grad)u = -grad(p)/rho + nu lap(u) - beta (T - Tref) g
///   div(u) = 0
///   dT/dt + u.grad(T) = alpha lap(T)
/// Finite volumes on a staggered (MAC) grid; Chorin projection: explicit momentum and energy steps
/// with a donor-cell / central blend for convection, then a pressure Poisson equation solved with
/// Jacobi-preconditioned conjugate gradients so the new velocity is divergence free.
/// Cells are numbered 1..nx, 1..ny with one ghost layer; u(i,j) sits on the right face of cell
/// (i,j) and v(i,j) on its top face (the layout of Griebel, Dornseifer & Neunhoeffer).
class NavierStokes2D {
public:
    explicit NavierStokes2D(const Params& p);

    /// Advance one adaptive time step; returns the step used.
    double step();
    /// Step until simulated time reaches t (or maxSteps).
    void runUntil(double t, int maxSteps = 10000000);
    /// Step until the largest relative change rate of u, v, T falls below tol (per unit of the
    /// reference time scale tScale), or maxSteps; returns true when converged.
    bool runToSteady(double tol, double tScale, int maxSteps = 2000000);
    void reset();

    Params& params() { return p_; }
    const Params& params() const { return p_; }
    /// Change boundary conditions or properties on the fly (live in the viewer).
    void setSide(int which, const Side& s); ///< 0 left, 1 right, 2 bottom, 3 top

    // ---- fields (cell centres, i in 1..nx, j in 1..ny) ----
    int nx() const { return p_.nx; }
    int ny() const { return p_.ny; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double uc(int i, int j) const { return 0.5 * (U(i, j) + U(i - 1, j)); }
    double vc(int i, int j) const { return 0.5 * (V(i, j) + V(i, j - 1)); }
    double pc(int i, int j) const { return P(i, j); }    ///< kinematic pressure p/rho, m^2/s^2
    double tc(int i, int j) const { return Tm(i, j); }
    double uFace(int i, int j) const { return U(i, j); } ///< i = 0..nx
    double vFace(int i, int j) const { return V(i, j); } ///< j = 0..ny
    /// Velocity at an arbitrary point (bilinear on the staggered grid), for tracer particles.
    void velocityAt(double x, double y, double& u, double& v) const;
    double temperatureAt(double x, double y) const;
    /// Stream function psi(x_i, y_j) at cell corners (i = 0..nx, j = 0..ny), psi = 0 on the bottom-left.
    std::vector<double> streamFunction() const;

    // ---- diagnostics ----
    double time() const { return time_; }
    long steps() const { return steps_; }
    double lastDt() const { return dt_; }
    double maxSpeed() const;
    double maxDivergence() const;          ///< 1/s, after the projection (should be ~0)
    double kineticEnergy() const;          ///< 0.5 * integral |u|^2 dA, m^4/s^2
    int lastPressureIterations() const { return cgIter_; }
    /// Mean heat flux INTO the fluid through a side, W/m^2 (uses k), second-order one-sided gradient.
    double wallHeatFlux(int which) const;
    /// Bulk (velocity-weighted) temperature over column i (flow in x).
    double bulkTemperatureColumn(int i) const;
    /// Mean kinematic pressure of column i.
    double meanPressureColumn(int i) const;

private:
    double& U(int i, int j) { return u_[i * (p_.ny + 2) + j]; }
    double U(int i, int j) const { return u_[i * (p_.ny + 2) + j]; }
    double& V(int i, int j) { return v_[i * (p_.ny + 1) + j]; }
    double V(int i, int j) const { return v_[i * (p_.ny + 1) + j]; }
    double& P(int i, int j) { return pr_[i * (p_.ny + 2) + j]; }
    double P(int i, int j) const { return pr_[i * (p_.ny + 2) + j]; }
    double& Tm(int i, int j) { return t_[i * (p_.ny + 2) + j]; }
    double Tm(int i, int j) const { return t_[i * (p_.ny + 2) + j]; }

    void applyVelocityBC();
    void applyTemperatureBC();
    void applyPressureBC(std::vector<double>& p) const;
    double computeDt() const;
    void solvePressure(const std::vector<double>& rhs);
    void laplacian(const std::vector<double>& x, std::vector<double>& out) const; ///< -lap with BCs
    double inflowProfile(const Side& s, double frac) const;
    bool anyOutflow() const;

    Params p_;
    double dx_, dy_;
    std::vector<double> u_, v_, pr_, t_, f_, g_;
    double time_ = 0.0, dt_ = 0.0;
    long steps_ = 0;
    int cgIter_ = 0;
};

} // namespace mf::cfd
