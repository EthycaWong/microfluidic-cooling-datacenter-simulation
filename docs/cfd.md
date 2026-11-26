# Navier–Stokes and energy equations: what the simulation solves

## Two levels of fluid physics

| Level | Where | Equations | Why |
|---|---|---|---|
| **Reduced-order** (runs the rack, plant and ML co-simulation) | `src/core/thermal`, `hydraulics`, `system` | Energy equation as a finite-volume network (silicon conduction + coolant advection, implicit, energy conserved to about 1e-8 W). Momentum as 1-D pressure balances, with Δp from laminar duct correlations (Shah & London f·Re, Hagenbach K∞) and heat transfer from Nusselt correlations. | Fast enough to run 24 chips, 8 servers and the cooling plant for a simulated day. |
| **First-principles CFD** (verifies the reduced model, and supports boundary-driven studies) | `src/core/cfd` | Full 2-D incompressible Navier–Stokes + energy with Boussinesq buoyancy. Fully developed Navier–Stokes + energy on the duct cross-section. | Shows that the correlations above are what the governing equations give for this geometry. |

The Shah & London correlations are not empirical guesses: they are tabulated solutions of the laminar Navier–Stokes and energy equations for rectangular ducts. `mfsim cfd --case duct` re-derives them for the actual 200 × 500 µm channel:

| Quantity | Navier–Stokes (64 × 160 cells) | Chip model (correlation) | Difference |
|---|---|---|---|
| Fanning f·Re | 16.359 | 16.377 | −0.1 % |
| Nu, three heated walls (H1) | 4.925 | 4.893 | +0.65 % |

## The equations (2-D solver)

```
∂u/∂t + (u·∇)u = −∇p/ρ + ν ∇²u − β (T − T_ref) g      momentum (Boussinesq)
∇·u = 0                                                mass (incompressible)
∂T/∂t + u·∇T = α ∇²T                                   energy
```

Here u = (u, v) is velocity, p pressure, T temperature, ν kinematic viscosity, α thermal diffusivity, β thermal expansion and g gravity.

**Method:**
- finite volumes on a **staggered (MAC) grid**: pressure and temperature at cell centres, u on vertical faces, v on horizontal faces;
- **Chorin projection**: an explicit momentum predictor with donor-cell / central blending for convection, then a pressure Poisson equation solved by Jacobi-preconditioned conjugate gradients, then a velocity correction, so ∇·u ≈ 0 to round-off (about 1e-7 to 1e-5 s⁻¹);
- **adaptive time step**: τ × min(CFL, viscous and thermal diffusion limits).

## Boundary conditions (all editable, live)

Each of the four walls has two independent settings.

**Velocity:**
- `Wall`: both components prescribed. The default is **u = v = 0** (no-slip, at rest). A tangential value gives a moving wall (lid-driven flow). A normal value gives an inflow, uniform or parabolic.
- `Outflow`: zero gradient, pressure fixed at 0. The outflow velocity comes from the pressure projection, so mass in = mass out exactly.

**Temperature:**
- **fixed temperature**: a hot or cold wall;
- **fixed heat flux**: q, where q = 0 is insulated.

How to change them:
- **Viewer** (CFD scene, key 5): pick Left, Right, Bottom or Top, then set *Wall u*, *Wall v*, *Hot / Cold / Insulated*, *Wall T* or *Outflow*. Changes act on the running flow.
- **Command line:** `--left-u`, `--left-v`, `--left-T`, `--left-q`, `--left-insulated`, `--left-outflow` (and the same for `right`, `bottom`, `top`).

## Ready-made cases (`mfsim cfd --case …`)

| Case | Boundaries | Benchmark / purpose | Result |
|---|---|---|---|
| `lid-cavity` | lid u = 1, other walls u = v = 0 | Ghia, Ghia & Shin (1982), Re = 100 | centreline velocity within 0.005 of the lid speed |
| `heated-cavity` | hot left (T = 1), cold right (T = 0), insulated top and bottom, all walls u = v = 0 | de Vahl Davis (1983) | Nu 2.252 vs 2.243 (Ra 1e4); 4.575 vs 4.519 (Ra 1e5), 64 × 64 grid |
| `rayleigh-benard` | hot floor, cold lid, insulated sides, all walls at rest | onset of convection near Ra ≈ 1708 | Nu = 1.00 at Ra 1e3; 1.46 at 3e3; 3.91 at 1e5 |
| `mixed-convection` | cold lid moving at u = 1 over a hot floor | shear plus buoyancy | — |
| `poiseuille` | uniform cold inflow, outflow, hot plates | analytic dp/dx = −12 μU/h², u_max = 1.5 U | pressure gradient within 0.3 % |
| `microchannel` | water at 1.67 m/s and 35 °C into the 200 µm gap between silicon ribs held at 65 °C (real water properties) | thermal entrance in the actual channel | the local Nu stays above the fully developed 7.54 over all 20 mm, as the chip model assumes |
| `duct` | fully developed flow in the 200 × 500 µm cross-section | Shah & London | see the table above |

```bash
./build/mfsim cfd --case heated-cavity --value 1e5 --n 64
```

```bash
./build/mfsim cfd --case lid-cavity --value 400 --top-u -1 --bottom-T 1 --top-T 0
```

```bash
./build/mfsim cfd --case microchannel --value 1.0 --bottom-insulated
```

Every run writes `fields.csv` (x, y, u, v, p, T) and `streamfunction.csv`. The microchannel and channel cases also write `along_flow.csv`: pressure, bulk temperature and local Nu along the flow.

## Limits

- **2-D and laminar.** The real channel is 3-D: the slice shows the gap between two ribs and leaves out the floor and the cover. The fully developed duct solve covers the 3-D cross-section, but not entrance effects together with it.
- **Explicit time stepping.** Steady-state runs on very thin channels take tens of seconds.
- **Constant properties** within a CFD run (Boussinesq): density varies only in the buoyancy term.
