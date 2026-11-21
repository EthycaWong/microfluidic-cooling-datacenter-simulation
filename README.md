# Microfluidic cooling in a data center: simulation

A C++ simulation with a 3D viewer of a liquid-cooled 42U rack. It is built up from a single chip so that each level can be checked before the next is added. The research question is whether **coolant viscosity** limits the energy savings of **silicon microchannel cooling** combined with a **closed-loop coolant system**. The two ideas follow the NVIDIA and Microsoft work.

![Rack view](docs/screenshots/rack.png)

| Server (2 microchannel GPUs + cold-plate CPU) | Chip (exploded: silicon heat map + microchannels) |
|---|---|
| ![Server view](docs/screenshots/server.png) | ![Chip view](docs/screenshots/chip.png) |

## What is modelled (this milestone)

| Level | Model |
|---|---|
| **Coolant** | Water from IAPWS-IF97, IAPWS 2008 viscosity and IAPWS 2011 conductivity, checked against the releases' test values. A viscosity-only multiplier is available for the 1×/2×/4× experiment. |
| **Chip** | 20 × 20 mm die with 50 microchannels (200 × 500 µm). A 10 × 10 thermal RC grid with backward-Euler time stepping; coolant warms along the channels. Laminar developing-flow pressure drop and heat transfer with rib fin efficiency. |
| **Server** | Two 300 W microchannel GPUs and one 200 W cold-plate CPU in parallel branches with hoses, quick disconnects and balancing valves. Branch flows are solved from pressure balance. A 200 W air-cooled group with fan power. |
| **Rack** | 8 × 4U servers between rear supply and return manifolds (direct or reverse return), a pump curve with affinity laws and an efficiency curve, and three pump modes: fixed speed, fixed flow, or PI control on the hottest chip. A placeholder CDU resistance, hydrostatics, per-server workloads and a mixed return temperature. |
| **Energy** | Heat-balance residuals at chip, server and rack level, pump and fan electricity, and cumulative energy. |

Parameter tables, boundary conditions, assumptions and the hand calculations are in **[docs/parameters.md](docs/parameters.md)**. Labelled schematics: [chip](docs/schematics/chip.svg), [server](docs/schematics/server.svg), [rack](docs/schematics/rack.svg).

**Not yet modelled:**
- the CDU heat exchanger, facility loop and dry cooler (the supply is held at 35 °C for now);
- glycol coolants;
- the CNN workload. The `Workload` interface is ready for it.

## Build

Requirements: CMake 3.24 or newer and a C++20 compiler. On macOS that means the Xcode command-line tools. On Linux, raylib also needs the X11/OpenGL development packages. The first configure downloads four dependencies, each pinned by version and SHA-256: raylib 5.5, raygui 4.0, doctest 2.4.11 and nlohmann/json 3.11.3.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

```bash
cmake --build build -j
```

```bash
ctest --test-dir build --output-on-failure
```

To build only the headless parts, add `-DMF_BUILD_VIEWER=OFF`.

## Run the 3D viewer

```bash
./build/mfsim_viewer --config configs/rack_baseline.json
```

| Action | How |
|---|---|
| Orbit / pan / zoom | left-drag / right-drag / wheel |
| Open a server | click a server, double-click to open it |
| Switch scene | keys 1, 2, 3 for rack / server / chip |
| Pause | Space |

The left panel controls:
- the per-server load, or a square-wave workload;
- coolant viscosity and supply temperature;
- the pump mode and its target;
- the manifold layout and the balancing valves.

The right panel shows live flows, pressures, temperatures, pump and fan power, and the energy residual. To render one frame to a file without interacting, add `--scene chip --screenshot out.png`.

## Run experiments (headless)

```bash
./build/mfsim info --config configs/chip_baseline.json
```

```bash
./build/mfsim steady --config configs/rack_baseline.json
```

```bash
./build/mfsim sweep-viscosity --config configs/rack_baseline.json
```

```bash
./build/mfsim resolution --config configs/chip_baseline.json
```

```bash
./build/mfsim uneven-load --config configs/rack_baseline.json
```

```bash
./build/mfsim run --config configs/rack_square_wave.json --duration 60 --dt 0.01
```

`--level chip|server|rack` changes the level of any config. Every command writes a CSV to `results/`.

## First findings

These are preliminary: the CDU and filter curves are placeholders.

- The model reproduces the proposal's hand numbers: Re ≈ 660, channel Δp ≈ 13.5 kPa, GPU outlet 43.7 °C, and a mixed rack return of 42.70 °C for 6.4 kW at 12 L/min.
- **At a given flow, viscosity does not change chip temperatures**, because Re·Pr contains no μ. It changes pumping power instead. Holding the hottest chip at 85 °C takes about 10.1 L/min at 1×, 2× and 4× viscosity, while pump power rises from 14.7 W to 21.8 W to 37.5 W.
- If the pump isn't re-adjusted, 4× viscosity cuts flow to 8.2 L/min and the hottest chip reaches 89.4 °C.
- Direct-return manifolds give the top server the least flow (2 % spread). Reverse return evens this out (0.7 %).

## Layout

```
src/core/      headless physics library (fluid, hydraulics, thermal, system, workload, config)
apps/mfsim/    experiments CLI
apps/viewer/   raylib + raygui 3D viewer
configs/       JSON configs (keys carry units)
tests/         doctest unit tests (IAPWS checks, analytic limits, energy and mass balance, controls)
docs/          parameters, schematics, screenshots
```

## License

MIT
