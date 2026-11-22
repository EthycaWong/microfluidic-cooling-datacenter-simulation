# Parameters, boundary conditions and assumptions

This page lists every input of the current model (chip → server → rack → cooling plant). Values tagged **ASSUMED** are study assumptions, not measurements; replace them as better data arrives. All of them live in [`configs/rack_baseline.json`](../configs/rack_baseline.json) (keys carry their units) and in the code defaults in [`src/core/config/Config.cpp`](../src/core/config/Config.cpp).

Schematics: [chip](schematics/chip.svg) · [server](schematics/server.svg) · [rack](schematics/rack.svg) · [cooling plant](schematics/plant.svg)

Reference numbers like [3] point to the reference list in the research proposal. Extra method references are marked A1, A2, and so on, and are listed at the end of this page.

## 1. One GPU die (silicon microchannels)

| Property | Value | Note |
|---|---|---|
| Heated silicon footprint | 20 mm × 20 mm | proposal |
| Silicon thickness before etching | 0.8 mm | proposal |
| Parallel microchannels | 50 | proposal; 50 × (200 + 200 µm) = 20.0 mm ✔ |
| Channel width / depth / length | 200 µm / 500 µm / 20 mm | proposal |
| Rib (fin) thickness | 200 µm | proposal |
| Floor under channels | 300 µm | = 0.8 mm − 500 µm ✔ |
| Hydraulic diameter | 285.7 µm | 2wd/(w+d) |
| Thermal grid | 10 × 10 (i along flow) | resolution-checked, see §6 |
| Silicon conductivity | 148 W/m·K × (300 K/T)^1.3 | temperature dependent |
| Silicon ρc_p | 2329 × 705 J/m³·K | |
| GPU power | P(u) = 60 + u·(300 − 60) W | **ASSUMED** idle 60 W; linear in activity u, needs calibration |
| GPU power map | central block, 50 % × 50 % of area, carries 50 % of power | **ASSUMED** (flux 2× average in the core) |

**Hydraulic model** ([1], A1). The fully developed friction uses Shah & London f·Re(α), which gives 16.38 (Fanning) at α = 0.4. Two further losses are added: the Hagenbach developing-flow excess K(∞) = 1.28, and header contraction and expansion losses of 0.5 + 1.0 velocity heads.

**Heat transfer.** The fully developed Nusselt number for three heated walls is 4.89, with the bonded cover adiabatic (A1). The thermal-entrance enhancement uses the Shah & London local solution for thermally developing flow (H boundary condition), applied as a ratio to its fully developed value. This is a documented approximation for rectangular ducts. Rib fin efficiency is η = tanh(mH)/(mH).

## 2. One CPU (conventional cold plate)

| Property | Value | Note |
|---|---|---|
| Die | 25 mm × 25 mm × 0.8 mm | **ASSUMED** |
| Power | P(u) = 50 + u·(200 − 50) W | 200 W full from proposal, idle **ASSUMED** |
| Power map | central 60 % × 60 % block carries 60 % | **ASSUMED** |
| TIM | 8 × 10⁻⁶ m²K/W (0.08 cm²K/W) | **ASSUMED**, typical |
| Copper base | 3 mm, k = 390 W/m·K | |
| Fin-side conductance | hA = 25 W/K × (Q / 0.5 L/min)^0.5 | **ASSUMED** |
| Pressure drop | 4 kPa at 0.5 L/min, 60 % laminar (∝ μ) | **ASSUMED** |

## 3. Air-cooled group (memory, storage, other electronics)

| Property | Value | Note |
|---|---|---|
| Power | 120–200 W (200 W at full load) | 200 W from proposal, idle **ASSUMED** |
| Heat path | to room air at 27 °C, never to the coolant | proposal: has its own heat-rejection path |
| Fan control | hold a 12 K air temperature rise | **ASSUMED** |
| Fan power | 150 Pa at 14.3 L/s, η 0.35, P_fan ∝ V̇³ | **ASSUMED** |
| Thermal mass | 2000 J/K, UA 10 W/K at reference airflow | **ASSUMED**; time constant ≈ 200 s |

The rack's air-cooled load is 8 × 200 W = 1.6 kW. Until its heat-rejection path is modelled, water-use statements apply only to the liquid loop.

## 4. Server (one 4U)

```
Server supply ─ QD ─┬─ hose ─ GPU 1 microchannels ─ valve ─┬─ QD ─ Server return
                    ├─ hose ─ GPU 2 microchannels ─ valve ─┤
                    └─ hose ─ CPU cold plate      ─ valve ─┘
```

| Property | Value | Note |
|---|---|---|
| Branch target | 0.5 L/min each, 1.5 L/min per server | proposal |
| Quick disconnects | 3 kPa each at 1.5 L/min, quadratic | **ASSUMED**, typical |
| Branch hoses | 6 mm ID, 0.6 m, K = 3 | **ASSUMED** |
| Balancing valves | sized so that every branch carries 0.5 L/min at 35 °C | GPU valves 0 kPa, CPU valve ≈ 9.5 kPa |
| Server rig boundary | 35 °C supply, 1.5 L/min total (fixed) | temporary, replaced by the rack pump |

Branch flows are **solved** from equal branch pressure drop ([6]). Without valves the CPU branch takes more than 1.3× the flow of a GPU branch (unit test).

## 5. Rack and loop

| Property | Value | Note |
|---|---|---|
| Rack | 42U; 8 × 4U servers in U7–U38; CDU + pump U1–U6; network U39–40; spare U41–42 | proposal: 10U for cooling/network/spare |
| Manifolds | 20 mm ID, vertical at the rear; Churchill friction (A2) + K = 0.5 per tap | **ASSUMED** size |
| Layout | direct return (default) or reverse return | |
| CDU heat exchanger, equipment-side Δp | 25 kPa at 12 L/min, 50 % laminar | **ASSUMED** (thermal model in §5b) |
| Filter | 8 kPa at 12 L/min, 30 % laminar | **ASSUMED** |
| Expansion vessel | fixes 150 kPa at the pump suction | reference pressure |
| Pump | H = s²·100 kPa − bQ², H = 0 at 20 L/min (s = 1); η_max 0.5 at 12 L/min; motor 0.9 | **ASSUMED** small in-rack pump |
| Pump control | fixed speed / fixed flow (12 L/min) / PI on hottest chip (target 85 °C) | |
| Supply temperature | `rack` level: fixed 35 °C. `plant` level: computed by the CDU model (§5b) | [8] |
| Hydrostatics | one loop density; ρgz cancels around the closed loop | [14]; unit test |
| Electrical pump power | P = Δp·V̇ / (η_pump·η_motor) | [10] |

## 5b. Cooling plant (CDU, facility loop, dry cooler, trim chiller)

Run it with `--level plant` or [`configs/plant_baseline.json`](../configs/plant_baseline.json). The equipment loop and the facility loop are sealed and never mix ([7]). Heat leaves through a dry cooler, so nothing evaporates in operation ([8]).

| Component | Value | Note |
|---|---|---|
| CDU plate HX | counterflow ε-NTU; each side hA = 2800 W/K at design (12 / 14 L/min) | **ASSUMED** sizing: UA ≈ 1440 W/K, effectiveness 0.66, ~4 K approach at 6.4 kW |
| Side conductance scaling | hA ∝ (ṁ/ṁ_ref)^0.6 (μ/μ_ref)^(1/3−0.6) (k/k_ref)^(2/3) (c_p/c_p,ref)^(1/3) | Dittus–Boelter form, so **thicker coolant lowers hA** |
| Facility loop | water, 14 L/min fixed; Δp = 20 (CDU) + 30 (coil) + 10 (piping) kPa; pump η 0.5 × motor 0.9 | **ASSUMED**; ≈ 31 W |
| Dry cooler | crossflow ε-NTU (air mixed); UA 2000 W/K at full fans; 75 % of resistance on the air side, hA_air ∝ V̇^0.6; 2.0 m³/s and 450 W at full speed, P ∝ speed³ | **ASSUMED**; full fans hold 35 °C supply up to ≈ 29 °C outdoor |
| Fan control | PI on equipment supply temperature, Kp 0.08 /K, Ki 0.0015 /(K·s) | tuned on a 30 % → 100 % load step (settles within ~0.03 K in 10 min) |
| Dry-cooler bypass | when outdoor air ≥ facility water, or when (while chilling) its fans cost more than the chiller power they save | economizer rule |
| Trim chiller | only with fans at 100 % and supply above setpoint; COP = 0.5 × Carnot (T_evap = target − 5 K, T_cond = outdoor + 10 K), ≤ 10; 10 kW capacity | **ASSUMED** |
| Thermal mass | 6 L equipment supply side, 8 L return side, 25 L facility loop | **ASSUMED** |
| Outdoor air | constant or daily sine (mean, amplitude, period) | |
| Supply setpoint | 35 °C | proposal |
| Reported | cooling electricity = rack pump + facility pump + dry-cooler fans + chiller; partial PUE = (IT + server fans + cooling)/(IT + server fans); water evaporated (0) and evaporative-tower equivalent (heat / 2.43 MJ/kg) | server fans count as IT, following PUE convention |

**Climate sweep** at full load (`mfsim climate-sweep --config configs/plant_baseline.json`):

| Outdoor °C | Dry-cooler fans | Chiller electric | Cooling electric | Partial PUE | Supply °C |
|---|---|---|---|---|---|
| 0 | 10 % (min) | 0 W | 56 W | 1.007 | 29.4 (below setpoint) |
| 20 | 23 % | 0 W | 60 W | 1.007 | 35.0 |
| 25 | 41 % | 0 W | 85 W | 1.011 | 35.0 |
| 28 | 79 % | 0 W | 277 W | 1.034 | 35.0 |
| 30 | 100 % | 94 W | 599 W | 1.074 | 35.0 |
| 35 | bypassed | 809 W | 864 W | 1.107 | 35.0 |
| 40 | bypassed | 1023 W | 1077 W | 1.134 | 35.0 |

With the chiller off, the supply temperature floats up instead: 41 °C at 35 °C outdoor, where the hottest chip reaches 87.8 °C. An evaporative tower rejecting the same 6.4 kW would evaporate about 9.5 L/h; the dry-cooler loop evaporates none. The cost of that water saving is the fan and chiller electricity on hot days ([9]).

## 6. Hand calculations and what the model gives

| Quantity | Hand calculation | Model |
|---|---|---|
| Channel velocity at 0.5 L/min | 1.667 m/s | 1.667 m/s |
| Reynolds number (35 °C) | ≈ 660 | 658 |
| Hydrodynamic / thermal entrance length | ≈ 9 mm / ≈ 45 mm | 9.4 mm / 45.4 mm (still developing over the 20 mm channel) |
| Channel Δp | ≈ 9.6 kPa friction, 12–15 kPa total | 9.6 + 1.8 + 2.1 = 13.5 kPa |
| Hydraulic power through one GPU's channels | ≈ 0.1 W | 0.11 W |
| GPU outlet at 300 W | 35 + 8.7 = 43.7 °C | 43.67 °C |
| Peak silicon, uniform 300 W | 65–70 °C | 65.9 °C |
| Peak silicon, assumed hotspot | — | 81.9 °C |
| Rack mixed return, 6.4 kW at 12 L/min | 35 + 7.7 = 42.7 °C | 42.70 °C |
| Silicon thermal time constant | ≈ 30 ms | implicit solver used, so no stability limit |

**Grid resolution check** for the hotspot map at 0.5 L/min (`mfsim resolution`):

| Grid | 5×5 | 10×10 | 20×20 | 40×40 |
|---|---|---|---|---|
| Peak °C | 80.12 | 81.85 | 82.64 | 82.49 |

The 10 × 10 grid is within about 0.7 K of the refined result. Use 20 × 20 when a comparison depends on differences below about 1 K.

## 7. Viscosity results

### 7a. Rack level (fixed 35 °C supply)

The experiment changes **viscosity only** (1×, 2×, 4× water); density, c_p and k stay those of water. All servers run at full load with the hotspot map. Balancing valves are re-commissioned for each coolant. Run it with `mfsim sweep-viscosity --config configs/rack_baseline.json`.

| μ | Case | Rack flow | Pump head | Pump electric | Peak °C |
|---|---|---|---|---|---|
| 1× | hold peak at 85 °C | 10.07 L/min | 39.4 kPa | 14.7 W | 85.00 |
| 2× | hold peak at 85 °C | 10.06 L/min | 58.0 kPa | 21.8 W | 85.00 |
| 4× | hold peak at 85 °C | 10.06 L/min | 95.1 kPa | 37.5 W | 85.00 |
| 4× | pump speed fixed at the 1× value | 8.15 L/min | 71.8 kPa | 23.5 W | 89.38 |

What this shows so far:

1. **At a given flow, viscosity does not change chip temperature.** The group Re·Pr = ρuD_h·c_p/k contains no μ, so the thermal entrance and h are unchanged. The flow needed to hold 85 °C is the same at every viscosity (unit-tested).
2. **Viscosity costs pumping energy instead.** Holding 85 °C at 4× viscosity needs 2.55× the pump power, because the laminar parts of the loop (microchannels, cold plates, part of the HX and filter) scale with μ.
3. **If the pump is not re-adjusted**, the thicker coolant reduces flow and the hottest chip exceeds the limit (89.4 °C at 4×).
4. **In absolute terms the pump is small** compared with the 8 kW IT load (0.18 % → 0.47 %) and with the server fans (≈ 49 W).

### 7b. Plant level (supply from the CDU, outdoor 25 °C and 30 °C)

`mfsim sweep-viscosity --config configs/plant_baseline.json [--outdoor 30]`, holding the hottest chip at 85 °C:

| μ | CDU UA | Rack pump | Dry-cooler fans | Chiller | **Cooling total** | Partial PUE |
|---|---|---|---|---|---|---|
| 1×, 25 °C | 1369 W/K | 14.7 W | 25.4 W | 0 W | **71.2 W** | 1.0088 |
| 2×, 25 °C | 1235 W/K | 21.8 W | 32.0 W | 0 W | **84.9 W** | 1.0106 |
| 4×, 25 °C | 1104 W/K | 37.5 W | 43.6 W | 0 W | **112.2 W** | 1.0139 |
| 1×, 30 °C | 1369 W/K | 14.7 W | 450 W | 64.8 W | **560.6 W** | 1.0697 |
| 4×, 30 °C | 1104 W/K | 37.5 W | 450 W | 142.5 W | **661.1 W** | 1.0821 |

The plant level adds a second viscosity penalty that the rack-only model could not show. A thicker coolant lowers the CDU's equipment-side conductance (UA −19 % at 4×), so the facility loop must supply colder water to hold the same 35 °C supply. On a mild day that shows up as extra dry-cooler fan power. On a warm day it shows up as extra chiller power, which is larger. At 4× viscosity the extra cooling electricity is +41 W at 25 °C (+23 W pump, +18 W fans) and +100 W at 30 °C (+23 W pump, +78 W chiller).

## 8. Not in this milestone

- Glycol in the facility loop (freeze protection), and room HVAC for the 1.6 kW air-cooled load. Water-use statements cover only the liquid loop.
- Part-load optimisation of fans plus chiller together (the model uses either full fans with trim chilling, or a full bypass). Humidity, and pump and fan curves beyond simple affinity laws.
- Expansion-vessel dynamics. Glycol coolants need mixture-specific property data; the IAPWS pure-water formulations are not substituted ([11], [12]).
- Qu & Mudawar benchmark ([3]). The geometry can be configured, but the comparison needs data points taken from the paper.
- The CNN workload driving u(t) ([13]); the `Workload` interface is ready for it.
- Header flow maldistribution across the 50 channels (flow is assumed equal per channel), the simultaneously developing Nusselt number, and thermosiphon effects.

## Method references (in addition to the proposal's list)

- **A1** R. K. Shah and A. L. London, *Laminar Flow Forced Convection in Ducts*, Academic Press, 1978. Rectangular-duct f·Re, K(∞), three-wall Nu and thermal-entrance solutions.
- **A2** S. W. Churchill, "Friction-factor equation spans all fluid-flow regimes," *Chemical Engineering* 84 (24), 1977.
- **A3** IAPWS R7-97(2012) (IF97), IAPWS R12-08 (viscosity), IAPWS R15-11 (thermal conductivity). Implemented and checked against the releases' verification tables (these are the proposal's [11] and [12]).
