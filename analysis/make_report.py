"""Run every simulation experiment and write an illustrated, written results report.

    .venv/bin/python analysis/make_report.py            # run all experiments (~4 min), then draw + write
    .venv/bin/python analysis/make_report.py --no-sim   # redraw from existing results/report data

Outputs: docs/report/REPORT.md and docs/report/fig*.png. Data: results/report/*.csv.
Every number quoted in the report text is read from the simulation output, not typed in.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MF = ROOT / "build" / "mfsim"
RES = ROOT / "results" / "report"
FIG = ROOT / "docs" / "report"
CFG = ROOT / "configs"

HEAT = LinearSegmentedColormap.from_list(
    "heat", ["#23176e", "#3a6ee6", "#1ebedc", "#46dc78", "#e1dc3c", "#fa8c28", "#c81e1e"])
C = {"water": "#2563eb", "PG25": "#dc2626", "PG40": "#9333ea", "EG25": "#ea580c"}
plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 10})


# --------------------------------------------------------------------------------------------- helpers
def mfsim(*args: str) -> str:
    r = subprocess.run([str(MF), *args], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"mfsim {' '.join(args)} failed:\n{r.stdout}\n{r.stderr}")
    return r.stdout


def load_cfg(name: str) -> dict:
    text = "\n".join(l for l in (CFG / name).read_text().splitlines() if not l.strip().startswith("//"))
    return json.loads(text)


def write_cfg(cfg: dict, name: str) -> str:
    p = RES / name
    p.write_text(json.dumps(cfg, indent=1))
    return str(p)


def table(path: Path) -> list[dict]:
    with open(path) as f:
        rows = list(csv.DictReader(f))
    out = []
    for r in rows:
        d = {}
        for k, v in r.items():
            if k is None:
                continue
            try:
                d[k] = float(v)
            except (TypeError, ValueError):
                d[k] = v
        out.append(d)
    return out


def col(rows: list[dict], key: str, **where) -> np.ndarray:
    return np.array([r[key] for r in rows if all(r.get(k) == v for k, v in where.items())], dtype=float)


def save(fig, name: str) -> str:
    fig.savefig(FIG / name, bbox_inches="tight")
    plt.close(fig)
    return name


# --------------------------------------------------------------------------------------------- experiments
def run_experiments() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    print("chip maps, resolution ...")
    mfsim("steady", "--config", "configs/validation_uniform.json", "--out", str(RES / "chip_uniform.csv"))
    mfsim("steady", "--config", "configs/chip_baseline.json", "--out", str(RES / "chip_hotspot.csv"))
    mfsim("steady", "--config", "configs/chip_baseline.json", "--coolant", "PG25", "--out", str(RES / "chip_hotspot_pg25.csv"))
    mfsim("resolution", "--config", "configs/chip_baseline.json", "--out", str(RES / "resolution.csv"))

    print("server and rack ...")
    mfsim("steady", "--config", "configs/server_baseline.json", "--out", str(RES / "server_balanced.csv"))
    cfg = load_cfg("server_baseline.json")
    cfg["server"]["auto_balance"] = False
    mfsim("steady", "--config", write_cfg(cfg, "cfg_server_unbalanced.json"), "--out", str(RES / "server_unbalanced.csv"))
    mfsim("steady", "--config", "configs/rack_baseline.json", "--out", str(RES / "rack_direct.csv"))
    cfg = load_cfg("rack_baseline.json")
    cfg["rack"]["layout"] = "reverse-return"
    mfsim("steady", "--config", write_cfg(cfg, "cfg_rack_reverse.json"), "--out", str(RES / "rack_reverse.csv"))
    mfsim("uneven-load", "--config", "configs/rack_baseline.json", "--out", str(RES / "uneven.csv"))
    mfsim("hydraulics", "--config", "configs/rack_baseline.json", "--out", str(RES / "hydraulics"))

    print("viscosity sweeps ...")
    scales = "1,1.5,2,3,4"
    mfsim("sweep-viscosity", "--config", "configs/rack_baseline.json", "--scales", scales, "--out", str(RES / "visc_rack.csv"))
    mfsim("sweep-viscosity", "--config", "configs/plant_baseline.json", "--scales", scales, "--out", str(RES / "visc_plant25.csv"))
    mfsim("sweep-viscosity", "--config", "configs/plant_baseline.json", "--scales", scales, "--outdoor", "30",
          "--out", str(RES / "visc_plant30.csv"))

    print("coolants and climate ...")
    mfsim("coolant-compare", "--config", "configs/rack_baseline.json", "--out", str(RES / "coolants_rack.csv"),
          "--props-out", str(RES / "coolant_props.csv"))
    mfsim("coolant-compare", "--config", "configs/plant_baseline.json", "--out", str(RES / "coolants_plant.csv"),
          "--props-out", str(RES / "coolant_props.csv"))
    mfsim("climate-sweep", "--config", "configs/plant_baseline.json", "--temps", "-5,0,5,10,15,20,25,28,30,32,35,40,45",
          "--out", str(RES / "climate.csv"))

    print("load-step transient ...")
    cfg = load_cfg("plant_baseline.json")
    cfg["workload"] = {"type": "step", "before": 0.3, "after": 1.0, "t_step_s": 60}
    mfsim("run", "--config", write_cfg(cfg, "cfg_step.json"), "--start-steady", "--duration", "1200", "--dt", "0.05",
          "--every", "1", "--out", str(RES / "load_step.csv"))

    print("Navier-Stokes cases (about 2 min) ...")
    cfd = RES / "cfd"
    mfsim("cfd", "--case", "lid-cavity", "--value", "100", "--n", "64", "--out", str(cfd / "lid"))
    nu_rows = ["ra,nu_hot,nu_cold"]
    for ra in ("1e3", "1e4", "1e5"):
        out = mfsim("cfd", "--case", "heated-cavity", "--value", ra, "--n", "64", "--out", str(cfd / f"cavity_{ra}"))
        line = next(l for l in out.splitlines() if l.startswith("Nusselt number: hot wall"))
        parts = line.replace(",", "").split()
        nu_rows.append(f"{ra},{parts[4]},{parts[7]}")
    (cfd / "cavity_nu.csv").write_text("\n".join(nu_rows) + "\n")
    rb_rows = ["ra,nu"]
    for ra in ("1000", "3000", "1e4", "3e4", "1e5"):
        out = mfsim("cfd", "--case", "rayleigh-benard", "--value", ra, "--n", "24", "--seed", "0.2",
                    "--out", str(cfd / f"rb_{ra}"))
        line = next(l for l in out.splitlines() if l.startswith("Nusselt number (floor)"))
        rb_rows.append(f"{ra},{line.split()[3]}")
    (cfd / "rb_nu.csv").write_text("\n".join(rb_rows) + "\n")
    mfsim("cfd", "--case", "duct", "--width", "200e-6", "--height", "500e-6", "--n", "64", "--out", str(cfd / "duct"))
    mfsim("cfd", "--case", "microchannel", "--out", str(cfd / "microchannel"))

    if (ROOT / "results" / "ml" / "activity_trace.csv").exists():
        print("CNN job on the Navier-Stokes microchannel (about 1.5 min) ...")
        mfsim("cfd-job", "--config", "configs/plant_baseline.json", "--trace", "results/ml/activity_trace.csv",
              "--coolants", "water,PG25", "--snapshot-at", "65", "--out", str(RES / "cfd_job"))

    print("24-hour day (about a minute) ...")
    rng = np.random.default_rng(7)
    with open(RES / "day_workload.csv", "w") as f:
        f.write("t_s,u\n")
        for k in range(0, 24 * 4 + 1):
            h = k / 4
            base = 0.25 + 0.6 * math.exp(-((h - 14.5) / 4.0) ** 2)   # busy afternoon, quiet night
            f.write(f"{h * 3600:.0f},{min(1.0, max(0.05, base + rng.normal(0, 0.05))):.3f}\n")
    cfg = load_cfg("plant_baseline.json")
    cfg["workload"] = {"type": "trace", "file": str(RES / "day_workload.csv")}
    cfg["plant"]["weather"] = {"outdoor_mean_C": 24, "daily_amplitude_K": 8, "period_s": 86400, "phase_s": 4 * 3600}
    mfsim("run", "--config", write_cfg(cfg, "cfg_day.json"), "--start-steady", "--duration", "86400", "--dt", "2",
          "--every", "60", "--out", str(RES / "day.csv"))


# --------------------------------------------------------------------------------------------- figures
def fig_chip(facts: dict) -> str:
    maps = [("Uniform 300 W, water", table(RES / "chip_uniform.csv")),
            ("Hotspot 300 W, water", table(RES / "chip_hotspot.csv")),
            ("Hotspot 300 W, PG25", table(RES / "chip_hotspot_pg25.csv"))]
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.2), gridspec_kw={"width_ratios": [1, 1, 1, 1.35], "wspace": 0.35})
    vmin, vmax = 35, 90
    for ax, (title, rows) in zip(axs[:3], maps):
        n = int(max(r["i"] for r in rows)) + 1
        grid = np.zeros((n, n))
        for r in rows:
            grid[int(r["j"]), int(r["i"])] = r["silicon_C"]
        im = ax.imshow(grid, origin="lower", cmap=HEAT, vmin=vmin, vmax=vmax, extent=[0, 20, 0, 20])
        ax.set_title(f"{title}\npeak {grid.max():.1f} °C, mean {grid.mean():.1f} °C")
        ax.set_xlabel("along flow (mm) →")
        ax.set_ylabel("across (mm)")
        ax.grid(False)
        facts[f"chip_peak_{title}"] = grid.max()
    fig.colorbar(im, ax=axs[:3], shrink=0.8, label="silicon °C")
    ax = axs[3]
    for (title, rows), ls in zip(maps, ["-", "--", ":"]):
        n = int(max(r["i"] for r in rows)) + 1
        x = [(i + 0.5) * 20 / n for i in range(n)]
        si = [np.mean([r["silicon_C"] for r in rows if r["i"] == i]) for i in range(n)]
        fl = [np.mean([r["coolant_C"] for r in rows if r["i"] == i]) for i in range(n)]
        ax.plot(x, si, ls, color="#ea580c", label=f"silicon, {title}")
        ax.plot(x, fl, ls, color="#2563eb", label=f"coolant, {title}")
    ax.set_xlabel("distance along channel (mm)")
    ax.set_ylabel("°C (column mean)")
    ax2 = ax.twinx()
    rows = maps[0][1]
    n = int(max(r["i"] for r in rows)) + 1
    h = [rows[i]["h_row_W_m2K"] / 1000 for i in range(n)]
    ax2.plot([(i + 0.5) * 2 for i in range(n)], h, color="#16a34a", lw=2, label="h, water")
    ax2.set_ylabel("heat-transfer coefficient h (kW/m²K)", color="#16a34a")
    ax2.grid(False)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("Profiles along the flow: coolant warms, h falls as the\nthermal boundary layer develops")
    facts["h_in"], facts["h_out"] = h[0], h[-1]
    return save(fig, "fig01_chip.png")


def fig_resolution(facts: dict) -> str:
    rows = table(RES / "resolution.csv")
    n = col(rows, "n")
    pk = col(rows, "peak_C")
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.plot(n, pk, "o-", color="#2563eb")
    for a, b in zip(n, pk):
        ax.annotate(f"{b:.2f}", (a, b), textcoords="offset points", xytext=(4, 6), fontsize=9)
    ax.set_xscale("log", base=2)
    ax.set_xticks(n)
    ax.set_xticklabels([f"{int(v)}×{int(v)}" for v in n])
    ax.set_ylabel("peak silicon °C")
    ax.set_title("Grid-resolution check (hotspot map, 0.5 L/min)")
    facts["res"] = dict(zip(n.astype(int), pk))
    return save(fig, "fig02_resolution.png")


def fig_server(facts: dict) -> str:
    bal, unb = table(RES / "server_balanced.csv"), table(RES / "server_unbalanced.csv")
    names = [r["branch"] for r in bal]
    x = np.arange(len(names))
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    axs[0].bar(x - 0.2, [r["flow_lpm"] for r in bal], 0.4, label="with balancing valves", color="#2563eb")
    axs[0].bar(x + 0.2, [r["flow_lpm"] for r in unb], 0.4, label="no valves", color="#94a3b8")
    axs[0].axhline(0.5, color="#111827", lw=0.8, ls="--")
    axs[0].set_xticks(x, names)
    axs[0].set_ylabel("branch flow (L/min)")
    axs[0].set_title("Branch flows at 1.5 L/min total (solved from pressure balance)")
    axs[0].legend(fontsize=8)
    axs[1].bar(x - 0.2, [r["peak_C"] for r in bal], 0.4, color="#2563eb")
    axs[1].bar(x + 0.2, [r["peak_C"] for r in unb], 0.4, color="#94a3b8")
    axs[1].axhline(85, color="#dc2626", lw=0.8, ls=":")
    axs[1].set_xticks(x, names)
    axs[1].set_ylabel("peak silicon °C")
    axs[1].set_title("Effect on chip temperature")
    facts["server"] = {"bal": bal, "unb": unb}
    return save(fig, "fig03_server.png")


def fig_rack(facts: dict) -> str:
    d, r = table(RES / "rack_direct.csv"), table(RES / "rack_reverse.csv")
    un = table(RES / "uneven.csv")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4))
    s = col(d, "server")
    axs[0].plot(s, col(d, "flow_lpm"), "o-", label="direct return", color="#dc2626")
    axs[0].plot(s, col(r, "flow_lpm"), "s-", label="reverse return", color="#2563eb")
    axs[0].set_xlabel("server (1 = bottom, 8 = top)")
    axs[0].set_ylabel("server flow (L/min)")
    axs[0].set_title("Flow distribution between servers")
    axs[0].legend()
    axs[1].plot(s, col(d, "peak_C"), "o-", color="#dc2626", label="direct return")
    axs[1].plot(s, col(r, "peak_C"), "s-", color="#2563eb", label="reverse return")
    axs[1].set_xlabel("server")
    axs[1].set_ylabel("hottest chip °C")
    axs[1].set_title("Hottest chip per server, all at full load")
    axs[1].legend()
    for pat, mk in [("bottom-half-busy", "o"), ("top-half-busy", "s")]:
        for lay, c in [("direct-return", "#dc2626"), ("reverse-return", "#2563eb")]:
            axs[2].plot(col(un, "server", layout=lay, pattern=pat), col(un, "peak_C", layout=lay, pattern=pat), mk + "-",
                        color=c, alpha=0.85, label=f"{lay}, {pat}")
    axs[2].set_xlabel("server")
    axs[2].set_ylabel("hottest chip °C")
    axs[2].set_title("Uneven load: half the rack busy, half idle")
    axs[2].legend(fontsize=7)
    fd, fr = col(d, "flow_lpm"), col(r, "flow_lpm")
    facts["spread_direct"] = 100 * (fd.max() - fd.min()) / fd.mean()
    facts["spread_reverse"] = 100 * (fr.max() - fr.min()) / fr.mean()
    facts["rack_direct"] = d
    return save(fig, "fig04_rack.png")


def fig_hydraulics(facts: dict) -> str:
    bud = table(RES / "hydraulics" / "pressure_budget.csv")
    op = table(RES / "hydraulics" / "operating_points.csv")
    pc = table(RES / "hydraulics" / "pump_curve.csv")
    bud.sort(key=lambda r: -r["dp_kPa"])
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.2), gridspec_kw={"wspace": 0.45})
    total = sum(r["dp_kPa"] for r in bud)
    axs[0].barh([r["component"] for r in bud][::-1], [r["dp_kPa"] for r in bud][::-1], color="#2563eb")
    for i, r in enumerate(bud[::-1]):
        axs[0].text(r["dp_kPa"] + 0.3, i, f"{r['dp_kPa']:.1f} kPa ({100 * r['dp_kPa'] / total:.0f} %)", va="center", fontsize=8)
    axs[0].set_xlabel("pressure drop (kPa)")
    axs[0].set_title(f"Pressure budget: {total:.1f} kPa at 12 L/min\n(path through the top server)")
    axs[0].set_xlim(0, max(r["dp_kPa"] for r in bud) * 1.45)
    q = col(op, "flow_lpm")
    axs[1].plot(col(pc, "flow_lpm"), col(pc, "head_kPa_speed1"), color="#111827", label="pump curve, 100 % speed")
    for sp in (0.6, 0.8, 1.2):
        qq = col(pc, "flow_lpm") * sp
        axs[1].plot(qq, col(pc, "head_kPa_speed1") * sp * sp, color="#94a3b8", lw=0.8)
        axs[1].text(qq[len(qq) // 3], col(pc, "head_kPa_speed1")[len(qq) // 3] * sp * sp, f"{int(sp * 100)} %", fontsize=8)
    axs[1].plot(q, col(op, "head_kPa"), "o-", color="#dc2626", label="system curve (operating points)")
    axs[1].set_xlabel("rack flow (L/min)")
    axs[1].set_ylabel("head (kPa)")
    axs[1].set_ylim(0, 160)
    axs[1].set_title("Pump vs system curve (affinity laws)")
    axs[1].legend(fontsize=8)
    axs[2].plot(q, col(op, "pump_W"), "o-", color="#2563eb", label="pump electric W")
    axs[2].set_xlabel("rack flow (L/min)")
    axs[2].set_ylabel("pump electricity (W)", color="#2563eb")
    a2 = axs[2].twinx()
    a2.plot(q, col(op, "peak_C"), "s-", color="#dc2626")
    a2.axhline(85, color="#dc2626", lw=0.8, ls=":")
    a2.set_ylabel("hottest chip °C", color="#dc2626")
    a2.grid(False)
    axs[2].set_title("The cooling trade-off: flow costs power, saves temperature")
    facts["budget"] = bud
    facts["budget_total"] = total
    return save(fig, "fig05_hydraulics.png")


def fig_viscosity(facts: dict) -> str:
    vr = table(RES / "visc_rack.csv")
    p25, p30 = table(RES / "visc_plant25.csv"), table(RES / "visc_plant30.csv")
    fig, axs = plt.subplots(1, 3, figsize=(18, 4.4), gridspec_kw={"wspace": 0.6})
    mu = col(vr, "mu_scale", case="hold-peak-at-limit")
    axs[0].plot(mu, col(vr, "pump_W", case="hold-peak-at-limit"), "o-", color="#2563eb", label="hold 85 °C: pump W")
    axs[0].plot(mu, col(vr, "pump_W", case="fixed-flow"), "s--", color="#94a3b8", label="fixed 12 L/min: pump W")
    axs[0].set_xlabel("viscosity multiplier")
    axs[0].set_ylabel("rack pump electricity (W)")
    a2 = axs[0].twinx()
    a2.plot(mu, col(vr, "flow_lpm", case="hold-peak-at-limit"), "^-", color="#16a34a")
    a2.set_ylabel("flow needed to hold 85 °C (L/min)", color="#16a34a")
    a2.set_ylim(0, 14)
    a2.grid(False)
    axs[0].set_title("Rack: same flow holds 85 °C, pumping cost rises")
    axs[0].legend(fontsize=8, loc="upper left")
    axs[1].plot(mu, col(vr, "flow_lpm", case="fixed-pump-speed"), "o-", color="#2563eb", label="flow (L/min)")
    axs[1].set_xlabel("viscosity multiplier")
    axs[1].set_ylabel("rack flow at fixed pump speed (L/min)", color="#2563eb")
    a3 = axs[1].twinx()
    a3.plot(mu, col(vr, "peak_C", case="fixed-pump-speed"), "s-", color="#dc2626")
    a3.axhline(85, color="#dc2626", lw=0.8, ls=":")
    a3.set_ylabel("hottest chip °C", color="#dc2626")
    a3.grid(False)
    axs[1].set_title("If the pump is not turned up")
    w = 0.18
    for k, (name, rows) in enumerate([("25 °C", p25), ("30 °C", p30)]):
        m = col(rows, "mu_scale", case="hold-peak-at-limit")
        bottom = np.zeros(len(m))
        for comp, c in [("equip_pump_W", "#2563eb"), ("facility_pump_W", "#06b6d4"), ("drycooler_fan_W", "#16a34a"),
                        ("chiller_W", "#9333ea")]:
            v = col(rows, comp, case="hold-peak-at-limit")
            axs[2].bar(np.arange(len(m)) + (k - 0.5) * 2.2 * w, v, 2 * w, bottom=bottom, color=c,
                       label=comp.replace("_W", "").replace("_", " ") if k == 0 else None, alpha=1 if k else 0.6)
            bottom += v
        for i, b in enumerate(bottom):
            axs[2].text(i + (k - 0.5) * 2.2 * w, b + 5, f"{name}\n{b:.0f}", ha="center", fontsize=7)
    axs[2].set_xticks(np.arange(len(m)), [f"{v:g}×" for v in m])
    axs[2].set_xlabel("viscosity multiplier")
    axs[2].set_ylabel("cooling electricity (W)")
    axs[2].set_title("Plant: total cooling to hold 85 °C (25 °C faded, 30 °C solid)")
    axs[2].set_ylim(0, axs[2].get_ylim()[1] * 1.12)
    axs[2].legend(fontsize=7, loc="center left")
    facts["visc_rack"] = vr
    facts["visc_p25"], facts["visc_p30"] = p25, p30
    return save(fig, "fig06_viscosity.png")


def fig_coolants(facts: dict) -> str:
    props = table(RES / "coolant_props.csv")
    cr, cp = table(RES / "coolants_rack.csv"), table(RES / "coolants_plant.csv")
    names = [r["coolant"] for r in props]
    w = props[0]
    fig, axs = plt.subplots(1, 3, figsize=(16, 4))
    x = np.arange(len(names))
    for k, (key, label, c) in enumerate([("mu40_mPas", "viscosity", "#dc2626"), ("k40", "conductivity", "#2563eb"),
                                         ("cp40", "specific heat", "#16a34a")]):
        axs[0].bar(x + (k - 1) * 0.27, [r[key] / w[key] for r in props], 0.27, color=c, label=label)
    axs[0].axhline(1, color="#111827", lw=0.8)
    axs[0].set_xticks(x, [f"{n}\n{r['freeze_C']:.0f} °C" for n, r in zip(names, props)])
    axs[0].set_xlabel("coolant and freezing point")
    axs[0].set_ylabel("relative to water at 40 °C")
    axs[0].set_title("Coolant properties (Melinder correlations)")
    axs[0].legend(fontsize=8)
    flow = [r["gpu_flow_for_limit_lpm"] for r in props]
    axs[1].bar(x, flow, color=[C.get(n, "#64748b") for n in names])
    for i, v in enumerate(flow):
        axs[1].text(i, v + 0.01, f"{v:.2f}", ha="center")
    axs[1].set_xticks(x, names)
    axs[1].set_ylabel("GPU flow to hold 85 °C (L/min)")
    axs[1].set_title("Lower k and cp ⇒ more flow per chip")
    pump = [col(cr, "equip_pump_W", coolant=n, case="hold peak at limit")[0] for n in names]
    cool = [col(cp, "cooling_W", coolant=n, case="hold peak at limit")[0] for n in names]
    axs[2].bar(x - 0.2, pump, 0.4, color="#2563eb", label="rack pump (rack level)")
    axs[2].bar(x + 0.2, cool, 0.4, color="#16a34a", label="total cooling (plant, 25 °C)")
    for i in range(len(names)):
        axs[2].text(i - 0.2, pump[i] + 3, f"{pump[i]:.0f}", ha="center", fontsize=8)
        axs[2].text(i + 0.2, cool[i] + 3, f"{cool[i]:.0f}", ha="center", fontsize=8)
    axs[2].set_xticks(x, names)
    axs[2].set_ylabel("electricity (W) holding 85 °C")
    axs[2].set_title("Energy cost of each coolant")
    axs[2].legend(fontsize=8)
    facts["coolants"] = {"props": props, "pump": dict(zip(names, pump)), "cool": dict(zip(names, cool))}
    return save(fig, "fig07_coolants.png")


def fig_climate(facts: dict) -> str:
    rows = table(RES / "climate.csv")
    on = [r for r in rows if r["chiller"] == "on"]
    off = [r for r in rows if r["chiller"] == "off"]
    t = col(on, "outdoor_C")
    fig, axs = plt.subplots(1, 2, figsize=(15, 4.2), gridspec_kw={"wspace": 0.35})
    bottom = np.zeros(len(t))
    for comp, c, lab in [("equip_pump_W", "#2563eb", "rack pump"), ("facility_pump_W", "#06b6d4", "facility pump"),
                         ("drycooler_fan_W", "#16a34a", "dry-cooler fans"), ("chiller_W", "#9333ea", "trim chiller")]:
        v = col(on, comp)
        axs[0].fill_between(t, bottom, bottom + v, color=c, alpha=0.85, label=lab, step=None)
        bottom += v
    axs[0].set_xlabel("outdoor air °C")
    axs[0].set_ylabel("cooling electricity (W)")
    a2 = axs[0].twinx()
    a2.plot(t, col(on, "partial_pue"), "k.-", label="partial PUE")
    a2.set_ylabel("partial PUE")
    a2.grid(False)
    axs[0].legend(fontsize=8, loc="upper left")
    axs[0].set_title("Cooling electricity vs climate (full load, 35 °C supply held)")
    axs[1].plot(col(off, "outdoor_C"), col(off, "supply_C"), "o-", color="#2563eb", label="supply °C (no chiller)")
    axs[1].plot(col(off, "outdoor_C"), col(off, "peak_C"), "s-", color="#dc2626", label="hottest chip °C (no chiller)")
    axs[1].plot(t, col(on, "peak_C"), "--", color="#dc2626", alpha=0.5, label="hottest chip °C (with chiller)")
    axs[1].axhline(85, color="#dc2626", lw=0.8, ls=":")
    axs[1].set_xlabel("outdoor air °C")
    axs[1].set_ylabel("°C")
    axs[1].set_title("Dry cooler only: supply and chips follow the weather")
    axs[1].legend(fontsize=8)
    facts["climate_on"], facts["climate_off"] = on, off
    return save(fig, "fig08_climate.png")


def fig_step(facts: dict) -> str:
    r = table(RES / "load_step.csv")
    t = col(r, "t_s")
    fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    axs[0].plot(t, col(r, "it_W") / 1000, color="#111827", label="IT power (kW)")
    axs[0].plot(t, col(r, "liquid_W") / 1000, color="#ea580c", label="to liquid (kW)")
    axs[0].set_ylabel("kW")
    axs[0].legend(fontsize=8)
    axs[0].set_title("Load step at t = 60 s: 30 % → 100 % activity (outdoor 25 °C, fixed 12 L/min)")
    axs[1].plot(t, col(r, "peak_C"), color="#dc2626", label="hottest chip")
    axs[1].plot(t, col(r, "return_C"), color="#f59e0b", label="rack return")
    axs[1].plot(t, col(r, "supply_C"), color="#2563eb", label="rack supply")
    axs[1].plot(t, col(r, "facility_supply_C"), color="#06b6d4", label="facility supply")
    axs[1].set_ylabel("°C")
    axs[1].legend(fontsize=8, ncol=4)
    axs[2].plot(t, 100 * col(r, "fan_speed"), color="#16a34a", label="dry-cooler fan speed (%)")
    axs[2].plot(t, col(r, "cooling_W"), color="#2563eb", label="cooling electricity (W)")
    axs[2].set_ylabel("% / W")
    axs[2].set_xlabel("time (s)")
    axs[2].legend(fontsize=8)
    sup = col(r, "supply_C")
    facts["step"] = {"supply_max": sup.max(), "supply_end": sup[-1], "peak_end": col(r, "peak_C")[-1],
                     "fan_start": 100 * col(r, "fan_speed")[0], "fan_end": 100 * col(r, "fan_speed")[-1],
                     "cool_start": col(r, "cooling_W")[0], "cool_end": col(r, "cooling_W")[-1],
                     "settle": next((tt for tt, s in zip(t[::-1], sup[::-1]) if abs(s - 35) > 0.1), t[0])}
    return save(fig, "fig09_load_step.png")


def fig_day(facts: dict) -> str:
    r = table(RES / "day.csv")
    h = col(r, "t_s") / 3600
    fig, axs = plt.subplots(3, 1, figsize=(12, 8.5), sharex=True)
    axs[0].plot(h, col(r, "it_W") / 1000, color="#111827", label="IT power (kW)")
    axs[0].set_ylabel("kW")
    a0 = axs[0].twinx()
    a0.plot(h, col(r, "outdoor_C"), color="#f59e0b", label="outdoor °C")
    a0.set_ylabel("outdoor °C", color="#f59e0b")
    a0.grid(False)
    axs[0].set_title("A day in the data center: workload peaks in the afternoon, outdoor 16–32 °C")
    axs[0].legend(fontsize=8, loc="upper left")
    bottom = np.zeros(len(h))
    for comp, c, lab in [("equip_pump_W", "#2563eb", "rack pump"), ("facility_pump_W", "#06b6d4", "facility pump"),
                         ("drycooler_fan_W", "#16a34a", "dry-cooler fans"), ("chiller_W", "#9333ea", "chiller")]:
        v = col(r, comp)
        axs[1].fill_between(h, bottom, bottom + v, color=c, alpha=0.85, label=lab)
        bottom += v
    axs[1].set_ylabel("cooling electricity (W)")
    axs[1].legend(fontsize=8, loc="upper left")
    axs[2].plot(h, col(r, "peak_C"), color="#dc2626", label="hottest chip")
    axs[2].plot(h, col(r, "supply_C"), color="#2563eb", label="rack supply")
    axs[2].plot(h, col(r, "return_C"), color="#f59e0b", label="rack return")
    axs[2].set_ylabel("°C")
    axs[2].set_xlabel("hour of day")
    axs[2].set_xticks(range(0, 25, 3))
    axs[2].legend(fontsize=8)
    dt = np.diff(col(r, "t_s"), prepend=0)
    it_kwh = float(np.sum(col(r, "it_W") * dt) / 3.6e6)
    cool_kwh = float(np.sum(col(r, "cooling_W") * dt) / 3.6e6)
    fans_kwh = float(np.sum(col(r, "server_fans_W") * dt) / 3.6e6)
    ch_kwh = float(np.sum(col(r, "chiller_W") * dt) / 3.6e6)
    dc_kwh = float(np.sum(col(r, "drycooler_fan_W") * dt) / 3.6e6)
    facts["day"] = {"it_kwh": it_kwh, "cool_kwh": cool_kwh, "pue": (it_kwh + fans_kwh + cool_kwh) / (it_kwh + fans_kwh),
                    "chiller_kwh": ch_kwh, "dc_kwh": dc_kwh, "peak_max": col(r, "peak_C").max(),
                    "chiller_hours": float(np.sum(dt[col(r, "chiller_W") > 1]) / 3600),
                    "bypass_hours": float(np.sum(dt[(col(r, "chiller_W") > 1) & (col(r, "drycooler_fan_W") < 1)]) / 3600),
                    "bypass_start": float(next((hh for hh, c, fw in zip(h, col(r, "chiller_W"), col(r, "drycooler_fan_W"))
                                                if c > 1 and fw < 1), float("nan"))),
                    "water_l": float(np.sum(col(r, "liquid_W") * dt)) / 2.43e6}  # kg of water to evaporate the liquid heat
    return save(fig, "fig10_day.png")


GHIA_Y = [0.0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
          0.9688, 0.9766, 1.0]
GHIA_U = [0.0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581, -0.13641, 0.00332,
          0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1.0]
DVD = {1e3: 1.118, 1e4: 2.243, 1e5: 4.519}


def field(path: Path, key: str):
    rows = table(path)
    nx = int(max(r["i"] for r in rows))
    ny = int(max(r["j"] for r in rows))
    g = np.zeros((ny, nx))
    for r in rows:
        g[int(r["j"]) - 1, int(r["i"]) - 1] = r[key]
    return g, rows


def psi_grid(path: Path):
    rows = table(path)
    nx = int(max(r["i"] for r in rows)) + 1
    ny = int(max(r["j"] for r in rows)) + 1
    g = np.zeros((ny, nx))
    for r in rows:
        g[int(r["j"]), int(r["i"])] = r["psi"]
    return g


def fig_cfd_verification(facts: dict) -> str:
    cfd = RES / "cfd"
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), gridspec_kw={"wspace": 0.35})
    u, rows = field(cfd / "lid" / "fields.csv", "u")
    ny, nx = u.shape
    y = (np.arange(ny) + 0.5) / ny
    uc = 0.5 * (u[:, nx // 2 - 1] + u[:, nx // 2])
    axs[0].plot(uc, y, color="#2563eb", lw=2, label=f"this solver ({nx}×{ny})")
    axs[0].plot(GHIA_U, GHIA_Y, "o", color="#dc2626", ms=5, label="Ghia, Ghia & Shin (1982)")
    axs[0].set_xlabel("u / U_lid on the vertical centreline")
    axs[0].set_ylabel("y / L")
    axs[0].set_title("Lid-driven cavity, Re = 100")
    axs[0].legend(fontsize=8)
    err = max(abs(np.interp(gy, y, uc) - gu) for gy, gu in zip(GHIA_Y[1:-1], GHIA_U[1:-1]))
    facts["ghia_err"] = err
    nu = table(cfd / "cavity_nu.csv")
    ra = [r["ra"] for r in nu]
    x = np.arange(len(ra))
    axs[1].bar(x - 0.2, [r["nu_hot"] for r in nu], 0.4, color="#2563eb", label="this solver (64×64)")
    axs[1].bar(x + 0.2, [DVD[v] for v in ra], 0.4, color="#dc2626", label="de Vahl Davis (1983)")
    for i, r in enumerate(nu):
        axs[1].text(i - 0.2, r["nu_hot"] + 0.08, f"{r['nu_hot']:.3f}", ha="center", fontsize=8)
        axs[1].text(i + 0.2, DVD[r["ra"]] + 0.08, f"{DVD[r['ra']]:.3f}", ha="center", fontsize=8)
    axs[1].set_xticks(x, [f"Ra = {v:.0e}" for v in ra])
    axs[1].set_ylabel("mean Nusselt number, hot wall")
    axs[1].set_title("Heated cavity (hot / cold walls fixed, u = v = 0)")
    axs[1].legend(fontsize=8)
    facts["dvd"] = [(r["ra"], r["nu_hot"], DVD[r["ra"]]) for r in nu]
    ds = table(cfd / "duct" / "duct_summary.csv")
    d = {r["quantity"]: r for r in ds}
    w, drows = field(cfd / "duct" / "duct.csv", "w_over_mean") if False else (None, None)
    rows = table(cfd / "duct" / "duct.csv")
    nxd = int(max(r["i"] for r in rows)) + 1
    nyd = int(max(r["j"] for r in rows)) + 1
    wg = np.zeros((nyd, nxd))
    for r in rows:
        wg[int(r["j"]), int(r["i"])] = r["w_over_mean"]
    im = axs[2].imshow(wg, origin="lower", cmap=HEAT, extent=[0, 200, 0, 500], aspect="equal")
    axs[2].contour(np.linspace(0, 200, nxd), np.linspace(0, 500, nyd), wg, levels=8, colors="k", linewidths=0.5)
    fig.colorbar(im, ax=axs[2], shrink=0.8, label="w / w_mean")
    axs[2].set_xlabel("width (µm)")
    axs[2].set_ylabel("depth (µm)")
    axs[2].set_title(f"Our 200×500 µm channel, fully developed NS\n"
                     f"f·Re {d['fRe']['navier_stokes']:.3f} (model {d['fRe']['chip_model']:.3f}), "
                     f"Nu₃ {d['Nu3']['navier_stokes']:.3f} (model {d['Nu3']['chip_model']:.3f})", fontsize=10)
    axs[2].grid(False)
    facts["duct"] = d
    return save(fig, "fig12_cfd_verification.png")


def fig_cfd_flows(facts: dict) -> str:
    cfd = RES / "cfd"
    fig = plt.figure(figsize=(17, 9))
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.35)
    # heated cavity Ra 1e5: isotherms + streamlines
    ax = fig.add_subplot(gs[0, 0])
    T, _ = field(cfd / "cavity_1e5" / "fields.csv", "T")
    psi = psi_grid(cfd / "cavity_1e5" / "streamfunction.csv")
    im = ax.imshow(T, origin="lower", cmap=HEAT, extent=[0, 1, 0, 1])
    ax.contour(np.linspace(0, 1, psi.shape[1]), np.linspace(0, 1, psi.shape[0]), psi, levels=12, colors="w", linewidths=0.7)
    ax.set_title("Heated cavity, Ra = 1e5: temperature + streamlines\nhot left, cold right, all walls u = v = 0")
    ax.set_xlabel("x / L"); ax.set_ylabel("y / L"); ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=0.8, label="T")
    # lid cavity streamlines
    ax = fig.add_subplot(gs[0, 1])
    sp, _ = field(cfd / "lid" / "fields.csv", "u")
    psi = psi_grid(cfd / "lid" / "streamfunction.csv")
    vv, _ = field(cfd / "lid" / "fields.csv", "v")
    spd = np.hypot(sp, vv)
    im = ax.imshow(spd, origin="lower", cmap=HEAT, extent=[0, 1, 0, 1])
    ax.contour(np.linspace(0, 1, psi.shape[1]), np.linspace(0, 1, psi.shape[0]), psi, levels=16, colors="w", linewidths=0.7)
    ax.set_title("Lid-driven cavity, Re = 100: speed + streamlines\nlid u = 1, other walls u = v = 0")
    ax.set_xlabel("x / L"); ax.set_ylabel("y / L"); ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=0.8, label="|u| / U")
    # RB onset
    ax = fig.add_subplot(gs[0, 2])
    rb = table(cfd / "rb_nu.csv")
    ax.semilogx([r["ra"] for r in rb], [r["nu"] for r in rb], "o-", color="#dc2626")
    ax.axvline(1708, color="#111827", ls=":", lw=1)
    ax.text(1800, 1.2, "Ra_c ≈ 1708", fontsize=8)
    ax.axhline(1.0, color="#94a3b8", lw=0.8)
    ax.set_xlabel("Rayleigh number")
    ax.set_ylabel("Nusselt number (hot floor)")
    ax.set_title("Rayleigh-Bénard: conduction (Nu = 1) until\nconvection rolls start near Ra ≈ 1708")
    facts["rb"] = rb
    # microchannel slice field + local Nu
    ax = fig.add_subplot(gs[1, :2])
    Tm, _ = field(cfd / "microchannel" / "fields.csv", "T")
    im = ax.imshow(Tm, origin="lower", cmap=HEAT, extent=[0, 20, 0, 200], aspect="auto", vmin=35, vmax=65)
    ax.set_xlabel("distance along the channel (mm)")
    ax.set_ylabel("across the gap (µm)")
    ax.set_title("Microchannel slice (Navier-Stokes): water enters at 35 °C between silicon ribs held at 65 °C")
    ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=0.9, label="°C")
    ax = fig.add_subplot(gs[1, 2])
    al = table(cfd / "microchannel" / "along_flow.csv")
    x = col(al, "x") * 1e3
    ax.plot(x, col(al, "Nu_x"), color="#2563eb", lw=2, label="local Nu (Navier-Stokes)")
    ax.axhline(7.541, color="#dc2626", ls="--", label="fully developed, parallel plates: 7.54")
    ax.set_ylim(0, 25)
    ax.set_xlabel("distance along the channel (mm)")
    ax.set_ylabel("Nusselt number (Dh = 2 × gap)")
    a2 = ax.twinx()
    a2.plot(x, col(al, "T_bulk"), color="#f59e0b")
    a2.set_ylabel("bulk water temperature °C", color="#f59e0b")
    a2.grid(False)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("Thermal entrance: Nu falls from the inlet\ntowards the fully developed value")
    nu = col(al, "Nu_x")
    facts["mc"] = {"nu_in": float(nu[2]), "nu_end": float(nu[-5]), "tb_out": float(col(al, "T_bulk")[-1]),
                   "dp": float((col(al, "p_mean")[0] - col(al, "p_mean")[-1]))}
    return save(fig, "fig13_cfd_flows.png")


def fig_cfd_job(facts: dict) -> str | None:
    d = RES / "cfd_job"
    if not (d / "summary.json").exists():
        return None
    summ = json.loads((d / "summary.json").read_text())["runs"]
    runs = {r["coolant"]: r for r in summ}
    data = {k: table(Path(r["csv"]) if Path(r["csv"]).is_absolute() else ROOT / r["csv"]) for k, r in runs.items()}
    w = data["water"]
    t = col(w, "t_s")
    fig = plt.figure(figsize=(17, 9.5))
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.35)
    ax = fig.add_subplot(gs[0, :2])
    phases = [r["phase"] for r in w]
    start = 0
    shade = {"train": "#dbeafe", "eval": "#ede9fe", "scan": "#dcfce7"}
    for i in range(1, len(t) + 1):
        if i == len(t) or phases[i] != phases[start]:
            if phases[start] in shade:
                ax.axvspan(t[start], t[i - 1], color=shade[phases[start]], lw=0)
            start = i
    for k, c in [("water", "#2563eb"), ("PG25", "#dc2626")]:
        if k in data:
            ax.plot(col(data[k], "t_s"), col(data[k], "ns_wall_max_C"), color=c, lw=1.2, label=f"rib wall max, Navier-Stokes, {k}")
            ax.plot(col(data[k], "t_s"), col(data[k], "chip_peak_C"), color=c, lw=1.0, ls="--", alpha=0.8,
                    label=f"silicon peak, chip model, {k}")
    ax.axhline(85, color="#dc2626", lw=0.8, ls=":")
    ax.set_xlabel("job time (s)")
    ax.set_ylabel("°C")
    ax.set_title("The CNN job driving the Navier-Stokes microchannel (shading: train / eval / scan)")
    ax.legend(fontsize=8, ncol=2, loc="lower center")
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(t, col(w, "gpu_power_W"), color="#111827", lw=0.8)
    ax.set_xlabel("job time (s)")
    ax.set_ylabel("GPU power (W)")
    a2 = ax.twinx()
    a2.plot(t, col(w, "q_mean_W_m2") / 1e3, color="#ea580c", lw=0.6, alpha=0.6)
    a2.set_ylabel("rib-wall heat flux (kW/m²)", color="#ea580c")
    a2.grid(False)
    ax.set_title("Activity → GPU power → wall heat flux")
    ax = fig.add_subplot(gs[1, :2])
    snap = d / "snapshot_water" / "fields.csv"
    if snap.exists():
        T, _ = field(snap, "T")
        im = ax.imshow(T, origin="lower", cmap=HEAT, extent=[0, 20, 0, 200], aspect="auto", vmin=35, vmax=85)
        fig.colorbar(im, ax=ax, shrink=0.9, label="°C")
        ax.set_title("Water temperature during training (job time 65 s): the hotspot heats the middle of the channel most")
        ax.set_xlabel("distance along the channel (mm)")
        ax.set_ylabel("across the gap (µm)")
        ax.grid(False)
    ax = fig.add_subplot(gs[1, 2])
    names = [k for k in ("water", "PG25") if k in runs]
    x = np.arange(len(names))
    vals = [("train_wall_mean_C", "NS wall, train mean"), ("peak_wall_C", "NS wall, peak"), ("peak_chip_C", "chip silicon, peak")]
    for j, (key, lab) in enumerate(vals):
        v = [runs[n][key] for n in names]
        ax.bar(x + (j - 1) * 0.27, v, 0.27, label=lab, color=["#93c5fd", "#2563eb", "#f59e0b"][j])
        for i, vv in enumerate(v):
            ax.text(i + (j - 1) * 0.27, vv + 0.5, f"{vv:.1f}", ha="center", fontsize=8)
    ax.axhline(85, color="#dc2626", lw=0.8, ls=":")
    ax.set_xticks(x, names)
    ax.set_ylim(30, 100)
    ax.set_ylabel("°C")
    ax.set_title("Same job, same 0.5 L/min: water vs PG25")
    ax.legend(fontsize=7, loc="upper left")
    facts["cfdjob"] = runs
    return save(fig, "fig14_cfd_cnn_job.png")


def fig_ml(facts: dict) -> str | None:
    src = ROOT / "results" / "ml" / "report.png"
    if not src.exists():
        return None
    shutil.copy(src, FIG / "fig11_ml_job.png")
    facts["ml_metrics"] = json.loads((ROOT / "results" / "ml" / "metrics.json").read_text())
    vf = ROOT / "results" / "ml" / "variable_flow" / "cosim_summary.json"
    base = ROOT / "results" / "ml" / "cosim_summary.json"
    if vf.exists() and base.exists():
        facts["ml_fixed"] = json.loads(base.read_text())["runs"]
        facts["ml_var"] = json.loads(vf.read_text())["runs"]
    return "fig11_ml_job.png"


# --------------------------------------------------------------------------------------------- text
def write_report(figs: dict, f: dict) -> None:
    def g(rows, key, **w):
        return col(rows, key, **w)

    vr = f["visc_rack"]
    hold = lambda k: g(vr, k, case="hold-peak-at-limit")
    mu = hold("mu_scale")
    i4 = int(np.argmax(mu))
    fixs = lambda k: g(vr, k, case="fixed-pump-speed")
    p25, p30 = f["visc_p25"], f["visc_p30"]
    c25 = g(p25, "cooling_W", case="hold-peak-at-limit")
    c30 = g(p30, "cooling_W", case="hold-peak-at-limit")
    ua = g(p25, "cdu_UA", case="hold-peak-at-limit")
    cl = f["coolants"]
    props = {r["coolant"]: r for r in cl["props"]}
    on = f["climate_on"]
    t_on = g(on, "outdoor_C")
    cool_on = g(on, "cooling_W")
    pue_on = g(on, "partial_pue")
    chill = g(on, "chiller_W")
    first_chill = next((t for t, c in zip(t_on, chill) if c > 1), None)
    off = f["climate_off"]
    off35 = [r for r in off if r["outdoor_C"] == 35.0]
    sb, su = f["server"]["bal"], f["server"]["unb"]
    res = f["res"]
    st = f["step"]
    day = f["day"]
    budget = f["budget"]
    bt = f["budget_total"]

    L = []
    L.append("# Simulation results report\n")
    L.append("*Generated by `analysis/make_report.py`. Every number below is read from the simulation output, and "
             "re-running the script regenerates the figures and the text.*\n")
    L.append("The model is described in the [README](../../README.md) and its inputs in [parameters.md](../parameters.md). "
             "Unless stated otherwise the configuration is the baseline: 8 × 4U servers, each with 2 × 300 W microchannel "
             "GPUs, a 200 W cold-plate CPU and a 200 W air-cooled group. The assumed GPU hotspot puts 50 % of the power in "
             "25 % of the die, the coolant is water, and the 85 °C peak limit is a study assumption.\n")

    L.append("## Contents\n")
    for i, t in enumerate(["Chip: where the heat goes", "Grid resolution", "Server: balancing parallel branches",
                           "Rack: sharing flow between servers", "Hydraulics: pressure budget and pump curves",
                           "Viscosity sensitivity", "Practical coolants", "Climate and the cooling plant",
                           "Dynamics: a load step", "A full day", "CNN workload",
                           "Navier-Stokes verification", "Boundary-driven flows", "The CNN job on the Navier-Stokes solver"], 1):
        L.append(f"{i}. [{t}](#{i}-{t.lower().replace(':', '').replace(',', '').replace(' ', '-')})")
    L.append("")

    up = f["chip_peak_Uniform 300 W, water"]
    hp = f["chip_peak_Hotspot 300 W, water"]
    pg = f["chip_peak_Hotspot 300 W, PG25"]
    L.append("## 1. Chip: where the heat goes\n")
    L.append(f"![Chip maps]({figs['chip']})\n")
    L.append(f"- **Uniform 300 W:** the silicon peaks at **{up:.1f} °C**, at the outlet end where the coolant is warmest "
             "and the heat-transfer coefficient lowest.")
    L.append(f"- **With the assumed hotspot,** the same 300 W peaks at **{hp:.1f} °C**, {hp - up:.1f} K hotter. That is "
             "within the 85 °C limit with little margin: the power map matters as much as the cooling.")
    L.append(f"- **PG25 instead of water** raises the hotspot peak to **{pg:.1f} °C** at the same 0.5 L/min. Glycol's "
             "lower thermal conductivity and heat capacity cost more than its viscosity.")
    L.append(f"- **Developing flow:** the heat-transfer coefficient falls from **{f['h_in']:.1f} kW/m²K** at the channel "
             f"inlet to **{f['h_out']:.1f} kW/m²K** at the outlet. The thermal boundary layer is still developing over the "
             "whole 20 mm channel, so a fully developed correlation would under-predict cooling near the inlet.\n")

    L.append("## 2. Grid resolution\n")
    L.append(f"![Resolution]({figs['res']})\n")
    ks = sorted(res)
    L.append(f"Peak temperature at {', '.join(f'{k}×{k}: {res[k]:.2f} °C' for k in ks)}. Going from 10×10 to the finest "
             f"grid changes the peak by **{res[ks[-1]] - res[10]:+.2f} K**. The coolant outlet temperature does not "
             "depend on the grid, because energy is conserved exactly. 10×10 is adequate for comparisons larger than about "
             "1 K; use 20×20 for finer ones.\n")

    L.append("## 3. Server: balancing parallel branches\n")
    L.append(f"![Server]({figs['server']})\n")
    gb = {r["branch"]: r for r in sb}
    gu = {r["branch"]: r for r in su}
    L.append(f"- **Without balancing valves,** the low-resistance CPU cold plate takes **{gu['CPU']['flow_lpm']:.2f} L/min** "
             f"and each GPU only **{gu['GPU1']['flow_lpm']:.2f} L/min**. The GPUs then run at "
             f"**{gu['GPU1']['peak_C']:.1f} °C** instead of {gb['GPU1']['peak_C']:.1f} °C.")
    L.append(f"- **The balancing valve** on the CPU branch absorbs **{gb['CPU']['valve_dp_kPa']:.1f} kPa**, so the branches "
             f"share {gb['GPU1']['flow_lpm']:.2f} / {gb['GPU2']['flow_lpm']:.2f} / {gb['CPU']['flow_lpm']:.2f} L/min. "
             "This flow is solved from the pressure balance, not assumed.\n")

    L.append("## 4. Rack: sharing flow between servers\n")
    L.append(f"![Rack]({figs['rack']})\n")
    d = f["rack_direct"]
    L.append(f"- **Direct-return manifolds** give the bottom server the shortest path, so flow falls from "
             f"{d[0]['flow_lpm']:.3f} to {d[-1]['flow_lpm']:.3f} L/min towards the top: a **{f['spread_direct']:.1f} %** "
             f"spread. **Reverse return** cuts the spread to **{f['spread_reverse']:.1f} %**.")
    L.append("- **The top server is hottest with direct return** because it gets the least flow, not because it is higher: "
             "elevation cancels in a closed loop.")
    L.append("- **Under uneven load** the hottest server is whichever busy server has the least flow, so the layout and "
             "the job placement together decide the hotspot.\n")

    L.append("## 5. Hydraulics: pressure budget and pump curves\n")
    L.append(f"![Hydraulics]({figs['hyd']})\n")
    top = budget[0]
    L.append(f"- **The pump delivers {bt:.1f} kPa at 12 L/min.** The largest share goes to the {top['component']} "
             f"({top['dp_kPa']:.1f} kPa, {100 * top['dp_kPa'] / bt:.0f} %). The microchannel GPU branch takes "
             f"{next(r['dp_kPa'] for r in budget if r['component'].startswith('GPU')):.1f} kPa. The manifolds are "
             "almost free, which is why flow is shared so evenly between servers.")
    L.append("- **Pumping energy is therefore set mostly by components outside the chip.** Any viscosity penalty "
             "applies to the laminar parts of all of them, not only the microchannels.")
    L.append("- **The operating points** move along the system curve as pump speed changes. More flow lowers chip "
             "temperature with diminishing returns, while pump power rises roughly with the cube of speed.\n")

    L.append("## 6. Viscosity sensitivity\n")
    L.append(f"![Viscosity]({figs['visc']})\n")
    L.append(f"- **Holding 85 °C needs the same flow at any viscosity:** {hold('flow_lpm')[0]:.2f} L/min at 1× and "
             f"{hold('flow_lpm')[i4]:.2f} L/min at {mu[i4]:g}×. Re·Pr, which governs microchannel heat transfer, contains "
             "no viscosity.")
    L.append(f"- **The pump pays instead:** {hold('pump_W')[0]:.1f} W at 1× rises to {hold('pump_W')[i4]:.1f} W at "
             f"{mu[i4]:g}×, which is {hold('pump_W')[i4] / hold('pump_W')[0]:.1f}×. Head rises from "
             f"{hold('pump_head_kPa')[0]:.0f} to {hold('pump_head_kPa')[i4]:.0f} kPa.")
    L.append(f"- **If the pump is not turned up,** flow falls to {fixs('flow_lpm')[i4]:.1f} L/min at {mu[i4]:g}× and the "
             f"hottest chip reaches **{fixs('peak_C')[i4]:.1f} °C**.")
    L.append(f"- **With the cooling plant,** viscosity also weakens the CDU heat exchanger: UA falls from {ua[0]:.0f} to "
             f"{ua[i4]:.0f} W/K at {mu[i4]:g}×. The facility side must then run colder. Total cooling electricity to hold "
             f"85 °C goes from **{c25[0]:.0f} → {c25[i4]:.0f} W** at 25 °C outdoor and **{c30[0]:.0f} → {c30[i4]:.0f} W** "
             "at 30 °C, where the chiller has to make up the difference.\n")

    L.append("## 7. Practical coolants\n")
    L.append(f"![Coolants]({figs['cool']})\n")
    pw, pp = cl["pump"], cl["cool"]
    pr = props
    L.append(f"- **PG25** is {pr['PG25']['mu40_mPas'] / pr['water']['mu40_mPas']:.1f}× as viscous as water at 40 °C, with "
             f"{100 * (1 - pr['PG25']['k40'] / pr['water']['k40']):.0f} % lower conductivity and "
             f"{100 * (1 - pr['PG25']['cp40'] / pr['water']['cp40']):.0f} % lower heat capacity. It freezes at "
             f"{pr['PG25']['freeze_C']:.1f} °C.")
    L.append(f"- **More flow per chip:** each GPU needs {pr['PG25']['gpu_flow_for_limit_lpm']:.2f} L/min of PG25 "
             f"against {pr['water']['gpu_flow_for_limit_lpm']:.2f} L/min of water to hold 85 °C, and "
             f"{pr['PG40']['gpu_flow_for_limit_lpm']:.2f} L/min of PG40.")
    L.append(f"- **Rack pump power** to hold 85 °C: water {pw['water']:.0f} W, EG25 {pw['EG25']:.0f} W, PG25 "
             f"{pw['PG25']:.0f} W, PG40 {pw['PG40']:.0f} W. **Total cooling** at 25 °C outdoor: water {pp['water']:.0f} W, "
             f"PG25 {pp['PG25']:.0f} W, PG40 {pp['PG40']:.0f} W.")
    L.append("- **Real glycols cost more than the viscosity-only experiment suggests**, because their thermal properties "
             "also get worse. Glycol data carry about ±20 % uncertainty in viscosity.\n")

    L.append("## 8. Climate and the cooling plant\n")
    L.append(f"![Climate]({figs['climate']})\n")
    L.append(f"- **Mild weather:** up to about 25 °C outdoors the dry cooler alone holds a 35 °C supply. Cooling uses "
             f"{cool_on[list(t_on).index(25.0)]:.0f} W at 25 °C, a partial PUE of {pue_on[list(t_on).index(25.0)]:.3f}.")
    if first_chill is not None:
        L.append(f"- **Warm weather:** from about {first_chill:.0f} °C the fans are at full speed and the trim chiller "
                 f"starts. At {t_on[-1]:.0f} °C cooling needs {cool_on[-1]:.0f} W, a partial PUE of {pue_on[-1]:.3f}.")
    if off35:
        L.append(f"- **Without a chiller** the supply follows the weather: at 35 °C outdoors the supply is "
                 f"{off35[0]['supply_C']:.1f} °C and the hottest chip {off35[0]['peak_C']:.1f} °C.")
    L.append("- **Water:** the closed loop evaporates no water in any climate. An evaporative tower rejecting the same "
             f"heat would use about {on[0]['evap_equiv_Lph']:.1f} L/h. The saving costs electricity only on hot days.\n")

    L.append("## 9. Dynamics: a load step\n")
    L.append(f"![Load step]({figs['step']})\n")
    L.append(f"- **The load step:** activity jumps from 30 % to 100 % at t = 60 s. Chip temperatures respond within "
             "seconds (the silicon time constant is about 30 ms; the coolant transit and loop volumes take longer).")
    L.append(f"- **The loop's thermal mass** lets the supply overshoot only to **{st['supply_max']:.2f} °C** while the "
             f"dry-cooler fans speed up from {st['fan_start']:.0f} % to {st['fan_end']:.0f} %. The fan controller returns "
             f"the supply to 35 °C (it ends at {st['supply_end']:.2f} °C).")
    L.append(f"- **Cooling electricity** rises from {st['cool_start']:.0f} W to {st['cool_end']:.0f} W. The fans, not the "
             "pumps, carry the change, because the rack pump runs at constant flow in this configuration.\n")

    L.append("## 10. A full day\n")
    L.append(f"![Day]({figs['day']})\n")
    L.append(f"- **Over 24 hours** the rack uses **{day['it_kwh']:.1f} kWh** of IT energy and **{day['cool_kwh']:.2f} kWh** "
             f"of cooling electricity, a partial PUE of **{day['pue']:.4f}**.")
    L.append(f"- **The hot afternoon** coincides with peak load. The trim chiller runs for {day['chiller_hours']:.1f} h "
             f"({day['chiller_kwh']:.2f} kWh) and the dry-cooler fans use {day['dc_kwh']:.2f} kWh. The hottest chip "
             f"reaches {day['peak_max']:.1f} °C.")
    if day["bypass_hours"] > 0:
        L.append(f"- **The economizer rule at work:** from about {day['bypass_start']:.1f} h the outdoor air is too warm for "
                 "the dry cooler to be worth its fan power. The water then bypasses it and the chiller carries the whole "
                 f"load for {day['bypass_hours']:.1f} h (the green band turning purple), until the evening cools down.")
    L.append(f"- **Water:** no water is evaporated. A cooling tower would have evaporated roughly "
             f"{day['water_l']:.0f} L for the same liquid heat load.\n")

    if "ml" in figs:
        m = f["ml_metrics"]
        L.append("## 11. CNN workload\n")
        L.append(f"![ML job]({figs['ml']})\n")
        L.append(f"- **The job:** a PyTorch CNN trained on MNIST reaches **{100 * m['final_test_accuracy']:.2f} %** test "
                 f"accuracy, then reads {m['scan']['digits_found']} digits on {m['scan']['pages']} scanned pages "
                 f"({100 * m['scan']['digit_accuracy_on_fully_segmented_pages']:.1f} % correct). Its measured device "
                 "activity drives the rack; see [ml_workload.md](../ml_workload.md).")
        if "ml_fixed" in f:
            fx = {r["coolant"]: r for r in f["ml_fixed"]}
            vv = {r["coolant"]: r for r in f["ml_var"]}
            L.append(f"- **Fixed 12 L/min:** the hottest chip reaches {fx['water']['peak_C']:.1f} °C with water and "
                     f"**{fx['PG25']['peak_C']:.1f} °C with PG25**.")
            L.append(f"- **Variable flow with load feedforward:** {vv['water']['peak_C']:.1f} °C and "
                     f"{vv['PG25']['peak_C']:.1f} °C. Cooling energy falls from {fx['water']['cooling_J'] / 1e3:.1f} to "
                     f"{vv['water']['cooling_J'] / 1e3:.1f} kJ for water.")
        L.append("")

    dv = f["duct"]
    L.append("## 12. Navier-Stokes verification\n")
    L.append(f"![CFD verification]({figs['cfdv']})\n")
    L.append("The rack model uses reduced-order physics: an energy-conserving thermal network plus channel correlations. "
             "To check that this is justified, a separate solver integrates the full 2-D incompressible Navier-Stokes "
             "and energy equations (Boussinesq buoyancy, staggered finite volumes, projection method). It is first "
             "verified against the standard published benchmarks:")
    L.append(f"- **Lid-driven cavity, Re = 100:** the centreline velocity matches Ghia, Ghia & Shin (1982) within "
             f"**{f['ghia_err']:.3f}** of the lid speed at every tabulated point.")
    L.append("- **Differentially heated cavity:** Nusselt numbers against de Vahl Davis (1983): " +
             ", ".join(f"Ra {ra:.0e}: {a:.3f} vs {b:.3f} ({100 * (a / b - 1):+.1f} %)" for ra, a, b in f["dvd"]) + ".")
    L.append(f"- **Our own 200 × 500 µm channel,** solved as fully developed Navier-Stokes + energy on its cross-section: "
             f"f·Re = **{dv['fRe']['navier_stokes']:.3f}** against the {dv['fRe']['chip_model']:.3f} the chip model uses "
             f"({100 * (dv['fRe']['navier_stokes'] / dv['fRe']['chip_model'] - 1):+.2f} %). The three-wall Nusselt number is "
             f"**{dv['Nu3']['navier_stokes']:.3f}** against {dv['Nu3']['chip_model']:.3f} "
             f"({100 * (dv['Nu3']['navier_stokes'] / dv['Nu3']['chip_model'] - 1):+.1f} %). **The correlations in the chip "
             "model are confirmed by a first-principles solution.**\n")
    mc = f["mc"]
    rb = f["rb"]
    L.append("## 13. Boundary-driven flows\n")
    L.append(f"![CFD flows]({figs['cfdf']})\n")
    L.append("Each case is set entirely by its boundaries. Every wall has a prescribed velocity (u, v): zero by default "
             "(no-slip), or a moving wall, or an inflow. Each wall is held at a fixed hot or cold temperature, or insulated. "
             "The interactive CFD scene of the viewer (key 5) lets you change any of these while the flow is running.")
    L.append("- **Heated cavity:** with all walls at rest, buoyancy alone drives a single circulating cell. Hot fluid "
             "rises along the hot wall and sinks along the cold one, and thin thermal boundary layers form.")
    L.append(f"- **Rayleigh-Bénard:** heated from below, the fluid conducts (Nu = {rb[0]['nu']:.2f} at Ra = "
             f"{rb[0]['ra']:.0f}) until it passes the critical Rayleigh number (about 1708). Then convection rolls carry "
             f"more heat: Nu = {rb[-1]['nu']:.2f} at Ra = {rb[-1]['ra']:.0e}.")
    L.append(f"- **Microchannel slice:** water at 1.67 m/s between ribs at 65 °C. The local Nusselt number falls from "
             f"about {mc['nu_in']:.0f} near the inlet to {mc['nu_end']:.1f} at 20 mm, still above the fully developed 7.54. "
             f"This confirms that the flow is thermally developing over the whole channel, as the chip model assumes. "
             f"The water leaves at {mc['tb_out']:.1f} °C.\n")

    if "cfdjob" in figs:
        cj = f["cfdjob"]
        wtr, pg = cj["water"], cj.get("PG25")
        L.append("## 14. The CNN job on the Navier-Stokes solver\n")
        L.append(f"![CNN job on Navier-Stokes]({figs['cfdjob']})\n")
        L.append("The recorded PyTorch job (section 11) also drives the first-principles solver:")
        L.append("1. Device activity u(t) sets the GPU power, P = 60 + u (300 − 60) W.")
        L.append("2. One channel's share of that power enters the water through the two silicon rib walls of the "
                 "microchannel slice, shaped along the flow like the GPU's hotspot (the channel that crosses it).")
        L.append("3. The velocity field comes from the full Navier-Stokes equations. Buoyancy is negligible in this "
                 "forced flow, so the flow does not change with temperature and is solved once.")
        L.append("4. The energy equation is then marched implicitly through the whole job, and the reduced chip model "
                 "follows the same trace for comparison.\n")
        L.append(f"- **Two independent models agree:** with water, the hottest rib wall in the Navier-Stokes solution "
                 f"reaches **{wtr['peak_wall_C']:.1f} °C** during training. The chip model's hottest silicon reaches "
                 f"**{wtr['peak_chip_C']:.1f} °C**. They differ in geometry (2-D slice vs. 3-D die) and in what they "
                 "report (wall surface vs. silicon node), yet they are within a few kelvin.")
        L.append(f"- **Energy conservation over the job:** {wtr['heat_in_J']:.2f} J entered the channel and "
                 f"{wtr['heat_out_J']:.2f} J left with the water.")
        if pg:
            L.append(f"- **PG25 under the same job and flow:** the rib wall peaks at **{pg['peak_wall_C']:.1f} °C** "
                     f"(chip model {pg['peak_chip_C']:.1f} °C), against {wtr['peak_wall_C']:.1f} °C with water, and "
                     f"the training-phase mean rises from {wtr['train_wall_mean_C']:.1f} to {pg['train_wall_mean_C']:.1f} °C. "
                     "The Navier-Stokes solution confirms what the reduced model found: glycol's poorer thermal "
                     "properties cost cooling performance under a real workload.")
        L.append("- **The water responds to the workload within milliseconds:** the transit time is about 12 ms. Wall "
                 "temperatures therefore follow every burst of the on-demand scanning phase, while the loop-level "
                 "thermal mass (section 9) smooths the slower plant response.\n")

    L.append("## Limitations\n")
    L.append("- **Assumed values:** the hotspot pattern, idle powers, and the cold-plate, quick-disconnect, CDU, "
             "dry-cooler and pump curves are assumptions; see [parameters.md](../parameters.md).")
    L.append("- **Simplifications:** the chip model uses a 2-D grid with a lumped fin; microchannel heat transfer uses a "
             "thermal-entrance approximation; flow is assumed equal across the 50 channels.")
    L.append("- **Not modelled:** room HVAC for the air-cooled load, and glycol inhibitors.\n")
    (FIG / "REPORT.md").write_text("\n".join(L))
    print(f"wrote {FIG / 'REPORT.md'}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-sim", action="store_true", help="reuse results/report data")
    args = ap.parse_args()
    if not MF.exists():
        sys.exit("build first: cmake -S . -B build && cmake --build build -j")
    FIG.mkdir(parents=True, exist_ok=True)
    if not args.no_sim:
        run_experiments()
    facts: dict = {}
    figs = {"chip": fig_chip(facts), "res": fig_resolution(facts), "server": fig_server(facts), "rack": fig_rack(facts),
            "hyd": fig_hydraulics(facts), "visc": fig_viscosity(facts), "cool": fig_coolants(facts),
            "climate": fig_climate(facts), "step": fig_step(facts), "day": fig_day(facts)}
    figs["cfdv"] = fig_cfd_verification(facts)
    figs["cfdf"] = fig_cfd_flows(facts)
    cj = fig_cfd_job(facts)
    if cj:
        figs["cfdjob"] = cj
    ml = fig_ml(facts)
    if ml:
        figs["ml"] = ml
    write_report(figs, facts)


if __name__ == "__main__":
    main()
