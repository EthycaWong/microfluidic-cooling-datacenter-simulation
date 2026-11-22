"""Plot the ML job and the simulated data center's response (after train_mnist.py and `mfsim cosim`).

  .venv/bin/python ml/report.py [--dir results/ml] [--out results/ml/report.png]
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PHASE_COLORS = {"setup": "#e5e7eb", "train": "#dbeafe", "eval": "#ede9fe", "scan": "#dcfce7", "idle": "#f3f4f6"}
COOLANT_COLORS = ["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#ea580c"]


def read_csv(path: Path) -> dict[str, list]:
    with open(path) as f:
        rows = list(csv.DictReader(f))
    out: dict[str, list] = {k: [] for k in rows[0]}
    for r in rows:
        for k, v in r.items():
            try:
                out[k].append(float(v))
            except ValueError:
                out[k].append(v)
    return out


def shade_phases(ax, t, phases):
    start = 0
    for i in range(1, len(t) + 1):
        if i == len(t) or phases[i] != phases[start]:
            ph = phases[start]
            if ph in PHASE_COLORS:
                ax.axvspan(t[start], t[i - 1], color=PHASE_COLORS[ph], lw=0, zorder=0)
            start = i


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/ml")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    d = Path(args.dir)
    metrics = json.loads((d / "metrics.json").read_text())
    summary = json.loads((d / "cosim_summary.json").read_text())
    runs = summary["runs"]
    sims = [(r["coolant"], read_csv(Path(r["csv"]))) for r in runs]

    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(3, 3, hspace=0.45, wspace=0.38)
    fig.suptitle(f"PyTorch CNN (MNIST) on the simulated microfluidic-cooled rack: {summary['gpus']} GPUs on "
                 f"{summary['job_servers']} servers, outdoor {summary['outdoor_C']:.0f} °C", fontsize=14, fontweight="bold")

    # 1. training curve
    ax = fig.add_subplot(gs[0, 0])
    ep = [e["epoch"] for e in metrics["epochs"]]
    ax.plot(ep, [e["train_loss"] for e in metrics["epochs"]], "o-", color="#2563eb", label="training loss")
    ax.set_xlabel("epoch"); ax.set_ylabel("loss", color="#2563eb")
    ax2 = ax.twinx()
    ax2.plot(ep, [100 * e["test_accuracy"] for e in metrics["epochs"]], "s-", color="#16a34a", label="test accuracy")
    ax2.set_ylabel("test accuracy (%)", color="#16a34a")
    ax.set_title(f"Training: {100 * metrics['final_test_accuracy']:.2f} % test accuracy ({metrics['device']})")

    # 2. activity trace (simulated time)
    t0, s0 = sims[0][1]["t_s"], sims[0][1]
    ax = fig.add_subplot(gs[0, 1:])
    shade_phases(ax, t0, s0["phase"])
    ax.plot(t0, s0["u"], color="#111827", lw=0.8)
    ax.set_ylim(0, 1.05); ax.set_ylabel("activity u"); ax.set_xlabel("simulated time (s)")
    ax.set_title("Measured device activity drives every job server (shading: setup / train / eval / scan / idle)")
    for ph, col in PHASE_COLORS.items():
        ax.fill_between([], [], color=col, label=ph)
    ax.legend(loc="upper right", ncol=5, fontsize=8)

    # 3. power
    ax = fig.add_subplot(gs[1, :2])
    shade_phases(ax, t0, s0["phase"])
    ax.plot(t0, [v / 1000 for v in s0["it_W"]], color="#111827", lw=1.2, label="IT power (kW)")
    ax.set_ylabel("IT power (kW)")
    axc = ax.twinx()
    for i, (name, s) in enumerate(sims):
        axc.plot(s["t_s"], s["cooling_W"], color=COOLANT_COLORS[i % 5], lw=1.2, label=f"cooling, {name} (W)")
    axc.set_ylabel("cooling electricity (W)")
    ax.set_xlabel("simulated time (s)")
    ax.set_title("Rack IT power follows the job; cooling electricity = pumps + dry-cooler fans + chiller")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axc.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    # 4. scan example
    ax = fig.add_subplot(gs[1, 2])
    img = d / "scan_example.png"
    if img.exists():
        ax.imshow(plt.imread(img)); ax.axis("off")
        sc = metrics["scan"]
        ax.set_title(f"Scan phase: {sc['pages']} pages, {100 * sc['digit_accuracy_on_fully_segmented_pages']:.1f} % digits correct",
                     fontsize=10)

    # 5. temperatures
    ax = fig.add_subplot(gs[2, :2])
    shade_phases(ax, t0, s0["phase"])
    for i, (name, s) in enumerate(sims):
        c = COOLANT_COLORS[i % 5]
        ax.plot(s["t_s"], s["peak_C"], color=c, lw=1.3, label=f"hottest chip, {name}")
        ax.plot(s["t_s"], s["return_C"], color=c, lw=0.9, ls="--", label=f"coolant return, {name}")
    ax.axhline(85, color="#dc2626", lw=0.8, ls=":")
    ax.set_ylabel("°C"); ax.set_xlabel("simulated time (s)")
    ax.set_title("Temperatures (85 °C = study limit)")
    ax.legend(loc="upper right", fontsize=8, ncol=2)

    # 6. cooling energy by phase (IT energy is the same for every coolant: same job)
    ax = fig.add_subplot(gs[2, 2])
    names = [r["coolant"] for r in runs]
    groups = [("before job", ["waiting"]), ("train + eval", ["setup", "train", "eval"]), ("scan", ["scan"]),
              ("after job", ["idle", "done"])]
    colors = ["#cbd5e1", "#2563eb", "#16a34a", "#94a3b8"]
    bottom = [0.0] * len(runs)
    for (label, phs), col in zip(groups, colors):
        vals = [sum(r["phases"].get(p, {}).get("cooling_J", 0.0) for p in phs) / 1000 for r in runs]
        ax.bar(names, vals, bottom=bottom, color=col, label=label)
        bottom = [b + v for b, v in zip(bottom, vals)]
    for i, r in enumerate(runs):
        ax.text(i, bottom[i] * 1.02, f"{bottom[i]:.1f} kJ\npPUE {r['partial_pue']:.4f}\n"
                f"{1000 * r['J_per_train_image']:.1f} mJ/img, {1000 * r['J_per_digit']:.0f} mJ/digit",
                ha="center", va="bottom", fontsize=7.5)
    ax.set_ylim(0, max(bottom) * 1.9)
    ax.set_ylabel("cooling electricity (kJ)")
    ax.set_title("Cooling energy over the run, by phase")
    ax.legend(fontsize=7, loc="upper center", ncol=2)

    out = Path(args.out) if args.out else d / "report.png"
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
