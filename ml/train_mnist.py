"""Train a digit-recognition CNN on MNIST while recording the compute activity it creates.

The recorded activity u(t) becomes the workload of the simulated data center:
  * live:  results/ml/activity_live.csv   (the 3D viewer follows it in real time, --live-trace)
  * final: results/ml/activity_trace.csv  (input to `mfsim cosim` for the energy analysis)

Job phases (each labelled in the trace):
  setup  - download/load data, build the model (device idle)
  train  - epochs of minibatch SGD with on-device augmentation
  eval   - test-set accuracy after each epoch
  scan   - on-demand inference: synthetic scanned pages arrive at random times; each page is segmented
           and every digit on it is read by the CNN (inference requests)
  idle   - job finished, servers wait for the next job

Usage:
  .venv/bin/python ml/train_mnist.py --epochs 8 --scan-seconds 30
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from torchvision import datasets

from activity import ActivityMonitor
from digits import MNIST_MEAN, MNIST_STD, DigitCNN, crop_to_mnist, make_sheet, segment_page, to_tensor_batch

ROOT = Path(__file__).resolve().parent.parent


def pick_device(name: str) -> torch.device:
    if name != "auto":
        return torch.device(name)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def sync_fn(device: torch.device):
    if device.type == "cuda":
        return torch.cuda.synchronize
    if device.type == "mps":
        return torch.mps.synchronize
    return None


def augment(x: torch.Tensor, gen: torch.Generator) -> torch.Tensor:
    """Small random rotations/scales/shifts on the device (keeps the device, not the CPU, busy)."""
    n = x.shape[0]
    ang = (torch.rand(n, generator=gen) - 0.5) * math.radians(20)
    sc = 1.0 + (torch.rand(n, generator=gen) - 0.5) * 0.2
    tx, ty = [(torch.rand(n, generator=gen) - 0.5) * 0.2 for _ in range(2)]
    theta = torch.stack([torch.stack([sc * torch.cos(ang), -sc * torch.sin(ang), tx], 1),
                         torch.stack([sc * torch.sin(ang), sc * torch.cos(ang), ty], 1)], 1).to(x.device)
    grid = F.affine_grid(theta, list(x.shape), align_corners=False)
    bg = (0.0 - MNIST_MEAN) / MNIST_STD  # pad with background, not zero
    return F.grid_sample(x - bg, grid, align_corners=False, padding_mode="zeros") + bg


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--lr", type=float, default=1.0, help="Adadelta learning rate")
    ap.add_argument("--device", default="auto", help="auto | cpu | mps | cuda")
    ap.add_argument("--data", default=str(ROOT / "data"), help="dataset download folder")
    ap.add_argument("--out", default=str(ROOT / "results" / "ml"))
    ap.add_argument("--window", type=float, default=0.1, help="activity sampling window, s")
    ap.add_argument("--scan-seconds", type=float, default=30.0, help="length of the on-demand scanning phase")
    ap.add_argument("--scan-rate", type=float, default=40.0,
                    help="peak page arrivals per second; demand ramps 15%% -> 100%% -> 15%% over the phase")
    ap.add_argument("--idle-tail", type=float, default=5.0, help="seconds of idle recorded after the job")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(args.seed)
    device = pick_device(args.device)
    sync = sync_fn(device)
    mon = ActivityMonitor(out / "activity_live.csv", args.window)
    print(f"device: {device}  (torch {torch.__version__})")
    print(f"live activity trace: {out / 'activity_live.csv'}")

    # ---- setup: dataset (downloaded once from the torchvision MNIST mirror) ----
    mon.set_phase("setup")
    train_ds = datasets.MNIST(args.data, train=True, download=True)
    test_ds = datasets.MNIST(args.data, train=False, download=True)
    norm = lambda u8: ((u8.float() / 255.0 - MNIST_MEAN) / MNIST_STD)[:, None]
    x_train, y_train = norm(train_ds.data), train_ds.targets
    x_test, y_test = norm(test_ds.data), test_ds.targets
    model = DigitCNN().to(device)
    opt = torch.optim.Adadelta(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.StepLR(opt, step_size=1, gamma=0.7)
    n_params = sum(p.numel() for p in model.parameters())
    print(f"MNIST: {len(x_train)} training / {len(x_test)} test images; CNN with {n_params:,} parameters")
    gen = torch.Generator().manual_seed(args.seed)

    history = []
    t_train0 = mon.now()
    for epoch in range(1, args.epochs + 1):
        # ---- train ----
        mon.set_phase("train")
        model.train()
        perm = torch.randperm(len(x_train), generator=gen)
        loss_sum, seen, t_ep = 0.0, 0, time.perf_counter()
        for i in range(0, len(perm), args.batch):
            idx = perm[i:i + args.batch]
            xb, yb = x_train[idx], y_train[idx]          # host-side batch assembly (not device work)
            with mon.busy(sync):
                xb, yb = xb.to(device), yb.to(device)
                xb = augment(xb, gen)
                opt.zero_grad(set_to_none=True)
                loss = F.cross_entropy(model(xb), yb)
                loss.backward()
                opt.step()
                lv = loss.item()
            loss_sum += lv * len(idx)
            seen += len(idx)
            mon.count("train_samples", len(idx))
        sched.step()
        train_s = time.perf_counter() - t_ep
        # ---- eval ----
        mon.set_phase("eval")
        model.eval()
        correct = 0
        with torch.no_grad():
            for i in range(0, len(x_test), 1000):
                with mon.busy(sync):
                    pred = model(x_test[i:i + 1000].to(device)).argmax(1).cpu()
                correct += int((pred == y_test[i:i + 1000]).sum())
        acc = correct / len(x_test)
        history.append({"epoch": epoch, "train_loss": loss_sum / seen, "test_accuracy": acc, "train_seconds": train_s})
        print(f"epoch {epoch}: loss {loss_sum / seen:.4f}  test accuracy {100 * acc:.2f}%  ({train_s:.1f} s)")
    t_train1 = mon.now()

    # ---- scan: on-demand inference on synthetic scanned pages ----
    mon.set_phase("scan")
    rng = np.random.default_rng(args.seed)
    test_np, test_lbl = test_ds.data.numpy(), test_ds.targets.numpy()
    pages = digits_total = digits_found = digits_correct = pages_segmented = 0
    latencies = []
    example = None
    t_scan0 = mon.now()
    next_arrival = time.perf_counter()
    while mon.now() - t_scan0 < args.scan_seconds:
        # Wait for the next request. Arrivals are a Poisson process whose rate follows a demand curve
        # (quiet -> busy -> quiet), generated by thinning; the device idles between requests.
        while True:
            next_arrival += rng.exponential(1.0 / args.scan_rate)
            tau = min(1.0, (next_arrival - mon.t0 - t_scan0) / args.scan_seconds)
            if rng.random() < 0.15 + 0.85 * math.sin(math.pi * tau) ** 2:
                break
        delay = next_arrival - time.perf_counter()
        if delay > 0:
            time.sleep(delay)
        t_req = time.perf_counter()
        page, truth = make_sheet(test_np, test_lbl, rng)
        boxes = segment_page(page)                       # host-side preprocessing
        crops = [crop_to_mnist(page, b) for b in boxes]
        if crops:
            with mon.busy(sync), torch.no_grad():
                pred = model(to_tensor_batch(crops).to(device)).argmax(1).cpu().numpy().tolist()
        else:
            pred = []
        latencies.append(time.perf_counter() - t_req)
        pages += 1
        digits_total += len(truth)
        digits_found += len(pred)
        if len(pred) == len(truth):
            pages_segmented += 1
            digits_correct += sum(int(a == b) for a, b in zip(pred, truth))
        mon.count("infer_samples", len(pred))
        if example is None:
            example = (page, boxes, pred, truth)
    t_scan1 = mon.now()

    # ---- idle tail ----
    mon.set_phase("idle")
    time.sleep(args.idle_tail)
    mon.close()
    mon.write_trace(out / "activity_trace.csv")
    torch.save(model.state_dict(), out / "mnist_cnn.pt")

    rows = mon.rows
    mean_u = lambda ph: float(np.mean([r["u"] for r in rows if r["phase"] == ph] or [0]))
    metrics = {
        "device": str(device),
        "platform": platform.platform(),
        "torch": torch.__version__,
        "parameters": n_params,
        "flops_per_image_forward": DigitCNN.flops_per_image(),
        "epochs": history,
        "final_test_accuracy": history[-1]["test_accuracy"],
        "train_wall_seconds": t_train1 - t_train0,
        "train_samples": rows[-1]["train_samples"],
        "train_samples_per_second": rows[-1]["train_samples"] / max(t_train1 - t_train0, 1e-9),
        "scan": {
            "pages": pages,
            "digits_on_pages": digits_total,
            "digits_found": digits_found,
            "pages_fully_segmented": pages_segmented,
            "digit_accuracy_on_fully_segmented_pages": digits_correct / max(1, 32 * pages_segmented),
            "segmentation_recall": digits_found / max(1, digits_total),
            "seconds": t_scan1 - t_scan0,
            "digits_per_second": digits_found / max(t_scan1 - t_scan0, 1e-9),
            "median_page_latency_ms": 1000 * float(np.median(latencies)) if latencies else None,
        },
        "mean_activity": {ph: mean_u(ph) for ph in ("setup", "train", "eval", "scan", "idle")},
        "trace": str(out / "activity_trace.csv"),
        "window_s": args.window,
    }
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2))
    print(f"scan phase: {pages} pages ({pages_segmented} fully segmented), {digits_found}/{digits_total} digits found, "
          f"accuracy {100 * metrics['scan']['digit_accuracy_on_fully_segmented_pages']:.2f}%, "
          f"{metrics['scan']['digits_per_second']:.0f} digits/s")
    print("mean activity u by phase: " + ", ".join(f"{k} {v:.2f}" for k, v in metrics["mean_activity"].items()))
    print(f"wrote {out / 'activity_trace.csv'}, {out / 'metrics.json'}, {out / 'mnist_cnn.pt'}")

    if example is not None:
        save_scan_figure(example, out / "scan_example.png")


def save_scan_figure(example, path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError:
        return
    page, boxes, pred, truth = example
    fig, ax = plt.subplots(figsize=(10, 5.4))
    ax.imshow(page, cmap="gray", vmin=0, vmax=1)
    for k, (y0, y1, x0, x1) in enumerate(boxes):
        ok = k < len(truth) and len(pred) == len(truth) and pred[k] == truth[k]
        col = "#16a34a" if ok else "#dc2626"
        ax.add_patch(Rectangle((x0 - 1, y0 - 1), x1 - x0 + 2, y1 - y0 + 2, fill=False, ec=col, lw=1.5))
        ax.text(x0, y0 - 3, str(pred[k]), color=col, fontsize=11, fontweight="bold")
    ax.set_title("Scanned page: digits found by projection segmentation, read by the CNN (green = correct)")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
