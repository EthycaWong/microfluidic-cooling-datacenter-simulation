"""Measure how busy the compute device is, as a workload trace for the data-center simulation.

The simulation drives every server with an activity u(t) in [0, 1] (power = P_idle + u (P_full - P_idle)).
Here u is measured from the real training job: in each fixed wall-clock window, u is the fraction of
the window the device spent executing work (forward/backward/optimizer or inference), with the device
synchronised so asynchronous GPU/MPS queues are counted when they actually run. Data loading, logging
and idle gaps count as not busy.

The trace is written two ways:
  * live:  appended and flushed every window, so the 3D viewer can follow the job in real time;
  * final: the complete trace, used by `mfsim cosim` for a reproducible energy analysis.
"""
from __future__ import annotations

import csv
import threading
import time
from contextlib import contextmanager
from pathlib import Path


class ActivityMonitor:
    def __init__(self, live_path: Path | None, window_s: float = 0.1):
        self.window = window_s
        self.t0 = time.perf_counter()
        self._lock = threading.Lock()
        self._busy_intervals: list[tuple[float, float]] = []  # (start, end), seconds since t0
        self._phase = "setup"
        self._counters = {"train_samples": 0, "infer_samples": 0}
        self.rows: list[dict] = []
        self._stop = threading.Event()
        self._live = None
        if live_path is not None:
            live_path.parent.mkdir(parents=True, exist_ok=True)
            self._live = open(live_path, "w", newline="")
            self._writer = csv.writer(self._live)
            self._writer.writerow(["t_s", "u", "phase", "train_samples", "infer_samples"])
            self._live.flush()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    # ---- instrumentation API ----
    def now(self) -> float:
        return time.perf_counter() - self.t0

    def set_phase(self, phase: str) -> None:
        with self._lock:
            self._phase = phase

    def count(self, key: str, n: int) -> None:
        with self._lock:
            self._counters[key] += n

    @contextmanager
    def busy(self, sync=None):
        """Mark the enclosed block as device work. `sync` is called before the block ends so queued
        accelerator work is included (torch.mps.synchronize / torch.cuda.synchronize)."""
        start = self.now()
        try:
            yield
        finally:
            if sync is not None:
                sync()
            end = self.now()
            with self._lock:
                self._busy_intervals.append((start, end))

    # ---- sampling ----
    def _emit(self, w_start: float, w_end: float) -> None:
        with self._lock:
            busy = 0.0
            keep = []
            for s, e in self._busy_intervals:
                overlap = min(e, w_end) - max(s, w_start)
                if overlap > 0:
                    busy += overlap
                if e > w_end:  # still relevant for the next window
                    keep.append((s, e))
            self._busy_intervals = keep
            row = {
                "t_s": round(w_end, 4),
                "u": round(min(1.0, busy / (w_end - w_start)), 4),
                "phase": self._phase,
                "train_samples": self._counters["train_samples"],
                "infer_samples": self._counters["infer_samples"],
            }
        self.rows.append(row)
        if self._live is not None:
            self._writer.writerow(row.values())
            self._live.flush()

    def _run(self) -> None:
        k = 1
        while not self._stop.is_set():
            target = self.t0 + k * self.window
            delay = target - time.perf_counter()
            if delay > 0:
                time.sleep(delay)
            self._emit((k - 1) * self.window, k * self.window)
            k += 1

    def close(self) -> None:
        self._stop.set()
        self._thread.join()
        if self._live is not None:
            self._live.close()

    def write_trace(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["t_s", "u", "phase", "train_samples", "infer_samples"])
            w.writeheader()
            w.writerows(self.rows)
