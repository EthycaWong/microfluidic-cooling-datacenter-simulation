"""CNN for handwritten digits (MNIST) and a simple page scanner that finds and reads digits."""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

MNIST_MEAN, MNIST_STD = 0.1307, 0.3081


class DigitCNN(nn.Module):
    """Two conv layers + two fully connected layers (~1.2 M parameters), ~99 % on MNIST in a few epochs."""

    def __init__(self) -> None:
        super().__init__()
        self.conv1 = nn.Conv2d(1, 32, 3)
        self.conv2 = nn.Conv2d(32, 64, 3)
        self.drop1 = nn.Dropout(0.25)
        self.fc1 = nn.Linear(9216, 128)
        self.drop2 = nn.Dropout(0.5)
        self.fc2 = nn.Linear(128, 10)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = F.relu(self.conv1(x))
        x = F.max_pool2d(F.relu(self.conv2(x)), 2)
        x = torch.flatten(self.drop1(x), 1)
        x = self.drop2(F.relu(self.fc1(x)))
        return self.fc2(x)

    @staticmethod
    def flops_per_image() -> int:
        """Multiply-accumulates x2 for one forward pass (for reporting)."""
        c1 = 26 * 26 * 32 * (3 * 3 * 1)
        c2 = 24 * 24 * 64 * (3 * 3 * 32)
        f1 = 9216 * 128
        f2 = 128 * 10
        return 2 * (c1 + c2 + f1 + f2)


# ---------------------------------------------------------------------------------------------
# Page scanning: compose "scanned sheets" of handwritten digits, then segment and classify them.
# ---------------------------------------------------------------------------------------------
def make_sheet(images: np.ndarray, labels: np.ndarray, rng: np.random.Generator, rows: int = 4, cols: int = 8,
               cell: int = 40) -> tuple[np.ndarray, list[int]]:
    """Paste rows*cols MNIST digits (28x28, ink = high values) onto a white page with random jitter and
    scanner noise. Returns the page (float 0..1, ink dark like a real scan) and the labels in reading order."""
    h, w = rows * cell + 20, cols * cell + 20
    page = np.zeros((h, w), dtype=np.float32)
    truth: list[int] = []
    idx = rng.choice(len(images), rows * cols, replace=False)
    for n, k in enumerate(idx):
        r, c = divmod(n, cols)
        y = 10 + r * cell + int(rng.integers(0, cell - 28))
        x = 10 + c * cell + int(rng.integers(0, cell - 28))
        page[y:y + 28, x:x + 28] = np.maximum(page[y:y + 28, x:x + 28], images[k] / 255.0)
        truth.append(int(labels[k]))
    page = 1.0 - page                                            # dark ink on white paper
    page += rng.normal(0.0, 0.04, page.shape).astype(np.float32)  # scanner noise
    return np.clip(page, 0.0, 1.0), truth


def _runs(profile: np.ndarray, thresh: float, min_len: int) -> list[tuple[int, int]]:
    on = profile > thresh
    out, start = [], None
    for i, v in enumerate(on):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if i - start >= min_len:
                out.append((start, i))
            start = None
    if start is not None and len(on) - start >= min_len:
        out.append((start, len(on)))
    return out


def segment_page(page: np.ndarray) -> list[tuple[int, int, int, int]]:
    """Projection-profile segmentation (as used for scanned forms): find text lines from the row ink
    profile, then digits within each line from the column profile. Returns boxes (y0, y1, x0, x1) in
    reading order."""
    ink = np.clip(0.75 - page, 0.0, None)  # remove paper + noise, keep ink
    boxes = []
    for y0, y1 in _runs(ink.sum(axis=1), thresh=1.0, min_len=6):
        line = ink[y0:y1]
        runs = _runs(line.sum(axis=0), thresh=0.6, min_len=1)
        # Merge fragments of one digit (strokes with a small vertical gap, e.g. a broken 4 or 5) as
        # long as the merged width still fits one digit; drop specks.
        merged: list[list[int]] = []
        for x0, x1 in runs:
            if merged and x0 - merged[-1][1] <= 4 and x1 - merged[-1][0] <= 24:
                merged[-1][1] = x1
            else:
                merged.append([x0, x1])
        for x0, x1 in merged:
            if x1 - x0 < 2 and line[:, x0:x1].sum() < 3.0:
                continue
            ys = np.nonzero(line[:, x0:x1].sum(axis=1) > 0.3)[0]
            if len(ys):
                boxes.append((y0 + ys[0], y0 + ys[-1] + 1, x0, x1))
    return boxes


def crop_to_mnist(page: np.ndarray, box: tuple[int, int, int, int]) -> np.ndarray:
    """MNIST-style normalisation: fit the digit in a 20x20 box, centre by mass in a 28x28 field."""
    y0, y1, x0, x1 = box
    d = np.clip(0.75 - page[y0:y1, x0:x1], 0.0, None) / 0.75
    h, w = d.shape
    s = 20.0 / max(h, w)
    nh, nw = max(1, round(h * s)), max(1, round(w * s))
    t = torch.from_numpy(d)[None, None]
    small = F.interpolate(t, size=(nh, nw), mode="bilinear", align_corners=False)[0, 0].numpy()
    out = np.zeros((28, 28), dtype=np.float32)
    oy, ox = (28 - nh) // 2, (28 - nw) // 2
    out[oy:oy + nh, ox:ox + nw] = small
    tot = out.sum()
    if tot > 0:
        cy = (out.sum(axis=1) * np.arange(28)).sum() / tot
        cx = (out.sum(axis=0) * np.arange(28)).sum() / tot
        out = np.roll(out, (int(round(14 - cy)), int(round(14 - cx))), axis=(0, 1))
    return out


def to_tensor_batch(crops: list[np.ndarray]) -> torch.Tensor:
    x = torch.from_numpy(np.stack(crops))[:, None]
    return (x - MNIST_MEAN) / MNIST_STD
