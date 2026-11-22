#!/usr/bin/env bash
# Train the CNN while recording its activity, replay that job on the simulated rack + cooling plant
# for each coolant, and draw the report. Run from the repository root:
#   ml/run_job.sh [epochs] [coolants]          e.g.  ml/run_job.sh 8 water,PG25
set -euo pipefail
EPOCHS="${1:-8}"
COOLANTS="${2:-water,PG25}"
PY=.venv/bin/python
[ -x "$PY" ] || { echo "create the venv first: python3 -m venv .venv && .venv/bin/pip install torch torchvision matplotlib"; exit 1; }
[ -x build/mfsim ] || { echo "build first: cmake -S . -B build && cmake --build build -j"; exit 1; }
$PY ml/train_mnist.py --epochs "$EPOCHS" --scan-seconds 30
./build/mfsim cosim --config configs/plant_baseline.json --trace results/ml/activity_trace.csv \
    --coolants "$COOLANTS" --job-servers 8 --lead 60 --tail 300
$PY ml/report.py --dir results/ml
