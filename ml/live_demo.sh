#!/usr/bin/env bash
# Watch the simulated data center respond to a real PyTorch training job in real time.
# Opens the 3D viewer following results/ml/activity_live.csv, then starts training.
#   ml/live_demo.sh [epochs] [job servers: 2|4|8]
set -euo pipefail
EPOCHS="${1:-8}"
SERVERS="${2:-8}"
mkdir -p results/ml
rm -f results/ml/activity_live.csv
./build/mfsim_viewer --config configs/plant_baseline.json --job live --job-servers "$SERVERS" &
VIEWER=$!
sleep 3
.venv/bin/python ml/train_mnist.py --epochs "$EPOCHS" --scan-seconds 30
echo "Training finished; close the viewer window when done."
wait $VIEWER
