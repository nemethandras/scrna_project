#!/usr/bin/env bash
# Run the pipeline sequentially for multiple samples.
# Detached from the terminal — safe to close your laptop.
# Usage: source .env && bash run_batch.sh

set -euo pipefail

LOG=logs/batch_run.log
mkdir -p logs
exec > >(tee -a "$LOG") 2>&1

echo "========================================"
echo "Batch run started: $(date)"
echo "========================================"

run() {
    local sample=$1
    local run_id=$2
    local extra=${3:-}
    echo ""
    echo "--- $sample ($run_id) started: $(date) ---"
    python run_pipeline.py --run-id "$run_id" --sample "$sample" --foreground $extra
    echo "--- $sample ($run_id) finished: $(date) ---"
}

run SRR5071662 SRR5071662_hg38
run SRR5071667 SRR5071667_hg38
run SRR5071672 SRR5071672_hg38
run SRR5071691 SRR5071691_hg38
run SRR5071692 SRR5071692_hg38 "--no-db"

echo ""
echo "========================================"
echo "Batch run finished: $(date)"
echo "========================================"
