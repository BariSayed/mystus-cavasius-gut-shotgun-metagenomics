#!/usr/bin/env bash
# Step 1 — Read quality control with FastQC (v0.12.1), summarised with MultiQC (v1.28).
set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$FASTQC_DIR"
fastqc --threads 4 --outdir "$FASTQC_DIR" *.fastq.gz

# MultiQC aggregation (published Methods; command reconstructed, not in the
# original committed script).
multiqc "$FASTQC_DIR" -o multiqc_summary
