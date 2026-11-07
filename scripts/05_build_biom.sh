#!/usr/bin/env bash
# Step 5 — Combine the per-sample Kraken2 reports into a single BIOM table
# (JSON format) for import into phyloseq.
set -euo pipefail
source "$(dirname "$0")/config.sh"

REPORTS=()
for S in "${SAMPLES[@]}"; do
  REPORTS+=("${S}.kreport")
done

kraken-biom "${REPORTS[@]}" --fmt json -o output.biom
echo "Wrote output.biom from: ${REPORTS[*]}"
