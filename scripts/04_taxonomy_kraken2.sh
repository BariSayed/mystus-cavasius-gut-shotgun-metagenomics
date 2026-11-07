#!/usr/bin/env bash
# Step 4 — Taxonomic classification of host-free reads with Kraken2 (v2.1.3).
#
# NOTE (database): the committed script used the pre-built 'k2_pluspf_20210517'
#   database (set in config.sh as KRAKEN2_DB). The published Methods instead
#   describe the "standard RefSeq" database. See docs/REPRODUCIBILITY.md.
# NOTE (abundance): the published Methods report abundance re-estimation with
#   Bracken (v2.6.2). No Bracken command was present in the committed script, and
#   the downstream relative-abundance calculation in R/taxonomic_composition.R
#   is performed directly from the Kraken2 report counts via phyloseq. A Bracken
#   step is therefore NOT included here; add it if you intend to match the paper.
set -euo pipefail
source "$(dirname "$0")/config.sh"

for S in "${SAMPLES[@]}"; do
  kraken2 \
    --threads "$THREADS" \
    --quick \
    --paired \
    --db "$KRAKEN2_DB" \
    --report "${S}.kreport" \
    --output "${S}.kraken" \
    "${S}_host_unmap_R1.clean.fastq.gz" \
    "${S}_host_unmap_R2.clean.fastq.gz"
done
