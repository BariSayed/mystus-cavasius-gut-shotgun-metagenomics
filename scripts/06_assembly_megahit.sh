#!/usr/bin/env bash
# Step 6 — De novo metagenome assembly with MEGAHIT (v1.2.9), one assembly per
# sample, discarding contigs shorter than 500 bp.
#
# NOTE: the published Methods also report assembly assessment with QUAST
#   (v5.2.0). No QUAST command was present in the committed script; see
#   docs/REPRODUCIBILITY.md. A template invocation is given (commented) below.
set -euo pipefail
source "$(dirname "$0")/config.sh"

for S in "${SAMPLES[@]}"; do
  megahit \
    -1 "${S}_host_unmap_R1.clean.fastq.gz" \
    -2 "${S}_host_unmap_R2.clean.fastq.gz" \
    --min-contig-len 500 \
    --tmp-dir ./ \
    --memory 12 \
    --num-cpu-threads "$THREADS" \
    --out-dir "${S}_megahit" \
    --out-prefix "$S"

  # Reconstructed, not in the original committed script:
  # quast.py "${S}_megahit/${S}.contigs.fa" -o "${S}_quast"
done
