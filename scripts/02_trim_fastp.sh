#!/usr/bin/env bash
# Step 2 — Adapter and quality trimming with fastp (one run per sample).
#
# NOTE 1 (tool): the committed script and the working-directory outputs
#   (*_fastp.json / *_fastp.html) show that fastp was the tool actually used for
#   trimming. The published Methods instead list "Trimmomatic v1.28"; v1.28
#   matches the cited MultiQC version and is therefore most likely a
#   typographical error. See docs/REPRODUCIBILITY.md.
# NOTE 2 (flag): the original script wrote "--cut_by_quality3", which is not a
#   valid fastp option. It is replaced here with the standard 3'-end quality
#   trimming flag ("--cut_tail") using the same window size and mean-quality
#   values. The exact behaviour of the original invocation cannot be reproduced
#   verbatim because the original flag would not have been accepted by fastp.
set -euo pipefail
source "$(dirname "$0")/config.sh"

for S in "${SAMPLES[@]}"; do
  fastp \
    -i "${S}_R1.fastq.gz" -I "${S}_R2.fastq.gz" \
    -o "${S}_R1.clean.fastq.gz" -O "${S}_R2.clean.fastq.gz" \
    --detect_adapter_for_pe \
    --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --length_required 50 \
    --thread 8 \
    --html "${S}_fastp.html" --json "${S}_fastp.json"
done
