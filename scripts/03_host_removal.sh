#!/usr/bin/env bash
# Step 3 — Remove host (channel catfish) reads (one pass per sample).
# Align trimmed reads to the host index with Bowtie2 (v2.5.4), keep read pairs
# where BOTH mates are unmapped (samtools -f 12), and export them back to FASTQ
# with SAMtools (v1.21).
#
# NOTE: the original committed script placed "# ..." comments on the same lines
#   as the backslash line-continuations of the bowtie2 command. In bash a
#   comment after "\" breaks the continuation and the command fails. The
#   comments have been moved to their own lines here so the command runs.
set -euo pipefail
source "$(dirname "$0")/config.sh"

for S in "${SAMPLES[@]}"; do
  bowtie2 --threads "$THREADS" \
    -x "./${HOST_INDEX}" \
    -1 "${S}_R1.clean.fastq.gz" \
    -2 "${S}_R2.clean.fastq.gz" \
    -S "${S}_host.sam"

  # Keep pairs with both mates unmapped (-f 12), output BAM.
  samtools view -f 12 -b "${S}_host.sam" > "${S}_host_unmap.bam"

  # Convert the host-free BAM back to paired FASTQ.
  samtools fastq \
    -1 "${S}_host_unmap_R1.clean.fastq.gz" \
    -2 "${S}_host_unmap_R2.clean.fastq.gz" \
    "${S}_host_unmap.bam"

  # The large intermediate SAM can be removed once the BAM exists.
  rm -f "${S}_host.sam"
done
