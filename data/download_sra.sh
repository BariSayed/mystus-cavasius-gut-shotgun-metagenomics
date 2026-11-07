#!/usr/bin/env bash
# Download the raw shotgun metagenomic reads from NCBI SRA (BioProject
# PRJNA1282882) using the SRA Toolkit, and convert them to gzipped paired FASTQ.
#
# Requires: sra-tools (prefetch, fasterq-dump) and pigz/gzip.
# The five runs below are all runs linked to PRJNA1282882 (two controls, three
# treated), identified from the SRA run titles shown in brackets.
set -euo pipefail

# run accession   # SRA run title
RUNS=(
  SRR34283204     # Tengra Microbiome Control1
  SRR34283203     # Tengra Microbiome Control2
  SRR34283202     # Tengra Microbiome Dose1
  SRR34283201     # Tengra Microbiome Dose2
  SRR34283200     # Tengra Microbiome Dose3
)

for RUN in "${RUNS[@]}"; do
  prefetch "$RUN"
  fasterq-dump "$RUN" --split-files --threads 8 --outdir .
  gzip -f "${RUN}_1.fastq" "${RUN}_2.fastq"
done

echo "Downloaded ${#RUNS[@]} runs. Rename to the <sample>_R1/_R2.fastq.gz"
echo "convention used by the scripts once the run<->sample mapping is confirmed."
