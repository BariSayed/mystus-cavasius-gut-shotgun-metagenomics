#!/usr/bin/env bash
# Shared configuration sourced by the pipeline scripts.
# Edit the values here, then run the numbered scripts in order (or run_all.sh).
set -euo pipefail

# Number of CPU threads (values kept from the original committed script where given).
THREADS=12

# Sample base names. Each must have <sample>_R1.fastq.gz and <sample>_R2.fastq.gz
# in the working directory. These are the five samples used to build the .biom
# table in the original script.
SAMPLES=(Control_S1 Control_S8 SingleAbx_S2 SingleAbx_S3 SingleAbx_S4)

# Host reference (channel catfish, Ictalurus punctatus, GCF_001660625.3_Coco_2.0).
HOST_FNA="genome.fna"          # decompressed reference FASTA
HOST_INDEX="genome_db"         # Bowtie2 index prefix

# Kraken2 database directory.
# NOTE: the committed script used the pre-built 'k2_pluspf_20210517' database.
# The published Methods instead describe the "standard RefSeq" database.
# See docs/REPRODUCIBILITY.md for this discrepancy.
KRAKEN2_DB="./k2_pluspf_20210517/"

# Output directories.
FASTQC_DIR="fastqc_results"
