#!/usr/bin/env bash
# Convenience driver that runs the shotgun pipeline end to end.
# Review config.sh and metadata/samples.tsv first. Raw FASTQ files and the
# Kraken2 database must already be present in the working directory.
set -euo pipefail
HERE="$(dirname "$0")"

bash "$HERE/00_download_reference.sh"
bash "$HERE/01_qc_fastqc.sh"
bash "$HERE/02_trim_fastp.sh"
bash "$HERE/03_host_removal.sh"
bash "$HERE/04_taxonomy_kraken2.sh"
bash "$HERE/05_build_biom.sh"
bash "$HERE/06_assembly_megahit.sh"
echo "Pipeline finished. Run R/taxonomic_composition.R for diversity/visualisation."
