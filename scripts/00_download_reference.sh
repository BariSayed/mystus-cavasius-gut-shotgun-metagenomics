#!/usr/bin/env bash
# Step 0 — Download the host reference genome and build the Bowtie2 index.
# Host: Ictalurus punctatus (channel catfish), assembly GCF_001660625.3_Coco_2.0
# as reported in the published Methods.
set -euo pipefail
source "$(dirname "$0")/config.sh"

URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/660/625/GCF_001660625.3_Coco_2.0/GCF_001660625.3_Coco_2.0_genomic.fna.gz"

wget -c "$URL" -O GCF_001660625.3_Coco_2.0_genomic.fna.gz
gunzip -kf GCF_001660625.3_Coco_2.0_genomic.fna.gz

# NOTE: the original committed script downloaded the file above but then referred
# to a file named "genome.fna" when building the index. We make that explicit
# here by copying the decompressed reference to the configured HOST_FNA name.
cp -f GCF_001660625.3_Coco_2.0_genomic.fna "$HOST_FNA"

bowtie2-build "$HOST_FNA" "$HOST_INDEX"
echo "Host reference index built: ${HOST_INDEX}.*.bt2"
