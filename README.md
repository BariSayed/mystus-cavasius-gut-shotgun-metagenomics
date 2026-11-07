# Shotgun metagenomic analysis pipeline — gut microbiome of antibiotic-treated *Mystus cavasius*

Code, figures and reproducibility documentation for:

> Begum M, Barsha KF, Rahman MM, Sarkar MMH, Chowdhury SF, Bhowmik S, Shormi AS,
> Bari SM. (2025). **Gut microbiome profiling of antibiotic-treated *Mystus
> cavasius* using culture-based and shotgun metagenomic approaches.** *Antonie
> van Leeuwenhoek*, 118(2): 183. https://doi.org/10.1007/s10482-025-02195-5

The study profiles the gut microbiome of *Mystus cavasius* (Gulsha Tengra) fed
oxytetracycline-supplemented feed for 90 days, combining culture-based 16S rRNA
analysis with Illumina shotgun metagenomics (read QC, host-read removal, Kraken2
classification, MEGAHIT assembly) and diversity/visualisation in R.

> **Honesty note.** This repository was reorganised from a single script file,
> with the authors' R analysis scripts and the published figures added later.
> Some published Methods details differ from the committed code, and a few steps
> have no saved script. Every such case is documented in
> [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md); reconstructed or
> unverified items are labelled as such.

## Repository layout

```
.
├── README.md  ·  LICENSE  ·  CITATION.cff  ·  environment.yml  ·  .gitignore
├── scripts/            # command-line pipeline (bash)
│   ├── config.sh  ·  run_all.sh
│   ├── 00_download_reference.sh   01_qc_fastqc.sh      02_trim_fastp.sh
│   ├── 03_host_removal.sh         04_taxonomy_kraken2.sh
│   └── 05_build_biom.sh           06_assembly_megahit.sh
├── R/                  # analysis + figure scripts (phyloseq / vegan / ggtree …)
│   ├── cfu.R                              phylogenetic_tree.R
│   ├── antibiotic_sensitivity.R          taxonomic_composition.R
│   ├── relative_abundance.R              taxonomic_composition_count_table.R
│   ├── alpha_diversity.R                 beta_diversity.R
│   └── heatmaps.R
├── figures/            # published figures (Figure_1 … Figure_8)
├── data/               # accessions + SRA download script
├── metadata/           # sample sheet + metadata template
├── supplementary/      # Supplementary Tables S1–S4 (+ CSV exports of S4)
├── docs/REPRODUCIBILITY.md
└── legacy/original_script.txt     # the original committed `script`, verbatim
```

## Quick start

```bash
conda env create -f environment.yml && conda activate shotgun-mg

# 1. Raw reads (NCBI SRA, BioProject PRJNA1282882)
bash data/download_sra.sh            # then rename to <sample>_R1/_R2.fastq.gz

# 2. Run the command-line pipeline (provide a Kraken2 DB; see docs/REPRODUCIBILITY.md)
bash scripts/run_all.sh              # -> output.biom + assemblies

# 3. Figures/diversity (need output.biom + metadata.csv in the working directory)
Rscript R/taxonomic_composition.R
Rscript R/alpha_diversity.R
Rscript R/beta_diversity.R
```

## Pipeline (command line)

| # | Step | Tool (Methods version) | Script |
|---|------|------------------------|--------|
| 0 | Host reference + index | Bowtie2 2.5.4 | `scripts/00_download_reference.sh` |
| 1 | Read QC | FastQC 0.12.1 / MultiQC 1.28 | `scripts/01_qc_fastqc.sh` |
| 2 | Trimming | fastp *(see note)* | `scripts/02_trim_fastp.sh` |
| 3 | Host-read removal | Bowtie2 2.5.4 / SAMtools 1.21 | `scripts/03_host_removal.sh` |
| 4 | Classification | Kraken2 2.1.3 | `scripts/04_taxonomy_kraken2.sh` |
| 5 | BIOM table | kraken-biom | `scripts/05_build_biom.sh` |
| 6 | Assembly | MEGAHIT 1.2.9 | `scripts/06_assembly_megahit.sh` |

## Figures and their scripts

| Figure | Content | Script |
|--------|---------|--------|
| 1 | Experimental workflow (schematic) | — |
| 2 | Total viable bacterial count (CFU/mL) | `R/cfu.R` |
| 3 | 16S phylogenetic tree | `R/phylogenetic_tree.R` *(see note)* |
| 4 | Antibiotic (OTC) sensitivity | `R/antibiotic_sensitivity.R` |
| 5 | Taxonomic composition (phylum→species) | `R/taxonomic_composition.R` |
| 6 | Top-50 genus/species heatmaps | `R/heatmaps.R` |
| 7 | Alpha diversity (Observed, Chao1, Shannon, Simpson) | `R/alpha_diversity.R` |
| 8 | Beta diversity (Bray-Curtis PCoA/NMDS) | `R/beta_diversity.R` |

Supporting: `R/relative_abundance.R`, `R/taxonomic_composition_count_table.R`.
Notes and Methods discrepancies: [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Data

Raw reads, reference genomes, Kraken2 databases and large intermediates are not
tracked in Git (see `.gitignore`). Accessions and sources:

- **Shotgun reads:** NCBI SRA **PRJNA1282882** — 5 runs: `SRR34283200`–`SRR34283204`.
  Samples **CR1, CR2** (control) and **TR1, TR2, TR3** (treated).
- **16S Sanger (12 isolates):** GenBank `PV810321`–`PV810323`, `PV809930`–`PV809938`.
- **Host genome:** *Ictalurus punctatus* `GCF_001660625.3_Coco_2.0`.

Details in [`data/README.md`](data/README.md); supplementary tables in
[`supplementary/`](supplementary/).

## Changes made during reorganisation

- Split the single committed `script` into numbered bash steps + separate R
  scripts; generalised single-sample commands into loops.
- **Bug fixes:** inline comments after `\` in the `bowtie2` command; a duplicated
  comma in an R theme; the invalid fastp flag `--cut_by_quality3` → `--cut_tail`.
- Removed hardcoded Windows absolute paths from all R scripts (now relative);
  normalised line endings.
- Added README, LICENSE, CITATION.cff, environment.yml, .gitignore, data-download
  script, supplementary tables, the published figures, and reproducibility docs.
- Preserved the original `script` verbatim at `legacy/original_script.txt`.

## License & citation

Code: [MIT](LICENSE). The article is © The Author(s) under exclusive licence to
Springer Nature and is not covered by this licence. Please cite both the
repository ([`CITATION.cff`](CITATION.cff)) and the article
(doi:10.1007/s10482-025-02195-5).
