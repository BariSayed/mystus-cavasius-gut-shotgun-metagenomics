# Reproducibility notes and known discrepancies

How the results in Begum et al. (2025, *Antonie van Leeuwenhoek* 118:183,
doi:10.1007/s10482-025-02195-5) relate to the code in this repository. It
distinguishes **committed code**, the **published Methods**, and
**working-directory / figure evidence**. Where they disagree, the disagreement is
stated rather than silently resolved; scripts not originally saved are marked
**reconstructed** or **not available**.

## 1. Command-line pipeline (`scripts/`)

| Step | Script | Tool / version (Methods) | Status / notes |
|------|--------|--------------------------|----------------|
| Read QC | `01_qc_fastqc.sh` | FastQC 0.12.1; MultiQC 1.28 | FastQC was in the committed script; the **MultiQC** command is **reconstructed** (a `multiqc_summary/` folder exists). |
| Trimming | `02_trim_fastp.sh` | "Trimmomatic v1.28" | **Discrepancy.** The committed script and `*_fastp.json/html` outputs show **fastp** was used. "v1.28" matches MultiQC and is likely a typo. The invalid flag `--cut_by_quality3` was replaced with `--cut_tail`. fastp **version not recorded.** |
| Host-read removal | `03_host_removal.sh` | Bowtie2 2.5.4; SAMtools 1.21; host `GCF_001660625.3_Coco_2.0` | Matches committed script. **Bug fixed:** inline `#` comments after `\` line-continuations in the `bowtie2` command (broke the command). |
| Classification | `04_taxonomy_kraken2.sh` | Kraken2 2.1.3, "standard RefSeq" | **Discrepancy (database):** committed script used `k2_pluspf_20210517`, not standard RefSeq. |
| Abundance | *not available* | Bracken 2.6.2 | **No committed script.** Relative abundance is computed in R from the Kraken2/BIOM counts; add Bracken to match the paper exactly. |
| BIOM table | `05_build_biom.sh` | kraken-biom (not named in Methods) | Matches committed script. Version not recorded. |
| Assembly | `06_assembly_megahit.sh` | MEGAHIT 1.2.9 (min 500 bp) | Matches committed script. |
| Assembly QC | *not available* | QUAST 5.2.0 | **No committed script** (commented template in `06`). |

## 2. R analyses (`R/`) and figures (`figures/`)

These scripts were provided by the authors after the initial reorganisation and
reproduce the paper's figures. **Edit applied to all of them:** hardcoded Windows
absolute paths (`setwd("C:/Users/.../M_cavasius")` and absolute output paths)
were replaced with relative paths so the scripts run from the repository working
directory; logic is unchanged (all scripts pass `parse()`).

| Figure | Script | Notes |
|--------|--------|-------|
| `Figure_1_workflow_overview.png` | *none* | Schematic diagram; no code. |
| `Figure_2_cfu_tvbc.png` | `R/cfu.R` | Total viable count (CFU/mL); Kruskal-Wallis p = 0.0031. |
| `Figure_3_phylogenetic_tree.png` | `R/phylogenetic_tree.R` | **Discrepancy:** this script builds a **maximum-likelihood** tree with `ape`/`phangorn`/`ggtree`; the Methods describe **MEGA 12** with the **K2P** model + FigTree. Requires an input alignment `metagenome_36_seq_maft-align.fasta` (not included). |
| `Figure_4_antibiotic_sensitivity.png` | `R/antibiotic_sensitivity.R` | Disc-diffusion S/I/R summary. |
| `Figure_5_taxonomic_composition.png` | `R/taxonomic_composition.R` | Phylum/Class/Genus/Species stacked bars (needs `output.biom` + `metadata.csv`). |
| `Figure_6_heatmaps_top50.png` | `R/heatmaps.R` | **Note:** uses `ggplot2` `geom_tile`; the Methods name **ComplexHeatmap** (not required here). |
| `Figure_7_alpha_diversity.png` | `R/alpha_diversity.R` | Observed, Chao1, Shannon, Simpson; Kruskal-Wallis. |
| `Figure_8_beta_diversity.png` | `R/beta_diversity.R` | Bray-Curtis PCoA and NMDS (`vegan`). |
| (supporting) | `R/relative_abundance.R` | Relative-abundance computations/plots. |
| (supporting) | `R/taxonomic_composition_count_table.R` | Exports the OTU counts/relative-abundance table (`openxlsx`). |

Inputs required by the R scripts (`output.biom`, `metadata.csv`, and the tree
alignment FASTA) are produced by the pipeline or supplied by the user; large
inputs are not tracked in Git.

## 3. Samples

The paper and **Supplementary Table S4** use canonical labels **CR1, CR2**
(controls) and **TR1, TR2, TR3** (treated). These correspond to the five SRA runs
of **PRJNA1282882** via the run titles (*Control1/2*, *Dose1/2/3*). The internal
labels in the original processing script (`Control_S1`, `Control_S8`,
`SingleAbx_S2`-`S4`) map onto the same five samples, but the exact one-to-one
assignment to CR/TR was **not recorded** and should be confirmed by the authors.
`MultiAbx_S5`-`S7` in the working directory are **not** among the deposited runs
and are excluded.

## 4. Items that could not be verified

- Exact `fastp`, `kraken-biom`, `krona` versions.
- Kraken2 database identity intended by the paper (RefSeq vs `k2_pluspf_20210517`).
- QUAST / Bracken command lines actually used.
- Internal-label to SRA-run (and to CR/TR) correspondence.
- Whether Fig. 3 in print used the MEGA 12 (K2P) tree or the ML tree in
  `R/phylogenetic_tree.R`.
