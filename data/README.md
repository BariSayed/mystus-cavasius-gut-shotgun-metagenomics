# Data availability

Raw and large derived files are **not** stored in this repository (see
`.gitignore`). Obtain them from the public archives below.

## Shotgun metagenomic raw reads
- **Archive:** NCBI Sequence Read Archive (SRA)
- **BioProject:** `PRJNA1282882`
- **Runs** (Illumina MiniSeq, 150 bp paired-end), with their SRA run titles:

  | SRA run      | SRA run title                 | Group (from title) |
  |--------------|-------------------------------|--------------------|
  | SRR34283204  | Tengra Microbiome Control1    | Control            |
  | SRR34283203  | Tengra Microbiome Control2    | Control            |
  | SRR34283202  | Tengra Microbiome Dose1       | Treatment          |
  | SRR34283201  | Tengra Microbiome Dose2       | Treatment          |
  | SRR34283200  | Tengra Microbiome Dose3       | Treatment          |

  Download with `data/download_sra.sh`.

> **Note on sample naming.** The analysis scripts and the original committed
> script use internal sample labels (`Control_S1`, `Control_S8`, `SingleAbx_S2`,
> `SingleAbx_S3`, `SingleAbx_S4`). The exact correspondence between these labels
> and the five SRA runs above is **not recorded** in the paper or repository and
> must be confirmed by the authors. See `docs/REPRODUCIBILITY.md`.

## 16S rRNA gene (Sanger) sequences from cultured isolates
- **Archive:** NCBI GenBank
- **Accessions (12):** `PV810321`, `PV810322`, `PV810323`, `PV809930`,
  `PV809931`, `PV809932`, `PV809933`, `PV809934`, `PV809935`, `PV809936`,
  `PV809937`, `PV809938`
- These support the culture-dependent phylogenetic analysis in the paper and are
  not used by the shotgun pipeline in this repository.

## Host reference genome
- *Ictalurus punctatus* (channel catfish), assembly
  `GCF_001660625.3_Coco_2.0`, downloaded by `scripts/00_download_reference.sh`.

## Kraken2 database
- The committed script used the pre-built `k2_pluspf_20210517` database
  (not redistributed here). The published Methods describe the "standard RefSeq"
  database. See `docs/REPRODUCIBILITY.md`.
