# GWAS2Mech

From GWAS summary statistics to candidate molecular mechanisms:
QC → ancestry-specific LD → clumping → SuSiE fine-mapping → VEP → GTEx / eQTL Catalogue
QTLs → SpliceAI / Pangolin → colocalisation → gene prioritisation (MAGMA, PoPS) →
Open Targets / FUMA interpretation.

Everything is installed and downloaded **inside the directory you run the scripts from**,
so the same commands work on a laptop, a workstation or any SLURM cluster.

## Setup

Requirements: Linux (or macOS), Python ≥ 3.7, and `mamba`, `conda` or `micromamba`
on `PATH` (on a cluster e.g. `module load miniforge`).

```bash
# 1. Software: conda environments under envs/, tools under external_tools/
python Step00_Installation.py              # install + verify everything
python Step00_Installation.py --verify-only

# 2. Shared reference data: one SLURM array task per file (max 1000 queued jobs)
python Step00_Download_Resources.py        # plans what is missing, writes the job script
sbatch Step00_Download_Resources.sh        # or without SLURM: bash Step00_Download_Resources.sh local 4
python Step00_Download_Resources.py --finalize   # re-check every file + report

# 3. Check every downloaded file is readable and print its header
envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py          # add --full for a complete read
```

## Software installed by `Step00_Installation.py`

| Location | Contents |
|---|---|
| `envs/pipeline` | Python/R stack, GWASLab, plink2, bcftools, samtools, htslib, bedtools, VEP 116, susieR, coloc, hyprcoloc, gseapy, pyranges, pyliftover, CrossMap |
| `envs/spliceai` | SpliceAI 1.3.1 (TensorFlow 2.15) |
| `envs/pangolin` | Pangolin (pinned commit) |
| `envs/ldsc` | LDSC / S-LDSC 1.0.1 (Python 2.7) |
| `envs/pops` | PoPS (pinned commit) |
| `external_tools/magma` | MAGMA v1.10 (official binary, SHA-256 verified) |

## Shared resources downloaded by `Step00_Download_Resources.py`

One SLURM array task per file; pick groups with `--groups` / `--skip-groups`.

| Group | Contents | Approx. download |
|---|---|---|
| `gwas_catalog` | GWAS Catalog studies + ancestry tables | 0.2 GB |
| `genome` | GENCODE 50 GRCh38 FASTA (+ `.fai`) and GTF | 0.9 GB |
| `1000g` | 1000 Genomes GRCh38 (2019) VCF + TBI + sample panel | 13 GB |
| `vep` | Ensembl VEP 116 GRCh38 cache (extracted) | 26 GB |
| `gtex` | GTEx v11 eQTL / sQTL / apaQTL + SuSiE (extracted), variant lookup, GENCODE 47 | 7 GB |
| `eqtl_catalogue` | eQTL Catalogue metadata tables | < 1 GB |
| `opentargets` | Open Targets Platform 26.09 study / credible sets / L2G / colocalisation | 19 GB |
| `qtl_catalogue` | eQTL Catalogue dense + `.tbi` for Step11: GTEx v8 eQTL / exonQTL / isoQTL / sQTL + Sun_2018 pQTL | 320 GB |
| `qtl_celltype` | Jerber_2021 cell-type eQTL, dense + `.tbi` + SuSiE (Setup_11B) | 62 GB |
| `molecular_qtl` | ADSP FunGen brain pQTL / sQTL / haQTL (NIAGADS NG00184) + Kumasaka caQTL, MD5-verified | 230 GB |
| `meqtl` | ADSP FunGen brain mQTL (NIAGADS NG00184), MD5-verified | 1.5 TB |
| `annotation` | AlphaMissense, ClinVar, ENCODE SCREEN cCREs v4, gnomAD v4.1 constraint, UCSC liftover chains, phyloP100way | 11 GB |
| `cadd` | CADD v1.7 GRCh38 SNVs + gnomAD4 indels (+ `.tbi`) | 89 GB |
| `qtl_catalogue_susie` *(opt-in)* | eQTL Catalogue SuSiE lbf / cs for the Step11 contexts | 316 GB |
| `adsp_eqtl` *(opt-in)* | ADSP FunGen brain bulk + single-nucleus eQTL | 720 GB |

Very large archives (molecular QTL) are deleted after a verified extraction. The brain mQTL
archive alone is 1.47 TB; it can be planned as its own job with a longer walltime:
`python Step00_Download_Resources.py --groups meqtl --time 7-00:00:00 --output Step00_Download_meQTL.sh`.
REVEL cannot be downloaded automatically (its servers block scripted downloads).

## Configuration

- SLURM settings: the `#SBATCH` header of `HPC_CPU.sh` (partition, qos, account, time);
  `Step00_Download_Resources.py` also accepts `--partition`, `--qos`, `--account`,
  `--cpus`, `--mem`, `--time`, `--max-jobs` and `--parallel`.
