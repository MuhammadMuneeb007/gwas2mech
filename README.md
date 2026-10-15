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

# 3. Derived resources (1000G PGEN panels, Pangolin DB) and a full status report
python Step00_Check_Resources.py --resources-only
python Step00_Check_Resources.py --inspect
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

GENCODE 50 GRCh38 genome + annotation, 1000 Genomes GRCh38 (2019) VCFs, Ensembl VEP 116 cache,
GTEx v11 eQTL/sQTL (+ SuSiE), eQTL Catalogue metadata, Open Targets Platform 26.09
(study, credible sets, L2G, colocalisation) and GWAS Catalog metadata (~65 GB download).

## Configuration

- SLURM settings: `config/slurm.yaml` (see `gwas2m_config.py`), else the `#SBATCH`
  header of `HPC_CPU.sh`; every script also accepts command-line overrides.
- Resource registry and validation rules: `gwas2m_resources.py`.
