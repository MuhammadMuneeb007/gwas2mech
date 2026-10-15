# GWAS2Mech

From GWAS summary statistics to candidate molecular mechanisms:
QC → ancestry-specific LD → clumping → SuSiE fine-mapping → VEP → GTEx / eQTL Catalogue
QTLs → SpliceAI / Pangolin → colocalisation → gene prioritisation (MAGMA, PoPS) →
Open Targets / FUMA interpretation.

Everything is installed and downloaded **inside the directory you run the scripts from**,
so the same commands work on a laptop, a workstation or any SLURM cluster.
Run every command from the repository directory.

## 1. Install the software

Requirements: Linux (or macOS), Python ≥ 3.7, and `mamba`, `conda` or `micromamba`
on `PATH` (on a cluster e.g. `module load miniforge`).

The installer is the only script started with your system `python`, because
`envs/pipeline` does not exist yet:

```bash
python Step00_Installation.py                          # install + verify everything
python Step00_Installation.py --verify-only            # check only
python Step00_Installation.py --only pipeline pangolin # repair selected environments
```

Every later command uses the project's own interpreter, **`envs/pipeline/bin/python`**.

| Location | Contents |
|---|---|
| `envs/pipeline` | Python/R stack, GWASLab, plink2, bcftools, samtools, htslib, bedtools, VEP 116, susieR, coloc, hyprcoloc, gseapy, pyranges, pyliftover, CrossMap |
| `envs/spliceai` | SpliceAI 1.3.1 (TensorFlow 2.15) |
| `envs/pangolin` | Pangolin (pinned commit) |
| `envs/ldsc` | LDSC / S-LDSC 1.0.1 (Python 2.7) |
| `envs/pops` | PoPS (pinned commit) |
| `external_tools/magma` | MAGMA v1.10 (official binary, SHA-256 verified) |

## 2. Download the shared resources

`Step00_Download_Resources.py` plans the files that are still missing and writes a SLURM
script as **one SLURM array of at most 999 tasks + 1 final check = never more than 1,000 jobs**.
Up to 999 files, each task downloads one file; with more files, each task downloads several
files one after another. Every job uses the `#SBATCH` settings of `HPC_CPU.sh`.

```bash
# how many files are downloaded / partial / missing, per group (read-only, any time)
envs/pipeline/bin/python Step00_Download_Resources.py --status
```

### Everything at once

```bash
envs/pipeline/bin/python Step00_Download_Resources.py --skip-groups meqtl
sbatch Step00_Download_Resources.sh

# the 1.47 TB brain meQTL archive as its own long job
envs/pipeline/bin/python Step00_Download_Resources.py --groups meqtl --time 7-00:00:00 --output Step00_Download_meQTL.sh
sbatch Step00_Download_meQTL.sh
```

### One resource at a time

Each command writes its own job script (`--output`), so you can plan several and submit
them independently. Then submit with `sbatch <that script>`.

**Reference genome and annotation (GENCODE 50, GRCh38)** — 0.9 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups genome --output Step00_Download_Genome.sh
sbatch Step00_Download_Genome.sh
```

**1000 Genomes reference panels (GRCh38, 2019 release)** — 13 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups 1000g --output Step00_Download_1000G.sh
sbatch Step00_Download_1000G.sh
```

**Ensembl VEP 116 cache** — 26 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups vep --output Step00_Download_VEP.sh
sbatch Step00_Download_VEP.sh
```

**GWAS Catalog metadata** — 0.2 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups gwas_catalog --output Step00_Download_GWASCatalog.sh
sbatch Step00_Download_GWASCatalog.sh
```

**GTEx v11 eQTL / sQTL / apaQTL (+ SuSiE)** — 7 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups gtex --output Step00_Download_GTEx.sh
sbatch Step00_Download_GTEx.sh
```

**eQTL Catalogue QTLs used by Step11** (GTEx v8 eQTL / exonQTL / isoQTL / sQTL + Sun_2018 pQTL) — 320 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups eqtl_catalogue qtl_catalogue --output Step00_Download_QTLCatalogue.sh
sbatch Step00_Download_QTLCatalogue.sh
```

Narrower or wider eQTL Catalogue selections:
```bash
# only splicing QTL (leafcutter sQTL)
envs/pipeline/bin/python Step00_Download_Resources.py --groups qtl_catalogue --catalogue-qtl-types sQTL --output Step00_Download_sQTL.sh
# only expression QTL
envs/pipeline/bin/python Step00_Download_Resources.py --groups qtl_catalogue --catalogue-qtl-types eQTL --output Step00_Download_eQTL.sh
# only the Sun_2018 plasma pQTL
envs/pipeline/bin/python Step00_Download_Resources.py --groups qtl_catalogue --catalogue-study-regex '^Sun_2018$' --output Step00_Download_pQTL_Sun2018.sh
# every eQTL Catalogue study (very large)
envs/pipeline/bin/python Step00_Download_Resources.py --groups qtl_catalogue --catalogue-all-studies --output Step00_Download_QTLCatalogueAll.sh
# SuSiE lbf / credible sets for the Step11 contexts (316 GB)
envs/pipeline/bin/python Step00_Download_Resources.py --groups qtl_catalogue_susie --output Step00_Download_QTLSuSiE.sh
```

**Cell-type eQTL (Jerber_2021, as Setup_11B)** — 62 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups qtl_celltype --output Step00_Download_CellTypeQTL.sh
sbatch Step00_Download_CellTypeQTL.sh
# another study:  --celltype-study-regex '^<study_label>$'   (QTL classes: --celltype-qtl-types eQTL,sQTL)
```

**Brain pQTL / sQTL / haQTL (ADSP FunGen, NIAGADS NG00184) + caQTL (Kumasaka)** — 230 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups molecular_qtl --output Step00_Download_MolecularQTL.sh
sbatch Step00_Download_MolecularQTL.sh
```

**Brain meQTL / mQTL (ADSP FunGen, NIAGADS NG00184)** — 1.47 TB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups meqtl --time 7-00:00:00 --output Step00_Download_meQTL.sh
sbatch Step00_Download_meQTL.sh
```

**Brain bulk + single-nucleus eQTL (ADSP FunGen)** — 720 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups adsp_eqtl --time 7-00:00:00 --output Step00_Download_ADSP_eQTL.sh
sbatch Step00_Download_ADSP_eQTL.sh
```

**Open Targets Platform 26.09** (study, credible sets, L2G, colocalisation) — 19 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups opentargets --output Step00_Download_OpenTargets.sh
sbatch Step00_Download_OpenTargets.sh
```

**Variant annotation** (AlphaMissense, ClinVar, ENCODE SCREEN cCREs v4, gnomAD v4.1 constraint, UCSC liftover chains, phyloP100way) — 11 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups annotation --output Step00_Download_Annotation.sh
sbatch Step00_Download_Annotation.sh
```

**CADD v1.7 (GRCh38 SNVs + gnomAD4 indels)** — 89 GB
```bash
envs/pipeline/bin/python Step00_Download_Resources.py --groups cadd --output Step00_Download_CADD.sh
sbatch Step00_Download_CADD.sh
```

**Specific eQTL Catalogue datasets by id** (e.g. to match an existing project)
```bash
printf 'QTD000584\nQTD000141\n' > setup_manifests/eqtl_catalogue_datasets.txt   # one dataset id per line
envs/pipeline/bin/python Step00_Download_Resources.py --groups eqtl_catalogue --output Step00_Download_Datasets.sh
sbatch Step00_Download_Datasets.sh
```

### Useful options

```bash
envs/pipeline/bin/python Step00_Download_Resources.py --help
#   --skip-groups G ...        everything except these groups
#   --no-url-check             plan without checking URLs / sizes (faster)
#   --time / --mem / --cpus    override HPC_CPU.sh for this plan
#   --max-jobs N --parallel N  jobs limit (tasks + final check) / concurrently running tasks
#   --status                   progress of every plan: done / partial / missing per group
bash Step00_Download_Resources.sh local 4    # without SLURM: 4 files at a time
```

Re-running a plan is safe: valid files are skipped and only missing, partial or corrupt files
are planned again. Very large archives (molecular QTL) are MD5-verified, extracted, then deleted.
REVEL cannot be downloaded automatically (its servers block scripted downloads).

## 3. Check the downloads

```bash
# quick progress of all plans (read-only)
envs/pipeline/bin/python Step00_Download_Resources.py --status

# per plan: verify every file, write the report (also run automatically by the last job)
envs/pipeline/bin/python Step00_Download_Resources.py --finalize
envs/pipeline/bin/python Step00_Download_Resources.py --finalize --output Step00_Download_meQTL.sh

# every planned file: readable? not corrupt? print its header
envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py
envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py --group qtl_catalogue
envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py --full     # read every file to the end (submit as a job)
```

## 4. Find, download and QC the GWAS for a phenotype

`Step02_Plan_QC_GWAS.py` searches the GWAS Catalog (full summary statistics, strict ancestry
match, harmonised file, ranked by sample size) and writes one SLURM array in which **each task
downloads one GWAS and runs GWASLab QC on it**.

```bash
# plan + submit
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR
sbatch Step02_GWAS_migraine_european.sh

# what is downloaded / pending, and QC done / failed / pending, per study
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --status

# publication table: authors, journal, year, PMID, trait (EFO), population and countries,
# total N, cases, controls, effective N, platform, variants raw / after QC, file columns,
# genome-wide significant hits, minimum P, lambda GC, available fields
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --stats

# run one task (download + QC) on this machine
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --index 1
```

Options: `--alias "migraine disorder"` (extra search term), `--max-studies 5`, `--refresh`
(re-download GWAS Catalog metadata), and `--partition --qos --account --cpus --mem --time
--parallel` to override `HPC_CPU.sh`. Outputs go to `01_gwas_catalog/<phenotype>/<ancestry>/`
(manifests, `study_selection.tsv`, `GWAS_publication_stats.tsv`) and
`02_summary_stats/<phenotype>/<ancestry>/{raw,qc}/`.

## Resource groups

| Group | Contents | Approx. download |
|---|---|---|
| `gwas_catalog` | GWAS Catalog studies + ancestry tables | 0.2 GB |
| `genome` | GENCODE 50 GRCh38 FASTA (+ `.fai`) and GTF | 0.9 GB |
| `1000g` | 1000 Genomes GRCh38 (2019) VCF + TBI + sample panel | 13 GB |
| `vep` | Ensembl VEP 116 GRCh38 cache (extracted) | 26 GB |
| `gtex` | GTEx v11 eQTL / sQTL / apaQTL + SuSiE (extracted), variant lookup, GENCODE 47 | 7 GB |
| `eqtl_catalogue` | eQTL Catalogue metadata tables (+ datasets listed by id) | < 1 GB |
| `opentargets` | Open Targets Platform 26.09 study / credible sets / L2G / colocalisation | 19 GB |
| `qtl_catalogue` | eQTL Catalogue dense + `.tbi` for Step11: GTEx v8 eQTL / exonQTL / isoQTL / sQTL + Sun_2018 pQTL | 320 GB |
| `qtl_celltype` | Jerber_2021 cell-type eQTL, dense + `.tbi` + SuSiE | 62 GB |
| `molecular_qtl` | ADSP FunGen brain pQTL / sQTL / haQTL + Kumasaka caQTL, MD5-verified | 230 GB |
| `meqtl` | ADSP FunGen brain mQTL, MD5-verified | 1.5 TB |
| `annotation` | AlphaMissense, ClinVar, ENCODE SCREEN cCREs v4, gnomAD v4.1 constraint, UCSC liftover chains, phyloP100way | 11 GB |
| `cadd` | CADD v1.7 GRCh38 SNVs + gnomAD4 indels (+ `.tbi`) | 89 GB |
| `qtl_catalogue_susie` *(only when named)* | eQTL Catalogue SuSiE lbf / cs for the Step11 contexts | 316 GB |
| `adsp_eqtl` *(only when named)* | ADSP FunGen brain bulk + single-nucleus eQTL | 720 GB |

Files are stored under `resources/` (dense eQTL Catalogue files in
`resources/coloc/eqtl_catalogue/dense/<dataset_id>/`, other molecular QTL in
`resources/molecular_qtl/<type>/<source>/`, annotation in `resources/annotation/<source>/`).

## Configuration

SLURM settings come from the `#SBATCH` header of `HPC_CPU.sh` (account, partition, qos,
nodes, ntasks, cpus-per-task, mem, time) and apply to every generated job. Edit that file
for your cluster, or override a single plan with `--account`, `--partition`, `--qos`,
`--cpus`, `--mem` or `--time`.
