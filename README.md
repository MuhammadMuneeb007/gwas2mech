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
| `envs/pipeline` | Python/R stack, GWASLab, plink2, bcftools, samtools, htslib, bedtools, VEP 116, susieR, coloc, hyprcoloc, gseapy, pyranges, pyliftover, CrossMap; figures: ggplot2, patchwork, ggrastr, ggrepel, svglite, ComplexHeatmap, karyoploteR |
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

`Step02_Plan_QC_GWAS.py` works for any phenotype. It searches the GWAS Catalog for the phenotype
as given (plus any `--alias`), keeps studies with full summary statistics, a strict ancestry
match and a harmonised file, and ranks them by sample size. For every study it also records
whether the phenotype is **Binary (cases/controls) or Continuous**, because fine-mapping and
colocalisation differ (SuSiE-RSS / coloc `type="cc"` with case fraction `s`, versus
`type="quant"`). The type is decided **per study**, not per disease: a depression diagnosis GWAS is
Binary, a depression symptom score GWAS is Continuous. The evidence, strongest first:

| Evidence | Source | Points to |
|---|---|---|
| `case_control_study` true / false | summary-statistics `-meta.yaml` (GWAS-SSF) | Binary / Continuous (strong) |
| separate numbers of cases **and** controls | `-meta.yaml`, GWAS Catalog API sample description | Binary (strong) |
| odds ratios | `-meta.yaml` `effect_statistic`, first rows of the harmonised file | Binary (moderate) |
| sample described without cases/controls | GWAS Catalog API / table | Continuous (weak) |
| beta only, and no cases/controls reported | file / `-meta.yaml` | Continuous (weak); beta alone can be log-odds |
| measurement wording (`measurement`, `levels`, `score`, `age at`, ...) | trait / EFO terms | Continuous (weak) |

Strong evidence, or two agreeing weak ones, decides. **Conflicting or insufficient evidence gives
`Unknown` with `NEEDS_REVIEW=True`; it is never guessed.** Flagged studies are listed in
`01_gwas_catalog/<phenotype>/<ancestry>/phenotype_type_review.tsv` with the trait, sample
description, PMID, DOI and evidence. Fill `VERIFIED_PHENOTYPE_TYPE` (`Binary` / `Continuous`;
optionally `VERIFIED_N_CASES`, `VERIFIED_N_CONTROLS`) after reading the paper and re-run with `--qc`;
your verdicts are kept on every re-run. `--phenotype-type binary|continuous` forces one type for
all studies.

Every study gets `PHENOTYPE_TYPE`, `N_CASES`, `N_CONTROLS`, `N_TOTAL`, `N_EFF`, `EFFECT_TYPE`,
`EFFECT_REPORTED`, `CLASSIFICATION_SOURCE`, `CLASSIFICATION_CONFIDENCE` (high / moderate /
conflicting / insufficient / verified / user), `NEEDS_REVIEW`, `CLASSIFICATION_EVIDENCE`, and the
publication `PUBMEDID` and `DOI` (PMID converted with Europe PMC, NCBI E-utilities as fallback;
cached in `reference_metadata/gwas_catalog/pmid_doi.json`).

Missing Python packages (pandas, numpy, requests, PyYAML; scipy and gwaslab 4.2.3 for QC) are
installed automatically into the Python that runs the script, i.e. `envs/pipeline`.

```bash
# 1. discover + download (one light array task per GWAS)
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR
sbatch Step02_Download_GWAS_migraine_european.sh

# 2. GWASLab QC of every GWAS downloaded so far (re-run as more downloads finish)
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --qc
sbatch Step02_QC_migraine_european.sh

# downloaded / partial / pending and QC done / failed / pending, per study (read-only)
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --status

# publication table: authors, journal, year, PMID, trait / EFO, Binary or Continuous,
# N, cases, controls, case fraction, N_eff, population, countries, cohorts, platform,
# imputed, variants raw / after QC, file columns, genome-wide hits, minimum P, lambda GC
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --stats

# PubMed ID -> DOI only
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --pmid-to-doi 34737426 33987465

# study N of the download manifest end to end on this machine: download (if needed) + QC
# (only the planner has to have run; on a login node wrap it in srun)
envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --index 1
```

Options: `--alias "migraine disorder"` (extra search term, repeatable), `--max-studies 5`,
`--refresh` (re-download GWAS Catalog metadata), and `--partition --qos --account --cpus --mem
--time --max-parallel` to override `HPC_CPU.sh`, `--phenotype-type binary|continuous`. Outputs go
to `01_gwas_catalog/<phenotype>/<ancestry>/` and `02_summary_stats/<phenotype>/<ancestry>/{raw,qc}/`.

### Fine-mapping inputs (SuSiE-RSS / FINEMAP)

`qc/<GCST>/<GCST>_GRCh38_QC.tsv.gz` keeps **every variant** (fine-mapping needs all variants in a
locus, not only the significant ones) and starts with these columns, followed by the other
GWASLab columns (`SNPID`, `rsID`, `OR`, `STATUS`, ...):

| Column | Meaning |
|---|---|
| `CHR` `POS` | GRCh38 position |
| `SNP` | rsID, otherwise `CHR:POS:NEA:EA` |
| `EA` `NEA` | effect / non-effect allele (Z and BETA are for `EA`) |
| `BETA` | effect size; **log(OR) for binary traits** (from `OR` when the file has no beta) |
| `SE` | standard error (from the OR 95% CI, or from BETA and P, when missing) |
| `P` | P-value |
| `EAF` `MAF` | effect allele frequency; MAF = min(EAF, 1 − EAF) |
| `N` | per-variant N; when missing, the study `N_TOTAL` (never N_eff) |
| `INFO` | imputation quality (blank when not reported) |
| `Z` | BETA / SE, otherwise signed from P and the effect direction |

How many values were derived (`N_BETA_FROM_OR`, `N_SE_FROM_OR_CI`, `N_Z_FROM_P`, `N_SOURCE`,
`N_MISSING_Z`, ...) is recorded in `<GCST>_QC_summary.tsv`.

The study metadata is in `qc/<GCST>/<GCST>_study_metadata.tsv` (and `.json`), and for all studies in
`01_gwas_catalog/<phenotype>/<ancestry>/GWAS_study_metadata.tsv` (updated by every Step02 run):

| Column | Meaning |
|---|---|
| `GCST_ID` | GWAS Catalog accession |
| `PHENOTYPE_TYPE` | `Binary` or `Continuous` (`Unknown` when the catalogue does not say) |
| `N_CASES` `N_CONTROLS` | binary studies only |
| `N_TOTAL` | total sample size |
| `N_EFF` | binary only: 4·N_cases·N_controls / (N_cases + N_controls) |
| `ANCESTRY` | ancestry code (EUR, EAS, ...) |
| `GENOME_BUILD` | `GRCh38` (harmonised files) |
| `EFFECT_TYPE` | `log(OR)` for binary, `BETA` for continuous |

Provenance columns follow (`CASE_FRACTION`, `EFFECT_REPORTED`, `CLASSIFICATION_SOURCE`,
`CLASSIFICATION_CONFIDENCE`, `NEEDS_REVIEW`, `CLASSIFICATION_EVIDENCE`, `PMID`, `DOI`, `N_SOURCE`,
`QC_STATUS`, `LD_REFERENCE`). Step06 skips studies with `PHENOTYPE_TYPE=Unknown` or
`NEEDS_REVIEW=True`. Both `N_TOTAL` and `N_EFF` are kept; each fine-mapping method chooses which one
it uses. The LD matrix comes from the ancestry-matched panel in `LD_REFERENCE`
(`resources/1000G/<ANCESTRY>`, Step04), with the same variants, order and allele orientation as
the GWAS (Step06). Studies whose sample includes **proxy cases** (family history, GWAX) are flagged
for review: their effects are attenuated and N_EFF is only approximate.

## 5. LD reference panels (Step04)

Ancestry-specific PLINK2 panels from the 1000 Genomes GRCh38 (2019) VCFs downloaded by
`--groups 1000g`. They are shared by all phenotypes and built once per ancestry into
`resources/1000G/<ANC>/chr<N>_<ANC>_GRCh38.pgen/.pvar/.psam`.

```bash
envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR       # or: EUR EAS / all
sbatch Step04_Prepare_1000G_EUR_GRCh38.sh
envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --status
envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --stats  # people per population, variants
```

## 6. LD clumping (Step05)

Independent lead variants per study with PLINK2 `--clump` (p1 5e-8, p2 1e-2, r2 0.1, 1000 kb)
on the ancestry-matched panel. GWAS variants are matched to the panel by chromosome, GRCh38
position **and both alleles**; a different variant at the same position is never used.

```bash
envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype migraine --ancestry EUR
sbatch Step05_Clump_GWAS_migraine_european.sh
envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype migraine --ancestry EUR --status
envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype migraine --ancestry EUR --stats   # leads, loci, mapping, lost GWS variants
envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype migraine --ancestry EUR --plots   # Manhattan, QQ, chromosome locus map
```

## 7. Fine-mapping (Step06, SuSiE-RSS)

Per study: leads ±1 Mb merged into loci; every GWAS variant of the locus; Z oriented to the 1000G
REF allele; a new signed LD matrix (`plink2 --r-unphased ref-based`) on exactly those variants in the
same order; `susieR::susie_rss(z, R, n)` with `estimate_residual_variance = FALSE`.

N comes from Step02 (`GWAS_study_metadata.tsv`): **Binary → N_EFF** (`--susie-n total` for N_TOTAL),
**Continuous → N_TOTAL**; the N used and why is recorded per study.

```bash
envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR
sbatch Step06_FineMap_GWAS_migraine_european.sh
envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR --status
envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR --stats   # loci, credible sets, PIPs, locus QC
envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR --plots   # regional + PIP + genes, LD heatmap, genome PIP
```

Every locus is checked and gets `LOCUS_QC = PASS / WARN` with the reasons (`LOCUS_QC_NOTES`):

| Check | What it detects |
|---|---|
| `ESTIMATE_S_RSS` (susieR) | z-scores inconsistent with the reference LD (0 = consistent; > 0.1 warn, > 0.3 strong) |
| `kriging_rss` | variants whose z disagrees with their LD neighbours (possible allele flips) |
| `CONVERGED` | SuSiE did not converge within `--max-iter` (500) |
| credible-set purity | min \|r\| inside a credible set < 0.5 |
| `N_HIGH_PIP_OUTSIDE_CS` | PIP > 0.5 outside every credible set (typical LD-reference artefact) |
| `TOP_AGREES_L1` | top variant also top in the L = 1 fit, which does not use the LD matrix |
| LD cap / N spread / leads not in LD | locus trimmed to `--max-ld-variants`, per-variant N varies > 2-fold, Step05 lead not in the LD set |

**LD reference size.** 1000G EUR has ~503 people. For very large GWAS (z-scores of 15–25) its
sampling noise alone gives `ESTIMATE_S_RSS` around 0.5 and can create spurious PIPs — the checks
above flag those loci. Interpret WARN loci with care and, where available, repeat with a larger
ancestry-matched panel or in-sample LD. Rare variants (MAF < 1% in the panel, e.g. many
high-penetrance coding variants) are not fine-mapped with this panel.

## 8. Publication figures (Step07)

Reads Step02/05/06 outputs (nothing is recomputed) and draws every figure in R; dense point
layers are rasterised, text and gene models stay vector.

```bash
envs/pipeline/bin/python Step07_Plots.py --phenotype migraine --ancestry EUR                    # pdf + png
envs/pipeline/bin/python Step07_Plots.py --phenotype migraine --ancestry EUR --formats pdf,svg,png
envs/pipeline/bin/python Step07_Plots.py --phenotype migraine --ancestry EUR --status
envs/pipeline/bin/python Step07_Plots.py --phenotype migraine --ancestry EUR --slurm            # then sbatch the .sh
```

| Figure (`07_plots/<phenotype>/<ancestry>/`) | Content | Package |
|---|---|---|
| `<GCST>/manhattan` | genome-wide association, independent leads labelled with the nearest gene, fine-mapping windows shaded, λGC | ggplot2 + ggrastr + ggrepel |
| `<GCST>/ideogram` | chromosomes 1–22 (hg38 cytobands) with every locus | karyoploteR |
| `<GCST>/loci/<L>_region` | one genomic axis: −log10 P coloured by r² with the index variant · SuSiE PIP with 95% credible sets · gene models (canonical GENCODE transcripts) · LD r² triangle; title carries the locus QC (s, allele flips, artefacts) | ggplot2 + patchwork |
| `<GCST>/all_loci_multipage.pdf` | every locus figure in one PDF | |
| `combined/ideogram_all_gwas` | all loci of all GWAS, one coloured band per GWAS | karyoploteR |
| `combined/locus_by_gwas` | loci × GWAS matrix: min P and max PIP (number of credible sets) | ComplexHeatmap |
| `combined/shared_loci/<chr>_<pos>_<gene>` | loci found in ≥ 2 GWAS: one association panel per GWAS on a shared axis, PIPs of all GWAS, genes | ggplot2 + patchwork |
| `combined/locus_summary.tsv` | the table behind the matrix | |

Options: `--max-loci N` (locus figures per study), `--ld-max 120` (variants in the LD triangle),
`--ld-window-kb 500` (largest pair distance shown), `--dpi 300`. Gene tracks use the local GENCODE
GTF from `--groups gencode` (offline, same release as the rest of the pipeline). Studies that were
not fine-mapped (e.g. NEEDS_REVIEW) still get a Manhattan plot and ideogram from their Step05 leads.

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
