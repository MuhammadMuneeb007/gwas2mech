#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2m STEP 00 - DOWNLOAD SHARED RESOURCES (SLURM: ONE ARRAY TASK PER FILE)
===============================================================================

    envs/pipeline/bin/python Step00_Download_Resources.py        # plan + write the SLURM script
    sbatch Step00_Download_Resources.sh        # download everything

    # 1000 Genomes phase 3 LD panel inputs (Step04 --reference-panel 1000G, the default)
    envs/pipeline/bin/python Step00_Download_Resources.py --groups 1000g --output Step00_Download_1000G.sh
    sbatch Step00_Download_1000G.sh

    # 1000 Genomes 30x LD panel inputs (Step04 --reference-panel 1000G30x; never part
    # of the default plan): 22 phased VCFs + TBI and the 3202-sample population file,
    # ~28 GB, one array task per file
    envs/pipeline/bin/python Step00_Download_Resources.py --groups 1000g30x --output Step00_Download_1000G30x.sh
    sbatch Step00_Download_1000G30x.sh
    envs/pipeline/bin/python Step00_Download_Resources.py --status

This script only PLANS; it never downloads data itself. It:

  1. lists every shared resource FILE the pipeline needs (paths and URLs
     are defined in this script);
  2. skips files that are already VALID - non-empty, no <file>.aria2 control
     file and, when the remote size is known, exactly that size - and archives
     already extracted / FASTA already indexed (re-running only re-plans what
     is missing, partial or corrupt);
  3. checks every URL and records its size (the size is enforced after
     download);
  4. writes setup_manifests/download_resources.tsv (one row = one file) and
     Step00_Download_Resources.sh.

SLURM settings: EVERY generated job (download tasks and the final check) uses
the #SBATCH header of
HPC_CPU.sh - account, partition, qos, nodes, ntasks, cpus-per-task, mem, time.
Only --job-name / --output / --error / --array are set per job. Highest
priority first:
     command-line options (--account --partition --qos --cpus --mem --time)
  >  the #SBATCH header of HPC_CPU.sh (or --template FILE)
  >  built-in defaults (2 CPUs, 16G, 24:00:00)

Jobs: never more than --max-jobs (default 1000) - ONE array of at most
(max-jobs - 1) tasks plus one final check (python ... --finalize). Up to 999
files: one task per file. More files: each task downloads several files one
after another (task k gets rows k, k+999, k+1998, ... so large files are spread
over different tasks).

Progress at any time (read-only): python ... --status

Per task: resumable download (aria2c, else curl/wget into <file>.part) ->
exact-size check -> integrity check (gzip -t / tar -t / parquet footer / MD5)
-> post-processing (FASTA decompress + samtools faidx; tar extraction for
GTEx, the VEP cache and the molecular-QTL archives - the very large ones are
deleted after a verified extraction).

A task does NOT end before its file is complete: after a dropped or stalled
connection it resumes the download again and again (waits growing to 10 min),
a file that fails its integrity check is downloaded again (3 attempts), and
when the job's walltime is nearly over the running download is stopped cleanly
(resumable) and the task submits a continuation job for the same file (at most
10 in a row). It only gives up when the server says the file does not exist
(HTTP 404 / 410). If the server's file changed size since planning, the new
size is used.

Resource groups (--groups / --skip-groups); sizes are approximate downloads:

  gwas_catalog        GWAS Catalog studies + ancestry tables                  ~0.2 GB
  genome              GENCODE 50 GRCh38 FASTA (+ .fai) and GTF                ~0.9 GB
  1000g               1000 Genomes GRCh38 VCF + TBI + panel (skipped once the
                      population PGEN panels under resources/1000G/ exist)    ~13 GB
                      (Step04 --reference-panel 1000G, the default)
  vep                 Ensembl VEP 116 GRCh38 indexed cache (extracted)        ~26 GB
  gtex                GTEx v11 eQTL / sQTL / apaQTL + SuSiE (extracted),
                      variant lookup, GENCODE 47                              ~7 GB
  eqtl_catalogue      eQTL Catalogue metadata tables + repair of dense / SuSiE
                      files referenced locally or listed in
                      setup_manifests/eqtl_catalogue_datasets.txt
  opentargets         Open Targets Platform study / credible_set /
                      l2g_prediction / colocalisation (one task per parquet)  ~19 GB
  qtl_catalogue       eQTL Catalogue dense .tsv.gz + .tbi for the Step09
                      contexts (GTEx v8 eQTL / exonQTL / isoQTL / sQTL
                      (leafcutter) + Sun_2018 pQTL)                           ~320 GB
  qtl_celltype        cell-type eQTL for Step10 (default Jerber_2021;
                      --celltype-study-regex for more), dense + .tbi + SuSiE
                      credible sets + lbf (lbf sorted + tabix-indexed)        ~62 GB
  molecular_qtl       ADSP FunGen xQTL brain pQTL / sQTL / haQTL (NIAGADS
                      NG00184, GRCh38) + Kumasaka caQTL (GRCh38), MD5-checked ~230 GB
  meqtl               ADSP FunGen xQTL brain mQTL / meQTL (NG00184, GRCh38)   ~1.5 TB
  annotation          AlphaMissense, ClinVar, ENCODE SCREEN cCREs v4,
                      gnomAD v4.1 constraint, UCSC liftover chains, phyloP,
                      RegulomeDB v2.2 scores (dbSNP153 common SNVs)           ~12 GB
  cadd                CADD v1.7 GRCh38 SNVs + gnomAD4 indels (+ .tbi)         ~89 GB
  Only when named in --groups:
  qtl_catalogue_susie eQTL Catalogue SuSiE credible sets + lbf for the Step09
                      contexts; every lbf file is position-sorted + tabix-
                      indexed after download (read per locus)                ~316 GB
  1000g30x            1000 Genomes 30x on GRCh38 (NYGC, 3202 samples): the
                      filtered phased SNV/INDEL/SV panel per chromosome + TBI
                      and the 3202-sample pedigree / population file, for
                      Step04 --reference-panel 1000G30x (skipped once its
                      population panels under resources/LD_panels/1000G30x/
                      exist)                                                  ~28 GB
  adsp_eqtl           ADSP FunGen brain bulk eQTL + single-nucleus eQTL       ~720 GB

Layouts match the pipeline: dense QTL files under
resources/coloc/eqtl_catalogue/dense/<dataset_id>/ (as Step09
expect), other molecular QTL under resources/molecular_qtl/<type>/<source>/,
annotation under resources/annotation/<source>/.

Step08 (local variant annotation) uses: vep (VEP 116 cache), genome (GRCh38
FASTA + GENCODE GTF), annotation (AlphaMissense - tabix-indexed here -, ClinVar,
ENCODE cCRE V4, gnomAD v4.1 constraint / LOEUF, phyloP 100-way, RegulomeDB v2.2
- sorted, bgzipped and tabix-indexed here -) and cadd (CADD v1.7 SNVs + indels):
    envs/pipeline/bin/python Step00_Download_Resources.py --groups vep genome annotation cadd --output Step00_Download_Step08.sh

NOT done here: derived resources (1000G PGEN panels, sample lists, Pangolin
DB), REVEL (its server blocks automated downloads), SpliceAI precomputed
scores (login required; SpliceAI / Pangolin run locally instead), and
phenotype GWAS summary statistics.

After the downloads: envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py

Usage:

    envs/pipeline/bin/python Step00_Download_Resources.py
    envs/pipeline/bin/python Step00_Download_Resources.py --groups genome vep gtex
    envs/pipeline/bin/python Step00_Download_Resources.py --skip-groups meqtl cadd
    envs/pipeline/bin/python Step00_Download_Resources.py --groups qtl_catalogue --catalogue-study-regex '^Sun_2018$'
    envs/pipeline/bin/python Step00_Download_Resources.py --groups meqtl --time 7-00:00:00 --output Step00_Download_meQTL.sh
    envs/pipeline/bin/python Step00_Download_Resources.py --max-jobs 500 --parallel 50
    envs/pipeline/bin/python Step00_Download_Resources.py --account a_myproject
    envs/pipeline/bin/python Step00_Download_Resources.py --finalize      # re-check + report

Without SLURM (laptop / workstation):

    bash Step00_Download_Resources.sh local 4           # 4 files at a time

===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


# =============================================================================
# ROOT
# =============================================================================

ROOT = Path.cwd().resolve()

MANIFESTS = ROOT / "setup_manifests"
LOGS = ROOT / "logs" / "resource_setup"
DEFAULT_OUTPUT = "Step00_Download_Resources.sh"


# =============================================================================
# VERSIONS / SOURCES
# =============================================================================

GENCODE_RELEASE = "50"
VEP_RELEASE = "116"
GTEX_RELEASE = "v11"
OPEN_TARGETS_RELEASE = "26.09"
LD_ANCESTRIES = ("EUR", "AFR", "EAS", "SAS", "AMR")
# Groups planned as one set: if any file is missing, all of them get a task.
COMPLETE_SET_GROUPS = {"1000g", "1000g30x"}

GENCODE_BASE = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}"
ONE_K_BASE = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
              "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL")
# 1000 Genomes 30x (NYGC; Byrska-Bishop et al. 2022, Cell): filtered phased panel.
ONE_K_30X_BASE = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
                  "1000G_2504_high_coverage")
ONE_K_30X_PHASED = f"{ONE_K_30X_BASE}/working/20220422_3202_phased_SNV_INDEL_SV"
PHASE3_PANEL_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
                    "integrated_call_samples_v3.20130502.ALL.panel")
GTEX_BASE = "https://storage.googleapis.com/adult-gtex"
GWAS_CATALOG_BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest"
EQTL_RAW = "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master"
OT_HTTPS = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"
OT_DATASETS = ("study", "credible_set", "l2g_prediction", "colocalisation")

EQTL_ROOT = Path("resources") / "coloc" / "eqtl_catalogue"
EQTL_METADATA = {
    "tabix_ftp_paths.tsv": f"{EQTL_RAW}/tabix/tabix_ftp_paths.tsv",
    "tabix_ftp_paths_imported.tsv": f"{EQTL_RAW}/tabix/tabix_ftp_paths_imported.tsv",
    "dataset_metadata_r8_beta.tsv": f"{EQTL_RAW}/data_tables/dataset_metadata_r8_beta.tsv",
    "population_assignments.tsv": f"{EQTL_RAW}/data_tables/population_assignments.tsv",
}
EXTRA_EQTL_DATASETS = MANIFESTS / "eqtl_catalogue_datasets.txt"
PARTIAL_SUFFIXES = (".aria2", ".part", ".tmp")

# --- Molecular QTL: ADSP FunGen xQTL Atlas (NIAGADS NG00184.v1, brain, GRCh38) ---
# NIAGADS file id -> file; MD5s from the NG00184.v1 file manifest.
NIAGADS = "https://st1.niagads.org/portal/v1/download-public/NG00184/fileset"
ADSP_ARCHIVES = [
    # qtl, group, (manifest id, md5), (archive id, md5)
    ("pQTL", "molecular_qtl", (57082, "04a815b631ed31b841a82209dd2bca4b"), (57083, "383d52607d05f007a74012f20e619a6f")),
    ("sQTL", "molecular_qtl", (57092, "ddaea5121cc43f1c3311324795ad7bec"), (57093, "041b7360ab415504f68a88c605e154cb")),
    ("haQTL", "molecular_qtl", (57062, "5b31f3d3a33794913a291603f7bf81b6"), (57063, "aab954ce983522608de12fa93376350c")),
    ("mQTL", "meqtl", (57072, "6dd59e7c712437b0995f6a6bc1012515"), (57073, "8c0b98d5fd35bb7c534f17dc52f4fc8f")),
    ("eQTL", "adsp_eqtl", (57052, "cf39475b1078b3151fe9f3da6438cecc"), (57053, "e410de967d14e9268a256d8c6f4185d8")),
    ("snuc-eQTL", "adsp_eqtl", (57102, "2088e99e091aa4a523316fd7679bf17c"), (57103, "ed248beca9f0b2c72b0f401ec77b8012")),
]
KUMASAKA = "https://zenodo.org/records/13848268/files"   # caQTL, GRCh38 re-analysis (QTD100018)

# --- Annotation resources ---
ALPHAMISSENSE = "https://storage.googleapis.com/dm_alphamissense"
CLINVAR = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
SCREEN_CCRES = "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed"
# RegulomeDB v2.2 precalculated rank / probability for dbSNP v153 common SNVs (MAF >= 0.01),
# GRCh38, ENCODE annotation ENCSR938LSO (Boyle lab). Plain TSV, 0-based BED-like coordinates.
REGULOMEDB = "https://www.encodeproject.org/files/ENCFF250UJY/@@download/ENCFF250UJY.tsv"
REGULOMEDB_MD5 = "a00d7cecb907dc60007d7b021e714f38"
GNOMAD_CONSTRAINT = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/"
                     "gnomad.v4.1.constraint_metrics.tsv")
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath"
CADD = "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38"

GROUPS = ["gwas_catalog", "genome", "1000g", "vep", "gtex", "eqtl_catalogue", "opentargets",
          "qtl_catalogue", "qtl_celltype", "molecular_qtl", "meqtl", "annotation", "cadd",
          "qtl_catalogue_susie", "adsp_eqtl", "1000g30x"]
OPTIONAL_GROUPS = {"qtl_catalogue_susie", "adsp_eqtl", "1000g30x"}   # only when named in --groups


# =============================================================================
# ARGUMENTS
# =============================================================================

parser = argparse.ArgumentParser(
    description="Plan shared-resource downloads and write a SLURM array script (one task per file)."
)
parser.add_argument("--groups", nargs="+", choices=GROUPS,
                    default=[g for g in GROUPS if g not in OPTIONAL_GROUPS],
                    help="Resource groups to include (default: all except qtl_catalogue_susie, adsp_eqtl, "
                         "1000g30x).")
parser.add_argument("--skip-groups", nargs="+", choices=GROUPS, default=[],
                    help="Resource groups to leave out.")
parser.add_argument("--template", default="HPC_CPU.sh",
                    help="SLURM template whose #SBATCH header is the base (default: HPC_CPU.sh).")
parser.add_argument("--output", default=DEFAULT_OUTPUT,
                    help=f"Generated SLURM script (default: {DEFAULT_OUTPUT}). A different name gets its own "
                         "task manifest, so e.g. a long meqtl-only plan can run next to the main one.")
parser.add_argument("--max-jobs", type=int, default=1000,
                    help="Most jobs queued at once: array tasks + 1 final check (default 1000). "
                         "With more files than tasks, each task downloads several files.")
parser.add_argument("--parallel", type=int, default=0,
                    help="Cap on simultaneously RUNNING tasks, array %%N (default 0 = no cap).")
parser.add_argument("--cpus", type=int, default=None, help="CPUs per download task.")
parser.add_argument("--mem", default=None, help="Memory per download task, e.g. 16G.")
parser.add_argument("--time", default=None, help="Walltime per download task, e.g. 24:00:00.")
parser.add_argument("--partition", default=None, help="SLURM partition.")
parser.add_argument("--qos", default=None, help="SLURM qos.")
parser.add_argument("--account", default=None, help="SLURM account, if your cluster requires one.")
parser.add_argument("--opentargets-release", default=OPEN_TARGETS_RELEASE,
                    help=f"Open Targets Platform release (default {OPEN_TARGETS_RELEASE}).")
parser.add_argument("--no-url-check", action="store_true",
                    help="Skip URL/size checks (faster; sizes are then not enforced).")
parser.add_argument("--force", action="store_true",
                    help="Re-plan every file, even valid ones.")
parser.add_argument("--finalize", action="store_true",
                    help="Verify the planned files, write the report and Open Targets markers; plan nothing.")
parser.add_argument("--status", action="store_true",
                    help="READ-ONLY: how many planned files are downloaded / partial / missing, per group.")
# eQTL Catalogue selections (same defaults as Step09)
parser.add_argument("--catalogue-study-regex", default="",
                    help="qtl_catalogue: study_label regex (default: Step09's GTEx-labelled + Sun_2018).")
parser.add_argument("--catalogue-tissue-regex", default="",
                    help="qtl_catalogue: only datasets whose tissue / sample group matches, e.g. 'brain'.")
parser.add_argument("--catalogue-all-studies", action="store_true",
                    help="qtl_catalogue: every catalogue study (very large).")
parser.add_argument("--catalogue-qtl-types", default="eQTL,sQTL,isoQTL,exonQTL,pQTL,apaQTL",
                    help="qtl_catalogue: QTL classes to include (default: all Step09 classes).")
parser.add_argument("--celltype-study-regex", default=r"^Jerber_2021$",
                    help="qtl_celltype: study_label regex (default: ^Jerber_2021$).")
parser.add_argument("--celltype-qtl-types", default="eQTL",
                    help="qtl_celltype: QTL classes to include (default eQTL).")
args = parser.parse_args()

SELECTED_GROUPS = [g for g in GROUPS if g in args.groups and g not in args.skip_groups]
OT_RELEASE = args.opentargets_release
OT_ROOT = Path("resources") / "opentargets" / OT_RELEASE
OT_EXPECTED = MANIFESTS / f"opentargets_{OT_RELEASE}_expected.tsv"   # expected file sizes

# Task manifest / plan / report belong to one generated script.
_SUFFIX = "" if args.output == DEFAULT_OUTPUT else "_" + re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(args.output).stem)
MANIFEST = MANIFESTS / f"download_resources{_SUFFIX}.tsv"
PLAN = MANIFESTS / f"download_resources{_SUFFIX}_plan.tsv"
REPORT = MANIFESTS / f"download_resources{_SUFFIX}_report.tsv"


# =============================================================================
# HELPERS
# =============================================================================

def banner(text: str) -> None:
    print()
    print("=" * 80)
    print(text)
    print("=" * 80)


def human(size) -> str:
    if not size:
        return "?"
    value = float(size)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024:
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} PB"


def http_open(url: str, method: str = "GET", headers: dict | None = None, timeout: int = 60):
    request = urllib.request.Request(
        url, method=method, headers={"User-Agent": "GWAS2m-Step00-Download/1.0", **(headers or {})}
    )
    return urllib.request.urlopen(request, timeout=timeout)


def https_url(url: str) -> str:
    """EBI FTP links -> HTTPS (same files; works through cluster firewalls)."""
    return re.sub(r"^ftp://ftp\.ebi\.ac\.uk/", "https://ftp.ebi.ac.uk/", url.strip())


def basename(url: str) -> str:
    return str(url).strip().split("?", 1)[0].rstrip("/").rsplit("/", 1)[-1]


def local_file_valid(path: Path, expected_size: int | None = None) -> bool:
    """Valid = non-empty, no .aria2 control file, exact size when known."""
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if Path(str(path) + ".aria2").exists():
        return False
    return expected_size is None or path.stat().st_size == expected_size


# =============================================================================
# FILE ITEMS  (one item = one file = one SLURM array task)
#   post  : - | fasta (decompress + faidx) | tar (extract into post_dir)
#           | tabix (index a bgzip TSV: chrom col 1, position col 2, '#' header lines)
#           | bed-tabix (plain BED-like TSV with one header line -> sorted <dest>.gz + .tbi;
#             the plain file is removed afterwards)
#           | lbf-tabix (eQTL Catalogue SuSiE lbf_variable file, ordered by molecular trait ->
#             position-sorted bgzip <name>.sorted.tsv.gz + .tbi, read by locus in Step09 / Step10;
#             the unsorted download is removed once the row counts match)
#   check : - | gzip | tar | parquet       (integrity test after download)
#   mutable: "latest"/"master" URLs whose size may change before download,
#            so their size is reported but not enforced.
# =============================================================================

def lbf_index_path(dest) -> Path:
    """<id>.lbf_variable.txt.gz -> <id>.lbf_variable.sorted.tsv.gz (same rule as Step09 / Step10)."""
    dest = Path(dest)
    name = dest.name
    stem = name[:-7] if name.endswith(".txt.gz") else name[:-3] if name.endswith(".gz") else name
    return dest.with_name(stem + ".sorted.tsv.gz")


def item(resource, group, url, dest, *, check="-", post="-", post_dir="-", post_marker="-",
         mutable=False, size=None):
    return dict(resource=resource, group=group, url=url, dest=str(dest), check=check, post=post,
                post_dir=str(post_dir), post_marker=post_marker, mutable=mutable, size=size)


def plan_preparation(g: Path | None = None) -> dict:
    """Which raw VCFs of a reference panel are still needed (raw files are only
    inputs to the Step04 PGEN panels under g/<ANC>/)."""
    g = g or ROOT / "resources" / "1000G"

    def panel_ok(prefix: Path) -> bool:
        return all(local_file_valid(Path(str(prefix) + ext)) for ext in (".pgen", ".pvar", ".psam"))

    population_missing = [(anc, c) for anc in LD_ANCESTRIES for c in range(1, 23)
                          if not panel_ok(g / anc / f"chr{c}_{anc}_GRCh38")]
    return {
        # Step04 builds every population panel from the raw VCF of its chromosome.
        "all_chromosomes": sorted({c for _, c in population_missing}),
        # The small sample panel: Step04 writes samples.txt from it on every plan.
        "samples": True,
        "panels_missing": len(population_missing),
    }


# =============================================================================
# 1000 GENOMES 30x LD REFERENCE PANEL (group 1000g30x; Step04 --reference-panel 1000G30x)
#
# NYGC 30x on GRCh38 (Byrska-Bishop et al. 2022, Cell): the FILTERED PHASED
# SNV/INDEL/SV panel per chromosome, not the raw "recalibrated_variants" genotype
# calls (unphased, unfiltered, > 140 GB for chr1 alone). Step04 keeps the 2504
# unrelated samples (those of the phase 3 sample panel), so related samples do
# not distort LD.
#   resources/LD_panels/1000G30x/raw/1kGP_high_coverage_Illumina.chr<N>.filtered.SNV_INDEL_SV_phased_panel.vcf.gz (+ .tbi)
#   resources/LD_panels/1000G30x/metadata/20130606_g1k_3202_samples_ped_population.txt
#   resources/1000G/metadata/integrated_call_samples_v3.20130502.ALL.panel (shared with 1000g)
# =============================================================================

def thousand_genomes_30x_items(notes: list[str]) -> list[dict]:
    items: list[dict] = []
    root = Path("resources/LD_panels/1000G30x")
    prep = plan_preparation(ROOT / root)

    for chrom in prep["all_chromosomes"]:
        name = f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        items.append(item("1000G30X_PHASED_VCF", "1000g30x", f"{ONE_K_30X_PHASED}/{name}", root / "raw" / name,
                          check="gzip"))
        items.append(item("1000G30X_PHASED_VCF", "1000g30x", f"{ONE_K_30X_PHASED}/{name}.tbi",
                          root / "raw" / f"{name}.tbi"))

    # Population of each of the 3202 samples.
    items.append(item("1000G30X_METADATA", "1000g30x",
                      f"{ONE_K_30X_BASE}/20130606_g1k_3202_samples_ped_population.txt",
                      root / "metadata" / "20130606_g1k_3202_samples_ped_population.txt"))
    # The 2504 unrelated samples (planned by the 1000g group when that is selected too).
    if "1000g" not in SELECTED_GROUPS:
        items.append(item("1000G30X_METADATA", "1000g30x", PHASE3_PANEL_URL,
                          Path("resources/1000G/metadata/integrated_call_samples_v3.20130502.ALL.panel")))

    if "1000g30x" in SELECTED_GROUPS and not prep["panels_missing"]:
        notes.append("1000g30x: all population PGEN panels are built, raw 30x files are not needed")

    return items


def static_items(notes: list[str]) -> list[dict]:
    items: list[dict] = []

    # --- GWAS Catalog metadata ---
    items += [
        item("GWAS_CATALOG_STUDIES", "gwas_catalog",
             f"{GWAS_CATALOG_BASE}/gwas-catalog-download-studies-v1.0.3.1.txt",
             "reference_metadata/gwas_catalog/gwas_catalog_studies.tsv", mutable=True),
        item("GWAS_CATALOG_ANCESTRY", "gwas_catalog",
             f"{GWAS_CATALOG_BASE}/gwas-catalog-download-ancestries-v1.0.3.1.txt",
             "reference_metadata/gwas_catalog/gwas_catalog_ancestry.tsv", mutable=True),
    ]

    # --- GENCODE genome + annotation ---
    gtf = f"gencode.v{GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"
    items += [
        item("GRCH38_FASTA", "genome", f"{GENCODE_BASE}/GRCh38.primary_assembly.genome.fa.gz",
             "resources/genome/GRCh38/GRCh38.primary_assembly.genome.fa.gz", check="gzip", post="fasta"),
        item("GENCODE_GTF", "genome", f"{GENCODE_BASE}/{gtf}",
             f"resources/gencode/release{GENCODE_RELEASE}/{gtf}", check="gzip"),
    ]

    # --- 1000 Genomes: raw files are only inputs to the PGEN panels ---
    prep = plan_preparation()
    raw, meta = Path("resources/1000G/raw"), Path("resources/1000G/metadata")
    for chrom in prep["all_chromosomes"]:
        name = f"ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
        items.append(item("1000G_RAW_VCF", "1000g", f"{ONE_K_BASE}/{name}", raw / name, check="gzip"))
        items.append(item("1000G_RAW_VCF", "1000g", f"{ONE_K_BASE}/{name}.tbi", raw / f"{name}.tbi"))
    if prep["samples"]:
        items.append(item("1000G_PANEL_METADATA", "1000g", PHASE3_PANEL_URL,
                          meta / "integrated_call_samples_v3.20130502.ALL.panel"))
    if prep["all_chromosomes"]:
        for name in ("20190312_biallelic_SNV_and_INDEL_README.txt", "20190312_biallelic_SNV_and_INDEL_MANIFEST.txt"):
            items.append(item("1000G_PANEL_METADATA", "1000g", f"{ONE_K_BASE}/{name}", meta / name))
    if "1000g" in SELECTED_GROUPS and not prep["panels_missing"]:
        notes.append("1000g: all population PGEN panels are built, raw 1000G files are not needed")

    # --- 1000 Genomes 30x (Step04 --reference-panel 1000G30x) ---
    items += thousand_genomes_30x_items(notes)

    # --- Ensembl VEP cache (extracted into resources/vep/cache) ---
    vep_url = (f"https://ftp.ensembl.org/pub/release-{VEP_RELEASE}/variation/indexed_vep_cache/"
               f"homo_sapiens_vep_{VEP_RELEASE}_GRCh38.tar.gz")
    items.append(item("VEP_CACHE", "vep", vep_url, f"resources/vep/{VEP_RELEASE}/{basename(vep_url)}",
                      check="tar", post="tar", post_dir="resources/vep/cache",
                      post_marker=f".vep{VEP_RELEASE}_GRCh38.complete"))

    # --- GTEx v11 (archives extracted next to themselves) ---
    for rid, path, sub, name in [
        ("GTEX_EQTL", "bulk-qtl/v11/single-tissue-cis-qtl/GTEx_Analysis_v11_eQTL.tar", "qtl", "eQTL"),
        ("GTEX_SQTL", "bulk-qtl/v11/single-tissue-cis-qtl/GTEx_Analysis_v11_sQTL.tar", "qtl", "sQTL"),
        ("GTEX_EQTL_SUSIE", "bulk-qtl/v11/susie-qtl/GTEx_Analysis_v11_eQTL_SuSiE.tar", "susie", "eQTL_SuSiE"),
        ("GTEX_SQTL_SUSIE", "bulk-qtl/v11/susie-qtl/GTEx_Analysis_v11_sQTL_SuSiE.tar", "susie", "sQTL_SuSiE"),
        # alternative polyadenylation QTL (read by Step09)
        ("GTEX_APAQTL", "bulk-qtl/v11/single-tissue-cis-qtl/GTEx_Analysis_v11_apaQTL.tar", "qtl", "apaQTL"),
        ("GTEX_APAQTL_SUSIE", "bulk-qtl/v11/susie-qtl/GTEx_Analysis_v11_apaQTL_SuSiE.tar", "susie", "apaQTL_SuSiE"),
    ]:
        directory = f"resources/gtex/{GTEX_RELEASE}/{sub}"
        items.append(item(rid, "gtex", f"{GTEX_BASE}/{path}", f"{directory}/{basename(path)}", check="tar",
                          post="tar", post_dir=f"{directory}/{name}", post_marker=".complete"))
    for rid, url in [
        ("GTEX_VARIANT_LOOKUP", f"{GTEX_BASE}/references/v11/reference-tables/"
                                "GTEx_Analysis_2021-02-11_v11_WholeGenomeSeq_953Indiv.lookup_table.txt.gz"),
        ("GTEX_GENCODE47", f"{GTEX_BASE}/references/v11/reference-tables/gencode.v47.genes.gtf"),
    ]:
        items.append(item(rid, "gtex", url, f"resources/gtex/{GTEX_RELEASE}/reference/{basename(url)}",
                          check="gzip" if url.endswith(".gz") else "-"))

    # --- eQTL Catalogue metadata tables ---
    for name, url in EQTL_METADATA.items():
        items.append(item("EQTL_CATALOGUE_METADATA", "eqtl_catalogue", url, EQTL_ROOT / name, mutable=True))

    return [i for i in items if i["group"] in SELECTED_GROUPS]


# -----------------------------------------------------------------------------
# eQTL Catalogue dense + provider SuSiE (repair of files already referenced locally)
# -----------------------------------------------------------------------------

_TABLE_CACHE: dict[str, str] = {}


def catalogue_table(name: str) -> str:
    """An eQTL Catalogue path table: the local copy if present, else fetched (not saved)."""
    if name not in _TABLE_CACHE:
        local = ROOT / EQTL_ROOT / name
        if local.is_file() and local.stat().st_size > 0:
            _TABLE_CACHE[name] = local.read_text(encoding="utf-8", errors="replace")
        else:
            with http_open(EQTL_METADATA[name], timeout=120) as response:
                _TABLE_CACHE[name] = response.read().decode("utf-8", errors="replace")
    return _TABLE_CACHE[name]


def read_catalogue() -> tuple[dict, dict]:
    aliases = {"dataset_id": ["dataset_id"], "ftp_path": ["ftp_path", "sumstats_path"],
               "ftp_cs_path": ["ftp_cs_path"], "ftp_lbf_path": ["ftp_lbf_path"]}
    by_basename, by_dataset = {}, {}
    for name in ("tabix_ftp_paths.tsv", "tabix_ftp_paths_imported.tsv"):  # later file wins
        for row in csv.DictReader(io.StringIO(catalogue_table(name)), delimiter="\t"):
            values = {k: next((row[a] for a in v if a in row and row[a]), "") for k, v in aliases.items()}
            dataset = values["dataset_id"].strip()
            if not dataset:
                continue
            by_dataset[dataset] = values
            if values["ftp_path"].strip():
                by_basename[basename(values["ftp_path"])] = (dataset, values["ftp_path"].strip())
    return by_basename, by_dataset


def usable_url(value: str) -> bool:
    return bool(value.strip()) and value.strip().lower() not in {"na", "none", "."}


def eqtl_catalogue_items(notes: list[str]) -> list[dict]:
    if "eqtl_catalogue" not in SELECTED_GROUPS:
        return []
    try:
        by_basename, by_dataset = read_catalogue()
    except (urllib.error.URLError, OSError) as exc:
        sys.exit(f"ERROR: cannot read the eQTL Catalogue path tables ({exc}).")
    extra = []
    if EXTRA_EQTL_DATASETS.is_file():
        extra = [x.split("#", 1)[0].strip() for x in EXTRA_EQTL_DATASETS.read_text(encoding="utf-8").splitlines()
                 if x.split("#", 1)[0].strip()]
        unknown_ids = [d for d in extra if d not in by_dataset]
        if unknown_ids:
            notes.append(f"eqtl_catalogue: {len(unknown_ids)} listed dataset ids are not in the catalogue tables "
                         f"(e.g. {', '.join(unknown_ids[:3])})")

    items = []
    base = ROOT / EQTL_ROOT
    dense = base / "dense"
    referenced: dict[str, Path] = {}
    if dense.is_dir():
        for path in dense.rglob("*"):
            if not path.is_file():
                continue
            name = path.name
            for suffix in PARTIAL_SUFFIXES:
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
                    break
            if name.endswith(".tbi"):
                name = name[:-4]
            if name.endswith(".tsv.gz"):
                referenced.setdefault(name, path.parent)
    for dataset in extra:
        url = by_dataset.get(dataset, {}).get("ftp_path", "")
        if usable_url(url):
            referenced.setdefault(basename(url), dense / dataset)

    unresolved = 0
    for name, directory in sorted(referenced.items()):
        if name not in by_basename:
            unresolved += 1
            continue
        url = https_url(by_basename[name][1])
        rel = directory.relative_to(ROOT)
        items.append(item("EQTL_CATALOGUE_DENSE", "eqtl_catalogue", url, rel / name, check="gzip"))
        items.append(item("EQTL_CATALOGUE_DENSE", "eqtl_catalogue", url + ".tbi", rel / f"{name}.tbi"))

    provider = base / "susie_provider"
    datasets = sorted({p.name for p in provider.iterdir() if p.is_dir()} if provider.is_dir() else set())
    datasets += [d for d in extra if not (provider / d).is_dir()]
    unknown = 0
    for dataset in datasets:
        record = by_dataset.get(dataset)
        if record is None:
            unknown += 1
            continue
        for column in ("ftp_lbf_path", "ftp_cs_path"):
            url = record.get(column, "")
            if usable_url(url):
                items.append(item("EQTL_CATALOGUE_SUSIE_PROVIDER", "eqtl_catalogue", https_url(url),
                                  (provider / dataset / basename(url)).relative_to(ROOT), check="gzip",
                                  post="lbf-tabix" if "lbf_variable" in basename(url) else "-"))
    if unresolved or unknown:
        notes.append(f"eqtl_catalogue: {unresolved + unknown} local datasets are not in the catalogue tables")
    return items


# -----------------------------------------------------------------------------
# Open Targets (HTTPS listing; sizes from the URL check)
# -----------------------------------------------------------------------------

def list_http_directory(url: str, depth: int = 0) -> list[str]:
    with http_open(url, timeout=120) as response:
        html = response.read().decode("utf-8", errors="replace")
    files = []
    for href in re.findall(r'href="([^"]+)"', html):
        if href.startswith(("?", "/", "#", "..", "http:", "https:", "mailto:")):
            continue
        href = urllib.parse.unquote(href)
        if href.endswith("/"):
            if depth < 3:
                files += [href + sub for sub in list_http_directory(url + href, depth + 1)]
        else:
            files.append(href)
    return files


def opentargets_items() -> list[dict]:
    if "opentargets" not in SELECTED_GROUPS:
        return []
    items = []
    for dataset in OT_DATASETS:
        base = f"{OT_HTTPS}/{OT_RELEASE}/output/{dataset}/"
        try:
            files = sorted(f for f in list_http_directory(base) if f.endswith(".parquet"))
        except (urllib.error.URLError, OSError) as exc:
            sys.exit(f"ERROR: cannot list Open Targets {OT_RELEASE}/{dataset} ({exc}).\n"
                     "       Check --opentargets-release, or use --skip-groups opentargets.")
        if not files:
            sys.exit(f"ERROR: no parquet files found under {base}")
        print(f"  Open Targets {OT_RELEASE}/{dataset:<16} {len(files):>4} parquet files")
        for rel in files:
            items.append(item("OPEN_TARGETS", f"opentargets/{dataset}", base + rel, OT_ROOT / dataset / rel,
                              check="parquet"))
    return items


# -----------------------------------------------------------------------------
# eQTL Catalogue dense summary statistics for colocalisation
#   Layout read by Step09:
#     resources/coloc/eqtl_catalogue/dense/<dataset_id>/<file>.tsv.gz (+ .tbi)
#     resources/coloc/eqtl_catalogue/susie_provider/<dataset_id>/<lbf|cs file>
# -----------------------------------------------------------------------------

def quant_method_to_qtl_type(value: str) -> str | None:
    """Same mapping as Step09."""
    q = re.sub(r"[^a-z0-9]+", " ", str(value).lower()).strip()
    if q in {"ge", "gene expression", "gene", "microarray"}:
        return "eQTL"
    if "leafcutter" in q or "majiq" in q or q in {"sqtl", "splicing"}:
        return "sQTL"
    if q in {"tx", "transcript", "transcript usage", "txrev", "txrevise"}:
        return "isoQTL"
    if q in {"exon", "exon expression"}:
        return "exonQTL"
    if "protein" in q or "aptamer" in q or q in {"prot", "pqtl"}:
        return "pQTL"
    if "apa" in q or "polyadenyl" in q:
        return "apaQTL"
    return None


def catalogue_rows() -> list[dict]:
    """Both eQTL Catalogue path tables (local copy if present, else fetched); later table wins."""
    rows: dict = {}
    for name in ("tabix_ftp_paths.tsv", "tabix_ftp_paths_imported.tsv"):
        for raw in csv.DictReader(io.StringIO(catalogue_table(name)), delimiter="\t"):
            r = {k: (v or "").strip() for k, v in raw.items() if k}
            row = {
                "dataset_id": r.get("dataset_id", ""),
                "study_label": r.get("study_label") or r.get("study", ""),
                "quant_method": r.get("quant_method") or r.get("quantification_method", ""),
                "ftp_path": r.get("ftp_path") or r.get("sumstats_path", ""),
                "ftp_lbf_path": r.get("ftp_lbf_path", ""),
                "ftp_cs_path": r.get("ftp_cs_path", ""),
                "tissue": " ".join(x for x in (r.get("tissue_label", ""), r.get("sample_group", ""),
                                               r.get("qtl_group", "")) if x),
            }
            row["qtl_type"] = quant_method_to_qtl_type(row["quant_method"])
            if row["ftp_path"] and row["qtl_type"]:
                rows[(row["dataset_id"], row["ftp_path"])] = row
    return list(rows.values())


def dense_items(row: dict, group: str) -> list[dict]:
    name = basename(row["ftp_path"])
    directory = EQTL_ROOT / "dense" / row["dataset_id"] if row["dataset_id"] else EQTL_ROOT / "dense"
    url = https_url(row["ftp_path"])
    return [item("EQTL_CATALOGUE_DENSE", group, url, directory / name, check="gzip"),
            item("EQTL_CATALOGUE_DENSE", group, url + ".tbi", directory / f"{name}.tbi")]


def provider_items(row: dict, group: str) -> list[dict]:
    out = []
    for column in ("ftp_lbf_path", "ftp_cs_path"):
        url = row.get(column, "")
        if row["dataset_id"] and usable_url(url):
            out.append(item("EQTL_CATALOGUE_SUSIE_PROVIDER", group, https_url(url),
                            EQTL_ROOT / "susie_provider" / row["dataset_id"] / basename(url), check="gzip",
                            post="lbf-tabix" if "lbf_variable" in basename(url) else "-"))
    return out


def summarise(label: str, rows: list[dict]) -> None:
    counts: dict = {}
    for r in rows:
        counts[(r["study_label"], r["qtl_type"])] = counts.get((r["study_label"], r["qtl_type"]), 0) + 1
    print(f"  {label}: {len(rows)} datasets")
    for (study, qtl), n in sorted(counts.items()):
        print(f"      {study:<24}{qtl:<10}{n:>4}")


def catalogue_qtl_items() -> list[dict]:
    wanted = {"qtl_catalogue", "qtl_catalogue_susie", "qtl_celltype"} & set(SELECTED_GROUPS)
    if not wanted:
        return []
    try:
        rows = catalogue_rows()
    except (urllib.error.URLError, OSError) as exc:
        sys.exit(f"ERROR: cannot read the eQTL Catalogue path tables ({exc}).")
    items: list[dict] = []

    if {"qtl_catalogue", "qtl_catalogue_susie"} & wanted:
        types = {t.strip().lower() for t in args.catalogue_qtl_types.split(",") if t.strip()}
        selected = [r for r in rows if r["qtl_type"].lower() in types]
        if args.catalogue_study_regex:
            pattern = re.compile(args.catalogue_study_regex, flags=re.I)
            selected = [r for r in selected if pattern.search(r["study_label"])]
        elif not args.catalogue_all_studies:      # Step09 default contexts
            selected = [r for r in selected
                        if "gtex" in r["study_label"].lower() or r["study_label"].lower() == "sun_2018"]
        if args.catalogue_tissue_regex:
            tissue = re.compile(args.catalogue_tissue_regex, flags=re.I)
            selected = [r for r in selected if tissue.search(r.get("tissue", ""))]
        summarise("eQTL Catalogue (Step09 contexts)", selected)
        for r in selected:
            if "qtl_catalogue" in wanted:
                items += dense_items(r, "qtl_catalogue")
            if "qtl_catalogue_susie" in wanted:
                items += provider_items(r, "qtl_catalogue_susie")

    if "qtl_celltype" in wanted:
        pattern = re.compile(args.celltype_study_regex, flags=re.I)
        types = {t.strip().lower() for t in args.celltype_qtl_types.split(",") if t.strip()}
        selected = [r for r in rows if pattern.search(r["study_label"]) and r["qtl_type"].lower() in types]
        summarise("eQTL Catalogue (cell-type contexts)", selected)
        for r in selected:
            items += dense_items(r, "qtl_celltype") + provider_items(r, "qtl_celltype")
    return items


# -----------------------------------------------------------------------------
# Other molecular QTL (pQTL / sQTL / haQTL / meQTL / caQTL) and annotation resources
#   Archives are MD5-verified, extracted, then deleted ("tar-rm") to save space.
# -----------------------------------------------------------------------------

def molecular_qtl_items() -> list[dict]:
    items = []
    for qtl, group, (manifest_id, manifest_md5), (archive_id, archive_md5) in ADSP_ARCHIVES:
        if group not in SELECTED_GROUPS:
            continue
        directory = Path("resources/molecular_qtl") / qtl.lower() / "adsp_fungen_brain"
        stem = f"ADSP_FunGen_xQTL.v1.{qtl}.all"
        items.append(item(f"ADSP_{qtl.upper()}_MANIFEST", group, f"{NIAGADS}/{manifest_id}",
                          directory / f"{stem}.manifest.txt", check=f"md5:{manifest_md5}"))
        items.append(item(f"ADSP_{qtl.upper()}", group, f"{NIAGADS}/{archive_id}", directory / f"{stem}.tar",
                          check=f"md5:{archive_md5}", post="tar-rm", post_dir=directory / "all",
                          post_marker=".complete"))
    if "molecular_qtl" in SELECTED_GROUPS:
        directory = Path("resources/molecular_qtl/caqtl/kumasaka_grch38")
        items += [
            item("KUMASAKA_CAQTL", "molecular_qtl", f"{KUMASAKA}/QTD100018.all.tsv.gz?download=1",
                 directory / "QTD100018.all.tsv.gz", check="md5:5f92d59ebcff47d151b10550cbfc67bd"),
            item("KUMASAKA_CAQTL_PEAKS", "molecular_qtl", f"{KUMASAKA}/QTD100018_peak_metadata.tsv.gz?download=1",
                 directory / "QTD100018_peak_metadata.tsv.gz", check="md5:4a77a273b0728a7b3ca35f0d914955b7"),
        ]
    return items


def annotation_items() -> list[dict]:
    items = []
    base = Path("resources/annotation")
    if "annotation" in SELECTED_GROUPS:
        items += [
            item("ALPHAMISSENSE", "annotation", f"{ALPHAMISSENSE}/AlphaMissense_hg38.tsv.gz",
                 base / "alphamissense/AlphaMissense_hg38.tsv.gz", check="gzip", post="tabix"),
            item("ALPHAMISSENSE_GENE", "annotation", f"{ALPHAMISSENSE}/AlphaMissense_gene_hg38.tsv.gz",
                 base / "alphamissense/AlphaMissense_gene_hg38.tsv.gz", check="gzip"),
            item("CLINVAR", "annotation", f"{CLINVAR}/clinvar.vcf.gz", base / "clinvar/clinvar.vcf.gz",
                 check="gzip", mutable=True),
            item("CLINVAR", "annotation", f"{CLINVAR}/clinvar.vcf.gz.tbi", base / "clinvar/clinvar.vcf.gz.tbi",
                 mutable=True),
            item("ENCODE_SCREEN_CCRES_V4", "annotation", SCREEN_CCRES, base / "encode_screen/GRCh38-cCREs.V4.bed"),
            item("GNOMAD_V4_1_CONSTRAINT", "annotation", GNOMAD_CONSTRAINT,
                 base / "gnomad/gnomad.v4.1.constraint_metrics.tsv"),
            item("LIFTOVER_HG19_TO_HG38", "annotation", f"{UCSC}/hg19/liftOver/hg19ToHg38.over.chain.gz",
                 base / "liftover/hg19ToHg38.over.chain.gz", check="gzip"),
            item("LIFTOVER_HG38_TO_HG19", "annotation", f"{UCSC}/hg38/liftOver/hg38ToHg19.over.chain.gz",
                 base / "liftover/hg38ToHg19.over.chain.gz", check="gzip"),
            item("PHYLOP_100WAY", "annotation", f"{UCSC}/hg38/phyloP100way/hg38.phyloP100way.bw",
                 base / "conservation/hg38.phyloP100way.bw"),
            item("REGULOMEDB_V2_2", "annotation", REGULOMEDB,
                 base / "regulomedb/RegulomeDB_v2.2_dbSNP153_common_GRCh38.tsv",
                 check=f"md5:{REGULOMEDB_MD5}", post="bed-tabix"),
        ]
    if "cadd" in SELECTED_GROUPS:
        for name in ("whole_genome_SNVs.tsv.gz", "gnomad.genomes.r4.0.indel.tsv.gz"):
            items += [
                item("CADD_V1_7", "cadd", f"{CADD}/{name}", base / "cadd/v1.7" / name, check="gzip"),
                item("CADD_V1_7", "cadd", f"{CADD}/{name}.tbi", base / "cadd/v1.7" / f"{name}.tbi"),
            ]
    return items


# =============================================================================
# STATUS
# =============================================================================

def post_done(entry: dict) -> bool:
    if entry["post"] == "fasta":
        fasta = ROOT / entry["dest"][: -len(".gz")]
        return local_file_valid(fasta) and local_file_valid(Path(str(fasta) + ".fai"))
    if entry["post"] in ("tar", "tar-rm"):   # marker written after a successful extraction
        return (ROOT / entry["post_dir"] / entry["post_marker"]).is_file()
    if entry["post"] == "tabix":
        return local_file_valid(ROOT / entry["dest"]) and local_file_valid(ROOT / (entry["dest"] + ".tbi"))
    if entry["post"] == "bed-tabix":
        return local_file_valid(ROOT / (entry["dest"] + ".gz")) and local_file_valid(ROOT / (entry["dest"] + ".gz.tbi"))
    if entry["post"] == "lbf-tabix":
        index = lbf_index_path(ROOT / entry["dest"])
        return local_file_valid(index) and local_file_valid(Path(str(index) + ".tbi"))
    return True


def is_complete(entry: dict) -> bool:
    if entry["post"] != "-":
        return post_done(entry)        # extracted / indexed: the archive itself may be gone
    size = None if entry["mutable"] else entry["size"]
    return local_file_valid(ROOT / entry["dest"], size)


def probe(entry: dict) -> tuple[dict, str | None]:
    """Reachability + remote size (HEAD, falling back to a 1-byte ranged GET)."""
    last = None
    for method, headers in (("HEAD", None), ("GET", {"Range": "bytes=0-0"})):
        for attempt in range(3):
            try:
                with http_open(entry["url"], method=method, headers=headers, timeout=60) as response:
                    if method == "HEAD":
                        size = int(response.headers.get("Content-Length") or 0)
                    else:
                        match = re.search(r"/(\d+)$", response.headers.get("Content-Range") or "")
                        size = int(match.group(1)) if match else 0
                    entry["size"] = size or None
                    return entry, None
            except urllib.error.HTTPError as exc:
                last = f"HTTP {exc.code}"
                if exc.code in (404, 410):
                    return entry, last
                if exc.code in (403, 405, 501):
                    break
            except (urllib.error.URLError, OSError) as exc:
                last = f"{type(exc).__name__}: {exc}"
            time.sleep(2 * (attempt + 1))
    return entry, last


# =============================================================================
# SLURM SETTINGS: CLI > HPC_CPU.sh > built-in
# =============================================================================

BUILTIN = {"cpus-per-task": "2", "mem": "16G", "time": "24:00:00", "nodes": "1", "ntasks": "1"}

# Set per generated job; everything else in the template applies to all jobs.
PER_JOB_KEYS = {"job-name", "output", "error", "array"}


def read_template(path: Path) -> tuple[str, dict]:
    shebang, values = "#!/bin/bash --login", {}
    if not path.is_file():
        return shebang, values
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("#!"):
            shebang = line
            continue
        match = re.match(r"^#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
        if match and (match.group(2) or "").strip():
            values[match.group(1)] = match.group(2).strip()
    return shebang, values


def slurm_settings() -> dict:
    shebang, template = read_template(ROOT / args.template)
    if not template:
        print(f"[WARN] No #SBATCH settings found in {args.template}; using built-in defaults.")
    values, source = dict(BUILTIN), {k: "built-in" for k in BUILTIN}

    # HPC_CPU.sh is the project's SLURM configuration.
    for key, value in template.items():
        if key not in PER_JOB_KEYS:
            values[key], source[key] = value, args.template

    for key, value in (("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time),
                       ("partition", args.partition), ("qos", args.qos), ("account", args.account)):
        if value not in (None, ""):
            values[key], source[key] = str(value), "command line"

    if args.max_jobs < 2:
        sys.exit("ERROR: --max-jobs must be at least 2 (one download task + the final check)")
    return dict(shebang=shebang, values=values, source=source,
                chunk=args.max_jobs - 1, max_jobs=args.max_jobs, parallel=args.parallel)


# =============================================================================
# GENERATED SLURM SCRIPT
# =============================================================================

SCRIPT_TEMPLATE = r'''@@SHEBANG@@
@@DIRECTIVES@@

# =============================================================================
# GWAS2m STEP 00 - SHARED RESOURCE DOWNLOAD (generated @@DATE@@)
# Generated by Step00_Download_Resources.py - re-run it instead of editing.
#
# ONE ARRAY TASK = ONE FILE (row of setup_manifests/download_resources.tsv)
# Files to download : @@TOTAL@@
# Array tasks       : @@NTASKS@@ (+ 1 final check = @@QUEUED@@ jobs; limit @@MAX_JOBS@@)
# Files per task    : up to @@PER_TASK@@ (task k downloads rows k, k+@@NTASKS@@, k+2x@@NTASKS@@, ...)
#
#   sbatch @@SCRIPT_NAME@@             # SLURM cluster
#   bash   @@SCRIPT_NAME@@ local 4     # no SLURM: 4 files at a time
#   bash   @@SCRIPT_NAME@@ task 12     # just row 12 (debugging)
#   bash   @@SCRIPT_NAME@@ finalize    # verify everything + report
#
# Check progress / readability (project interpreter):
#   envs/pipeline/bin/python Step00_Download_Resources.py --finalize@@OUTPUT_ARG@@
#   envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py
# =============================================================================

set -uo pipefail

ROOT="@@ROOT@@"
SCRIPT="$ROOT/@@SCRIPT_NAME@@"
PLANNER="@@PLANNER@@"
MANIFEST="@@MANIFEST@@"
LOGS="@@LOGS@@"
TOTAL=@@TOTAL@@
NTASKS=@@NTASKS@@
THROTTLE="@@THROTTLE@@"
SBATCH_COMMON=(@@SBATCH_COMMON@@)

cd "$ROOT" || exit 1
mkdir -p "$LOGS"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

# Prefer the project's own tools (envs/pipeline), else whatever is on PATH.
tool() {
    if [[ -x "$ROOT/envs/pipeline/bin/$1" ]]; then echo "$ROOT/envs/pipeline/bin/$1"; else command -v "$1" || true; fi
}
ARIA2="$(tool aria2c)"
CURL="$(tool curl)"
WGET="$(tool wget)"
PIGZ="$(tool pigz)"
SAMTOOLS="$(tool samtools)"
TABIX="$(tool tabix)"
BGZIP="$(tool bgzip)"
# The final check runs with the project's interpreter, envs/pipeline/bin/python.
PYTHON="$ROOT/envs/pipeline/bin/python"
if [[ ! -x "$PYTHON" ]]; then
    PYTHON="$(command -v python3 || command -v python || true)"
    echo "WARNING: $ROOT/envs/pipeline/bin/python not found (run: python Step00_Installation.py); using ${PYTHON:-none}" >&2
fi

log() { echo "[$(date '+%F %T')] $*"; }
file_size() { stat -c %s "$1" 2>/dev/null || stat -f %z "$1" 2>/dev/null || echo 0; }

# -----------------------------------------------------------------------------
# validity (gwas2m rule): non-empty, no aria2 control file, exact size if known
# -----------------------------------------------------------------------------
valid() {
    [[ -s "$OUT" && ! -e "$OUT.aria2" ]] || return 1
    [[ "$SIZE" == "-" ]] && return 0
    [[ "$(file_size "$OUT")" == "$SIZE" ]]
}

post_done() {
    case "$POST" in
        fasta) [[ -s "${OUT%.gz}" && -s "${OUT%.gz}.fai" ]] ;;
        tar|tar-rm) [[ -e "$ROOT/$PDIR/$PMARKER" ]] ;;
        tabix) [[ -s "$OUT" && -s "$OUT.tbi" ]] ;;
        bed-tabix) [[ -s "$OUT.gz" && -s "$OUT.gz.tbi" ]] ;;
        lbf-tabix) [[ -s "$(lbf_index "$OUT")" && -s "$(lbf_index "$OUT").tbi" ]] ;;
        *)     return 0 ;;
    esac
}

# <id>.lbf_variable.txt.gz -> <id>.lbf_variable.sorted.tsv.gz (same rule as Step09 / Step10)
lbf_index() {
    local stem="$1"
    case "$stem" in *.txt.gz) stem="${stem%.txt.gz}" ;; *.gz) stem="${stem%.gz}" ;; esac
    echo "$stem.sorted.tsv.gz"
}

md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi
}

# -----------------------------------------------------------------------------
# download: aria2c (resumable via .aria2), else curl/wget into <file>.part
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# walltime: seconds left in this SLURM job (huge outside SLURM). Every download
# attempt is stopped RESERVE_S before the end, so it can resume in a new job.
# -----------------------------------------------------------------------------
RESERVE_S=900
MAX_CONTINUATIONS=10
MAX_LOCAL_ROUNDS=50

time_left() {
    local left rest d=0 h=0 m=0 s=0 parts
    if [[ -z "${SLURM_JOB_ID:-}" ]] || ! command -v squeue >/dev/null 2>&1; then echo 999999999; return; fi
    left="$(squeue -h -j "$SLURM_JOB_ID" -o %L 2>/dev/null | head -1)"
    if [[ -z "$left" || "$left" == UNLIMITED || "$left" == INVALID || "$left" == NOT_SET ]]; then
        echo 999999999; return
    fi
    rest="$left"
    if [[ "$rest" == *-* ]]; then d="${rest%%-*}"; rest="${rest#*-}"; fi
    IFS=: read -r -a parts <<< "$rest"
    case ${#parts[@]} in
        3) h="${parts[0]}"; m="${parts[1]}"; s="${parts[2]}" ;;
        2) m="${parts[0]}"; s="${parts[1]}" ;;
        1) s="${parts[0]}" ;;
    esac
    echo $(( 10#${d:-0} * 86400 + 10#${h:-0} * 3600 + 10#${m:-0} * 60 + 10#${s:-0} ))
}

# Run one download attempt, stopped (TERM: aria2c / curl keep what they have)
# RESERVE_S before the walltime ends.
run_limited() {
    local budget=$(( $(time_left) - RESERVE_S ))
    (( budget < 60 )) && return 124
    if command -v timeout >/dev/null 2>&1 && (( budget < 999999999 - RESERVE_S )); then
        timeout --signal=TERM "$budget" "$@"
    else
        "$@"
    fi
}

http_status() {
    [[ -n "$CURL" ]] || { echo 0; return; }
    "$CURL" -sIL -o /dev/null -w '%{http_code}' --connect-timeout 60 --max-time 120 "$URL" 2>/dev/null || echo 0
}

remote_size() {
    [[ -n "$CURL" ]] || return 0
    "$CURL" -sIL --connect-timeout 60 --max-time 120 "$URL" 2>/dev/null | tr -d '\r' \
        | awk 'tolower($1) == "content-length:" { n = $2 } END { if (n != "") print n }'
}

aria2() {
    run_limited "$ARIA2" --continue=true --max-connection-per-server=4 --split=4 --min-split-size=16M \
        --file-allocation=none --auto-file-renaming=false --allow-overwrite=true \
        --max-tries=10 --retry-wait=10 --connect-timeout=60 --timeout=120 \
        --console-log-level=warn --summary-interval=60 \
        --dir "$(dirname "$OUT")" --out "$(basename "$OUT")" "$URL"
}

to_part() {
    if [[ -n "$CURL" ]]; then
        # --speed-limit/--speed-time: a stalled connection (< 1 KB/s for 10 min) is
        # dropped and resumed instead of hanging until the walltime.
        run_limited "$CURL" --fail --location --retry 10 --retry-delay 10 \
            --connect-timeout 60 --speed-limit 1024 --speed-time 600 \
            --continue-at - --output "$OUT.part" "$URL" || return 1
    elif [[ -n "$WGET" ]]; then
        run_limited "$WGET" --continue --tries=10 --waitretry=10 --timeout=120 -O "$OUT.part" "$URL" || return 1
    else
        log "ERROR: none of aria2c / curl / wget is available"; return 1
    fi
    mv "$OUT.part" "$OUT"
}

# One download attempt (resumes whatever is there).
fetch_once() {
    if [[ -e "$OUT.aria2" && -n "$ARIA2" ]]; then
        log "resuming aria2 download ($(file_size "$OUT") bytes so far)"; aria2
    elif [[ -e "$OUT.part" || -z "$ARIA2" ]]; then
        [[ -e "$OUT.part" ]] && log "resuming $OUT.part ($(file_size "$OUT.part") bytes so far)"
        to_part
    else
        rm -f "$OUT"            # only reached when OUT is empty / wrong size
        aria2 && return 0
        log "aria2c failed; trying curl/wget"
        to_part
    fi
}

# Keep resuming until the file is complete. Returns 0 complete, 2 the file does
# not exist on the server, 3 walltime nearly over (continue in a new job).
fetch() {
    local round=0 wait code remote left
    while true; do
        round=$(( round + 1 ))
        fetch_once
        valid && return 0
        # Complete download, but not the planned size: the server's file changed?
        if [[ -s "$OUT" && ! -e "$OUT.aria2" && "$SIZE" != "-" ]]; then
            remote="$(remote_size)"
            if [[ -n "$remote" && "$remote" == "$(file_size "$OUT")" ]]; then
                log "WARNING: the server's file is $remote bytes (planned: $SIZE); using the server's size"
                SIZE="$remote"; return 0
            fi
        fi
        code="$(http_status)"
        if [[ "$code" == 404 || "$code" == 410 ]]; then
            log "FAILED: the server returns HTTP $code for $URL (file not available)"; return 2
        fi
        left="$(time_left)"
        if (( left < RESERVE_S + 300 )); then
            log "walltime nearly over after $round attempt(s); $(file_size "$OUT") / $(file_size "$OUT.part") bytes kept"
            return 3
        fi
        if [[ -z "${SLURM_JOB_ID:-}" ]] && (( round >= MAX_LOCAL_ROUNDS )); then
            log "FAILED: $round attempts without SLURM"; return 1
        fi
        wait=$(( 60 * round < 600 ? 60 * round : 600 ))
        log "attempt $round incomplete (HTTP $code); resuming in ${wait}s"
        sleep "$wait"
    done
}

# Submit a job that continues this array task (same rows) after the walltime.
continue_in_new_job() {
    local n=$(( ${G2M_CONTINUATION:-0} + 1 )) job
    if (( n > MAX_CONTINUATIONS )); then
        log "FAILED: still incomplete after $MAX_CONTINUATIONS continuation jobs"; return 1
    fi
    job="$(sbatch --parsable --array="$SLURM_ARRAY_TASK_ID" --export=ALL,G2M_CONTINUATION="$n" "$SCRIPT")" || {
        log "FAILED: could not submit a continuation job"; return 1; }
    log "continuing the unfinished download(s) of task $SLURM_ARRAY_TASK_ID in job $job ($n/$MAX_CONTINUATIONS)"
}

integrity() {
    case "$CHECK" in
        gzip)    if [[ -n "$PIGZ" ]]; then "$PIGZ" -t "$OUT"; else gzip -t "$OUT"; fi ;;
        tar)     tar -tf "$OUT" >/dev/null ;;
        parquet) [[ "$(tail -c 4 "$OUT")" == "PAR1" ]] ;;
        md5:*)   log "computing MD5 ..."; [[ "$(md5_of "$OUT")" == "${CHECK#md5:}" ]] ;;
        *)       return 0 ;;
    esac
}

post_process() {
    local fasta
    case "$POST" in
        fasta)
            fasta="${OUT%.gz}"
            if [[ ! -s "$fasta" ]]; then
                log "decompressing -> ${fasta#$ROOT/}"
                if [[ -n "$PIGZ" ]]; then "$PIGZ" -dc "$OUT" > "$fasta.part"; else gzip -dc "$OUT" > "$fasta.part"; fi \
                    && mv "$fasta.part" "$fasta" || return 1
            fi
            [[ -n "$SAMTOOLS" ]] || { log "ERROR: samtools not found (run Step00_Installation.py first)"; return 1; }
            [[ -s "$fasta.fai" ]] || "$SAMTOOLS" faidx "$fasta" || return 1
            log "indexed: ${fasta#$ROOT/}.fai" ;;
        tar|tar-rm)
            mkdir -p "$ROOT/$PDIR"
            log "extracting -> $PDIR"
            if [[ "$OUT" == *.gz ]]; then tar -xzf "$OUT" -C "$ROOT/$PDIR"; else tar -xf "$OUT" -C "$ROOT/$PDIR"; fi || return 1
            touch "$ROOT/$PDIR/$PMARKER"
            log "extracted: $PDIR/$PMARKER"
            if [[ "$POST" == "tar-rm" ]]; then
                rm -f "$OUT"; log "removed the verified archive to save space: ${OUT#$ROOT/}"
            fi ;;
        tabix)
            [[ -n "$TABIX" ]] || { log "ERROR: tabix not found (run Step00_Installation.py first)"; return 1; }
            "$TABIX" -f -s 1 -b 2 -e 2 "$OUT" || return 1
            log "indexed: ${OUT#$ROOT/}.tbi" ;;
        bed-tabix)
            [[ -n "$TABIX" && -n "$BGZIP" ]] || { log "ERROR: tabix / bgzip not found (run Step00_Installation.py first)"; return 1; }
            log "sorting + bgzip -> ${OUT#$ROOT/}.gz"
            { head -n 1 "$OUT"
              tail -n +2 "$OUT" | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -S 1G -T "$(dirname "$OUT")"
            } | "$BGZIP" -c > "$OUT.gz.part" && mv "$OUT.gz.part" "$OUT.gz" || return 1
            "$TABIX" -f -0 -s 1 -b 2 -e 3 -S 1 "$OUT.gz" || return 1
            rm -f "$OUT"
            log "indexed: ${OUT#$ROOT/}.gz.tbi (removed the plain TSV)" ;;
        lbf-tabix)
            # Provider SuSiE LBFs are ordered by molecular trait; Step09 / Step10 read them by locus, so
            # sort by chromosome + position once, bgzip and tabix-index (header line kept as '#...').
            [[ -n "$TABIX" && -n "$BGZIP" ]] || { log "ERROR: tabix / bgzip not found (run Step00_Installation.py first)"; return 1; }
            local idx header n_in n_out threads="${SLURM_CPUS_PER_TASK:-2}"
            idx="$(lbf_index "$OUT")"
            header="$(gzip -dc "$OUT" 2>/dev/null | head -n 1 || true)"
            [[ "$header" == *molecular_trait_id*position* ]] || { log "ERROR: unexpected LBF header: ${header:0:120}"; return 1; }
            log "sorting by chromosome + position + bgzip -> ${idx#$ROOT/}"
            { printf '#%s\n' "$header"
              gzip -dc "$OUT" | tail -n +2 | LC_ALL=C sort -t $'\t' -k4,4 -k5,5n -S 2G --parallel="$threads" -T "$(dirname "$OUT")"
            } | "$BGZIP" -@ "$threads" -c > "$idx.part" || { rm -f "$idx.part"; return 1; }
            mv "$idx.part" "$idx"
            "$TABIX" -f -s 4 -b 5 -e 5 "$idx" || return 1
            n_in="$(gzip -dc "$OUT" | tail -n +2 | wc -l)"
            n_out="$(gzip -dc "$idx" | grep -vc '^#')"
            if [[ "$n_in" != "$n_out" ]]; then
                log "ERROR: sorted LBF has $n_out rows, download has $n_in - keeping the download"
                rm -f "$idx" "$idx.tbi"; return 1
            fi
            rm -f "$OUT"
            log "indexed: ${idx#$ROOT/}.tbi ($n_out rows; removed the unsorted download)" ;;
    esac
}

# -----------------------------------------------------------------------------
# one task = one manifest row
# -----------------------------------------------------------------------------
run_task() {
    # Always release the per-file lock, whether the task succeeded or failed.
    OUT=""
    _run_task "$@"
    local rc=$?
    [[ -n "$OUT" ]] && rm -f "$OUT.lock"
    return $rc
}

_run_task() {
    local idx="$1" line I RESOURCE GROUP URL DEST SIZE CHECK POST PDIR PMARKER
    line="$(awk -F'\t' -v i="$idx" 'NR > 1 && $1 == i { print; exit }' "$MANIFEST")"
    [[ -n "$line" ]] || { log "ERROR: no manifest row $idx"; return 2; }
    IFS=$'\t' read -r I RESOURCE GROUP URL DEST SIZE CHECK POST PDIR PMARKER <<< "$line"
    OUT="$ROOT/$DEST"

    echo "=============================================================================="
    echo "Task     : $idx / $TOTAL   (SLURM ${SLURM_ARRAY_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-$idx})"
    echo "Resource : $RESOURCE   [$GROUP]"
    echo "Source   : $URL"
    echo "Target   : $DEST   (expected size: $SIZE)"
    echo "Node     : $(hostname)    Start: $(date)"
    echo "=============================================================================="

    if [[ "$POST" != "-" ]] && post_done; then log "[SKIP] already extracted / indexed"; return 0; fi
    mkdir -p "$(dirname "$OUT")"
    if command -v flock >/dev/null 2>&1; then exec 9>"$OUT.lock"; flock 9; fi

    local attempt rc
    for attempt in 1 2 3; do
        if valid; then
            (( attempt == 1 )) && log "[SKIP] target already valid"
        else
            fetch; rc=$?
            if (( rc == 3 )); then CONTINUE_LATER=1; return 3; fi
            (( rc == 0 )) || { log "FAILED: download"; return 1; }
        fi
        integrity && break
        log "integrity check ($CHECK) failed (attempt $attempt/3): deleting the corrupt file and downloading again"
        rm -f "$OUT" "$OUT.part" "$OUT.aria2"
        (( attempt == 3 )) && { log "FAILED: integrity check ($CHECK) after 3 downloads"; return 1; }
    done
    log "verified: $(file_size "$OUT") bytes ($CHECK)"

    post_done || post_process || { log "FAILED: post-processing ($POST)"; return 1; }
    rm -f "$OUT.lock"
    log "[OK] $RESOURCE: $DEST"
}

# -----------------------------------------------------------------------------
# final check: task 1 queues ONE job that runs after the whole array has ended
# (NTASKS array tasks + 1 final job = never more than the job limit)
# -----------------------------------------------------------------------------
submit_finalize() {
    local name=Step00_Download_finalize
    mkdir "$LOGS/.final_${SLURM_ARRAY_JOB_ID}" 2>/dev/null || return 0   # once per array
    sbatch --parsable ${SBATCH_COMMON[@]+"${SBATCH_COMMON[@]}"} --dependency="afterany:${SLURM_ARRAY_JOB_ID}" \
        --job-name="$name" \
        --output="$LOGS/${name}_%j.log" --error="$LOGS/${name}_%j.err" \
        --export=ALL,G2M_MODE=finalize --wrap="bash '$SCRIPT'" \
        && log "queued the final check (runs after array ${SLURM_ARRAY_JOB_ID})"
}

finalize() {
    [[ -n "$PYTHON" ]] || { log "ERROR: no python to run the final check"; return 1; }
    "$PYTHON" "$PLANNER" --finalize @@PLANNER_ARGS@@
}

# -----------------------------------------------------------------------------
# dispatch
# -----------------------------------------------------------------------------
case "${1:-}" in
    task)     run_task "${2:?row number}"; exit $? ;;
    finalize) finalize; exit $? ;;
    local)
        log "local mode: $TOTAL files, ${2:-4} at a time"
        seq 1 "$TOTAL" | xargs -P "${2:-4}" -I{} bash "$SCRIPT" task {}
        finalize; exit $? ;;
esac

if [[ "${G2M_MODE:-}" == "finalize" ]]; then finalize; exit $?; fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Submit with:  sbatch @@SCRIPT_NAME@@      (or without SLURM: bash @@SCRIPT_NAME@@ local 4)"
    exit 2
fi

[[ "$SLURM_ARRAY_TASK_ID" == "1" && -z "${G2M_CONTINUATION:-}" ]] && submit_finalize
[[ -n "${G2M_CONTINUATION:-}" ]] && log "continuation job ${G2M_CONTINUATION}/$MAX_CONTINUATIONS of task $SLURM_ARRAY_TASK_ID"

# Task k downloads rows k, k+NTASKS, k+2*NTASKS, ... (one row when TOTAL <= NTASKS).
# Interleaving spreads the large files over different tasks.
failed=0
CONTINUE_LATER=0
for (( row = SLURM_ARRAY_TASK_ID; row <= TOTAL; row += NTASKS )); do
    (( CONTINUE_LATER )) && break          # out of walltime: the continuation job does the rest
    run_task "$row"; rc=$?
    (( rc == 0 || rc == 3 )) || failed=$(( failed + 1 ))
done
if (( CONTINUE_LATER )); then
    continue_in_new_job || failed=$(( failed + 1 ))
fi
if (( failed > 0 )); then
    log "$failed file(s) of this task failed (re-plan and resubmit to retry them)"
    exit 1
fi
'''


def render_script(total: int, settings: dict) -> str:
    values = settings["values"]
    ntasks = min(total, settings["chunk"])        # at most max_jobs - 1 tasks (+ 1 final check)
    throttle = f"%{settings['parallel']}" if settings["parallel"] and settings["parallel"] > 0 else ""
    directives = {**values, "job-name": Path(args.output).stem, "array": f"1-{ntasks}{throttle}",
                  "output": f"{LOGS}/Step00_Download_%A_%a.log", "error": f"{LOGS}/Step00_Download_%A_%a.err"}
    order = ["job-name", "account", "partition", "qos", "reservation", "constraint", "nodes", "ntasks",
             "cpus-per-task", "mem", "time", "array", "output", "error"]
    lines = [f"#SBATCH --{k}={directives[k]}" for k in order if directives.get(k)]
    lines += [f"#SBATCH --{k}={v}" for k, v in directives.items() if k not in order and v]
    # The final-check job gets exactly the same settings as the download tasks.
    common = [f"--{k}={v}" for k, v in directives.items() if k not in PER_JOB_KEYS and v]

    planner_args = []
    if args.output != DEFAULT_OUTPUT:
        planner_args.append(f"--output {args.output}")
    if args.opentargets_release != OPEN_TARGETS_RELEASE:
        planner_args.append(f"--opentargets-release {args.opentargets_release}")

    replacements = {
        "@@SHEBANG@@": settings["shebang"],
        "@@DIRECTIVES@@": "\n".join(lines),
        "@@DATE@@": time.strftime("%Y-%m-%d %H:%M"),
        "@@TOTAL@@": str(total),
        "@@NTASKS@@": str(ntasks),
        "@@QUEUED@@": str(ntasks + 1),
        "@@MAX_JOBS@@": str(settings["max_jobs"]),
        "@@PER_TASK@@": str(-(-total // ntasks)),
        "@@THROTTLE@@": throttle,
        "@@SBATCH_COMMON@@": " ".join(f'"{c}"' for c in common),
        "@@ROOT@@": str(ROOT),
        "@@LOGS@@": str(LOGS),
        "@@SCRIPT_NAME@@": args.output,
        "@@PLANNER@@": str(Path(__file__).resolve()),
        "@@PLANNER_ARGS@@": " ".join(planner_args),
        "@@MANIFEST@@": str(MANIFEST),
        "@@OUTPUT_ARG@@": "" if args.output == DEFAULT_OUTPUT else f" --output {args.output}",
    }
    text = SCRIPT_TEMPLATE
    for key, value in replacements.items():
        text = text.replace(key, value)
    return text


# =============================================================================
# MANIFESTS
# =============================================================================

MANIFEST_COLUMNS = ["INDEX", "RESOURCE", "GROUP", "SOURCE", "DEST", "EXPECTED_SIZE", "CHECK", "POST",
                    "POST_DIR", "POST_MARKER"]


def write_tsv(path: Path, header: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with open(temporary, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    temporary.replace(path)


def read_manifest() -> list[dict]:
    if not MANIFEST.is_file():
        return []
    with open(MANIFEST, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [item(r["RESOURCE"], r["GROUP"], r["SOURCE"], r["DEST"], check=r["CHECK"], post=r["POST"],
                 post_dir=r["POST_DIR"], post_marker=r["POST_MARKER"],
                 size=None if r["EXPECTED_SIZE"] == "-" else int(r["EXPECTED_SIZE"])) for r in rows]


def finalize_opentargets() -> list[str]:
    """Write <dataset>/.complete only when every expected Open Targets file has its exact size."""
    if not OT_EXPECTED.is_file():
        return []
    expected: dict[str, list[tuple[str, int]]] = {}
    with open(OT_EXPECTED, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expected.setdefault(row["DATASET"], []).append((row["RELATIVE_PATH"], int(row["SIZE_BYTES"])))
    messages = []
    for dataset, files in expected.items():
        directory = ROOT / OT_ROOT / dataset
        bad = [rel for rel, size in files if not local_file_valid(directory / rel, size)]
        if bad:
            messages.append(f"Open Targets {dataset}: {len(bad)}/{len(files)} files missing or wrong size; "
                            ".complete NOT written")
        else:
            (directory / ".complete").write_text(f"files={len(files)}\nverified_against={OT_EXPECTED.name}\n",
                                                 encoding="utf-8")
            messages.append(f"Open Targets {dataset}: {len(files)}/{len(files)} files verified; .complete written")
    return messages


# =============================================================================
# FINALIZE (run by the last SLURM job, or by hand)
# =============================================================================

def finalize() -> int:
    banner("GWAS2m STEP 00 - VERIFY SHARED RESOURCE DOWNLOADS")
    entries = read_manifest()
    if not entries:
        print(f"No task manifest ({MANIFEST}); nothing to verify.")
        return 0
    rows, missing = [], []
    for index, entry in enumerate(entries, start=1):
        path = ROOT / entry["dest"]
        ok = is_complete(entry) if entry["post"] != "-" else local_file_valid(path, entry["size"])
        stamp = (datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()
                 if path.is_file() else "")
        rows.append([index, entry["resource"], entry["group"], "OK" if ok else "MISSING",
                     path.stat().st_size if path.is_file() else "", stamp, entry["dest"]])
        if not ok:
            missing.append(entry)
    write_tsv(REPORT, ["INDEX", "RESOURCE", "GROUP", "STATUS", "SIZE_BYTES", "MTIME_UTC", "FILE"], rows)

    for message in finalize_opentargets():
        print(message)

    print(f"\nFiles verified : {len(entries) - len(missing)} / {len(entries)}")
    print(f"Report         : {REPORT}")
    if missing:
        print(f"\n{len(missing)} file(s) missing or failed:")
        for entry in missing[:50]:
            print(f"  - {entry['resource']:<30} {entry['dest']}")
        print("\nRe-plan only what is missing and resubmit:")
        print(f"  envs/pipeline/bin/python {Path(__file__).name} && sbatch {args.output}")
        return 1
    print("\nAll planned files downloaded and verified.")
    print("Next: envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py   (readability + headers)")
    return 0


# =============================================================================
# STATUS (read-only, offline: from the saved plans + the files on disk)
# =============================================================================

def file_state(row: dict) -> tuple[str, int]:
    """(DONE | PARTIAL | MISSING, bytes on disk) for one planned file."""
    dest = ROOT / row["DEST"]
    post_dir, marker = row.get("POST_DIR") or "-", row.get("POST_MARKER") or "-"
    expected = int(row["SIZE_BYTES"]) if str(row.get("SIZE_BYTES", "")).isdigit() else None
    if post_dir != "-" and (ROOT / post_dir / marker).is_file():
        return "DONE", expected or (dest.stat().st_size if dest.is_file() else 0)
    if local_file_valid(dest, expected):
        return "DONE", dest.stat().st_size
    partial = [p for p in (dest, Path(str(dest) + ".part")) if p.is_file() and p.stat().st_size > 0]
    if partial or Path(str(dest) + ".aria2").exists():
        return "PARTIAL", max((p.stat().st_size for p in partial), default=0)
    return "MISSING", 0


def slurm_jobs() -> str:
    if not shutil.which("squeue"):
        return ""
    try:
        out = subprocess.run(["squeue", "--me", "-h", "-o", "%j|%T"], capture_output=True, text=True,
                             timeout=60).stdout
    except (OSError, subprocess.SubprocessError):
        return ""
    counts: dict = {}
    for line in out.splitlines():
        name, _, state = line.partition("|")
        if name.startswith("Step00_Download"):
            counts[state] = counts.get(state, 0) + 1
    return ", ".join(f"{state.lower()} {n}" for state, n in sorted(counts.items())) or "none in the queue"


def plan_script(plan: Path) -> str:
    """download_resources_<X>_plan.tsv was written with --output <X>.sh."""
    stem = plan.name[len("download_resources"):-len("_plan.tsv")].lstrip("_")
    return f"{stem or Path(DEFAULT_OUTPUT).stem}.sh"


def status() -> int:
    # Oldest first: a file planned again by a newer plan takes that plan's row.
    plans = sorted(MANIFESTS.glob("download_resources*_plan.tsv"), key=lambda p: p.stat().st_mtime)
    if not plans:
        print("No plan yet. Run: envs/pipeline/bin/python Step00_Download_Resources.py")
        return 1
    rows: dict = {}
    for plan in plans:
        with open(plan, encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row["_PLAN"] = plan
                rows[row["DEST"]] = row

    banner("GWAS2m STEP 00 - DOWNLOAD STATUS (read-only)")
    jobs = slurm_jobs()
    if jobs:
        print(f"SLURM jobs     : {jobs}")

    # One line per plan: when it was made, its script, and what of it is unfinished.
    print(f"\n{'PLAN (setup_manifests/)':<56}{'MADE':<18}{'FILES':>7}{'NOT DONE':>10}  SCRIPT")
    states = {dest: file_state(row)[0] for dest, row in rows.items()}
    for plan in plans:
        own = [d for d, r in rows.items() if r["_PLAN"] == plan]
        made = datetime.fromtimestamp(plan.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        print(f"{plan.name:<56}{made:<18}{len(own):>7}{sum(states[d] != 'DONE' for d in own):>10}  "
              f"{plan_script(plan)}")

    groups: dict = {}
    pending = []
    for row in rows.values():
        state, size = file_state(row)
        group = row["GROUP"].split("/", 1)[0]
        g = groups.setdefault(group, {"files": 0, "DONE": 0, "PARTIAL": 0, "MISSING": 0, "bytes": 0, "planned": 0})
        g["files"] += 1
        g[state] += 1
        g["bytes"] += size
        g["planned"] += int(row["SIZE_BYTES"]) if str(row.get("SIZE_BYTES", "")).isdigit() else 0
        if state != "DONE":
            pending.append((state, size, row))

    print(f"\n{'GROUP':<22}{'FILES':>7}{'DONE':>7}{'PARTIAL':>9}{'MISSING':>9}{'DOWNLOADED':>14}{'PLANNED':>12}{'%':>7}")
    total = {"files": 0, "DONE": 0, "PARTIAL": 0, "MISSING": 0, "bytes": 0, "planned": 0}
    for group, g in groups.items():
        for key in total:
            total[key] += g[key]
        share = f"{100 * g['DONE'] / g['files']:.0f}%" if g["files"] else "-"
        print(f"{group:<22}{g['files']:>7}{g['DONE']:>7}{g['PARTIAL']:>9}{g['MISSING']:>9}"
              f"{human(g['bytes']) if g['bytes'] else '0 B':>14}{human(g['planned']):>12}{share:>7}")
    print("-" * 87)
    print(f"{'TOTAL':<22}{total['files']:>7}{total['DONE']:>7}{total['PARTIAL']:>9}{total['MISSING']:>9}"
          f"{human(total['bytes']) if total['bytes'] else '0 B':>14}{human(total['planned']):>12}"
          f"{100 * total['DONE'] / max(total['files'], 1):>6.0f}%")

    if pending:
        print(f"\nNot finished ({len(pending)}; first 40) - with the script of the plan that lists it:")
        for state, size, row in sorted(pending, key=lambda x: (x[0] != "PARTIAL", x[2]["DEST"]))[:40]:
            print(f"  {state:<8}{human(size) if size else '':>10}  {row['DEST']}   [{plan_script(row['_PLAN'])}]")
        print("\nPartial files resume automatically. To (re)download what is missing, re-plan with the")
        print("SAME --groups / --output as that plan, then sbatch its script, e.g.:")
        print("  envs/pipeline/bin/python Step00_Download_Resources.py --groups 1000g30x --output Step00_Download_1000G30x.sh")
        print("  sbatch Step00_Download_1000G30x.sh")
        print("Re-planning replaces the plan's rows with what the current Step00 needs; files only an older")
        print("Step00 version planned (e.g. other folder layouts) then drop out of this status.")
    else:
        print("\nEvery planned file is downloaded.")
        print("Readability check: envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py")
    return 0


# =============================================================================
# MAIN (plan)
# =============================================================================

def main() -> int:
    if args.status:
        return status()
    if args.finalize:
        return finalize()

    banner("GWAS2m STEP 00 - PLAN SHARED RESOURCE DOWNLOADS")
    print(f"Project root : {ROOT}")
    print(f"Groups       : {', '.join(SELECTED_GROUPS) or '(none)'}")
    settings = slurm_settings()

    print("\nCollecting file list ...")
    notes: list[str] = []
    entries = (static_items(notes) + eqtl_catalogue_items(notes) + opentargets_items()
               + catalogue_qtl_items() + molecular_qtl_items() + annotation_items())
    # One task per FILE: the same file can be wanted by two groups (e.g. a dense file
    # that the repair group and qtl_catalogue both list) - keep the first.
    entries = list({e["dest"]: e for e in reversed(entries)}.values())[::-1]

    if not args.no_url_check and entries:
        print(f"\nChecking {len(entries)} URLs and sizes ...")
        failed = []
        with ThreadPoolExecutor(max_workers=16) as pool:
            for entry, error in pool.map(probe, entries):
                if error:
                    failed.append((entry, error))
        if failed:
            print("\nERROR: these URLs are not reachable:")
            for entry, error in failed:
                print(f"  - {entry['resource']:<30} {error}\n      {entry['url']}")
            print("\nNothing was written. Fix the network/URLs, or leave the group out with --skip-groups.")
            return 1
        print("[OK] all URLs reachable")

    group_of = lambda e: e["group"].split("/", 1)[0]  # noqa: E731  (opentargets/x -> opentargets)
    complete = {id(e) for e in entries if is_complete(e)}
    missing = [e for e in entries if id(e) not in complete]
    pending = entries if args.force else missing

    # Groups used as one set (1000g: Step04 needs every VCF + the sample panel):
    # when any file of the set is missing, EVERY file of it gets a task, so the
    # array always covers the whole set; a task whose file is valid only verifies it.
    incomplete_sets = {group_of(e) for e in missing if group_of(e) in COMPLETE_SET_GROUPS}
    if incomplete_sets:
        tasks = {id(e) for e in pending}
        pending = [e for e in entries if id(e) in tasks or group_of(e) in incomplete_sets]

    print(f"\n{'RESOURCE':<32}{'EXPECTED':>9}{'VALID':>8}{'TASKS':>8}")
    for resource in dict.fromkeys(e["resource"] for e in entries):
        n_all = sum(e["resource"] == resource for e in entries)
        n_valid = sum(e["resource"] == resource and id(e) in complete for e in entries)
        n_tasks = sum(e["resource"] == resource for e in pending)
        print(f"{resource:<32}{n_all:>9}{n_valid:>8}{n_tasks:>8}")
    for note in notes:
        print(f"NOTE: {note}")
    for group in sorted(incomplete_sets):
        print(f"NOTE: {group}: files are missing, so all {sum(group_of(e) == group for e in entries)} files of "
              "the set get a task (valid ones are only verified)")

    print(f"\n{'GROUP':<22}{'FILES':>7}{'TO DOWNLOAD':>13}{'SIZE':>12}{'TASKS':>8}")
    for group in dict.fromkeys(group_of(e) for e in entries):
        todo = [e for e in missing if group_of(e) == group] if not args.force \
            else [e for e in pending if group_of(e) == group]
        print(f"{group:<22}{sum(group_of(e) == group for e in entries):>7}{len(todo):>13}"
              f"{human(sum(e['size'] or 0 for e in todo)):>12}{sum(group_of(e) == group for e in pending):>8}")

    MANIFESTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    write_tsv(PLAN, ["RESOURCE", "GROUP", "STATUS", "SIZE_BYTES", "SOURCE", "DEST", "POST_DIR", "POST_MARKER"],
              [[e["resource"], e["group"], "PLANNED" if e in pending else "VALID", e["size"] or "",
                e["url"], e["dest"], e["post_dir"], e["post_marker"]] for e in entries])

    ot = [e for e in entries if e["group"].startswith("opentargets/")]
    if ot and all(e["size"] for e in ot):
        write_tsv(OT_EXPECTED, ["DATASET", "RELATIVE_PATH", "SIZE_BYTES"],
                  [[e["group"].split("/", 1)[1], e["url"].split(f"/output/{e['group'].split('/', 1)[1]}/", 1)[1],
                    e["size"]] for e in ot])
        for message in finalize_opentargets():
            if "written" in message and ".complete NOT" not in message:
                print(message)

    if not pending:
        if MANIFEST.exists():
            MANIFEST.unlink()
        banner("NOTHING TO DOWNLOAD")
        print("Every planned file is already valid.")
        print("Next: envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py   (readability + headers)")
        return 0

    write_tsv(MANIFEST, MANIFEST_COLUMNS, [
        [index, e["resource"], e["group"], e["url"], e["dest"],
         "-" if e["mutable"] or not e["size"] else e["size"], e["check"], e["post"], e["post_dir"], e["post_marker"]]
        for index, e in enumerate(pending, start=1)
    ])

    total_bytes = sum(e["size"] or 0 for e in pending)
    # tar-rm archives are deleted after extraction, so they need no extra space at the end.
    extract_bytes = sum(e["size"] or 0 for e in pending if e["post"] in ("tar", "fasta"))
    free = shutil.disk_usage(ROOT).free

    # A single file this big may not finish within one walltime; re-running resumes it.
    huge = [e for e in pending if (e["size"] or 0) > 300e9]
    walltime = str(settings["values"].get("time", ""))
    long_walltime = "-" in walltime and walltime.split("-", 1)[0].isdigit() and int(walltime.split("-", 1)[0]) >= 2
    if huge and not long_walltime:
        print(f"\n[NOTE] {len(huge)} file(s) larger than 300 GB (largest {human(max(e['size'] for e in huge))}).")
        print(f"       Each needs download + MD5 + extraction in one task; with --time "
              f"{settings['values'].get('time')} it may time out.")
        print("       A timed-out download resumes when you re-run this planner, or give it longer, e.g.:")
        print("         envs/pipeline/bin/python Step00_Download_Resources.py --skip-groups meqtl          # everything else")
        print("         envs/pipeline/bin/python Step00_Download_Resources.py --groups meqtl --time 7-00:00:00 "
              "--output Step00_Download_meQTL.sh")

    script = ROOT / args.output
    with open(script, "w", newline="\n", encoding="utf-8") as handle:
        handle.write(render_script(len(pending), settings))
    script.chmod(0o755)

    ntasks = min(len(pending), settings["chunk"])
    per_task = -(-len(pending) // ntasks)
    banner("SLURM SCRIPT READY")
    print(f"Script         : {script}")
    print(f"Task manifest  : {MANIFEST}  ({len(pending)} files)")
    print(f"Download size  : {human(total_bytes)}  (+ ~{human(extract_bytes)} after extraction)")
    print(f"Free disk      : {human(free)}" + ("   [WARN] not enough for everything planned"
                                                if total_bytes + extract_bytes > free else ""))
    print(f"Jobs           : 1 array of {ntasks} tasks + 1 final check = {ntasks + 1} jobs "
          f"(limit {settings['max_jobs']})")
    print(f"Files per task : " + ("1" if per_task == 1 else
                                  f"up to {per_task} ({len(pending)} files > {ntasks} tasks, so each task "
                                  f"downloads several files one after another)"))
    print("\nSBATCH per download task:")
    for key in ("partition", "qos", "account", "cpus-per-task", "mem", "time"):
        if settings["values"].get(key):
            print(f"  --{key:<14}{settings['values'][key]:<14}<- {settings['source'][key]}")
    if settings["parallel"]:
        print(f"  {'throttle':<16}%{settings['parallel']}")
    print(f"\nSubmit         : sbatch {args.output}")
    print("Monitor        : squeue --me")
    print(f"Without SLURM  : bash {args.output} local 4")
    return 0


if __name__ == "__main__":
    sys.exit(main())
