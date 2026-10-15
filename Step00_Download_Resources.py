#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2m STEP 00 - DOWNLOAD SHARED RESOURCES (SLURM: ONE ARRAY TASK PER FILE)
===============================================================================

    python Step00_Download_Resources.py        # plan + write the SLURM script
    sbatch Step00_Download_Resources.sh        # download everything

This script only PLANS; it never downloads data itself. It:

  1. lists every shared resource FILE the pipeline needs, with the paths,
     URLs and validation rules of gwas2m_resources.py and
     gwas2m_setup_downloads.py (used directly when they sit next to this
     script; an identical built-in copy is used otherwise);
  2. skips files that are already VALID - non-empty, no <file>.aria2 control
     file and, when the remote size is known, exactly that size - and archives
     already extracted / FASTA already indexed (re-running only re-plans what
     is missing, partial or corrupt);
  3. checks every URL and records its size (the size is enforced after
     download);
  4. writes setup_manifests/download_resources.tsv (one row = one file) and
     Step00_Download_Resources.sh.

SLURM settings, highest priority first:
     command-line options
  >  config/slurm.yaml, stage "setup_download" (via gwas2m_config.py)
  >  the #SBATCH header of HPC_CPU.sh
  >  built-in defaults (2 CPUs, 16G, 24:00:00)

In the generated script ONE ARRAY TASK DOWNLOADS EXACTLY ONE FILE. No more
than --max-jobs (default: array_limit = 1000) jobs are ever queued: larger
plans are split into chained arrays of (max-jobs - 1) tasks; a 1-task
controller job submits the next array only after the previous one finished,
and a last job verifies everything (python ... --finalize).

Per task: resumable download (aria2c, else curl/wget into <file>.part) ->
exact-size check -> integrity check (gzip -t / tar -t / parquet footer) ->
post-processing (FASTA decompress + samtools faidx; tar extraction for GTEx
and the VEP cache). A corrupt file is deleted so the next run fetches it again.

Resource groups (--groups / --skip-groups):

  gwas_catalog    GWAS Catalog studies + ancestry tables (Step01)
  genome          GENCODE 50 GRCh38 FASTA (+ .fai) and GTF
  1000g           1000 Genomes GRCh38 VCF + TBI + panel (only while the
                  PGEN panels still have to be built, as gwas2m plans it)
  vep             Ensembl VEP 116 GRCh38 indexed cache (extracted)
  gtex            GTEx v11 eQTL / sQTL / SuSiE archives (extracted) + references
  eqtl_catalogue  eQTL Catalogue metadata tables + repair of dense / SuSiE
                  files referenced locally or listed in
                  setup_manifests/eqtl_catalogue_datasets.txt
  opentargets     Open Targets Platform study / credible_set / l2g_prediction /
                  colocalisation (one task per parquet file)

NOT done here: derived resources (1000G PGEN panels, sample lists, Pangolin
DB) -> python Step00_Check_Resources.py --resources-only (it detects these
downloads and skips them); molecular QTL -> Setup_10_*; ENCODE -> Setup_07D_*;
phenotype GWAS summary statistics -> Step01 / Step02.

Usage:

    python Step00_Download_Resources.py
    python Step00_Download_Resources.py --groups genome vep gtex
    python Step00_Download_Resources.py --skip-groups opentargets
    python Step00_Download_Resources.py --max-jobs 500 --parallel 50
    python Step00_Download_Resources.py --account a_myproject
    python Step00_Download_Resources.py --finalize      # re-check + report

Without SLURM (laptop / workstation):

    bash Step00_Download_Resources.sh local 4           # 4 files at a time

===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


# =============================================================================
# ROOT / SHARED MODULES
# =============================================================================

ROOT = Path.cwd().resolve()
HERE = Path(__file__).resolve().parent

sys.path.insert(0, str(HERE))
try:
    import gwas2m_resources as gr  # paths, URLs, validators (single source of truth)
except Exception:
    gr = None
try:
    import gwas2m_config as gc     # config/slurm.yaml
except Exception:
    gc = None

MANIFESTS = ROOT / "setup_manifests"
LOGS = ROOT / "logs" / "resource_setup"
MANIFEST = MANIFESTS / "download_resources.tsv"
PLAN = MANIFESTS / "download_resources_plan.tsv"
REPORT = MANIFESTS / "download_resources_report.tsv"


# =============================================================================
# VERSIONS / SOURCES (taken from gwas2m_resources when available)
# =============================================================================

GENCODE_RELEASE = getattr(gr, "GENCODE_RELEASE", "50")
VEP_RELEASE = getattr(gr, "VEP_RELEASE", "116")
GTEX_RELEASE = getattr(gr, "GTEX_RELEASE", "v11")
OPEN_TARGETS_RELEASE = getattr(gr, "OPEN_TARGETS_RELEASE", "26.09")
LD_ANCESTRIES = tuple(getattr(gr, "LD_ANCESTRIES", ("EUR", "AFR", "EAS", "SAS", "AMR")))

GENCODE_BASE = getattr(gr, "GENCODE_BASE",
                       f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}")
ONE_K_BASE = getattr(gr, "ONE_K_BASE",
                     "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
                     "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL")
GTEX_BASE = getattr(gr, "GTEX_BASE", "https://storage.googleapis.com/adult-gtex")
EQTL_RAW = "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master"
OT_HTTPS = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"
OT_DATASETS = ("study", "credible_set", "l2g_prediction", "colocalisation")

EQTL_ROOT = Path("resources") / "coloc" / "eqtl_catalogue"
EQTL_METADATA = {  # same files/locations as gwas2m_setup_downloads / Step11
    "tabix_ftp_paths.tsv": f"{EQTL_RAW}/tabix/tabix_ftp_paths.tsv",
    "tabix_ftp_paths_imported.tsv": f"{EQTL_RAW}/tabix/tabix_ftp_paths_imported.tsv",
    "dataset_metadata_r8_beta.tsv": f"{EQTL_RAW}/data_tables/dataset_metadata_r8_beta.tsv",
    "population_assignments.tsv": f"{EQTL_RAW}/data_tables/population_assignments.tsv",
}
EXTRA_EQTL_DATASETS = MANIFESTS / "eqtl_catalogue_datasets.txt"
PARTIAL_SUFFIXES = (".aria2", ".part", ".tmp")

GROUPS = ["gwas_catalog", "genome", "1000g", "vep", "gtex", "eqtl_catalogue", "opentargets"]


# =============================================================================
# ARGUMENTS
# =============================================================================

parser = argparse.ArgumentParser(
    description="Plan shared-resource downloads and write a SLURM array script (one task per file)."
)
parser.add_argument("--groups", nargs="+", choices=GROUPS, default=GROUPS,
                    help="Resource groups to include (default: all).")
parser.add_argument("--skip-groups", nargs="+", choices=GROUPS, default=[],
                    help="Resource groups to leave out.")
parser.add_argument("--template", default="HPC_CPU.sh",
                    help="SLURM template whose #SBATCH header is the base (default: HPC_CPU.sh).")
parser.add_argument("--slurm-config", default=None,
                    help="SLURM config (default: config/slurm.yaml, or $GWAS2M_SLURM_CONFIG).")
parser.add_argument("--output", default="Step00_Download_Resources.sh",
                    help="Generated SLURM script (default: Step00_Download_Resources.sh).")
parser.add_argument("--max-jobs", type=int, default=None,
                    help="Most jobs ever queued at once incl. 1 controller (default: array_limit, else 1000).")
parser.add_argument("--parallel", type=int, default=None,
                    help="Cap on simultaneously RUNNING tasks, array %%N (default: max_parallel, else none).")
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
args = parser.parse_args()

SELECTED_GROUPS = [g for g in GROUPS if g in args.groups and g not in args.skip_groups]
OT_RELEASE = args.opentargets_release
OT_ROOT = Path("resources") / "opentargets" / OT_RELEASE
OT_EXPECTED = MANIFESTS / f"opentargets_{OT_RELEASE}_expected.tsv"   # gwas2m_setup_downloads format


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
    """gwas2m_setup_downloads rule: non-empty, no .aria2, exact size when known."""
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if Path(str(path) + ".aria2").exists():
        return False
    return expected_size is None or path.stat().st_size == expected_size


# =============================================================================
# FILE ITEMS  (one item = one file = one SLURM array task)
#   post  : - | fasta (decompress + faidx) | tar (extract into post_dir)
#   check : - | gzip | tar | parquet       (integrity test after download)
#   mutable: "latest"/"master" URLs whose size may change before download,
#            so their size is reported but not enforced.
# =============================================================================

def item(resource, group, url, dest, *, check="-", post="-", post_dir="-", post_marker="-",
         mutable=False, size=None):
    return dict(resource=resource, group=group, url=url, dest=str(dest), check=check, post=post,
                post_dir=str(post_dir), post_marker=post_marker, mutable=mutable, size=size)


def registry(resource_id: str, key: str, default: str) -> str:
    if gr is not None and resource_id in getattr(gr, "RESOURCE_BY_ID", {}):
        return gr.RESOURCE_BY_ID[resource_id].get(key) or default
    return default


def plan_preparation() -> dict:
    """Same rule as gwas2m_setup_downloads.plan_preparation(): which 1000G inputs are still needed."""
    g = ROOT / "resources" / "1000G"

    def panel_ok(prefix: Path) -> bool:
        return all(local_file_valid(Path(str(prefix) + ext)) for ext in (".pgen", ".pvar", ".psam"))

    all_missing = [c for c in range(1, 23) if not panel_ok(g / "ALL" / f"chr{c}")]
    population_missing = [(anc, c) for anc in LD_ANCESTRIES for c in range(1, 23)
                          if not panel_ok(g / anc / f"chr{c}_{anc}_GRCh38")]
    samples_needed = any(not local_file_valid(g / anc / "samples.txt") for anc in LD_ANCESTRIES)
    return {
        "all_chromosomes": all_missing if population_missing else [],
        "samples": bool(population_missing) and samples_needed,
        "panels_missing": len(population_missing),
    }


def static_items(notes: list[str]) -> list[dict]:
    items: list[dict] = []

    # --- GWAS Catalog metadata (Step01 discovery) ---
    items += [
        item("GWAS_CATALOG_STUDIES", "gwas_catalog",
             registry("GWAS_CATALOG_STUDIES", "url",
                      "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
                      "gwas-catalog-download-studies-v1.0.3.1.txt"),
             registry("GWAS_CATALOG_STUDIES", "path", "reference_metadata/gwas_catalog/gwas_catalog_studies.tsv"),
             mutable=True),
        item("GWAS_CATALOG_ANCESTRY", "gwas_catalog",
             registry("GWAS_CATALOG_ANCESTRY", "url",
                      "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
                      "gwas-catalog-download-ancestries-v1.0.3.1.txt"),
             registry("GWAS_CATALOG_ANCESTRY", "path", "reference_metadata/gwas_catalog/gwas_catalog_ancestry.tsv"),
             mutable=True),
    ]

    # --- GENCODE genome + annotation ---
    fasta = "resources/genome/GRCh38/GRCh38.primary_assembly.genome.fa"
    gtf = registry("GENCODE_GTF", "path", f"resources/gencode/release{GENCODE_RELEASE}/"
                                          f"gencode.v{GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz")
    items += [
        item("GRCH38_FASTA", "genome", f"{GENCODE_BASE}/GRCh38.primary_assembly.genome.fa.gz",
             fasta + ".gz", check="gzip", post="fasta"),
        item("GENCODE_GTF", "genome",
             registry("GENCODE_GTF", "url",
                      f"{GENCODE_BASE}/gencode.v{GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"),
             gtf, check="gzip"),
    ]

    # --- 1000 Genomes: raw files are only inputs to the PGEN panels ---
    prep = plan_preparation()
    raw, meta = Path("resources/1000G/raw"), Path("resources/1000G/metadata")
    for chrom in prep["all_chromosomes"]:
        name = f"ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
        items.append(item("1000G_RAW_VCF", "1000g", f"{ONE_K_BASE}/{name}", raw / name, check="gzip"))
        items.append(item("1000G_RAW_VCF", "1000g", f"{ONE_K_BASE}/{name}.tbi", raw / f"{name}.tbi"))
    if prep["samples"]:
        items.append(item("1000G_PANEL_METADATA", "1000g",
                          "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
                          "integrated_call_samples_v3.20130502.ALL.panel",
                          meta / "integrated_call_samples_v3.20130502.ALL.panel"))
    if prep["all_chromosomes"]:
        for name in ("20190312_biallelic_SNV_and_INDEL_README.txt", "20190312_biallelic_SNV_and_INDEL_MANIFEST.txt"):
            items.append(item("1000G_PANEL_METADATA", "1000g", f"{ONE_K_BASE}/{name}", meta / name))
    if "1000g" in SELECTED_GROUPS and not prep["panels_missing"]:
        notes.append("1000g: all population PGEN panels are built, raw 1000G files are not needed")

    # --- Ensembl VEP cache (extracted into resources/vep/cache) ---
    vep_url = registry("VEP_CACHE", "url",
                       f"https://ftp.ensembl.org/pub/release-{VEP_RELEASE}/variation/indexed_vep_cache/"
                       f"homo_sapiens_vep_{VEP_RELEASE}_GRCh38.tar.gz")
    items.append(item("VEP_CACHE", "vep", vep_url, f"resources/vep/{VEP_RELEASE}/{basename(vep_url)}",
                      check="tar", post="tar", post_dir="resources/vep/cache",
                      post_marker=f".vep{VEP_RELEASE}_GRCh38.complete"))

    # --- GTEx v11 (same items/locations as gwas2m_setup_downloads.GTEX_ITEMS) ---
    for rid, path, sub, name in [
        ("GTEX_EQTL", "bulk-qtl/v11/single-tissue-cis-qtl/GTEx_Analysis_v11_eQTL.tar", "qtl", "eQTL"),
        ("GTEX_SQTL", "bulk-qtl/v11/single-tissue-cis-qtl/GTEx_Analysis_v11_sQTL.tar", "qtl", "sQTL"),
        ("GTEX_EQTL_SUSIE", "bulk-qtl/v11/susie-qtl/GTEx_Analysis_v11_eQTL_SuSiE.tar", "susie", "eQTL_SuSiE"),
        ("GTEX_SQTL_SUSIE", "bulk-qtl/v11/susie-qtl/GTEx_Analysis_v11_sQTL_SuSiE.tar", "susie", "sQTL_SuSiE"),
    ]:
        directory = f"resources/gtex/{GTEX_RELEASE}/{sub}"
        items.append(item(rid, "gtex", f"{GTEX_BASE}/{path}", f"{directory}/{basename(path)}", check="tar",
                          post="tar", post_dir=f"{directory}/{name}", post_marker=".complete"))
    for rid, default in [
        ("GTEX_VARIANT_LOOKUP", f"{GTEX_BASE}/references/v11/reference-tables/"
                                "GTEx_Analysis_2021-02-11_v11_WholeGenomeSeq_953Indiv.lookup_table.txt.gz"),
        ("GTEX_GENCODE47", f"{GTEX_BASE}/references/v11/reference-tables/gencode.v47.genes.gtf"),
    ]:
        url = registry(rid, "url", default)
        items.append(item(rid, "gtex", url, f"resources/gtex/{GTEX_RELEASE}/reference/{basename(url)}",
                          check="gzip" if url.endswith(".gz") else "-"))

    # --- eQTL Catalogue metadata tables ---
    for name, url in EQTL_METADATA.items():
        items.append(item("EQTL_CATALOGUE_METADATA", "eqtl_catalogue", url, EQTL_ROOT / name, mutable=True))

    return [i for i in items if i["group"] in SELECTED_GROUPS]


# -----------------------------------------------------------------------------
# eQTL Catalogue dense + provider SuSiE (repair; same rules as gwas2m_setup_downloads)
# -----------------------------------------------------------------------------

def read_catalogue() -> tuple[dict, dict]:
    aliases = {"dataset_id": ["dataset_id"], "ftp_path": ["ftp_path", "sumstats_path"],
               "ftp_cs_path": ["ftp_cs_path"], "ftp_lbf_path": ["ftp_lbf_path"]}
    by_basename, by_dataset = {}, {}
    for name in ("tabix_ftp_paths.tsv", "tabix_ftp_paths_imported.tsv"):  # later file wins (as Step11)
        path = ROOT / EQTL_ROOT / name
        if not path.is_file():
            continue
        with open(path, encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
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
    by_basename, by_dataset = read_catalogue()
    extra = []
    if EXTRA_EQTL_DATASETS.is_file():
        extra = [x.split("#", 1)[0].strip() for x in EXTRA_EQTL_DATASETS.read_text(encoding="utf-8").splitlines()
                 if x.split("#", 1)[0].strip()]
    if not by_dataset:
        if extra:
            notes.append("eqtl_catalogue: datasets are listed but the metadata tables are not downloaded yet; "
                         "re-run this planner after this download to plan them")
        return []

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
                                  (provider / dataset / basename(url)).relative_to(ROOT), check="gzip"))
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


# =============================================================================
# STATUS
# =============================================================================

def post_done(entry: dict) -> bool:
    if entry["post"] == "fasta":
        fasta = ROOT / entry["dest"][: -len(".gz")]
        return local_file_valid(fasta) and local_file_valid(Path(str(fasta) + ".fai"))
    if entry["post"] == "tar":
        validators = getattr(gr, "VALIDATORS", {}) if gr is not None else {}
        if entry["resource"] in validators:
            return validators[entry["resource"]](ROOT)["status"] == gr.AVAILABLE
        if entry["resource"] == "VEP_CACHE" and gr is not None:
            return gr.check_resource(ROOT, gr.RESOURCE_BY_ID["VEP_CACHE"])["status"] == gr.AVAILABLE
        return (ROOT / entry["post_dir"] / entry["post_marker"]).is_file()
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
# SLURM SETTINGS: CLI > config/slurm.yaml (setup_download) > HPC_CPU.sh > built-in
# =============================================================================

BUILTIN = {"cpus-per-task": "2", "mem": "16G", "time": "24:00:00", "nodes": "1", "ntasks": "1"}

# HPC_CPU.sh is a compute template (8 CPUs / 64G). A download needs the
# setup_download stage's CPUs/memory, not that: only its scheduling keys are used.
TEMPLATE_SKIP = {"cpus-per-task", "mem", "job-name", "output", "error", "array"}


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
    values, source = dict(BUILTIN), {k: "built-in" for k in BUILTIN}
    for key, value in template.items():
        if key not in TEMPLATE_SKIP:
            values[key], source[key] = value, args.template

    extra, array_limit, max_parallel, config_name = [], 1000, 0, None
    if gc is not None:
        try:
            config = gc.load_slurm_config(ROOT, args.slurm_config)
            if config.get("_source", "").startswith("built-in"):
                raise FileNotFoundError
            r = gc.get_stage_resources("setup_download", config)
            config_name = Path(config["_path"]).name
            for key, value in (("cpus-per-task", r["cpus"]), ("mem", r["memory"]), ("time", r["time"]),
                               ("nodes", r["nodes"]), ("ntasks", r["ntasks"]), ("partition", r["partition"]),
                               ("account", r["account"]), ("qos", r["qos"]),
                               ("reservation", r["reservation"]), ("constraint", r["constraint"])):
                if value not in (None, ""):
                    values[key], source[key] = str(value), f"{config_name} (setup_download)"
            extra = list(r.get("extra_sbatch_directives") or [])
            array_limit, max_parallel = int(r["array_limit"]), int(r["max_parallel"] or 0)
        except FileNotFoundError:
            pass
        except gc.SlurmConfigError as exc:
            sys.exit(f"SLURM CONFIGURATION : INVALID\n{exc}")
    elif args.slurm_config:
        sys.exit("ERROR: --slurm-config needs gwas2m_config.py next to this script.")

    for key, value in (("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time),
                       ("partition", args.partition), ("qos", args.qos), ("account", args.account)):
        if value not in (None, ""):
            values[key], source[key] = str(value), "command line"

    max_jobs = args.max_jobs or array_limit
    if max_jobs < 2:
        sys.exit("ERROR: --max-jobs must be at least 2 (one download task + one controller job)")
    parallel = args.parallel if args.parallel is not None else max_parallel
    return dict(shebang=shebang, values=values, source=source, extra=extra, config=config_name,
                chunk=min(max_jobs - 1, array_limit), max_jobs=max_jobs, parallel=parallel)


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
# Tasks per array   : @@CHUNK@@  (never more than @@QUEUED@@ jobs queued; limit @@MAX_JOBS@@)
# Arrays            : @@N_ARRAYS@@ (each next array is submitted after the previous one ends)
#
#   sbatch @@SCRIPT_NAME@@             # SLURM cluster
#   bash   @@SCRIPT_NAME@@ local 4     # no SLURM: 4 files at a time
#   bash   @@SCRIPT_NAME@@ task 12     # just row 12 (debugging)
#   bash   @@SCRIPT_NAME@@ finalize    # verify everything + report
# =============================================================================

set -uo pipefail

ROOT="@@ROOT@@"
SCRIPT="$ROOT/@@SCRIPT_NAME@@"
PLANNER="@@PLANNER@@"
MANIFEST="$ROOT/setup_manifests/download_resources.tsv"
LOGS="@@LOGS@@"
TOTAL=@@TOTAL@@
CHUNK=@@CHUNK@@
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
PYTHON="$(tool python)"; [[ -n "$PYTHON" ]] || PYTHON="$(command -v python3 || command -v python || true)"

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
        tar)   [[ -e "$ROOT/$PDIR/$PMARKER" ]] ;;
        *)     return 0 ;;
    esac
}

# -----------------------------------------------------------------------------
# download: aria2c (resumable via .aria2), else curl/wget into <file>.part
# -----------------------------------------------------------------------------
aria2() {
    "$ARIA2" --continue=true --max-connection-per-server=4 --split=4 --min-split-size=16M \
        --file-allocation=none --auto-file-renaming=false --allow-overwrite=true \
        --max-tries=10 --retry-wait=10 --connect-timeout=60 --timeout=120 \
        --console-log-level=warn --summary-interval=60 \
        --dir "$(dirname "$OUT")" --out "$(basename "$OUT")" "$URL"
}

to_part() {
    if [[ -n "$CURL" ]]; then
        "$CURL" --fail --location --retry 10 --retry-delay 10 --retry-all-errors \
            --connect-timeout 60 --continue-at - --output "$OUT.part" "$URL" || return 1
    elif [[ -n "$WGET" ]]; then
        "$WGET" --continue --tries=10 --waitretry=10 --timeout=120 -O "$OUT.part" "$URL" || return 1
    else
        log "ERROR: none of aria2c / curl / wget is available"; return 1
    fi
    mv "$OUT.part" "$OUT"
}

fetch() {
    local round
    for round in 1 2 3; do
        if [[ -e "$OUT.aria2" && -n "$ARIA2" ]]; then
            log "resuming aria2 download"; aria2 && return 0
        elif [[ -e "$OUT.part" || -z "$ARIA2" ]]; then
            [[ -e "$OUT.part" ]] && log "resuming $OUT.part"
            to_part && return 0
        else
            rm -f "$OUT"            # only reached when OUT is empty / wrong size
            aria2 && return 0
            log "aria2c failed; trying curl/wget"
            to_part && return 0
        fi
        log "download round $round failed; retrying in $((60 * round))s"
        sleep $((60 * round))
    done
    return 1
}

integrity() {
    case "$CHECK" in
        gzip)    if [[ -n "$PIGZ" ]]; then "$PIGZ" -t "$OUT"; else gzip -t "$OUT"; fi ;;
        tar)     tar -tf "$OUT" >/dev/null ;;
        parquet) [[ "$(tail -c 4 "$OUT")" == "PAR1" ]] ;;
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
        tar)
            mkdir -p "$ROOT/$PDIR"
            log "extracting -> $PDIR"
            if [[ "$OUT" == *.gz ]]; then tar -xzf "$OUT" -C "$ROOT/$PDIR"; else tar -xf "$OUT" -C "$ROOT/$PDIR"; fi || return 1
            touch "$ROOT/$PDIR/$PMARKER"
            log "extracted: $PDIR/$PMARKER" ;;
    esac
}

# -----------------------------------------------------------------------------
# one task = one manifest row
# -----------------------------------------------------------------------------
run_task() {
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

    if valid; then
        log "[SKIP] target already valid"
    else
        fetch || { log "FAILED: download"; return 1; }
        valid || { log "FAILED: target missing, empty, still partial or wrong size ($(file_size "$OUT") bytes)"; return 1; }
    fi
    if ! integrity; then
        log "FAILED: integrity check ($CHECK); deleting the corrupt file so the next run re-downloads it"
        rm -f "$OUT"; return 1
    fi
    log "verified: $(file_size "$OUT") bytes ($CHECK)"

    post_done || post_process || { log "FAILED: post-processing ($POST)"; return 1; }
    rm -f "$OUT.lock"
    log "[OK] $RESOURCE: $DEST"
}

# -----------------------------------------------------------------------------
# chaining: never more than CHUNK + 1 jobs queued
# -----------------------------------------------------------------------------
submit_array() {
    local offset="$1" n=$(( TOTAL - $1 ))
    (( n > CHUNK )) && n=$CHUNK
    sbatch --parsable --array="1-${n}${THROTTLE}" --export=ALL,G2M_MODE=download,G2M_OFFSET="$offset" "$SCRIPT"
}

submit_followup() {
    # Task 1 of every array queues ONE small job that starts after the whole
    # array has ended: the next array's controller, or the final check.
    local next=$(( $1 + CHUNK )) mode name
    mkdir "$LOGS/.chain_${SLURM_ARRAY_JOB_ID}" 2>/dev/null || return 0   # once per array
    if (( next < TOTAL )); then mode=next; name=Step00_Download_next; else mode=finalize; name=Step00_Download_finalize; fi
    sbatch --parsable ${SBATCH_COMMON[@]+"${SBATCH_COMMON[@]}"} --dependency="afterany:${SLURM_ARRAY_JOB_ID}" \
        --job-name="$name" --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=4G --time=02:00:00 \
        --output="$LOGS/${name}_%j.log" --error="$LOGS/${name}_%j.err" \
        --export=ALL,G2M_MODE="$mode",G2M_OFFSET="$next" --wrap="bash '$SCRIPT'" \
        && log "queued $mode job (runs after array ${SLURM_ARRAY_JOB_ID})"
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

case "${G2M_MODE:-download}" in
    next)     log "submitting the next array (rows $(( G2M_OFFSET + 1 ))..)"; submit_array "$G2M_OFFSET"; exit $? ;;
    finalize) finalize; exit $? ;;
esac

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Submit with:  sbatch @@SCRIPT_NAME@@      (or without SLURM: bash @@SCRIPT_NAME@@ local 4)"
    exit 2
fi

OFFSET="${G2M_OFFSET:-0}"
[[ "$SLURM_ARRAY_TASK_ID" == "1" ]] && submit_followup "$OFFSET"
run_task $(( OFFSET + SLURM_ARRAY_TASK_ID ))
'''


def render_script(total: int, settings: dict) -> str:
    values = settings["values"]
    chunk = min(total, settings["chunk"])
    throttle = f"%{settings['parallel']}" if settings["parallel"] and settings["parallel"] > 0 else ""
    directives = {**values, "job-name": "Step00_Download_Resources", "array": f"1-{chunk}{throttle}",
                  "output": f"{LOGS}/Step00_Download_%A_%a.log", "error": f"{LOGS}/Step00_Download_%A_%a.err"}
    order = ["job-name", "account", "partition", "qos", "reservation", "constraint", "nodes", "ntasks",
             "cpus-per-task", "mem", "time", "array", "output", "error"]
    lines = [f"#SBATCH --{k}={directives[k]}" for k in order if directives.get(k)]
    lines += [f"#SBATCH --{k}={v}" for k, v in directives.items() if k not in order and v]
    lines += [f"#SBATCH {x.strip() if x.strip().startswith('-') else '--' + x.strip()}" for x in settings["extra"]]
    common = [f"--{k}={values[k]}" for k in ("account", "partition", "qos", "reservation", "constraint")
              if values.get(k)]

    planner_args = []
    if args.opentargets_release != OPEN_TARGETS_RELEASE:
        planner_args.append(f"--opentargets-release {args.opentargets_release}")

    replacements = {
        "@@SHEBANG@@": settings["shebang"],
        "@@DIRECTIVES@@": "\n".join(lines),
        "@@DATE@@": time.strftime("%Y-%m-%d %H:%M"),
        "@@TOTAL@@": str(total),
        "@@CHUNK@@": str(chunk),
        "@@QUEUED@@": str(chunk + 1),
        "@@MAX_JOBS@@": str(settings["max_jobs"]),
        "@@N_ARRAYS@@": str(-(-total // chunk)),
        "@@THROTTLE@@": throttle,
        "@@SBATCH_COMMON@@": " ".join(f'"{c}"' for c in common),
        "@@ROOT@@": str(ROOT),
        "@@LOGS@@": str(LOGS),
        "@@SCRIPT_NAME@@": args.output,
        "@@PLANNER@@": str(Path(__file__).resolve()),
        "@@PLANNER_ARGS@@": " ".join(planner_args),
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
    """gwas2m_setup_downloads.finalize_opentargets(): .complete only when every expected file has its size."""
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
    if gr is not None:
        print("\nShared resources (gwas2m_resources validators):")
        for resource in gr.RESOURCES:
            result = gr.check_resource(ROOT, resource)
            print(f"  {resource['id']:<32}{result['status']:<18}{(result.get('detail') or '')[:70]}")
    if missing:
        print(f"\n{len(missing)} file(s) missing or failed:")
        for entry in missing[:50]:
            print(f"  - {entry['resource']:<30} {entry['dest']}")
        print("\nRe-plan only what is missing and resubmit:")
        print(f"  python {Path(__file__).name} && sbatch {args.output}")
        return 1
    print("\nAll planned files downloaded and verified.")
    print("Next (derived resources: 1000G PGEN panels, sample lists, Pangolin DB):")
    print("  python Step00_Check_Resources.py --resources-only")
    return 0


# =============================================================================
# MAIN (plan)
# =============================================================================

def main() -> int:
    if args.finalize:
        return finalize()

    banner("GWAS2m STEP 00 - PLAN SHARED RESOURCE DOWNLOADS")
    print(f"Project root : {ROOT}")
    print(f"Registry     : {'gwas2m_resources.py' if gr is not None else 'built-in copy (gwas2m_resources.py not found)'}")
    print(f"Groups       : {', '.join(SELECTED_GROUPS) or '(none)'}")
    settings = slurm_settings()

    print("\nCollecting file list ...")
    notes: list[str] = []
    entries = static_items(notes) + eqtl_catalogue_items(notes) + opentargets_items()

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

    pending = entries if args.force else [e for e in entries if not is_complete(e)]

    print(f"\n{'RESOURCE':<32}{'EXPECTED':>9}{'VALID':>8}{'TASKS':>8}")
    for resource in dict.fromkeys(e["resource"] for e in entries):
        n_all = sum(e["resource"] == resource for e in entries)
        n_todo = sum(e["resource"] == resource for e in pending)
        print(f"{resource:<32}{n_all:>9}{n_all - n_todo:>8}{n_todo:>8}")
    for note in notes:
        print(f"NOTE: {note}")

    MANIFESTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    write_tsv(PLAN, ["RESOURCE", "GROUP", "STATUS", "SIZE_BYTES", "SOURCE", "DEST"],
              [[e["resource"], e["group"], "PLANNED" if e in pending else "VALID", e["size"] or "",
                e["url"], e["dest"]] for e in entries])

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
        print("Next: python Step00_Check_Resources.py --resources-only   (derived resources)")
        return 0

    write_tsv(MANIFEST, MANIFEST_COLUMNS, [
        [index, e["resource"], e["group"], e["url"], e["dest"],
         "-" if e["mutable"] or not e["size"] else e["size"], e["check"], e["post"], e["post_dir"], e["post_marker"]]
        for index, e in enumerate(pending, start=1)
    ])

    total_bytes = sum(e["size"] or 0 for e in pending)
    extract_bytes = sum(e["size"] or 0 for e in pending if e["post"] != "-")
    free = shutil.disk_usage(ROOT).free

    script = ROOT / args.output
    with open(script, "w", newline="\n", encoding="utf-8") as handle:
        handle.write(render_script(len(pending), settings))
    script.chmod(0o755)

    chunk = min(len(pending), settings["chunk"])
    banner("SLURM SCRIPT READY")
    print(f"Script         : {script}")
    print(f"Task manifest  : {MANIFEST}  ({len(pending)} files)")
    print(f"Download size  : {human(total_bytes)}  (+ ~{human(extract_bytes)} after extraction)")
    print(f"Free disk      : {human(free)}" + ("   [WARN] not enough for everything planned"
                                                if total_bytes + extract_bytes > free else ""))
    print(f"Arrays         : {-(-len(pending) // chunk)} x up to {chunk} tasks; never more than {chunk + 1} "
          f"jobs queued (limit {settings['max_jobs']})")
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
