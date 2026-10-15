#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2m STEP 00 - ONE-COMMAND SOFTWARE INSTALLATION
===============================================================================

Installs, repairs and verifies ALL software used by:

GWAS
 -> QC                    (gwaslab)
 -> ancestry-specific LD  (plink2)
 -> clumping              (plink2)
 -> SuSiE                 (R susieR)
 -> figures               (R ggplot2, patchwork, ggrastr, ggrepel, svglite,
                           ComplexHeatmap, karyoploteR; matplotlib)
 -> VEP                   (ensembl-vep 116)
 -> GTEx / eQTL Catalogue (tabix, pandas, aria2c)
 -> SpliceAI              (isolated env: TensorFlow 2.15 / Keras 2.15)
 -> Pangolin              (isolated env: PyTorch, pinned git commit)
 -> colocalisation        (R coloc, data.table; hyprcoloc multi-trait)
 -> gene prioritisation   (MAGMA v1.10, PoPS)
 -> heritability enrich.  (LDSC / S-LDSC)
 -> pathways / intervals  (gseapy, pyranges)
 -> liftover              (pyliftover, CrossMap)
 -> Open Targets / FUMA   (duckdb, pyarrow, polars, rsync)
 -> resource downloads    (aria2c, curl, wget, pigz, samtools, tabix)

Everything is created beneath the CURRENT WORKING DIRECTORY:

    envs/pipeline        main pipeline stack (from environment.yml)
    envs/spliceai        SpliceAI 1.3.1 (incompatible with main env)
    envs/pangolin        Pangolin (incompatible with main env)
    envs/ldsc            LDSC 1.0.1 (Python 2.7; incompatible with main env)
    envs/pops            PoPS (pins numpy 1.19 / scikit-learn 0.24)
    external_tools/Pangolin
    external_tools/pops  PoPS source (pinned commit)
    external_tools/magma/v1.10/   MAGMA binary (official CTG build, SHA-256 checked)
    setup_manifests/     Pangolin/PoPS commit locks, installed-software manifest,
                         per-environment package lists (reproducibility)
    logs/installation/   full installation log
    .gwas2m_cache/       conda package cache, pip cache and temp files
                         (nothing is written to $HOME or /tmp; delete with
                         --clean-cache once installation succeeds)

Works on any Linux machine or HPC (and macOS) with mamba, conda or micromamba
on PATH - e.g. `module load miniforge` on a cluster, or Miniforge on a laptop.
No path is hard-coded: the project root is wherever you run the command from.
Conda environments cannot be moved afterwards; if the project folder is moved
or copied, rebuild with --force-recreate (the script detects this).

Usage (from the project directory):

    python Step00_Installation.py                       # install + verify all
    python Step00_Installation.py --only pipeline       # one environment
    python Step00_Installation.py --only spliceai pangolin
    python Step00_Installation.py --only ldsc pops magma
    python Step00_Installation.py --verify-only         # read-only check
    python Step00_Installation.py --dry-run             # print commands only
    python Step00_Installation.py --force-recreate --only spliceai

Re-running is safe: existing environments are updated/repaired, never rebuilt
unless --force-recreate is given. Only SOFTWARE is installed. Reference DATA
(genome, GTEx, 1000G, VEP cache, MAGMA gene locations, LDSC baseline / weights,
PoPS features, liftover chain files, ...) is NOT downloaded here; use
Step00_Download_Resources.py for shared resources.

Launch with any Python >= 3.7 (e.g. the one that comes with conda/Miniforge).

===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import shlex
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path


# =============================================================================
# ROOT
# =============================================================================

ROOT = Path.cwd().resolve()

ENVS = ROOT / "envs"
EXTERNAL = ROOT / "external_tools"
MANIFESTS = ROOT / "setup_manifests"
LOGS = ROOT / "logs" / "installation"

CORE_ENV = ENVS / "pipeline"
SPLICEAI_ENV = ENVS / "spliceai"
PANGOLIN_ENV = ENVS / "pangolin"

PANGOLIN_URL = "https://github.com/tkzeng/Pangolin.git"
PANGOLIN_REPO = EXTERNAL / "Pangolin"
PANGOLIN_LOCK = MANIFESTS / "Pangolin.commit.txt"

ENVIRONMENT_YML = ROOT / "environment.yml"

# Every cache / temp file also stays inside the project (nothing in $HOME or
# /tmp), so the install behaves the same on a laptop, a workstation or any HPC.
CACHE = ROOT / ".gwas2m_cache"
CONDA_PKGS = CACHE / "conda_pkgs"
PIP_CACHE = CACHE / "pip"
TMP = CACHE / "tmp"

# Conda environments hard-code their absolute path; this records where they
# were built so a moved/copied project is detected instead of failing oddly.
ROOT_MARKER = ENVS / ".installed_root"


# =============================================================================
# PLATFORM
# =============================================================================

SYSTEM = platform.system()
IS_LINUX = SYSTEM == "Linux"
IS_MAC = SYSTEM == "Darwin"
# Many bioconda tools have no osx-arm64 build; Apple Silicon uses osx-64 (Rosetta).
CONDA_SUBDIR = "osx-64" if IS_MAC and platform.machine() == "arm64" else None


# =============================================================================
# VERSIONS
# =============================================================================

VEP_RELEASE = "116"

GWASLAB_SPEC = "gwaslab==4.2.3"
POLARS_SPEC = "polars[rtcompat]"

SPLICEAI_CONDA = [
    "python=3.10",
    "pip",
    "setuptools<81",
    "numpy<2",
]

SPLICEAI_PIP = [
    "setuptools<81",
    "numpy<2",
    # tensorflow-cpu has no macOS wheels; plain tensorflow is CPU-only there.
    "tensorflow==2.15.1" if IS_MAC else "tensorflow-cpu==2.15.1",
    "keras==2.15.0",
    "spliceai==1.3.1",
]

PANGOLIN_CONDA = [
    "python=3.10",
    "pip",
    "numpy<2",
    "pytorch",
    "gffutils",
    "biopython",
    "pandas",
    "pyfastx",
]

PANGOLIN_PIP = [
    "PyVCF3",
]

# --- hyprcoloc (multi-trait colocalisation; R, GitHub only, commit pinned) ---
HYPRCOLOC_REPO = "jrs95/hyprcoloc"
HYPRCOLOC_COMMIT = "0348bbdd977be731d82e4efda115a9bed63cd44f"

# --- LDSC / S-LDSC (Python 2.7 only -> isolated env; bioconda package) ---
LDSC_ENV = ENVS / "ldsc"
LDSC_CONDA = ["ldsc=1.0.1"]

# --- PoPS (pins numpy 1.19 / scikit-learn 0.24 and patches sklearn internals
#     -> isolated env; GitHub only, commit pinned) ---
POPS_ENV = ENVS / "pops"
POPS_URL = "https://github.com/FinucaneLab/pops.git"
POPS_COMMIT = "76eb86cba10254490003c8f4dc7ff5ce492d3667"
POPS_REPO = EXTERNAL / "pops"
POPS_CONDA = ["python=3.8", "pip"]

# --- MAGMA v1.10 (CTG lab; official binaries, SHA-256 verified) ---
#     NOTE: the conda package called "magma" is an unrelated GPU linear-algebra
#     library, so MAGMA is always taken from the official CTG download.
MAGMA_VERSION = "1.10"
MAGMA_DIR = EXTERNAL / "magma" / f"v{MAGMA_VERSION}"
MAGMA_BIN = MAGMA_DIR / "magma"
MAGMA_DOWNLOADS = {
    # Statically linked build: no dependency on the host's glibc / Intel libs.
    "Linux-x86_64": ("magma_v1.10_static.zip",
                     "https://vu.data.surfsara.nl/index.php/s/lxDgt2dNdNr6DYt/download",
                     "d8b20778b773f47b4fb0f1020baefebc92b11ee65946e708a618d494d6819e39"),
    "Darwin": ("magma_v1.10_mac.zip",
               "https://vu.data.surfsara.nl/index.php/s/1M1d9vHtVidEwvU/download",
               "5ebfc624a8d693be22bf1f0830b4fdc5810119384c60f98044aa168494eeba71"),
}

# Rough on-disk size of all environments + package cache.
EXPECTED_DISK_GB = 30


# =============================================================================
# MAIN ENVIRONMENT SPECIFICATION (written to environment.yml)
# =============================================================================

ENVIRONMENT_YML_TEXT = f"""name: pipeline

channels:
  - conda-forge
  - bioconda

dependencies:

  # Python
  - python=3.11

  # Core Python scientific stack
  - pandas
  - numpy
  - scipy
  - pyarrow
  - duckdb
  - pysam
  - pybigwig
  - cyvcf2
  - scikit-learn
  - statsmodels
  - matplotlib
  - requests
  - urllib3
  - httpx
  - tenacity
  - tqdm
  - orjson
  - pyyaml
  - networkx
  - gprofiler-official

  # Genomics
  - plink2
  - bcftools
  - htslib
  - samtools
  - bedtools

  # Download / compression
  - aria2
  - wget
  - curl
  - pigz
  - zstd
  - parallel
  - rsync

  # VEP
  - ensembl-vep={VEP_RELEASE}

  # Gene-set enrichment / genomic intervals / liftover / tables + figures
  - gseapy
  - pyranges
  - pyliftover
  - crossmap
  - seaborn
  - openpyxl
  - xlsxwriter

  # R
  - r-base
  - r-susier
  - r-coloc
  - r-data.table

  # R dependencies of hyprcoloc (hyprcoloc itself is installed from GitHub)
  - r-remotes
  - r-iterpc
  - r-rcpp
  - r-rcppeigen
  - r-rmpfr
  - r-matrix
  - r-knitr
  - r-rmarkdown
  - r-pheatmap
  - r-rcolorbrewer

  # Publication figures (Step07_Plots.py)
  - r-ggplot2
  - r-patchwork
  - r-ggrastr
  - r-ggrepel
  - r-scales
  - r-svglite
  - bioconductor-complexheatmap
  - bioconductor-karyoploter

  # Toolchain for compiling R packages from source inside the environment
  - compilers
  - make

  # General
  - git
  - pip

  # pip-only packages
  - pip:
      - {GWASLAB_SPEC}
      - {POLARS_SPEC}
"""


# =============================================================================
# SOFTWARE CHECKED AFTER INSTALLATION
# =============================================================================

# name, env, version command, required, step, purpose, install spec
TOOLS = [
    ("python", "pipeline", ["--version"], True, 1, "pipeline interpreter", "python=3.11"),
    ("Rscript", "pipeline", ["--version"], True, 6, "SuSiE-RSS / coloc", "r-base"),
    ("plink2", "pipeline", ["--version"], True, 5, "clumping / LD", "plink2"),
    ("bcftools", "pipeline", ["--version"], True, 9, "VCF normalisation", "bcftools"),
    ("samtools", "pipeline", ["--version"], True, 9, "FASTA index", "samtools"),
    ("tabix", "pipeline", ["--version"], True, 11, "eQTL Catalogue region queries", "htslib"),
    ("bgzip", "pipeline", ["--version"], False, 11, "compression", "htslib"),
    ("bedtools", "pipeline", ["--version"], False, None, "declared; not used by Steps 01-11", "bedtools"),
    ("vep", "pipeline", ["--help"], True, 7, "Ensembl VEP annotation", f"ensembl-vep={VEP_RELEASE}"),
    ("aria2c", "pipeline", ["--version"], True, None, "shared-resource downloads (setup)", "aria2"),
    ("pigz", "pipeline", ["--version"], True, None, "FASTA decompression (setup)", "pigz"),
    ("wget", "pipeline", ["--version"], False, 2, "GWAS download (first choice)", "wget"),
    ("curl", "pipeline", ["--version"], False, 2, "GWAS download (fallback)", "curl"),
    ("rsync", "pipeline", ["--version"], False, None, "Open Targets download (setup, optional)", "rsync"),
    ("git", "pipeline", ["--version"], True, None, "Pangolin checkout (setup)", "git"),
    ("spliceai", "spliceai", None, True, 9, "SpliceAI executable", "spliceai==1.3.1"),
    ("pangolin", "pangolin", None, True, 10, "Pangolin executable", "git+tkzeng/Pangolin"),
]

# import name, distribution name, env, required, step, install method, spec, used by
PYTHON_PACKAGES = [
    ("pandas", "pandas", "pipeline", True, 1, "conda", "pandas", "all steps"),
    ("numpy", "numpy", "pipeline", True, 1, "conda", "numpy", "all steps"),
    ("requests", "requests", "pipeline", True, 1, "conda", "requests", "Step01, Step06"),
    ("urllib3", "urllib3", "pipeline", True, 1, "conda", "urllib3", "Step01, Step06"),
    ("gwaslab", "gwaslab", "pipeline", True, 3, "pip", GWASLAB_SPEC, "Step03 QC"),
    ("yaml", "PyYAML", "pipeline", False, 6, "conda", "pyyaml", "Step06 metadata"),
    ("scipy", "scipy", "pipeline", False, None, "conda", "scipy", "Pipeline_Report, Step13, Step14"),
    ("polars", "polars", "pipeline", False, None, "pip", POLARS_SPEC, "Diagnose_Step11, Pipeline_Report"),
    ("orjson", "orjson", "pipeline", False, None, "conda", "orjson", "Diagnose_Step11"),
    ("pyarrow", "pyarrow", "pipeline", False, None, "conda", "pyarrow", "parquet (Open Targets)"),
    ("duckdb", "duckdb", "pipeline", False, None, "conda", "duckdb", "Step13"),
    ("pysam", "pysam", "pipeline", True, 8, "conda", "pysam", "Step08 CADD / AlphaMissense / ClinVar lookups"),
    ("pyBigWig", "pybigwig", "pipeline", True, 8, "conda", "pybigwig", "Step08 phyloP conservation"),
    ("cyvcf2", "cyvcf2", "pipeline", False, None, "conda", "cyvcf2", "declared"),
    ("sklearn", "scikit-learn", "pipeline", False, None, "conda", "scikit-learn", "declared"),
    ("statsmodels", "statsmodels", "pipeline", False, None, "conda", "statsmodels", "Step14"),
    ("httpx", "httpx", "pipeline", False, None, "conda", "httpx", "declared"),
    ("tenacity", "tenacity", "pipeline", False, None, "conda", "tenacity", "declared"),
    ("tqdm", "tqdm", "pipeline", False, None, "conda", "tqdm", "Diagnose_Step11"),
    ("matplotlib", "matplotlib", "pipeline", False, None, "conda", "matplotlib", "declared"),
    ("networkx", "networkx", "pipeline", False, None, "conda", "networkx", "Step15"),
    ("gprofiler", "gprofiler-official", "pipeline", False, None, "conda", "gprofiler-official", "Step14"),
    ("spliceai", "spliceai", "spliceai", True, 9, "pip", "spliceai==1.3.1", "Step09"),
    ("tensorflow", "tensorflow-cpu", "spliceai", True, 9, "pip", "tensorflow-cpu==2.15.1", "Step09"),
    ("keras", "keras", "spliceai", True, 9, "pip", "keras==2.15.0", "Step09"),
    ("pangolin", "pangolin", "pangolin", True, 10, "git+pip -e", "tkzeng/Pangolin (commit pinned)", "Step10"),
    ("gffutils", "gffutils", "pangolin", True, 10, "conda", "gffutils", "Step10 / Pangolin DB"),
    ("pyfastx", "pyfastx", "pangolin", True, 10, "conda", "pyfastx", "Step10"),
    ("Bio", "biopython", "pangolin", True, 10, "conda", "biopython", "Step10"),
]

# name, required, step, spec, used by
R_PACKAGES = [
    ("susieR", True, 6, "r-susier", "Step06 SuSiE-RSS"),
    ("coloc", True, 11, "r-coloc", "Step11 coloc.bf_bf"),
    ("data.table", True, 11, "r-data.table", "coloc dependency"),
    ("ggplot2", True, 7, "r-ggplot2", "Step07 figures"),
    ("patchwork", True, 7, "r-patchwork", "Step07 multi-panel locus figures"),
    ("ggrastr", True, 7, "r-ggrastr", "Step07 rasterised dense point layers"),
    ("ggrepel", True, 7, "r-ggrepel", "Step07 variant / gene labels"),
    ("scales", True, 7, "r-scales", "Step07 axes"),
    ("svglite", True, 7, "r-svglite", "Step07 SVG output"),
    ("ComplexHeatmap", True, 7, "bioconductor-complexheatmap", "Step07 locus x GWAS matrix"),
    ("karyoploteR", True, 7, "bioconductor-karyoploter", "Step07 chromosome ideograms"),
]

# -----------------------------------------------------------------------------
# GWAS -> mechanism add-ons. Same tuple layouts as above.
# -----------------------------------------------------------------------------

EXTRA_TOOLS = [
    ("CrossMap", "pipeline", None, True, None, "coordinate liftover (VCF/BED/GFF)", "crossmap"),
    ("ldsc.py", "ldsc", None, True, None, "LD score regression / S-LDSC", "ldsc=1.0.1"),
    ("munge_sumstats.py", "ldsc", None, True, None, "LDSC summary-statistic munging", "ldsc=1.0.1"),
    ("pops", "pops", None, True, None, "PoPS gene prioritisation (wrapper)", f"FinucaneLab/pops@{POPS_COMMIT[:7]}"),
    ("magma", "magma", ["--version"], True, None, "MAGMA gene / gene-set analysis", f"MAGMA v{MAGMA_VERSION}"),
]

EXTRA_PYTHON_PACKAGES = [
    ("gseapy", "gseapy", "pipeline", True, None, "conda", "gseapy", "GSEA / Enrichr pathway enrichment"),
    ("pyranges", "pyranges", "pipeline", True, None, "conda", "pyranges", "variant x ENCODE interval overlap"),
    ("pyliftover", "pyliftover", "pipeline", True, None, "conda", "pyliftover", "GRCh37 <-> GRCh38 positions"),
    ("cmmodule", "CrossMap", "pipeline", True, None, "conda", "crossmap", "CrossMap liftover library"),
    ("seaborn", "seaborn", "pipeline", True, None, "conda", "seaborn", "figures"),
    ("openpyxl", "openpyxl", "pipeline", True, None, "conda", "openpyxl", "xlsx supplementary tables"),
    ("xlsxwriter", "XlsxWriter", "pipeline", True, None, "conda", "xlsxwriter", "xlsx supplementary tables"),
    ("ldscore", "ldsc", "ldsc", True, None, "conda", "ldsc=1.0.1", "LDSC core library"),
    ("bitarray", "bitarray", "ldsc", True, None, "conda", "bitarray", "LDSC dependency"),
    ("numpy", "numpy", "pops", True, None, "pip", "numpy==1.19.5", "PoPS"),
    ("scipy", "scipy", "pops", True, None, "pip", "scipy==1.5.2", "PoPS"),
    ("pandas", "pandas", "pops", True, None, "pip", "pandas==1.0.5", "PoPS"),
    ("sklearn", "scikit-learn", "pops", True, None, "pip", "scikit_learn==0.24.2", "PoPS"),
]

EXTRA_R_PACKAGES = [
    ("hyprcoloc", True, None, f"github {HYPRCOLOC_REPO}@{HYPRCOLOC_COMMIT[:7]}", "multi-trait colocalisation"),
]

# Executables that have an alternative name in some package versions.
TOOL_ALIASES = {"CrossMap": ["CrossMap", "CrossMap.py"]}

# Executables that do not live in envs/<env>/bin.
TOOL_PATHS = {"magma": MAGMA_BIN}

# Extra runtime checks that are not plain imports (env, label, python code).
RUNTIME_CHECKS = [
    ("spliceai", "SpliceAI Annotator + pkg_resources",
     "import pkg_resources; from spliceai.utils import Annotator"),
    ("pangolin", "Pangolin PyTorch backend", "import torch"),
]

# Runtime checks that run a command (component, label, command).
COMMAND_CHECKS = [
    ("pangolin", "pangolin --help", [ENVS / "pangolin" / "bin" / "pangolin", "--help"]),
    ("ldsc", "ldsc.py -h", [LDSC_ENV / "bin" / "python", LDSC_ENV / "bin" / "ldsc.py", "-h"]),
    ("pops", "pops --help", [POPS_ENV / "bin" / "pops", "--help"]),
    ("magma", "magma --version", [MAGMA_BIN, "--version"]),
]

# Host tools used by the SLURM job scripts; never installed here.
HOST_TOOLS = [
    ("sbatch", "submits the Step00_Download_Resources.sh job arrays"),
    ("squeue", "job monitoring"),
]

TOOLS = TOOLS + EXTRA_TOOLS
PYTHON_PACKAGES = PYTHON_PACKAGES + EXTRA_PYTHON_PACKAGES
R_PACKAGES = R_PACKAGES + EXTRA_R_PACKAGES

# Conda environments under envs/, plus MAGMA (a binary under external_tools/).
CONDA_COMPONENTS = ["pipeline", "spliceai", "pangolin", "ldsc", "pops"]
COMPONENTS = CONDA_COMPONENTS + ["magma"]


# =============================================================================
# ARGUMENTS
# =============================================================================

parser = argparse.ArgumentParser(
    description="Install, repair and verify all GWAS2m software environments."
)
parser.add_argument("--only", nargs="+", choices=COMPONENTS, default=COMPONENTS,
                    help="Install only these environments (default: all three).")
parser.add_argument("--verify-only", action="store_true",
                    help="READ-ONLY: verify what is installed; install nothing.")
parser.add_argument("--dry-run", action="store_true",
                    help="Print every command without running it.")
parser.add_argument("--skip-verify", action="store_true",
                    help="Skip the verification stage after installation.")
parser.add_argument("--force-recreate", action="store_true",
                    help="Delete and rebuild the selected environments from scratch.")
parser.add_argument("--retries", type=int, default=2,
                    help="Retries for each conda/pip/git command on failure (default 2).")
parser.add_argument("--pkgs-dir", default=None,
                    help=f"Conda package cache directory (default: {CONDA_PKGS.relative_to(ROOT)} "
                         "inside the project).")
parser.add_argument("--clean-cache", action="store_true",
                    help="Delete the project-local download cache after a successful install.")
parser.add_argument("--min-free-gb", type=float, default=EXPECTED_DISK_GB,
                    help=f"Warn when free disk under the project is below this (default {EXPECTED_DISK_GB}).")
args = parser.parse_args()

if args.verify_only and (args.force_recreate or args.dry_run):
    parser.error("--verify-only cannot be combined with --force-recreate or --dry-run")

SELECTED = [c for c in COMPONENTS if c in args.only]


# =============================================================================
# LOGGING
# =============================================================================

LOGS.mkdir(parents=True, exist_ok=True)
LOG_FILE = LOGS / f"Step00_Installation_{time.strftime('%Y%m%d_%H%M%S')}.log"


def log(message: str = "", end: str = "\n") -> None:
    print(message, end=end, flush=True)
    with open(LOG_FILE, "a", encoding="utf-8") as handle:
        handle.write(message + end)


def banner(text: str) -> None:
    log()
    log("=" * 80)
    log(text)
    log("=" * 80)


def section(text: str) -> None:
    log()
    log("-" * 80)
    log(text)
    log("-" * 80)


# =============================================================================
# COMMAND RUNNER (live output + heartbeat + log + retries)
# =============================================================================

PKGS_DIR = Path(args.pkgs_dir).resolve() if args.pkgs_dir else CONDA_PKGS


def child_env() -> dict:
    """Environment for every child process.

    Keeps all downloads, caches and temp files under ROOT, and stops anything
    in the user's home (~/.local Python, ~/R library, PYTHONPATH, pip --user)
    from leaking into the installed or verified environments.
    """
    env = dict(os.environ)

    # Python isolation
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    env["PYTHONNOUSERSITE"] = "1"
    env["PIP_USER"] = "false"
    env["PIP_NO_INPUT"] = "1"
    env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"

    # R isolation (non-existent user library is ignored by R)
    env.pop("R_LIBS", None)
    env["R_LIBS_USER"] = str(CACHE / "no_user_R_library")

    # Project-local caches / temp
    env["CONDA_PKGS_DIRS"] = str(PKGS_DIR)
    env["PIP_CACHE_DIR"] = str(PIP_CACHE)
    env["XDG_CACHE_HOME"] = str(CACHE)
    if TMP.is_dir():  # created before installing; verify-only leaves the system default
        for name in ("TMPDIR", "TMP", "TEMP"):
            env[name] = str(TMP)

    env["CONDA_CHANNEL_PRIORITY"] = "strict"
    # Answer "yes" without -y: older conda rejects "-y" on "conda env create/update".
    env["CONDA_ALWAYS_YES"] = "true"
    if CONDA_SUBDIR:
        env["CONDA_SUBDIR"] = CONDA_SUBDIR
    return env


# Explicit channels for "create" commands: never mix in Anaconda "defaults" from
# the user's .condarc (gffutils / pyfastx / ldsc only exist on bioconda).
CONDA_FORGE = ["--override-channels", "-c", "conda-forge"]
FORGE_BIOCONDA = ["--override-channels", "-c", "conda-forge", "-c", "bioconda"]


def _stream(command: list[str], cwd: Path, label: str) -> int:
    import selectors

    started = time.time()
    process = subprocess.Popen(
        command,
        cwd=cwd,
        env=child_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ)

    with open(LOG_FILE, "a", encoding="utf-8") as handle:
        while True:
            events = selector.select(timeout=30.0)
            if events:
                line = process.stdout.readline()
                if line:
                    print(line, end="", flush=True)
                    handle.write(line)
                    continue
            if process.poll() is not None:
                for line in process.stdout:
                    print(line, end="", flush=True)
                    handle.write(line)
                break
            if not events:
                elapsed = int(time.time() - started)
                print(f"[RUNNING] {label} | elapsed {elapsed // 60}m{elapsed % 60:02d}s", flush=True)

    return int(process.wait())


def run(command: list, *, label: str | None = None, cwd: Path | None = None,
        check: bool = True, retries: int | None = None) -> int:
    command = [str(x) for x in command]
    label = label or Path(command[0]).name
    attempts = 1 + (args.retries if retries is None else retries)

    log()
    log(f">>> {label}")
    log("$ " + " ".join(shlex.quote(x) for x in command))

    if args.dry_run:
        return 0

    rc = 1
    for attempt in range(1, attempts + 1):
        started = time.time()
        rc = _stream(command, cwd or ROOT, label)
        log(f"[DONE] exit={rc} elapsed={(time.time() - started) / 60:.1f} min"
            + (f" (attempt {attempt}/{attempts})" if attempts > 1 else ""))
        if rc == 0:
            break
        if attempt < attempts:
            wait = 20 * attempt
            log(f"[RETRY] {label} failed; retrying in {wait}s ...")
            time.sleep(wait)

    if check and rc != 0:
        raise subprocess.CalledProcessError(rc, command)
    return rc


def capture(command: list, timeout: int = 300) -> tuple[int, str]:
    """Quiet probe; never raises."""
    try:
        result = subprocess.run([str(x) for x in command], capture_output=True, text=True,
                                timeout=timeout, env=child_env())
        return result.returncode, (result.stdout or "") + (result.stderr or "")
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 127, f"{type(exc).__name__}: {exc}"


# =============================================================================
# PACKAGE MANAGER
# =============================================================================

def find_package_managers() -> list[str]:
    """mamba first (fast solver), conda as fallback; micromamba if nothing else."""
    found = []
    for name in ["mamba", "conda", "micromamba"]:
        path = shutil.which(name)
        if path:
            found.append(str(Path(path).resolve()))
    return found


def is_micromamba(pm: str) -> bool:
    return Path(pm).name.startswith("micromamba")


def env_exists(env: Path) -> bool:
    return (env / "conda-meta").is_dir()


def env_bin(env: Path, name: str) -> Path:
    return env / "bin" / name


def with_fallback(pms: list[str], build_command, label: str) -> None:
    """Run a package-manager command with each available manager until one succeeds."""
    last: Exception | None = None
    for pm in pms:
        try:
            run(build_command(pm), label=f"{label} [{Path(pm).name}]")
            return
        except subprocess.CalledProcessError as exc:
            last = exc
            if pm != pms[-1]:
                log(f"[FALLBACK] {Path(pm).name} failed; trying {Path(pms[pms.index(pm) + 1]).name}")
    if last:
        raise last


def remove_env(env: Path) -> None:
    if not env.exists():
        return
    log(f"[RECREATE] Removing {env}")
    if not args.dry_run:
        shutil.rmtree(env)


# =============================================================================
# PREFLIGHT
# =============================================================================

def preflight(pms: list[str]) -> bool:
    banner("PREFLIGHT CHECKS")
    ok = True

    log(f"Project root     : {ROOT}")
    log(f"Host             : {socket.gethostname()}")
    log(f"Platform         : {platform.system()} {platform.release()} ({platform.machine()})")
    log(f"Launcher python  : {sys.executable} ({platform.python_version()})")
    log(f"Components       : {', '.join(SELECTED)}")
    log(f"Log file         : {LOG_FILE}")

    if IS_MAC:
        log("[WARN] macOS: supported, but plink2/VEP/SpliceAI are primarily tested on Linux.")
        if CONDA_SUBDIR:
            log(f"[NOTE] Apple Silicon: environments use {CONDA_SUBDIR} packages (Rosetta 2 required).")
    elif not IS_LINUX:
        log(f"[ERROR] {SYSTEM} is not supported (bioconda tools are Linux/macOS only).")
        log("        On Windows, run this inside WSL2 (Ubuntu) instead.")
        if not args.dry_run and not args.verify_only:
            ok = False

    if not os.access(ROOT, os.W_OK):
        log(f"[ERROR] Project root is not writable: {ROOT}")
        ok = False

    if " " in str(ROOT):
        log("[WARN] The project path contains spaces; conda/Perl (VEP) scripts can break. "
            "Prefer a path without spaces.")

    if ROOT_MARKER.is_file():
        built_at = ROOT_MARKER.read_text(encoding="utf-8").strip()
        if built_at and built_at != str(ROOT):
            log(f"[ERROR] envs/ was built at {built_at} and has since been moved/copied to {ROOT}.")
            log("        Conda environments cannot be relocated. Rebuild them here with:")
            log("          python Step00_Installation.py --force-recreate")
            if not args.force_recreate:
                ok = False

    if not args.verify_only:
        if pms:
            log(f"Package managers : {', '.join(pms)}")
        else:
            log("[ERROR] mamba / conda / micromamba not found on PATH.")
            log("        HPC:      module load miniforge   (or anaconda3 / miniconda3)")
            log("        Laptop:   install Miniforge from https://github.com/conda-forge/miniforge")
            ok = False

    if os.environ.get("CONDA_PREFIX") and Path(os.environ["CONDA_PREFIX"]).resolve() in \
            [(ENVS / c).resolve() for c in SELECTED] and args.force_recreate:
        log(f"[ERROR] Cannot recreate the currently active environment {os.environ['CONDA_PREFIX']}."
            " Run `conda deactivate` first.")
        ok = False

    if os.environ.get("PYTHONPATH"):
        log("[NOTE] PYTHONPATH is set in your shell; it is removed for every install/verify command.")

    if shutil.which("sbatch") and not os.environ.get("SLURM_JOB_ID") and not args.verify_only:
        log("[NOTE] SLURM detected but this is not a batch job (probably a login node).")
        log("       A full install takes 30-90 min; consider submitting it with sbatch.")

    try:
        free_gb = shutil.disk_usage(ROOT).free / 1024 ** 3
        log(f"Free disk        : {free_gb:.1f} GB under {ROOT}")
        if free_gb < args.min_free_gb and not args.verify_only:
            log(f"[WARN] Less than {args.min_free_gb:.0f} GB free; environments need ~{EXPECTED_DISK_GB} GB.")
    except OSError:
        pass

    log(f"Conda pkg cache  : {PKGS_DIR}")
    log(f"pip cache / tmp  : {PIP_CACHE} , {TMP}")

    if not args.verify_only and not args.dry_run:
        for name, url in [
            ("conda-forge", "https://conda.anaconda.org/conda-forge/noarch/repodata.json"),
            ("bioconda", "https://conda.anaconda.org/bioconda/noarch/repodata.json"),
            ("PyPI", "https://pypi.org/simple/pip/"),
            ("GitHub", "https://github.com/tkzeng/Pangolin"),
            ("CRAN", "https://cloud.r-project.org/"),
            ("MAGMA (CTG)", "https://vu.data.surfsara.nl/"),
        ]:
            try:
                request = urllib.request.Request(url, method="HEAD")
                with urllib.request.urlopen(request, timeout=15):
                    pass
                log(f"[OK]      network: {name}")
            except Exception as exc:
                log(f"[WARN]    network: {name} unreachable ({type(exc).__name__}); installs may fail")

    log()
    log("Host tools (SLURM; not installed here):")
    for name, purpose in HOST_TOOLS:
        path = shutil.which(name)
        log(f"  {'[OK]     ' if path else '[ABSENT] '} {name:<8} {path or ''}  - {purpose}")

    return ok


# =============================================================================
# INSTALLERS
# =============================================================================

def install_pipeline(pms: list[str]) -> None:
    banner("MAIN PIPELINE ENVIRONMENT  ->  envs/pipeline")

    if args.force_recreate:
        remove_env(CORE_ENV)

    log(f"Specification : {ENVIRONMENT_YML}")
    if not args.dry_run:
        ENVIRONMENT_YML.write_text(ENVIRONMENT_YML_TEXT, encoding="utf-8")

    update_error: Exception | None = None
    if env_exists(CORE_ENV):
        log("[CACHE] Environment exists; updating / repairing it.")

        def command(pm):
            if is_micromamba(pm):
                return [pm, "install", "-y", "-p", CORE_ENV, "-f", ENVIRONMENT_YML]
            return [pm, "env", "update", "-p", CORE_ENV, "-f", ENVIRONMENT_YML, "--prune"]

        try:
            with_fallback(pms, command, "Update envs/pipeline")
        except subprocess.CalledProcessError as exc:
            # The environment still exists: carry on with the steps below, report at the end.
            update_error = exc
            log("[WARN] update failed; continuing with pip / hyprcoloc / VEP on the existing environment")
    else:
        def command(pm):
            if is_micromamba(pm):
                return [pm, "create", "-y", "-p", CORE_ENV, "-f", ENVIRONMENT_YML]
            return [pm, "env", "create", "-p", CORE_ENV, "-f", ENVIRONMENT_YML]

        with_fallback(pms, command, "Create envs/pipeline")

    section("GWASLab + CPU-compatible Polars runtime (pip)")
    run([env_bin(CORE_ENV, "python"), "-m", "pip", "install", "--upgrade",
         GWASLAB_SPEC, POLARS_SPEC], label="pip: gwaslab + polars")

    pin_vep_perl()
    install_hyprcoloc()

    if update_error is not None:
        raise update_error


VEP_SCRIPTS = ("vep", "filter_vep", "variant_recoder", "haplo", "vep_install", "vep_convert_cache")
LAUNCHER_TAG = "# GWAS2m launcher"


def pin_vep_perl() -> None:
    """Make VEP use envs/pipeline's Perl even when the environment is not activated.

    VEP's scripts start with '#!/usr/bin/env perl', so running envs/pipeline/bin/vep
    directly (as Step07 does) picks up the first perl on PATH - often the system perl,
    which lacks DBI ("Can't locate DBI.pm"). Each bin/ link is replaced by a tiny
    launcher that runs the real script with this environment's perl.
    """
    section("Ensembl VEP: always use this environment's Perl")
    perl = env_bin(CORE_ENV, "perl")
    if args.dry_run:
        log(f"Would point {', '.join(VEP_SCRIPTS)} at {perl}")
        return
    if not perl.exists():
        raise RuntimeError(f"{perl} missing; the ensembl-vep package did not install correctly")

    for name in VEP_SCRIPTS:
        entry = env_bin(CORE_ENV, name)
        if entry.is_symlink():
            real = entry.resolve()
        else:
            hits = sorted((CORE_ENV / "share").glob(f"ensembl-vep*/{name}"))
            real = hits[-1] if hits else None
        if real is None or not real.is_file():
            continue
        with open(real, "rb") as handle:
            if b"perl" not in handle.readline():
                continue
        if entry.is_file() and not entry.is_symlink() and LAUNCHER_TAG not in entry.read_text(errors="replace"):
            continue  # a real script living in bin/ (not the bioconda layout): leave it alone
        if entry.is_symlink() or entry.exists():
            entry.unlink()
        entry.write_text(
            "#!/bin/sh\n"
            f"{LAUNCHER_TAG}: run Ensembl VEP with envs/pipeline's perl, whatever is on PATH\n"
            "unset PERL5LIB PERL_LOCAL_LIB_ROOT\n"
            f'exec "{perl}" "{real}" "$@"\n',
            encoding="utf-8",
        )
        entry.chmod(0o755)
        log(f"[OK] {name} -> {perl} {real.relative_to(CORE_ENV)}")


def install_hyprcoloc() -> None:
    """hyprcoloc into envs/pipeline's own R library (never ~/R), commit pinned."""
    section(f"hyprcoloc (R, GitHub {HYPRCOLOC_REPO}@{HYPRCOLOC_COMMIT[:7]})")
    code = f"""
options(repos = c(CRAN = "https://cloud.r-project.org"))
lib  <- file.path(R.home(), "library")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib)
want <- "{HYPRCOLOC_COMMIT}"
desc <- suppressWarnings(utils::packageDescription("hyprcoloc", lib.loc = lib))
have <- if (is.list(desc) && !is.null(desc$RemoteSha)) desc$RemoteSha else ""
if (!identical(have, want)) {{
  remotes::install_github("{HYPRCOLOC_REPO}", ref = want, lib = lib,
                          dependencies = FALSE, upgrade = "never",
                          build_vignettes = FALSE, force = TRUE)
}} else {{
  cat("[CACHE] hyprcoloc already at pinned commit\\n")
}}
suppressPackageStartupMessages(library(hyprcoloc, lib.loc = lib))
cat("[OK] hyprcoloc", as.character(packageVersion("hyprcoloc", lib.loc = lib)), "in", lib, "\\n")
"""
    run([env_bin(CORE_ENV, "Rscript"), "-e", code], label="R: install hyprcoloc")


def install_spliceai(pms: list[str]) -> None:
    banner("SPLICEAI ISOLATED ENVIRONMENT  ->  envs/spliceai")

    if args.force_recreate:
        remove_env(SPLICEAI_ENV)

    if env_exists(SPLICEAI_ENV):
        log(f"[CACHE] SpliceAI environment exists: {SPLICEAI_ENV}")
    else:
        with_fallback(
            pms,
            lambda pm: [pm, "create", "-y", "-p", SPLICEAI_ENV, *CONDA_FORGE, *SPLICEAI_CONDA],
            "Create envs/spliceai",
        )

    section("Pinned SpliceAI runtime (pip)")
    run([env_bin(SPLICEAI_ENV, "python"), "-m", "pip", "install", "--upgrade", *SPLICEAI_PIP],
        label="pip: tensorflow-cpu + keras + spliceai")


def git_executable() -> str:
    return shutil.which("git") or str(env_bin(CORE_ENV, "git"))


def require_git() -> str:
    git = git_executable()
    if not args.dry_run:
        EXTERNAL.mkdir(parents=True, exist_ok=True)
        MANIFESTS.mkdir(parents=True, exist_ok=True)
        if not Path(git).exists() and shutil.which(git) is None:
            raise RuntimeError("git not found on PATH or in envs/pipeline; install the pipeline env first.")
    return git


def git_checkout(git: str, url: str, repo: Path, commit: str, name: str) -> None:
    """Clone (once) and detach at an exact commit."""
    if not (repo / ".git").is_dir():
        if repo.exists() and any(repo.iterdir()):
            raise RuntimeError(f"{repo} exists but is not a git checkout; move it aside.")
        run([git, "clone", url, repo], label=f"git clone {name}")
    if not args.dry_run:
        rc, out = capture([git, "-C", repo, "rev-parse", "HEAD"])
        if rc == 0 and out.strip() == commit:
            log(f"[CACHE] {name} already at pinned commit {commit}")
            return
        run([git, "fetch", "--all", "--tags", "--prune"], cwd=repo, label=f"git fetch {name}")
    run([git, "checkout", "--detach", commit], cwd=repo, label=f"git checkout {name} {commit[:7]}")


def install_pangolin(pms: list[str]) -> None:
    banner("PANGOLIN ISOLATED ENVIRONMENT  ->  envs/pangolin")

    git = require_git()

    if args.force_recreate:
        remove_env(PANGOLIN_ENV)

    section("Pangolin source checkout")
    if not (PANGOLIN_REPO / ".git").is_dir():
        if PANGOLIN_REPO.exists() and any(PANGOLIN_REPO.iterdir()):
            raise RuntimeError(f"{PANGOLIN_REPO} exists but is not a git checkout; move it aside.")
        run([git, "clone", PANGOLIN_URL, PANGOLIN_REPO], label="git clone Pangolin")

    lock = PANGOLIN_LOCK.read_text().strip() if PANGOLIN_LOCK.is_file() else ""
    if lock:
        log(f"Pinned Pangolin commit: {lock}")
        run([git, "fetch", "--all", "--tags", "--prune"], cwd=PANGOLIN_REPO, label="git fetch")
        run([git, "checkout", "--detach", lock], cwd=PANGOLIN_REPO, label="git checkout pinned commit")
    elif not args.dry_run:
        rc, out = capture([git, "-C", PANGOLIN_REPO, "rev-parse", "HEAD"])
        commit = out.strip()
        if rc != 0 or not commit:
            raise RuntimeError(f"Could not read Pangolin commit: {out}")
        PANGOLIN_LOCK.write_text(commit + "\n", encoding="utf-8")
        log(f"Locked Pangolin commit: {commit}  ->  {PANGOLIN_LOCK}")

    section("Pangolin conda environment")
    if env_exists(PANGOLIN_ENV):
        log(f"[CACHE] Pangolin environment exists: {PANGOLIN_ENV}")
    else:
        with_fallback(
            pms,
            lambda pm: [pm, "create", "-y", "-p", PANGOLIN_ENV, *FORGE_BIOCONDA, *PANGOLIN_CONDA],
            "Create envs/pangolin",
        )

    py = env_bin(PANGOLIN_ENV, "python")
    run([py, "-m", "pip", "install", "--upgrade", *PANGOLIN_PIP], label="pip: PyVCF3")
    run([py, "-m", "pip", "install", "-e", PANGOLIN_REPO, "--no-deps"], label="pip -e Pangolin")


def install_ldsc(pms: list[str]) -> None:
    banner("LDSC / S-LDSC ISOLATED ENVIRONMENT (Python 2.7)  ->  envs/ldsc")

    if args.force_recreate:
        remove_env(LDSC_ENV)

    if env_exists(LDSC_ENV):
        log(f"[CACHE] LDSC environment exists: {LDSC_ENV}")
        return
    with_fallback(
        pms,
        lambda pm: [pm, "create", "-y", "-p", LDSC_ENV, *FORGE_BIOCONDA, *LDSC_CONDA],
        "Create envs/ldsc",
    )


def write_wrapper(path: Path, target: Path) -> None:
    """Small launcher so a git-checkout script runs as a normal command."""
    log(f"Wrapper : {path}  ->  {target}")
    if args.dry_run:
        return
    path.write_text(
        "#!/bin/sh\n"
        f'exec "{env_bin(POPS_ENV, "python")}" "{target}" "$@"\n',
        encoding="utf-8",
    )
    path.chmod(0o755)


def install_pops(pms: list[str]) -> None:
    banner(f"PoPS ISOLATED ENVIRONMENT  ->  envs/pops  (FinucaneLab/pops@{POPS_COMMIT[:7]})")

    git = require_git()

    if args.force_recreate:
        remove_env(POPS_ENV)

    section("PoPS source checkout")
    git_checkout(git, POPS_URL, POPS_REPO, POPS_COMMIT, "PoPS")
    if not args.dry_run:
        (MANIFESTS / "PoPS.commit.txt").write_text(POPS_COMMIT + "\n", encoding="utf-8")

    section("PoPS conda environment")
    if env_exists(POPS_ENV):
        log(f"[CACHE] PoPS environment exists: {POPS_ENV}")
    else:
        with_fallback(
            pms,
            lambda pm: [pm, "create", "-y", "-p", POPS_ENV, *CONDA_FORGE, *POPS_CONDA],
            "Create envs/pops",
        )

    # The authors' exact pins (numpy 1.19.5, scipy 1.5.2, pandas 1.0.5, sklearn 0.24.2).
    run([env_bin(POPS_ENV, "python"), "-m", "pip", "install", "-r", POPS_REPO / "requirements.txt"],
        label="pip: PoPS requirements.txt")

    write_wrapper(env_bin(POPS_ENV, "pops"), POPS_REPO / "pops.py")
    write_wrapper(env_bin(POPS_ENV, "pops_munge_features"), POPS_REPO / "munge_feature_directory.py")


def _sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def install_magma(pms: list[str]) -> None:
    banner(f"MAGMA v{MAGMA_VERSION}  ->  {MAGMA_DIR.relative_to(ROOT)}")

    import zipfile

    if IS_LINUX and platform.machine() in ("x86_64", "AMD64"):
        key = "Linux-x86_64"
    elif IS_MAC:
        key = "Darwin"
    else:
        raise RuntimeError(f"No official MAGMA binary for {SYSTEM} {platform.machine()}; "
                           "build from source: https://cncr.nl/research/magma/")
    filename, url, sha256 = MAGMA_DOWNLOADS[key]
    marker = MAGMA_DIR / ".sha256"

    if args.force_recreate and MAGMA_DIR.exists():
        log(f"[RECREATE] Removing {MAGMA_DIR}")
        if not args.dry_run:
            shutil.rmtree(MAGMA_DIR)

    if MAGMA_BIN.is_file() and marker.is_file() and marker.read_text().strip() == sha256:
        log(f"[CACHE] MAGMA already installed: {MAGMA_BIN}")
    else:
        archive = CACHE / "downloads" / filename
        log(f"Download : {url}")
        log(f"Archive  : {archive}")
        log(f"SHA-256  : {sha256}")
        if args.dry_run:
            return

        archive.parent.mkdir(parents=True, exist_ok=True)
        if not (archive.is_file() and _sha256(archive) == sha256):
            for attempt in range(1, args.retries + 2):
                try:
                    partial = archive.with_suffix(".part")
                    request = urllib.request.Request(url, headers={"User-Agent": "GWAS2m-Step00"})
                    with urllib.request.urlopen(request, timeout=120) as response, open(partial, "wb") as out:
                        shutil.copyfileobj(response, out)
                    partial.replace(archive)
                    break
                except OSError as exc:
                    log(f"[RETRY] MAGMA download failed ({exc}); attempt {attempt}")
                    if attempt > args.retries:
                        raise
                    time.sleep(20 * attempt)

        actual = _sha256(archive)
        if actual != sha256:
            archive.unlink()
            raise RuntimeError(f"MAGMA checksum mismatch: expected {sha256}, got {actual}")
        log("[OK] checksum verified")

        MAGMA_DIR.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive) as bundle:
            for member in bundle.namelist():
                target = (MAGMA_DIR / member).resolve()
                if MAGMA_DIR.resolve() not in target.parents:
                    raise RuntimeError(f"Unsafe path in MAGMA archive: {member}")
            bundle.extractall(MAGMA_DIR)
        MAGMA_BIN.chmod(0o755)
        marker.write_text(sha256 + "\n", encoding="utf-8")
        log(f"[OK] MAGMA installed: {MAGMA_BIN}")

    # Also on PATH whenever envs/pipeline is activated.
    link = env_bin(CORE_ENV, "magma")
    if env_exists(CORE_ENV) and not args.dry_run:
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(MAGMA_BIN)
        log(f"Linked   : {link}  ->  {MAGMA_BIN}")


INSTALLERS = {
    "pipeline": install_pipeline,
    "spliceai": install_spliceai,
    "pangolin": install_pangolin,
    "ldsc": install_ldsc,
    "pops": install_pops,
    "magma": install_magma,
}


def location(component: str) -> str:
    if component == "magma":
        return str(MAGMA_DIR.relative_to(ROOT))
    return f"envs/{component}"


# =============================================================================
# VERIFICATION (read-only)
# =============================================================================

def _first_line(name: str, output: str) -> str:
    """The informative version line of a tool's --version / --help output."""
    import re

    lines = [x.strip() for x in output.splitlines() if x.strip()]
    if name == "vep":
        hits = [x for x in lines if re.match(r"^ensembl(-vep)?\s*:", x)]
        return " | ".join(hits[:2])[:80] if hits else (lines[0][:80] if lines else "")
    if name == "Rscript":
        hits = [x for x in lines if "version" in x.lower()]
        return (hits[0] if hits else (lines[0] if lines else ""))[:80]
    return lines[0][:80] if lines else ""


def verify_tools(env_ids: list[str]) -> list[dict]:
    rows = []
    for name, env, version_args, required, step, purpose, spec in TOOLS:
        if env not in env_ids:
            continue
        candidates = [TOOL_PATHS[name]] if name in TOOL_PATHS else \
            [env_bin(ENVS / env, alias) for alias in TOOL_ALIASES.get(name, [name])]
        path = next((p for p in candidates if p.exists()), candidates[0])
        status, version = "MISSING", ""
        if path.exists():
            status = "OK"
            if version_args:
                rc, out = capture([path, *version_args], timeout=120)
                version = _first_line(name, out)
                # A load error is not a version line (e.g. VEP: "Can't locate DBI.pm").
                error = next((line.strip() for line in out.splitlines() if any(
                    marker in line for marker in ("Can't locate", "error while loading shared libraries",
                                                  "command not found", "Traceback", "bad interpreter"))), "")
                if error or (rc != 0 and not version):
                    status, version = "BROKEN", (error or out.strip()[-80:])[:80]
        rows.append(dict(env=env, type="tool", name=name, status=status, version=version,
                         required=required, step=step, used_by=purpose, spec=spec))

    # wget and curl are each optional, but downloads need at least one of them.
    if "pipeline" in env_ids:
        either = any(r["status"] == "OK" for r in rows if r["env"] == "pipeline" and r["name"] in ("wget", "curl"))
        rows.append(dict(env="pipeline", type="tool", name="GWAS downloader (wget or curl)",
                         status="OK" if either else "MISSING", version="", required=True, step=2,
                         used_by="Step02 runtime GWAS download", spec="wget / curl"))
    return rows


PROBE = r"""
import importlib, json, sys
try:
    from importlib import metadata
except ImportError:
    metadata = None
result = {}
for module, dist in json.loads(sys.argv[1]):
    try:
        mod = importlib.import_module(module)
        version = ""
        if metadata is not None:
            try:
                version = metadata.version(dist)
            except Exception:
                version = ""
        version = version or str(getattr(mod, "__version__", ""))
        result[module] = ["OK", version]
    except Exception as exc:
        result[module] = ["MISSING", type(exc).__name__ + ": " + str(exc)[:120]]
print("@@PROBE@@" + json.dumps(result))
"""


def verify_python(env_ids: list[str]) -> list[dict]:
    rows = []
    for env in env_ids:
        packages = [p for p in PYTHON_PACKAGES if p[2] == env]
        if not packages:
            continue
        python = env_bin(ENVS / env, "python")
        probed: dict = {}
        if python.exists():
            rc, out = capture([python, "-c", PROBE, json.dumps([[p[0], p[1]] for p in packages])],
                              timeout=600)
            marker = [line for line in out.splitlines() if line.startswith("@@PROBE@@")]
            if marker:
                probed = json.loads(marker[-1][len("@@PROBE@@"):])
        for module, dist, _, required, step, method, spec, used in packages:
            status, version = probed.get(module, ["MISSING", "interpreter missing" if not python.exists() else ""])
            rows.append(dict(env=env, type="python", name=dist, status=status,
                             version=version if status == "OK" else "", detail="" if status == "OK" else version,
                             required=required, step=step, used_by=used, spec=spec))
    return rows


def verify_runtime(env_ids: list[str]) -> list[dict]:
    rows = []
    for env, label, code in RUNTIME_CHECKS:
        if env not in env_ids:
            continue
        python = env_bin(ENVS / env, "python")
        rc, out = capture([python, "-c", code], timeout=600) if python.exists() else (1, "interpreter missing")
        last_line = (out.strip().splitlines() or [""])[-1]
        rows.append(dict(env=env, type="runtime", name=label, status="OK" if rc == 0 else "BROKEN",
                         version="", detail="" if rc == 0 else last_line,
                         required=True, step=None, used_by="", spec=""))
    for component, label, command in COMMAND_CHECKS:
        if component not in env_ids:
            continue
        exe = Path(command[0])
        rc, out = capture(command, timeout=300) if exe.exists() else (1, "executable missing")
        # MAGMA prints its version banner but may exit non-zero without a full command.
        ok = rc == 0 or (component == "magma" and "MAGMA" in out)
        last_line = (out.strip().splitlines() or [""])[-1]
        rows.append(dict(env=component, type="runtime", name=label,
                         status="OK" if ok else "BROKEN", version="", detail="" if ok else last_line[:200],
                         required=True, step=None, used_by="", spec=""))
    return rows


def verify_r(env_ids: list[str]) -> list[dict]:
    if "pipeline" not in env_ids:
        return []
    rscript = env_bin(CORE_ENV, "Rscript")
    versions: dict[str, str] = {}
    if rscript.exists():
        names = ",".join(f'"{p[0]}"' for p in R_PACKAGES)
        code = ("cat('@@R@@', 'R', paste0(R.version$major, '.', R.version$minor), '\\n'); "
                f"for (p in c({names})) "
                "cat('@@R@@', p, if (requireNamespace(p, quietly=TRUE)) "
                "as.character(packageVersion(p)) else 'MISSING', '\\n')")
        rc, out = capture([rscript, "-e", code], timeout=600)
        for line in out.splitlines():
            parts = line.split()
            if len(parts) == 3 and parts[0] == "@@R@@":
                versions[parts[1]] = parts[2]
    rows = [dict(env="pipeline", type="R", name="R runtime", status="OK" if "R" in versions else "MISSING",
                 version=versions.get("R", ""), required=True, step=6, used_by="R runtime", spec="r-base")]
    for name, required, step, spec, used in R_PACKAGES:
        version = versions.get(name, "MISSING")
        ok = version != "MISSING"
        rows.append(dict(env="pipeline", type="R", name=name, status="OK" if ok else "MISSING",
                         version=version if ok else "", required=required, step=step, used_by=used, spec=spec))
    return rows


def verify(env_ids: list[str]) -> list[dict]:
    banner("VERIFICATION (read-only)")

    rows: list[dict] = []
    for env in env_ids:
        if env == "magma":
            status = "OK" if MAGMA_BIN.is_file() else "MISSING"
        elif not (ENVS / env).exists():
            status = "MISSING"
        elif not env_exists(ENVS / env) or not env_bin(ENVS / env, "python").exists():
            status = "BROKEN"   # a usable env needs conda-meta/ and bin/python
        else:
            status = "OK"
        rows.append(dict(env=env, type="environment", name=location(env), status=status, version="",
                         required=True, step=None, used_by="", spec=str(ROOT / location(env))))

    rows += verify_tools(env_ids)
    rows += verify_python(env_ids)
    rows += verify_r(env_ids)
    rows += verify_runtime(env_ids)

    for env in env_ids:
        section(location(env))
        log(f"{'TYPE':<12}{'NAME':<36}{'STATUS':<10}{'REQUIRED':<10}VERSION / DETAIL")
        for r in rows:
            if r["env"] != env:
                continue
            detail = r.get("version") or r.get("detail") or ""
            log(f"{r['type']:<12}{r['name'][:35]:<36}{r['status']:<10}"
                f"{'YES' if r['required'] else 'no':<10}{str(detail)[:70]}")

    return rows


# =============================================================================
# REPRODUCIBILITY MANIFESTS
# =============================================================================

def write_manifests(rows: list[dict], pms: list[str]) -> None:
    MANIFESTS.mkdir(parents=True, exist_ok=True)

    path = MANIFESTS / "installed_software.tsv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["ENV", "TYPE", "NAME", "STATUS", "VERSION", "REQUIRED", "FIRST_STEP",
                         "SPECIFICATION", "USED_BY", "DETAIL"])
        for r in rows:
            writer.writerow([r["env"], r["type"], r["name"], r["status"], r.get("version", ""),
                             "YES" if r["required"] else "NO", r.get("step") or "setup",
                             r.get("spec", ""), r.get("used_by", ""), r.get("detail", "")])
    log(f"Software manifest : {path}")

    lister = next((pm for pm in pms if not is_micromamba(pm)), pms[0] if pms else None)
    for env in SELECTED:
        prefix = ENVS / env
        if not env_exists(prefix):
            continue
        if lister:
            rc, out = capture([lister, "list", "-p", prefix], timeout=300)
            if rc == 0:
                (MANIFESTS / f"env_{env}.conda_list.txt").write_text(out, encoding="utf-8")
        rc, out = capture([env_bin(prefix, "python"), "-m", "pip", "freeze"], timeout=300)
        if rc == 0:
            (MANIFESTS / f"env_{env}.pip_freeze.txt").write_text(out, encoding="utf-8")
    log(f"Package lists     : {MANIFESTS}/env_<name>.conda_list.txt / .pip_freeze.txt")


# =============================================================================
# MAIN
# =============================================================================

def main() -> int:
    banner("GWAS2m STEP 00 - SOFTWARE INSTALLATION")
    started = time.time()

    pms = find_package_managers()
    if not preflight(pms):
        log("\nPreflight failed; nothing was installed.")
        return 2

    failed: dict[str, str] = {}

    if not args.verify_only:
        if not args.dry_run:
            for directory in (ENVS, PKGS_DIR, PIP_CACHE, TMP):
                directory.mkdir(parents=True, exist_ok=True)

        # Pangolin needs git; when the host has none it comes from envs/pipeline.
        for component in SELECTED:
            try:
                INSTALLERS[component](pms)
            except (subprocess.CalledProcessError, OSError, RuntimeError) as exc:
                log(f"\n[FAILED] {component}: {exc}")
                failed[component] = str(exc)

        if not args.dry_run:
            ROOT_MARKER.write_text(str(ROOT) + "\n", encoding="utf-8")
            shutil.rmtree(TMP, ignore_errors=True)

    if args.dry_run:
        banner("DRY RUN COMPLETE")
        log("No command was executed. Remove --dry-run to install.")
        return 0

    if args.skip_verify:
        banner("INSTALLATION FINISHED (verification skipped)")
        for component in SELECTED:
            log(f"{component:<10} {'FAILED' if component in failed else 'done'}")
        return 1 if failed else 0

    rows = verify(SELECTED)
    write_manifests(rows, pms)

    banner("SUMMARY")
    missing_envs = {r["env"] for r in rows if r["type"] == "environment" and r["status"] != "OK"}
    # A missing environment is reported once, not once per package inside it.
    reportable = [r for r in rows if r["env"] not in missing_envs or r["type"] == "environment"]
    blocking = [r for r in reportable if r["required"] and r["status"] != "OK"]
    optional = [r for r in reportable if not r["required"] and r["status"] != "OK"]

    for env in SELECTED:
        env_rows = [r for r in rows if r["env"] == env]
        bad = [r for r in env_rows if r["required"] and r["status"] != "OK"]
        status = "FAILED" if env in failed else ("READY" if not bad else "INCOMPLETE")
        ok = sum(r["status"] == "OK" for r in env_rows)
        log(f"{location(env):<30}{status:<12}{ok}/{len(env_rows)} checks OK")

    if failed:
        log("\nInstall errors:")
        for component, message in failed.items():
            log(f"  - {component}: {message}")
    if blocking:
        log("\nBLOCKING (required):")
        for r in blocking:
            log(f"  - {location(r['env'])}: {r['type']} {r['name']} ({r['status']})")
    if optional:
        log("\nOptional / not blocking:")
        for r in optional:
            log(f"  - {location(r['env'])}: {r['type']} {r['name']} ({r['status']})")

    log(f"\nElapsed  : {(time.time() - started) / 60:.1f} min")
    log(f"Full log : {LOG_FILE}")

    if failed or blocking:
        log("\nSOFTWARE INSTALLATION INCOMPLETE. Re-run to resume/repair:")
        log(f"  python Step00_Installation.py --only {' '.join(sorted(set(failed) | {r['env'] for r in blocking}))}")
        return 1

    if args.clean_cache and not args.verify_only:
        log(f"\n[CLEAN] Removing download cache {CACHE}")
        shutil.rmtree(CACHE, ignore_errors=True)

    log("\nSOFTWARE INSTALLATION COMPLETE.")
    log("Use the software with:")
    log(f"  conda activate {CORE_ENV}        # Python/R stack, VEP, plink2, hyprcoloc, CrossMap, magma")
    log(f"  {env_bin(SPLICEAI_ENV, 'spliceai')}")
    log(f"  {env_bin(PANGOLIN_ENV, 'pangolin')}")
    log(f"  {env_bin(LDSC_ENV, 'python')} {env_bin(LDSC_ENV, 'ldsc.py')}")
    log(f"  {env_bin(POPS_ENV, 'pops')}")
    log(f"  {MAGMA_BIN}")
    log("Reference data (MAGMA gene locations, LDSC baseline/weights, PoPS features,")
    log("liftover chain files) is NOT downloaded by this script.")
    log("Next (shared resources, one SLURM task per file):")
    log("  envs/pipeline/bin/python Step00_Download_Resources.py && sbatch Step00_Download_Resources.sh")
    return 0


if __name__ == "__main__":
    sys.exit(main())
