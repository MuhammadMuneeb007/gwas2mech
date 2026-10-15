#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2m STEP 02 - GWAS DISCOVERY, DOWNLOAD AND GWASLAB QC (generic phenotype)
===============================================================================

One file with the two original stages (Step01_Plan_GWAS.py and
Step03_QC_One_GWAS.py, logic unchanged) plus status / statistics:

    # 1. Discover + plan downloads (GWAS Catalog: phenotype match, full summary
    #    statistics, STRICT ancestry, harmonised file, ranked by sample size).
    #    Also records per study whether the phenotype is BINARY (cases/controls)
    #    or CONTINUOUS - from the summary-statistics metadata (-meta.yaml,
    #    "case_control_study", effect_statistic), the sample description (GWAS
    #    Catalog REST API v2), the effect columns of the file and the trait
    #    wording - because fine-mapping / colocalisation differ (SuSiE-RSS /
    #    coloc type "cc" with case fraction s, versus type "quant"). Ambiguous
    #    studies are Unknown + NEEDS_REVIEW (phenotype_type_review.tsv).
    #    PMID and DOI per study (Europe PMC / NCBI).
    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype "Parkinson's disease" --ancestry EUR
    sbatch Step02_Download_GWAS_parkinson_s_disease_european.sh

    # 2. QC is automatic: each array task downloads ONE GWAS and then runs the
    #    GWASLab QC on it. --qc is only needed to re-run QC as its own array
    #    (e.g. after a QC failure):
    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype "Parkinson's disease" --ancestry EUR --qc

    # Progress: downloaded / partial / pending, QC done / failed / pending
    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype "Parkinson's disease" --ancestry EUR --status

    # Publication table per study (authors, journal, PMID, trait/EFO, binary or
    # continuous, N, cases, controls, N_eff, population, countries,
    # cohorts, platform, imputed, variants raw/QC, columns, genome-wide hits,
    # minimum P, lambda GC, available fields) -> GWAS_publication_stats.tsv
    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype "Parkinson's disease" --ancestry EUR --stats

    # Study 2 of the download manifest end to end on this machine: download it
    # (if needed) and QC it - only the planner above has to have run
    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype "Parkinson's disease" --ancestry EUR --index 2

    # PubMed ID -> DOI
    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --pmid-to-doi 34737426 33987465

Missing Python packages are pip-installed into the running interpreter.

Generic: the phenotype is searched as given, plus any --alias terms
(no built-in phenotype list).

SLURM settings for every generated script come from the #SBATCH header of
HPC_CPU.sh (or --template); --partition --qos --account --cpus --mem --time
--max-parallel override them. QC workers run envs/pipeline/bin/python.

Outputs (read by Step04 / Step05):

    01_gwas_catalog/<phenotype>/<ancestry>/
        candidate_studies.tsv, study_selection.tsv, discovery_summary.json
        GWAS_download_manifest.tsv   (+ PUBMEDID, DOI, PHENOTYPE_TYPE, N_CASES, N_CONTROLS, N_TOTAL,
                                      N_EFF, EFFECT_TYPE, CLASSIFICATION_SOURCE / _CONFIDENCE /
                                      _EVIDENCE, NEEDS_REVIEW)
        phenotype_type_review.tsv    studies to verify by hand (VERIFIED_PHENOTYPE_TYPE)
        GWAS_study_metadata.tsv      GCST_ID PHENOTYPE_TYPE N_CASES N_CONTROLS N_TOTAL N_EFF
                                     ANCESTRY GENOME_BUILD EFFECT_TYPE (+ provenance, LD_REFERENCE)
        GWAS_QC_manifest.tsv, GWAS_QC_summary.tsv, GWAS_QC_missing_downloads.tsv
        GWAS_publication_stats.tsv   (--stats)
    02_summary_stats/<phenotype>/<ancestry>/
        raw/<GCST>/<harmonised file>
        qc/<GCST>/<GCST>_GRCh38_QC.tsv.gz      ALL variants; first columns
                      CHR POS SNP EA NEA BETA SE P EAF MAF N INFO Z
                      (BETA = log OR for binary traits; Z = BETA/SE, else signed from P;
                       MAF = min(EAF, 1-EAF); missing N = study N_TOTAL, never N_EFF)
                  <GCST>_GRCh38_significant.tsv.gz
                  <GCST>_study_metadata.tsv/.json
                  <GCST>_QC_summary.tsv/.json (QC_STATUS, CLUMPING_READY, SUSIE_READY,
                      PHENOTYPE_TYPE, N_CASES, N_CONTROLS, N_EFF, EFFECT_TYPE, derivation counts)

Binary vs Continuous is decided per STUDY (the same disease can be a diagnosis
or a score); conflicting or insufficient evidence is never guessed.
--phenotype-type binary|continuous overrides it. N_TOTAL and N_EFF are both kept: each
fine-mapping method chooses its own N.
===============================================================================
"""

from __future__ import annotations

import argparse
import gzip
import html
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote, urljoin


def ensure_packages(packages):
    """pip-install missing packages into THIS interpreter (envs/pipeline when run
    as envs/pipeline/bin/python), with the pip cache inside the project."""

    import importlib
    import importlib.util

    missing = [spec for module, spec in packages.items() if importlib.util.find_spec(module) is None]

    if not missing:
        return

    print(f"Installing missing Python packages into {sys.executable}: {' '.join(missing)}", flush=True)

    env = dict(os.environ)
    env.setdefault("PIP_CACHE_DIR", str(Path.cwd() / ".gwas2m_cache" / "pip"))
    env["PYTHONNOUSERSITE"] = "1"

    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "--no-input", *missing],
        env=env,
    )

    if result.returncode != 0:
        raise SystemExit(
            f"Could not install {' '.join(missing)} (pip exit code {result.returncode}).\n"
            "Run: python Step00_Installation.py"
        )

    importlib.invalidate_caches()


ensure_packages({
    "numpy": "numpy",
    "pandas": "pandas",
    "requests": "requests",
    "urllib3": "urllib3",
    "yaml": "pyyaml",
})

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


PROJECT_ROOT = Path.cwd().resolve()

# Phenotype type per study (download manifest -> QC manifest -> QC summary
# -> study metadata). PHENOTYPE_TYPE is Binary, Continuous or Unknown.
TRAIT_FIELDS = [
    "PHENOTYPE_TYPE",
    "N_CASES",
    "N_CONTROLS",
    "CASE_FRACTION",
    "N_TOTAL",
    "N_EFF",
    "EFFECT_TYPE",
    "EFFECT_REPORTED",
    "CLASSIFICATION_SOURCE",
    "CLASSIFICATION_CONFIDENCE",
    "NEEDS_REVIEW",
    "CLASSIFICATION_EVIDENCE",
]

# Manual verdicts for studies flagged NEEDS_REVIEW (phenotype_type_review.tsv),
# keyed by accession; loaded by main().
REVIEWS = {}

REVIEW_COLUMNS = [
    "STUDY_ACCESSION", "DISEASE_TRAIT", "MAPPED_TRAIT", "INITIAL_SAMPLE_SIZE", "PUBMEDID", "DOI",
    "PHENOTYPE_TYPE", "CLASSIFICATION_CONFIDENCE", "CLASSIFICATION_EVIDENCE",
    "VERIFIED_PHENOTYPE_TYPE", "VERIFIED_N_CASES", "VERIFIED_N_CONTROLS", "REVIEW_NOTE",
]

# Standardised fine-mapping columns (SuSiE-RSS / FINEMAP), first and in this order
# in every <GCST>_GRCh38_QC.tsv.gz. All variants are kept.
FINEMAP_COLUMNS = ["CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P", "EAF", "MAF", "N", "INFO", "Z"]

# GWASLab range-checks these columns in basic_check and REMOVES a variant when
# the value is missing. EAF, INFO and the OR / HR confidence limits are left out
# (default list minus those): a variant without them still has a valid Z, so it
# is kept, and out-of-range EAF / INFO are set to missing instead.
GWASLAB_RANGE_CHECK = [
    "P", "MLOG10P", "Z", "BETA", "SE", "CHISQ", "F", "N", "N_CASE", "N_CONTROL",
    "OR", "HR", "STATUS",
]

# Per-study metadata for fine-mapping (<GCST>_study_metadata.tsv, GWAS_study_metadata.tsv).
STUDY_METADATA_COLUMNS = [
    "GCST_ID", "PHENOTYPE_TYPE", "N_CASES", "N_CONTROLS", "N_TOTAL", "N_EFF",
    "ANCESTRY", "GENOME_BUILD", "EFFECT_TYPE",
]


def effective_n(cases, controls):
    """N_eff = 4 * cases * controls / (cases + controls) (binary GWAS only)."""

    cases, controls = parse_number(cases), parse_number(controls)

    return round(4 * cases * controls / (cases + controls)) if cases and controls else ""


def phenotype_type_value(value):

    value = str(value or "").strip().lower()

    return {"binary": "Binary", "case-control": "Binary", "quantitative": "Continuous",
            "continuous": "Continuous"}.get(value, "Unknown")


def normalise_trait(trait, accession=""):
    """Binary / Continuous / Unknown (older manifests used binary / quantitative),
    then the manual review verdict, then --phenotype-type; N_EFF for binary
    studies and EFFECT_TYPE."""

    trait = dict(trait)

    # older manifests: TRAIT_TYPE / TRAIT_TYPE_SOURCE / PHENOTYPE_TYPE_SOURCE
    if not str(trait.get("PHENOTYPE_TYPE", "") or "").strip() and "TRAIT_TYPE" in trait:
        trait["PHENOTYPE_TYPE"] = trait["TRAIT_TYPE"]
    trait.pop("TRAIT_TYPE", None)
    for old in ("PHENOTYPE_TYPE_SOURCE", "TRAIT_TYPE_SOURCE"):
        if old in trait:
            value = trait.pop(old)
            if not str(trait.get("CLASSIFICATION_SOURCE", "") or "").strip():
                trait["CLASSIFICATION_SOURCE"] = value

    value = phenotype_type_value(trait.get("PHENOTYPE_TYPE"))
    source = str(trait.get("CLASSIFICATION_SOURCE", "") or "")
    review = str(trait.get("NEEDS_REVIEW", "")).strip().lower() in ("true", "1", "yes")

    verdict = REVIEWS.get(str(accession).strip(), {})
    verified = phenotype_type_value(verdict.get("VERIFIED_PHENOTYPE_TYPE"))

    if verified != "Unknown":
        value, review = verified, False
        source = "+".join(x for x in (source, "manual review") if x)
        trait["CLASSIFICATION_CONFIDENCE"] = "verified"
        for field in ("N_CASES", "N_CONTROLS"):
            if parse_number(verdict.get(f"VERIFIED_{field}")):
                trait[field] = parse_number(verdict.get(f"VERIFIED_{field}"))

    forced = getattr(ARGS, "phenotype_type", "auto") if ARGS is not None else "auto"
    if forced in ("binary", "continuous"):
        value, review = forced.capitalize(), False
        if "--phenotype-type" not in source:
            source = "+".join(x for x in (source, "--phenotype-type") if x)
        trait["CLASSIFICATION_CONFIDENCE"] = "user"

    trait["PHENOTYPE_TYPE"] = value
    trait["CLASSIFICATION_SOURCE"] = source
    trait["NEEDS_REVIEW"] = review
    trait["EFFECT_TYPE"] = effect_type(value)

    if value == "Binary":
        cases, controls = parse_number(trait.get("N_CASES")), parse_number(trait.get("N_CONTROLS"))
        trait["N_EFF"] = effective_n(cases, controls)
        trait["CASE_FRACTION"] = round(cases / (cases + controls), 6) if cases and controls else ""
    else:
        trait.update({"N_CASES": "", "N_CONTROLS": "", "CASE_FRACTION": "", "N_EFF": ""})

    return trait


def effect_type(phenotype_type):
    """Effect scale of BETA in the standardised output."""

    return {"Binary": "log(OR)", "Continuous": "BETA"}.get(phenotype_type, "unknown (review)")

GWAS_API_V2 = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"

# Command-line arguments (set in main; read by the script generators).
ARGS = None


# =============================================================================
# GWAS CATALOG
# =============================================================================

STUDIES_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-download-studies-v1.0.3.1.txt"
)

ANCESTRY_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-download-ancestries-v1.0.3.1.txt"
)

SUMMARY_ROOT = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
)

# =============================================================================
# ANCESTRY ALIASES
# =============================================================================

ANCESTRY_ALIASES = {

    "eur": "European",
    "european": "European",

    "afr": "African",
    "african": "African",

    "eas": "East Asian",
    "east asian": "East Asian",

    "sas": "South Asian",
    "south asian": "South Asian",

    "amr": "Hispanic or Latin American",
    "hispanic": "Hispanic or Latin American",
    "latino": "Hispanic or Latin American",
    "hispanic or latin american": "Hispanic or Latin American",

    "mid": "Greater Middle Eastern",
    "middle eastern": "Greater Middle Eastern",
    "greater middle eastern": "Greater Middle Eastern",
}

ANCESTRY_CODES = {
    "European": "EUR",
    "African": "AFR",
    "East Asian": "EAS",
    "South Asian": "SAS",
    "Hispanic or Latin American": "AMR",
    "Greater Middle Eastern": "MID",
}


# =============================================================================
# HELPERS
# =============================================================================

def banner(text):

    print()
    print("=" * 80)
    print(text)
    print("=" * 80)


def normalize(value):

    if value is None:
        return ""

    value = str(value).lower()

    value = re.sub(
        r"[^a-z0-9]+",
        " ",
        value
    )

    return re.sub(
        r"\s+",
        " ",
        value
    ).strip()


def slugify(value):

    value = normalize(value)

    return value.replace(
        " ",
        "_"
    )


def canonical_ancestry(value):

    key = normalize(value)

    return ANCESTRY_ALIASES.get(
        key,
        value
    )


def ancestry_code(value):

    label = canonical_ancestry(value)

    return ANCESTRY_CODES.get(
        label,
        slugify(label).upper()
    )


def build_terms(
    phenotype,
    user_aliases
):

    terms = [
        phenotype
    ]

    terms.extend(
        user_aliases
    )

    final = []
    seen = set()

    for term in terms:

        n = normalize(
            term
        )

        if n and n not in seen:

            seen.add(n)

            final.append(
                term
            )

    return final


def make_session():

    session = requests.Session()

    retries = Retry(
        total=8,
        connect=8,
        read=8,
        status=8,
        backoff_factor=1.5,
        status_forcelist=[
            429,
            500,
            502,
            503,
            504
        ]
    )

    adapter = HTTPAdapter(
        max_retries=retries
    )

    session.mount(
        "https://",
        adapter
    )

    session.headers.update(
        {
            "User-Agent":
                "Generic-GWAS-Pipeline/1.0"
        }
    )

    return session


# =============================================================================
# DOWNLOAD SMALL METADATA FILES
# =============================================================================

def download_metadata(
    session,
    url,
    path,
    refresh=False
):

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if (
        path.exists()
        and path.stat().st_size > 0
        and not refresh
    ):

        print(
            f"Using cached metadata: {path}"
        )

        return

    print(
        f"Downloading:\n  {url}"
    )

    response = session.get(
        url,
        timeout=(30, 300)
    )

    response.raise_for_status()

    temporary = Path(
        str(path) + ".part"
    )

    with open(
        temporary,
        "wb"
    ) as handle:

        handle.write(
            response.content
        )

    os.replace(
        temporary,
        path
    )


# =============================================================================
# SAMPLE SIZE
# =============================================================================

def estimate_sample_size(text):

    if text is None:
        return 0

    numbers = re.findall(
        r"\d[\d,]*",
        str(text)
    )

    values = []

    for number in numbers:

        try:

            x = int(
                number.replace(
                    ",",
                    ""
                )
            )

            if x >= 100:

                values.append(
                    x
                )

        except ValueError:

            continue

    return sum(values)


def parse_number(value):

    if value is None:
        return 0

    numbers = re.findall(
        r"\d[\d,]*",
        str(value)
    )

    if not numbers:
        return 0

    try:

        return int(
            numbers[0].replace(
                ",",
                ""
            )
        )

    except ValueError:

        return 0


# =============================================================================
# SUMMARY STATISTICS AVAILABILITY
# =============================================================================

def has_full_summary_stats(value):

    if value is None:
        return False

    value = str(
        value
    ).strip()

    if not value:
        return False

    bad = {
        "",
        "na",
        "n/a",
        "nan",
        "none",
        "no",
        "false",
        "0",
        "not available",
        "nr",
    }

    return (
        value.lower()
        not in bad
    )


# =============================================================================
# PHENOTYPE MATCHING
# =============================================================================

def contains_term(
    series,
    term
):

    term = normalize(
        term
    )

    if not term:

        return pd.Series(
            False,
            index=series.index
        )

    if (
        " " not in term
        and len(term) <= 5
    ):

        pattern = (
            rf"(?<![a-z0-9])"
            rf"{re.escape(term)}"
            rf"(?![a-z0-9])"
        )

    else:

        pattern = re.escape(
            term
        )

    return series.str.contains(
        pattern,
        regex=True,
        na=False
    )


def find_phenotype_studies(
    studies,
    terms
):

    disease = (
        studies["DISEASE/TRAIT"]
        .fillna("")
        .map(normalize)
    )

    mapped = (
        studies["MAPPED_TRAIT"]
        .fillna("")
        .map(normalize)
        if "MAPPED_TRAIT" in studies.columns
        else pd.Series(
            "",
            index=studies.index
        )
    )

    title = (
        studies["STUDY"]
        .fillna("")
        .map(normalize)
        if "STUDY" in studies.columns
        else pd.Series(
            "",
            index=studies.index
        )
    )

    scores = pd.Series(
        0,
        index=studies.index,
        dtype=int
    )

    for term in terms:

        t = normalize(
            term
        )

        # Exact reported trait
        exact_disease = (
            disease == t
        )

        scores.loc[
            exact_disease
        ] = np.maximum(
            scores.loc[
                exact_disease
            ],
            100
        )

        # Exact ontology mapped trait
        exact_mapped = (
            mapped == t
        )

        scores.loc[
            exact_mapped
        ] = np.maximum(
            scores.loc[
                exact_mapped
            ],
            95
        )

        # Contains reported trait
        disease_hit = contains_term(
            disease,
            term
        )

        scores.loc[
            disease_hit
        ] = np.maximum(
            scores.loc[
                disease_hit
            ],
            80
        )

        # Contains mapped trait
        mapped_hit = contains_term(
            mapped,
            term
        )

        scores.loc[
            mapped_hit
        ] = np.maximum(
            scores.loc[
                mapped_hit
            ],
            70
        )

        # Study title
        title_hit = contains_term(
            title,
            term
        )

        scores.loc[
            title_hit
        ] = np.maximum(
            scores.loc[
                title_hit
            ],
            40
        )

    result = studies[
        scores > 0
    ].copy()

    result[
        "MATCH_SCORE"
    ] = scores.loc[
        result.index
    ]

    result = result.drop_duplicates(
        subset=[
            "STUDY ACCESSION"
        ]
    )

    return result


# =============================================================================
# ANCESTRY FILTER
# =============================================================================

def get_ancestry_information(
    ancestry,
    target
):

    target_norm = normalize(
        target
    )

    # Use initial/discovery stage where available.
    if "STAGE" in ancestry.columns:

        stage = (
            ancestry["STAGE"]
            .fillna("")
            .astype(str)
            .str.lower()
        )

        initial = ancestry[
            stage.str.contains(
                "initial",
                regex=False
            )
        ].copy()

        if not initial.empty:

            ancestry = initial

    records = {}

    for accession, group in ancestry.groupby(
        "STUDY ACCESSION"
    ):

        categories = {

            normalize(x)

            for x in group[
                "BROAD ANCESTRAL CATEGORY"
            ].dropna()

            if normalize(x)
        }

        # STRICT ancestry:
        # every discovery ancestry category must equal requested ancestry.
        strict_match = (
            categories
            == {
                target_norm
            }
        )

        total_n = 0

        if (
            "NUMBER OF INDIVIDUALS"
            in group.columns
        ):

            total_n = sum(
                parse_number(x)
                for x in group[
                    "NUMBER OF INDIVIDUALS"
                ]
            )

        records[
            str(accession)
        ] = {

            "STRICT_MATCH":
                strict_match,

            "ANCESTRY_CATEGORIES":
                " | ".join(
                    sorted(categories)
                ),

            "ANCESTRY_N":
                total_n,
        }

    return records


# =============================================================================
# GWAS CATALOG FTP PATH
# =============================================================================

def gwas_bucket(
    accession
):

    numeric = accession.replace(
        "GCST",
        ""
    )

    number = int(
        numeric
    )

    width = max(
        len(numeric),
        8
    )

    start = (
        ((number - 1) // 1000)
        * 1000
        + 1
    )

    end = (
        start
        + 999
    )

    return (
        f"GCST{start:0{width}d}"
        "-"
        f"GCST{end:0{width}d}"
    )


# =============================================================================
# FIND HARMONISED FILE
# =============================================================================

def find_harmonised_file(
    session,
    accession
):

    bucket = gwas_bucket(
        accession
    )

    directory = (
        f"{SUMMARY_ROOT}/"
        f"{bucket}/"
        f"{accession}/"
        f"harmonised/"
    )

    try:

        response = session.get(
            directory,
            timeout=(20, 60)
        )

    except requests.RequestException:

        return None

    if response.status_code != 200:

        return None

    hrefs = re.findall(
        r'href=["\']([^"\']+)["\']',
        response.text,
        flags=re.I
    )

    files = []

    for href in hrefs:

        href = html.unescape(
            href
        )

        name = href.split("/")[-1]

        lower = name.lower()

        if not name:
            continue

        if lower.endswith(
            (
                ".tbi",
                ".log",
                ".yaml",
                ".yml",
                ".md5",
                ".json",
            )
        ):
            continue

        # Strong preference for harmonised GWAS files.
        if (
            ".h.tsv" in lower
            or lower.endswith(".tsv.gz")
            or lower.endswith(".tsv")
        ):

            files.append(
                name
            )

    if not files:

        return None

    def priority(name):

        lower = name.lower()

        if lower.endswith(
            ".h.tsv.gz"
        ):
            return 0

        if lower.endswith(
            ".h.tsv"
        ):
            return 1

        if lower.endswith(
            ".tsv.gz"
        ):
            return 2

        if lower.endswith(
            ".tsv"
        ):
            return 3

        return 10

    files.sort(
        key=priority
    )

    filename = files[0]

    return {

        "URL":
            urljoin(
                directory,
                filename
            ),

        "FILENAME":
            filename,

        "DIRECTORY":
            directory,
    }


# =============================================================================
# GENERATE SLURM DOWNLOAD SCRIPT
# =============================================================================

def generate_slurm_script(
    manifest,
    script_path,
    output_dir,
    phenotype_slug,
    ancestry_slug,
    n_jobs,
    max_parallel,
    partition,
    time_limit,
):

    manifest = Path(manifest).resolve()
    output_dir = Path(output_dir).resolve()
    script_path = Path(script_path).resolve()

    logs = (
        script_path.parent
        / "logs"
        / "step02_download"
        / phenotype_slug
        / ancestry_slug
    ).resolve()

    logs.mkdir(
        parents=True,
        exist_ok=True
    )

    job_name = (
        f"GWAS_{phenotype_slug}_{ancestry_slug}"
    )[:100]

    # SLURM settings from HPC_CPU.sh (command-line overrides already applied).
    sbatch_header = sbatch_header_block(
        ARGS,
        job_name=job_name,
        output=f"{logs}/step02_%A_%a.out",
        error=f"{logs}/step02_%A_%a.err",
        n_tasks=n_jobs,
        max_parallel=max_parallel,
    )

    text = f"""#!/bin/bash
{sbatch_header}

set -euo pipefail

ROOT="{PROJECT_ROOT}"
cd "$ROOT"
export PATH="$ROOT/envs/pipeline/bin:$PATH"

MANIFEST="{manifest}"
OUTDIR="{output_dir}"
PYTHON="$ROOT/envs/pipeline/bin/python"
[[ -x "$PYTHON" ]] || PYTHON="$(command -v python3)"

ROW=$(awk -F '\\t' -v id="$SLURM_ARRAY_TASK_ID" '
    NR > 1 && $1 == id {{print; exit}}
' "$MANIFEST")

if [[ -z "$ROW" ]]; then
    echo "Could not find array task $SLURM_ARRAY_TASK_ID in $MANIFEST"
    exit 1
fi

IFS=$'\\t' read -r ARRAY_ID ACCESSION DOWNLOAD_URL FILE_NAME REST <<< "$ROW"

echo "============================================================"
echo "GWAS download + QC"
echo "============================================================"
echo "Array task : $ARRAY_ID"
echo "Accession  : $ACCESSION"
echo "URL        : $DOWNLOAD_URL"
echo "Output     : $OUTDIR/$ACCESSION/$FILE_NAME"
echo

# One job per GWAS: download (resumes a .part file; skipped when complete),
# then GWASLab QC (skipped when already complete) - no separate QC submission.
"$PYTHON" {shlex.quote(str(Path(__file__).resolve()))} \\
    --phenotype {shlex.quote(ARGS.phenotype.strip())} \\
    --ancestry {shlex.quote(ARGS.ancestry.strip())} \\
    --index "$ARRAY_ID"

echo
echo "DOWNLOAD + QC COMPLETE: $ACCESSION"
"""

    script_path.write_text(
        text
    )

    script_path.chmod(
        0o755
    )


# =============================================================================
# SLURM SETTINGS (#SBATCH header of HPC_CPU.sh; command line overrides)
# =============================================================================

PER_JOB_KEYS = {"job-name", "output", "error", "array"}


def read_slurm_template(path):

    values = {}

    if not Path(path).is_file():
        return values

    for line in Path(path).read_text(encoding="utf-8").splitlines():

        match = re.match(
            r"^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$",
            line
        )

        if match and (match.group(2) or "").strip():
            values[match.group(1)] = match.group(2).strip()

    return values


def validate_slurm_time(value):
    """SLURM walltime: minutes, MM:SS, HH:MM:SS or D-HH:MM:SS."""

    text = str(value).strip()

    if not re.fullmatch(r"(\d+-)?\d{1,3}(:\d{2}){0,2}", text):
        raise argparse.ArgumentTypeError(
            f"walltime {value!r} must look like 12:00:00 or 2-00:00:00"
        )

    return text


def resolve_slurm(args):
    """Fill args.partition / time / slurm_resources from HPC_CPU.sh + overrides."""

    template = read_slurm_template(PROJECT_ROOT / args.template)

    values = {
        "nodes": "1",
        "ntasks": "1",
        "cpus-per-task": "4",
        "mem": "16G",
        "time": "24:00:00",
    }

    values.update(
        {k: v for k, v in template.items() if k not in PER_JOB_KEYS}
    )

    for key, value in (
        ("partition", args.partition),
        ("qos", args.qos),
        ("account", args.account),
        ("cpus-per-task", args.cpus),
        ("mem", args.mem),
        ("time", args.time),
    ):
        if value not in (None, ""):
            values[key] = str(value)

    args.partition = values.get("partition")
    args.time = values.get("time")
    args.slurm_values = values
    args.slurm_resources = {
        "memory": values.get("mem"),
        "cpus": values.get("cpus-per-task"),
        "max_parallel": int(args.max_parallel or 0),
    }

    return args


def sbatch_header_block(
    args,
    job_name,
    output,
    error,
    n_tasks,
    max_parallel=0,
):

    if n_tasks > 1000:
        raise SystemExit(
            f"{n_tasks} array tasks exceed the 1000-job limit; "
            "use --max-studies."
        )

    values = dict(args.slurm_values)

    array = f"1-{n_tasks}" if n_tasks > 1 else "1"

    if max_parallel and int(max_parallel) > 0:
        array += f"%{int(max_parallel)}"

    directives = {
        **values,
        "job-name": str(job_name)[:100],
        "array": array,
        "output": output,
        "error": error,
    }

    order = [
        "job-name", "account", "partition", "qos", "nodes", "ntasks",
        "cpus-per-task", "mem", "time", "array", "output", "error",
    ]

    lines = [
        f"#SBATCH --{key}={directives[key]}"
        for key in order
        if directives.get(key)
    ]

    lines += [
        f"#SBATCH --{key}={value}"
        for key, value in directives.items()
        if key not in order and value
    ]

    return "\n".join(lines)


# =============================================================================
# DIRECT SINGLE-STUDY DOWNLOAD (same tools / flags as the generated array)
# =============================================================================

def download_study_file(
    download_url,
    target
):
    """wget (else curl) into <target>.part, renamed when complete.

    Returns "SKIPPED_ALREADY_COMPLETE" or "COMPLETE".
    """

    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)

    if target.exists() and target.stat().st_size > 0:
        return "SKIPPED_ALREADY_COMPLETE"

    temporary = Path(str(target) + ".part")
    tools = PROJECT_ROOT / "envs" / "pipeline" / "bin"
    wget = str(tools / "wget") if (tools / "wget").exists() else shutil.which("wget")
    curl = str(tools / "curl") if (tools / "curl").exists() else shutil.which("curl")

    if wget:
        command = [wget, "--continue", "--tries=10", "--timeout=60",
                   "-O", str(temporary), str(download_url)]
    elif curl:
        command = [curl, "--location", "--fail", "--retry", "10", "--retry-delay", "5",
                   "--continue-at", "-", "--output", str(temporary), str(download_url)]
    else:
        raise RuntimeError("neither wget nor curl is available.")

    subprocess.run(command, check=True)
    os.replace(temporary, target)

    return "COMPLETE"


# =============================================================================
# BINARY vs CONTINUOUS PHENOTYPE (study metadata, sample description, effects)
#
# Fine-mapping and colocalisation need it: SuSiE-RSS / coloc use type "cc"
# with the case fraction s for case-control GWAS and type "quant" otherwise.
# Sources, strongest first:
#   1. <harmonised file>-meta.yaml (GWAS-SSF): samples[].case_control_study,
#      sample_size, case_count / control_count when submitted
#   2. GWAS Catalog REST API v2 /studies/<GCST>: initial_sample_size text,
#      e.g. "4,110 European ancestry cases, 13,569 European ancestry controls"
#   3. the cached GWAS Catalog studies table (same text) when the API is down
# =============================================================================

# =============================================================================
# PMID -> DOI (Europe PMC, then NCBI E-utilities; cached)
# =============================================================================

def clean_pmid(value):

    match = re.search(r"\d+", str(value or ""))

    return match.group(0) if match and str(value).strip() not in ("nan", "") else ""


def pmid_to_doi(session, pmids):
    """{pmid: doi} for PubMed IDs; '' when the article has no DOI."""

    cache_file = PROJECT_ROOT / "reference_metadata" / "gwas_catalog" / "pmid_doi.json"

    try:
        cache = json.loads(cache_file.read_text(encoding="utf-8")) if cache_file.is_file() else {}
    except (OSError, ValueError):
        cache = {}

    wanted = sorted({clean_pmid(x) for x in pmids} - {""})
    todo = [x for x in wanted if x not in cache]

    for start in range(0, len(todo), 100):

        batch = todo[start:start + 100]
        query = " OR ".join(f"EXT_ID:{x}" for x in batch) + " AND SRC:MED"
        data = fetch_json(
            session,
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
            f"?query={quote(query)}&format=json&resultType=lite&pageSize=1000",
        ) or {}

        for hit in (data.get("resultList") or {}).get("result") or []:
            if str(hit.get("pmid", "")) in batch and hit.get("doi"):
                cache[str(hit["pmid"])] = hit["doi"]

        missing = [x for x in batch if not cache.get(x)]

        if missing:
            data = fetch_json(
                session,
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=pubmed&retmode=json&id={','.join(missing)}",
            ) or {}
            for pmid in missing:
                ids = ((data.get("result") or {}).get(pmid) or {}).get("articleids") or []
                doi = next((x.get("value", "") for x in ids if x.get("idtype") == "doi"), "")
                if doi or data:
                    cache[pmid] = doi

    try:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(cache, indent=1, sort_keys=True), encoding="utf-8")
    except OSError:
        pass

    return {x: cache.get(x, "") for x in wanted}


def parse_cases_controls(text):

    text = str(text or "")

    def total(word):
        return sum(
            int(x.replace(",", ""))
            for x in re.findall(
                rf"(\d[\d,]*)\s+[^,;\d]*?\b{word}s?\b",
                text,
                flags=re.I
            )
        )

    return total("case"), total("control")


def fetch_json(session, url):

    try:
        response = session.get(url, timeout=(20, 60))
    except requests.RequestException:
        return None

    if response.status_code != 200:
        return None

    try:
        return response.json()
    except ValueError:
        return None


def read_sumstats_metadata(session, harmonised_url):

    try:
        response = session.get(str(harmonised_url) + "-meta.yaml", timeout=(20, 60))
    except requests.RequestException:
        return None

    if response.status_code != 200:
        return None

    try:
        import yaml
        return yaml.safe_load(response.text) or None
    except Exception:
        return None


def read_effect_columns(session, url, max_bytes=262144, max_rows=500, max_seconds=30):
    """Which effect columns carry values in the first rows of the harmonised
    file (streamed, only the first ~256 KB, at most max_seconds): odds_ratio /
    beta / z."""

    import time
    import zlib

    deadline = time.monotonic() + max_seconds

    try:
        response = session.get(str(url), stream=True, timeout=(15, 20))
    except requests.RequestException:
        return set()

    if response.status_code != 200:
        return set()

    text = b""
    decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)

    try:
        for chunk in response.iter_content(65536):
            text += decompressor.decompress(chunk) if str(url).endswith(".gz") else chunk
            if len(text) >= max_bytes or text.count(b"\n") > max_rows or time.monotonic() > deadline:
                break
    except (requests.RequestException, zlib.error):
        pass
    finally:
        response.close()

    lines = text.decode("utf-8", errors="replace").splitlines()[:max_rows]

    if len(lines) < 2:
        return set()

    header = [x.strip().lower() for x in lines[0].split("\t")]
    rows = [x.split("\t") for x in lines[1:-1] or lines[1:]]
    groups = {
        "odds_ratio": {"odds_ratio", "hm_odds_ratio", "or"},
        "beta": {"beta", "hm_beta", "b"},
        "z": {"z", "zscore", "z_score"},
    }
    found = set()

    for name, columns in groups.items():
        for index, column in enumerate(header):
            if column in columns and any(
                index < len(row) and row[index].strip().lower() not in ("", "na", "nan", "null", "none", "#na")
                for row in rows
            ):
                found.add(name)

    return found


# Generic wording of a MEASURED trait (no phenotype names): EFO "... measurement",
# levels, scores, volumes, age at onset, ...
MEASUREMENT_WORDS = re.compile(
    r"\b(measurement|levels?|concentrations?|scores?|volumes?|thickness|density|"
    r"quantitative|age at|duration|counts?|index|height|weight|area|amount|"
    r"progression rate|rate of)\b",
    flags=re.I,
)


def study_trait_information(
    session,
    accession,
    download,
    catalog_row
):
    """Binary / Continuous / Unknown from the study's own metadata.

    Evidence (strongest first):
      strong   sumstats -meta.yaml case_control_study; separate case AND control
               counts (meta.yaml or the sample description)
      moderate odds ratios in the summary statistics (meta.yaml effect_statistic
               or the file's first rows) -> Binary
      weak     a sample description without cases/controls, beta effects,
               measurement wording in the trait -> Continuous
    Conflicting or too little evidence -> Unknown + NEEDS_REVIEW (never guessed);
    the verdict goes in phenotype_type_review.tsv.
    """

    cases = 0
    controls = 0
    n_total = 0
    case_control = None
    sources = []
    evidence = []
    effect_reported = set()

    meta = read_sumstats_metadata(session, download["URL"])

    if isinstance(meta, dict):

        sources.append("sumstats-meta.yaml")

        for sample in meta.get("samples") or []:

            if not isinstance(sample, dict):
                continue

            n_total += parse_number(sample.get("sample_size"))
            cases += parse_number(sample.get("case_count"))
            controls += parse_number(sample.get("control_count"))

            flag = sample.get("case_control_study")

            if isinstance(flag, bool):
                case_control = flag if case_control is None else (case_control or flag)

        statistic = str(meta.get("effect_statistic") or "").strip().lower()

        if "odds" in statistic or statistic in ("or", "log_odds"):
            effect_reported.add("odds_ratio")
        elif statistic in ("beta", "b"):
            effect_reported.add("beta")

    api = fetch_json(session, f"{GWAS_API_V2}/{accession}")

    if isinstance(api, dict) and api.get("initial_sample_size"):
        text = api["initial_sample_size"]
        sources.append("GWAS Catalog API v2")
    else:
        text = catalog_row.get("INITIAL SAMPLE SIZE", "")
        sources.append("GWAS Catalog table")

    text_cases, text_controls = parse_cases_controls(text)

    if not (cases and controls):
        cases, controls = text_cases, text_controls

    if not n_total:
        n_total = cases + controls if cases and controls else estimate_sample_size(text)

    # The file itself is read only when -meta.yaml does not state the effect.
    file_effects = set() if effect_reported else read_effect_columns(session, download["URL"])

    if file_effects:
        sources.append("file header")

    effect_reported |= file_effects

    trait_text = " ".join(
        str(x) for x in (
            catalog_row.get("DISEASE/TRAIT", ""),
            catalog_row.get("MAPPED_TRAIT", ""),
            " ".join(t.get("efo_trait", "") for t in (api or {}).get("efo_traits") or [] if isinstance(t, dict))
            if isinstance(api, dict) else "",
        )
        if str(x) not in ("", "nan")
    )
    case_words = re.search(r"\bcases?\b|\bcontrols?\b", str(text), flags=re.I)

    # ---- evidence ----------------------------------------------------------
    binary, continuous = [], []   # (strength, description)

    if case_control is True:
        binary.append(("strong", "meta.yaml case_control_study=true"))
    elif case_control is False:
        continuous.append(("strong", "meta.yaml case_control_study=false"))

    if cases and controls:
        binary.append(("strong", f"{cases:,} cases / {controls:,} controls"))
    elif cases and not controls:
        evidence.append(f"{cases:,} cases but no controls (case-only / within-case analysis?)")

    if "odds_ratio" in effect_reported:
        binary.append(("moderate", "odds ratios reported"))

    if str(text).strip() and not case_words:
        continuous.append(("weak", "sample described without cases/controls"))

    # Beta alone is compatible with both (binary GWAS often report log-odds as
    # beta), so it only counts when no cases/controls are reported.
    if effect_reported and "odds_ratio" not in effect_reported and case_control is not True \
            and not (cases and controls):
        continuous.append(("weak", "beta effects (no odds ratios)"))

    measurement = MEASUREMENT_WORDS.search(trait_text)
    if measurement:
        continuous.append(("weak", f"trait wording '{measurement.group(0)}'"))

    evidence += [f"Binary: {d} [{s}]" for s, d in binary]
    evidence += [f"Continuous: {d} [{s}]" for s, d in continuous]

    if not evidence:
        evidence.append("no case/control, sample or effect-size information found")

    strengths = {"strong": 3, "moderate": 2, "weak": 1}
    best_binary = max((strengths[s] for s, _ in binary), default=0)
    best_continuous = max((strengths[s] for s, _ in continuous), default=0)
    review = False

    # ---- decision ----------------------------------------------------------
    if best_binary and best_continuous:
        # Strong evidence wins only over WEAK contrary evidence.
        if best_binary == 3 and best_continuous == 1:
            phenotype_type, confidence = "Binary", "moderate"
        elif best_continuous == 3 and best_binary == 1:
            phenotype_type, confidence = "Continuous", "moderate"
        else:
            phenotype_type, confidence, review = "Unknown", "conflicting", True
    elif best_binary:
        phenotype_type = "Binary"
        confidence = "high" if best_binary == 3 else "moderate"
    elif best_continuous == 3 or len(continuous) >= 2:
        phenotype_type = "Continuous"
        confidence = "high" if best_continuous == 3 else "moderate"
    else:
        phenotype_type, confidence, review = "Unknown", "insufficient", True

    if cases and not controls:
        review = True

    # Proxy cases (family history, "GWAX") are not diagnosed cases: their effect
    # sizes are attenuated and 4*cases*controls/N overstates the information.
    if re.search(r"\bproxy\b|by[- ]proxy|family history", str(text), flags=re.I):
        review = True
        evidence.append("sample includes PROXY cases (GWAX / family history): N_EFF is approximate")

    forced = getattr(ARGS, "phenotype_type", "auto") if ARGS is not None else "auto"
    if forced in ("binary", "continuous"):
        phenotype_type, confidence, review = forced.capitalize(), "user", False
        sources.append("--phenotype-type")

    binary_type = phenotype_type == "Binary"

    # Counts are kept as reported (a later --phenotype-type binary can use them);
    # normalise_trait() blanks them for continuous studies.
    return {
        "PHENOTYPE_TYPE": phenotype_type,
        "N_CASES": cases or "",
        "N_CONTROLS": controls or "",
        "CASE_FRACTION": round(cases / (cases + controls), 6) if binary_type and cases and controls else "",
        "N_TOTAL": n_total or "",
        "N_EFF": effective_n(cases, controls) if binary_type else "",
        "EFFECT_TYPE": effect_type(phenotype_type),
        "EFFECT_REPORTED": "+".join(sorted(effect_reported)),
        "CLASSIFICATION_SOURCE": "+".join(sources),
        "CLASSIFICATION_CONFIDENCE": confidence,
        "NEEDS_REVIEW": review,
        "CLASSIFICATION_EVIDENCE": "; ".join(evidence),
    }


# =============================================================================
# STUDY SELECTION AUDIT TABLE
#
# Records EVERY phenotype-matched study and why it was or was not selected.
# Bookkeeping only: it never changes which studies are selected.
# =============================================================================

def new_selection_table(
    candidates,
    phenotype,
    ancestry_name,
    ancestry_info
):

    rows = []

    for _, row in candidates.iterrows():

        accession = str(
            row["STUDY ACCESSION"]
        )

        info = ancestry_info.get(
            accession,
            {}
        )

        rows.append(
            {
                "STUDY_ACCESSION": accession,
                "REQUESTED_PHENOTYPE": phenotype,
                "DISEASE_TRAIT": row.get("DISEASE/TRAIT", ""),
                "MAPPED_TRAIT": row.get("MAPPED_TRAIT", ""),
                "MATCH_SCORE": int(row["MATCH_SCORE"]),
                "REQUESTED_ANCESTRY": ancestry_name,
                "ANCESTRY_CATEGORIES": info.get("ANCESTRY_CATEGORIES", ""),
                "ANCESTRY_MATCH": bool(info.get("STRICT_MATCH", False)),
                "INITIAL_SAMPLE_SIZE": row.get("INITIAL SAMPLE SIZE", ""),
                "SUMMARY_STATS_AVAILABLE": has_full_summary_stats(
                    row.get("FULL SUMMARY STATISTICS")
                ),
                "HARMONISED_FILE_AVAILABLE": "NOT_CHECKED",
                "DOWNLOAD_URL": "",
                "ELIGIBLE": False,
                "EXCLUSION_REASON": "",
                "RANK": "",
            }
        )

    columns = [
        "STUDY_ACCESSION", "REQUESTED_PHENOTYPE", "DISEASE_TRAIT",
        "MAPPED_TRAIT", "MATCH_SCORE", "REQUESTED_ANCESTRY",
        "ANCESTRY_CATEGORIES", "ANCESTRY_MATCH", "INITIAL_SAMPLE_SIZE",
        "SUMMARY_STATS_AVAILABLE", "HARMONISED_FILE_AVAILABLE",
        "DOWNLOAD_URL", "ELIGIBLE", "EXCLUSION_REASON", "RANK",
    ]

    selection = pd.DataFrame(
        rows,
        columns=columns
    )

    selection.index = selection[
        "STUDY_ACCESSION"
    ].tolist()

    # Reasons follow the planner's filter order:
    # full summary statistics -> strict ancestry -> harmonised file.
    no_sumstats = ~selection[
        "SUMMARY_STATS_AVAILABLE"
    ].astype(bool)

    selection.loc[
        no_sumstats,
        "EXCLUSION_REASON"
    ] = "NO_FULL_SUMMARY_STATISTICS"

    no_ancestry = (
        (selection["EXCLUSION_REASON"] == "")
        & ~selection["ANCESTRY_MATCH"].astype(bool)
    )

    selection.loc[
        no_ancestry
        & (selection["ANCESTRY_CATEGORIES"] == ""),
        "EXCLUSION_REASON"
    ] = "ANCESTRY_NOT_REPORTED"

    selection.loc[
        no_ancestry
        & (selection["ANCESTRY_CATEGORIES"] != ""),
        "EXCLUSION_REASON"
    ] = "ANCESTRY_MISMATCH"

    return selection


def write_study_selection(
    state
):

    analysis_dir = state.get(
        "analysis_dir"
    )

    if analysis_dir is None:
        return

    selection = state.get(
        "selection"
    )

    if selection is None:

        selection = pd.DataFrame(
            columns=[
                "STUDY_ACCESSION",
                "REQUESTED_PHENOTYPE",
                "ELIGIBLE",
                "EXCLUSION_REASON",
            ]
        )

    selection = selection.copy()

    # Studies that were never checked for a harmonised file because
    # --max-studies was reached are explicitly labelled, not dropped.
    unchecked = (
        (selection["EXCLUSION_REASON"] == "")
        & ~selection["ELIGIBLE"].astype(bool)
    )

    selection.loc[
        unchecked,
        "EXCLUSION_REASON"
    ] = "NOT_CHECKED_MAX_STUDIES_LIMIT"

    selection.to_csv(
        analysis_dir / "study_selection.tsv",
        sep="\t",
        index=False
    )

    summary = {
        key: value
        for key, value in state.items()
        if key not in {"analysis_dir", "selection"}
    }

    (analysis_dir / "discovery_summary.json").write_text(
        json.dumps(
            summary,
            indent=2,
            default=str
        ),
        encoding="utf-8"
    )


# =============================================================================
# MAIN
# =============================================================================

def plan(args):
    """Run Step01 planning and always write study_selection.tsv.

    Returns a dict with the manifest/candidate paths and counts.
    Raises RuntimeError exactly where the planner always has
    (e.g. no eligible studies); the selection table is still written.
    """

    state = {}

    try:

        return _plan(
            args,
            state
        )

    finally:

        write_study_selection(
            state
        )


def _plan(args, state):

    phenotype = args.phenotype.strip()

    if not phenotype:

        raise RuntimeError(
            "Phenotype cannot be empty."
        )

    if args.max_studies < 0:

        raise RuntimeError(
            "--max-studies must be >= 0."
        )

    if args.max_parallel < 0:

        raise RuntimeError(
            "--max-parallel must be >= 0 (0 = no throttle)."
        )

    ancestry_name = canonical_ancestry(
        args.ancestry
    )

    ancestry_code_value = ancestry_code(
        ancestry_name
    )

    phenotype_slug = slugify(
        phenotype
    )

    ancestry_slug = slugify(
        ancestry_name
    )

    terms = build_terms(
        phenotype,
        args.alias
    )

    # -------------------------------------------------------------------------
    # Directories
    # -------------------------------------------------------------------------

    cache = Path(
        "reference_metadata/gwas_catalog"
    )

    analysis_dir = (
        Path("01_gwas_catalog")
        / phenotype_slug
        / ancestry_slug
    )

    download_dir = (
        Path("02_summary_stats")
        / phenotype_slug
        / ancestry_slug
        / "raw"
    )

    cache.mkdir(
        parents=True,
        exist_ok=True
    )

    analysis_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    download_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    state["analysis_dir"] = analysis_dir
    state["REQUESTED_PHENOTYPE"] = phenotype
    state["REQUESTED_ANCESTRY"] = ancestry_name
    state["ANCESTRY_CODE"] = ancestry_code_value
    state["SEARCH_TERMS"] = terms
    state["MAX_STUDIES"] = args.max_studies

    studies_file = (
        cache
        / "gwas_catalog_studies.tsv"
    )

    ancestry_file = (
        cache
        / "gwas_catalog_ancestry.tsv"
    )

    session = make_session()

    # -------------------------------------------------------------------------
    # Configuration
    # -------------------------------------------------------------------------

    banner(
        "GWAS DOWNLOAD PLANNER"
    )

    print(
        f"Phenotype : {phenotype}"
    )

    print(
        f"Ancestry  : {ancestry_name}"
    )

    print(
        f"Code      : {ancestry_code_value}"
    )

    print(
        f"Partition : {args.partition or 'scheduler default'}"
    )

    print(
        f"Walltime  : {args.time}"
    )

    print(
        "\nSearch terms:"
    )

    for term in terms:

        print(
            f"  - {term}"
        )

    # -------------------------------------------------------------------------
    # Metadata
    # -------------------------------------------------------------------------

    banner(
        "GWAS CATALOG METADATA"
    )

    download_metadata(
        session,
        STUDIES_URL,
        studies_file,
        args.refresh
    )

    download_metadata(
        session,
        ANCESTRY_URL,
        ancestry_file,
        args.refresh
    )

    studies = pd.read_csv(
        studies_file,
        sep="\t",
        dtype=str,
        low_memory=False
    )

    ancestry = pd.read_csv(
        ancestry_file,
        sep="\t",
        dtype=str,
        low_memory=False
    )

    state["GWAS_CATALOG_STUDIES_URL"] = STUDIES_URL
    state["GWAS_CATALOG_ANCESTRY_URL"] = ANCESTRY_URL
    state["N_CATALOG_STUDIES"] = int(len(studies))

    ancestry_info = get_ancestry_information(
        ancestry,
        ancestry_name
    )

    # -------------------------------------------------------------------------
    # Phenotype
    # -------------------------------------------------------------------------

    banner(
        "PHENOTYPE FILTER"
    )

    candidates = find_phenotype_studies(
        studies,
        terms
    )

    print(
        f"Phenotype matches: {len(candidates):,}"
    )

    state["selection"] = new_selection_table(
        candidates,
        phenotype,
        ancestry_name,
        ancestry_info
    )

    state["N_PHENOTYPE_MATCHED"] = int(len(candidates))

    if candidates.empty:

        raise RuntimeError(
            f"No GWAS studies found for phenotype: {phenotype}"
        )

    # -------------------------------------------------------------------------
    # Full summary statistics
    # -------------------------------------------------------------------------

    if (
        "FULL SUMMARY STATISTICS"
        not in candidates.columns
    ):

        raise RuntimeError(
            "GWAS Catalog studies file does not contain "
            "'FULL SUMMARY STATISTICS'."
        )

    candidates = candidates[
        candidates[
            "FULL SUMMARY STATISTICS"
        ].map(
            has_full_summary_stats
        )
    ].copy()

    print(
        f"With full summary statistics: {len(candidates):,}"
    )

    state["N_WITH_FULL_SUMMARY_STATISTICS"] = int(len(candidates))

    if candidates.empty:

        raise RuntimeError(
            "No phenotype-matching studies have full summary statistics."
        )

    # -------------------------------------------------------------------------
    # Ancestry
    # -------------------------------------------------------------------------

    banner(
        "ANCESTRY FILTER"
    )

    strict_accessions = {

        accession

        for accession, info
        in ancestry_info.items()

        if info[
            "STRICT_MATCH"
        ]
    }

    candidates = candidates[
        candidates[
            "STUDY ACCESSION"
        ].astype(str).isin(
            strict_accessions
        )
    ].copy()

    print(
        f"Strict {ancestry_name} studies: "
        f"{len(candidates):,}"
    )

    state["N_FULL_SUMMARY_STATISTICS_AND_ANCESTRY_MATCHED"] = int(len(candidates))

    if candidates.empty:

        raise RuntimeError(
            f"No strict {ancestry_name} GWAS studies found."
        )

    # -------------------------------------------------------------------------
    # Sample size ranking
    # -------------------------------------------------------------------------

    ancestry_n = []

    categories = []

    for accession in candidates[
        "STUDY ACCESSION"
    ].astype(str):

        info = ancestry_info.get(
            accession,
            {}
        )

        ancestry_n.append(
            int(
                info.get(
                    "ANCESTRY_N",
                    0
                )
            )
        )

        categories.append(
            info.get(
                "ANCESTRY_CATEGORIES",
                ""
            )
        )

    candidates[
        "ANCESTRY_N"
    ] = ancestry_n

    candidates[
        "ANCESTRY_CATEGORIES"
    ] = categories

    if (
        "INITIAL SAMPLE SIZE"
        in candidates.columns
    ):

        candidates[
            "ESTIMATED_INITIAL_N"
        ] = candidates[
            "INITIAL SAMPLE SIZE"
        ].map(
            estimate_sample_size
        )

    else:

        candidates[
            "ESTIMATED_INITIAL_N"
        ] = 0

    candidates[
        "RANKING_N"
    ] = candidates[
        [
            "ANCESTRY_N",
            "ESTIMATED_INITIAL_N"
        ]
    ].max(
        axis=1
    )

    candidates = candidates.sort_values(
        by=[
            "MATCH_SCORE",
            "RANKING_N"
        ],
        ascending=[
            False,
            False
        ]
    )

    for rank, accession in enumerate(
        candidates["STUDY ACCESSION"].astype(str),
        start=1
    ):

        state["selection"].at[accession, "RANK"] = rank

    candidate_file = (
        analysis_dir
        / "candidate_studies.tsv"
    )

    candidates.to_csv(
        candidate_file,
        sep="\t",
        index=False
    )

    # -------------------------------------------------------------------------
    # Check harmonised data
    # -------------------------------------------------------------------------

    banner(
        "CHECKING HARMONISED GWAS FILES"
    )

    eligible = []

    for _, row in candidates.iterrows():

        accession = str(
            row[
                "STUDY ACCESSION"
            ]
        )

        print(
            f"{accession}: ",
            end="",
            flush=True
        )

        result = find_harmonised_file(
            session,
            accession
        )

        if result is None:

            print(
                "no harmonised file"
            )

            state["selection"].at[accession, "HARMONISED_FILE_AVAILABLE"] = False
            state["selection"].at[accession, "EXCLUSION_REASON"] = "NO_HARMONISED_FILE"

            continue

        print(
            result[
                "FILENAME"
            ]
        )

        state["selection"].at[accession, "HARMONISED_FILE_AVAILABLE"] = True
        state["selection"].at[accession, "DOWNLOAD_URL"] = result["URL"]
        state["selection"].at[accession, "ELIGIBLE"] = True

        trait = study_trait_information(
            session,
            accession,
            result,
            row
        )

        print(
            f"    phenotype: {trait['PHENOTYPE_TYPE']} ({trait['CLASSIFICATION_CONFIDENCE']})"
            f" | cases {trait['N_CASES'] or '-'}"
            f" | controls {trait['N_CONTROLS'] or '-'}"
            f" | N {trait['N_TOTAL'] or '-'}"
            f" | N_eff {trait['N_EFF'] or '-'}"
            f" | effect {trait['EFFECT_REPORTED'] or '-'}"
            + ("  ** NEEDS REVIEW **" if trait["NEEDS_REVIEW"] else "")
        )
        print(f"      evidence: {trait['CLASSIFICATION_EVIDENCE'] or '-'}")

        eligible.append(
            {
                "ROW":
                    row,

                "DOWNLOAD":
                    result,

                "TRAIT":
                    trait,
            }
        )

        if (
            args.max_studies > 0
            and len(eligible)
            >= args.max_studies
        ):

            break

    if not eligible:

        raise RuntimeError(
            "No matching harmonised GWAS files were found."
        )

    # -------------------------------------------------------------------------
    # Manifest
    # -------------------------------------------------------------------------

    manifest_rows = []

    dois = pmid_to_doi(
        session,
        [item["ROW"].get("PUBMEDID", "") for item in eligible]
    )

    for i, item in enumerate(
        eligible,
        start=1
    ):

        row = item[
            "ROW"
        ]

        download = item[
            "DOWNLOAD"
        ]

        manifest_rows.append(
            {

                "ARRAY_ID":
                    i,

                "STUDY_ACCESSION":
                    str(
                        row[
                            "STUDY ACCESSION"
                        ]
                    ),

                "DOWNLOAD_URL":
                    download[
                        "URL"
                    ],

                "FILE_NAME":
                    download[
                        "FILENAME"
                    ],

                "PHENOTYPE":
                    phenotype,

                # Backward-compatible human-readable ancestry label.
                "ANCESTRY":
                    ancestry_name,

                # Explicit machine-readable ancestry code for downstream
                # resource selection (EUR/AFR/EAS/SAS/AMR/MID).
                "ANCESTRY_CODE":
                    ancestry_code_value,

                "ANCESTRY_LABEL":
                    ancestry_name,

                "RANKING_N":
                    int(
                        row[
                            "RANKING_N"
                        ]
                    ),

                "MATCH_SCORE":
                    int(
                        row[
                            "MATCH_SCORE"
                        ]
                    ),

                "DISEASE_TRAIT":
                    row.get(
                        "DISEASE/TRAIT",
                        ""
                    ),

                "MAPPED_TRAIT":
                    row.get(
                        "MAPPED_TRAIT",
                        ""
                    ),

                "FIRST_AUTHOR":
                    row.get(
                        "FIRST AUTHOR",
                        ""
                    ),

                "PUBMEDID":
                    clean_pmid(
                        row.get(
                            "PUBMEDID",
                            ""
                        )
                    ),

                "DOI":
                    dois.get(
                        clean_pmid(row.get("PUBMEDID", "")),
                        ""
                    ),

                "INITIAL_SAMPLE_SIZE":
                    row.get(
                        "INITIAL SAMPLE SIZE",
                        ""
                    ),

                # Binary vs continuous phenotype (fine-mapping / coloc need it).
                **item[
                    "TRAIT"
                ],

            }
        )

    manifest = pd.DataFrame(
        manifest_rows
    )

    manifest_file = (
        analysis_dir
        / "GWAS_download_manifest.tsv"
    )

    manifest.to_csv(
        manifest_file,
        sep="\t",
        index=False
    )

    state["N_ELIGIBLE"] = int(len(manifest))
    state["CANDIDATE_FILE"] = str(candidate_file)
    state["DOWNLOAD_MANIFEST"] = str(manifest_file)

    # -------------------------------------------------------------------------
    # Generate Step02
    # -------------------------------------------------------------------------

    script_file = Path(
        f"Step02_Download_GWAS_"
        f"{phenotype_slug}_"
        f"{ancestry_slug}.sh"
    )

    generate_slurm_script(

        manifest=manifest_file,

        script_path=script_file,

        output_dir=download_dir,

        phenotype_slug=phenotype_slug,

        ancestry_slug=ancestry_slug,

        n_jobs=len(manifest),

        max_parallel=args.max_parallel,

        partition=args.partition,

        time_limit=args.time,
    )

    # -------------------------------------------------------------------------
    # Finish
    # -------------------------------------------------------------------------

    banner(
        "PLANNING COMPLETE"
    )

    print(
        f"Phenotype           : {phenotype}"
    )

    print(
        f"Ancestry            : {ancestry_name}"
    )

    print(
        f"Ancestry code       : {ancestry_code_value}"
    )

    print(
        f"SLURM partition     : {args.partition or 'scheduler default'}"
    )

    print(
        f"SLURM walltime      : {args.time}"
    )

    print(
        f"Eligible GWAS files : {len(manifest):,}"
    )

    print(
        f"\nCandidate table:\n  {candidate_file}"
    )

    print(
        f"\nDownload manifest:\n  {manifest_file}"
    )

    print(
        f"\nGenerated SLURM script:\n  {script_file}"
    )

    print()

    show = [
        column
        for column in [
            "ARRAY_ID",
            "STUDY_ACCESSION",
            "RANKING_N",
            "PHENOTYPE_TYPE",
            "N_CASES",
            "N_CONTROLS",
            "FILE_NAME",
        ]
        if column in manifest.columns
    ]

    print(
        manifest[
            show
        ].to_string(
            index=False
        )
    )

    print()

    print(
        "NEXT COMMAND:"
    )

    print()

    print(
        f"  sbatch {script_file}"
    )

    print()

    print(
        "Each array task downloads one GWAS and then QCs it (no separate QC step)."
    )

    print()

    print(
        "THEN:"
    )

    print(
        f"  envs/pipeline/bin/python {Path(__file__).name} "
        f"--phenotype \"{phenotype}\" --ancestry {ancestry_code_value} --status"
    )

    print(
        "  (--qc only re-runs QC as its own array, e.g. after a QC failure)"
    )

    print()

    return {
        "phenotype": phenotype,
        "phenotype_slug": phenotype_slug,
        "ancestry_label": ancestry_name,
        "ancestry_code": ancestry_code_value,
        "ancestry_slug": ancestry_slug,
        "analysis_dir": analysis_dir,
        "download_dir": download_dir,
        "candidate_file": candidate_file,
        "manifest_file": manifest_file,
        "n_eligible": len(manifest),
    }


# =============================================================================
# CONFIGURATION
# =============================================================================

# SLURM partition/time/memory/CPUs come from HPC_CPU.sh (or the command line);
# CPUS_PER_TASK is only the GWASLab thread fallback outside SLURM.
CPUS_PER_TASK = 4

P_THRESHOLD = 5e-8

# =============================================================================
# TEXT HELPERS
# =============================================================================

def normalize_text(value):

    value = str(value).lower()

    value = re.sub(
        r"[^a-z0-9]+",
        " ",
        value,
    )

    value = re.sub(
        r"\s+",
        " ",
        value,
    )

    return value.strip()

# =============================================================================
# ANCESTRY NORMALIZATION
# =============================================================================

QC_ANCESTRY_ALIASES = {
    "eur": ("EUR", "European"),
    "european": ("EUR", "European"),
    "afr": ("AFR", "African"),
    "african": ("AFR", "African"),
    "eas": ("EAS", "East Asian"),
    "east asian": ("EAS", "East Asian"),
    "sas": ("SAS", "South Asian"),
    "south asian": ("SAS", "South Asian"),
    "amr": ("AMR", "Hispanic or Latin American"),
    "hispanic": ("AMR", "Hispanic or Latin American"),
    "latino": ("AMR", "Hispanic or Latin American"),
    "hispanic or latin american": ("AMR", "Hispanic or Latin American"),
    "mid": ("MID", "Greater Middle Eastern"),
    "middle eastern": ("MID", "Greater Middle Eastern"),
    "greater middle eastern": ("MID", "Greater Middle Eastern"),
}


def canonical_ancestry_pair(value):

    key = normalize_text(value)

    if key not in QC_ANCESTRY_ALIASES:

        supported = (
            "EUR, AFR, EAS, SAS, AMR, MID"
        )

        raise ValueError(
            f"Unsupported ancestry: {value!r}. "
            f"Supported values: {supported}"
        )

    return QC_ANCESTRY_ALIASES[key]

# =============================================================================
# ONE QC MANIFEST ROW (shared by planner_mode and Step01_10_Run.py)
# =============================================================================

def build_qc_row(
    task_id,
    accession,
    input_file,
    output_directory,
    phenotype,
    ancestry_code,
    ancestry_label,
    extra=None,
):

    return {
        "QC_TASK_ID":
            task_id,

        "STUDY_ACCESSION":
            accession,

        "INPUT_FILE":
            str(input_file),

        "OUTPUT_DIR":
            str(
                output_directory.resolve()
            ),

        "PHENOTYPE":
            phenotype,

        "ANCESTRY":
            ancestry_label,

        "ANCESTRY_CODE":
            ancestry_code,

        "ANCESTRY_LABEL":
            ancestry_label,

        "INPUT_SIZE_BYTES":
            input_file.stat().st_size,

        # Binary vs continuous phenotype from the download manifest.
        **(extra or {}),
    }


# =============================================================================
# COLUMN DETECTION
# =============================================================================

def find_column(
    dataframe,
    candidates,
):

    lookup = {

        str(column).lower():
            column

        for column in dataframe.columns
    }

    for candidate in candidates:

        candidate_lower = (
            candidate.lower()
        )

        if candidate_lower in lookup:

            return lookup[
                candidate_lower
            ]

    return None


# =============================================================================
# STANDARDISED FINE-MAPPING COLUMNS
#
# CHR POS SNP EA NEA BETA SE P EAF MAF N INFO Z first, in that order, then the
# remaining GWASLab columns. No variant is removed here: SuSiE-RSS / FINEMAP
# need every variant in a locus, not only the significant ones.
# =============================================================================

def standardize_finemap_columns(
    data,
    trait,
):

    from scipy.stats import norm

    info = {}

    def numeric(names):

        column = find_column(data, names)

        if column is None:
            return pd.Series(np.nan, index=data.index, dtype="float64"), None

        return pd.to_numeric(data[column], errors="coerce").astype("float64"), column

    out = pd.DataFrame(index=data.index)

    out["CHR"] = data[find_column(data, ["CHR"])]
    out["POS"] = data[find_column(data, ["POS"])]

    ea_col = find_column(data, ["EA"])
    nea_col = find_column(data, ["NEA"])
    ea = data[ea_col].astype("string") if ea_col else pd.Series(pd.NA, index=data.index, dtype="string")
    nea = data[nea_col].astype("string") if nea_col else pd.Series(pd.NA, index=data.index, dtype="string")

    # SNP: rsID when valid, otherwise CHR:POS:NEA:EA (CHR:POS without alleles).
    rsid_col = find_column(data, ["rsID", "RSID"])
    rsid = (
        data[rsid_col].astype("string").str.strip()
        if rsid_col else pd.Series(pd.NA, index=data.index, dtype="string")
    )
    valid_rsid = rsid.str.match(r"^rs\d+$", na=False)

    position_id = out["CHR"].astype("string") + ":" + out["POS"].astype("string")
    allele_id = position_id + ":" + nea + ":" + ea
    out["SNP"] = rsid.where(valid_rsid, allele_id.fillna(position_id))
    info["N_SNP_FROM_RSID"] = int(valid_rsid.sum())
    info["N_SNP_FROM_POSITION"] = int((~valid_rsid).sum())

    out["EA"] = ea
    out["NEA"] = nea

    # BETA (log OR for binary traits) and SE.
    beta, beta_col = numeric(["BETA"])
    se, _ = numeric(["SE"])
    odds, _ = numeric(["OR"])
    or_low, _ = numeric(["OR_95L"])
    or_high, _ = numeric(["OR_95U"])
    p, _ = numeric(["P"])

    with np.errstate(divide="ignore", invalid="ignore"):

        beta_from_or = beta.isna() & (odds > 0)
        beta = beta.where(~beta_from_or, np.log(odds))

        se_from_ci = se.isna() & (or_low > 0) & (or_high > or_low)
        se = se.where(~se_from_ci, (np.log(or_high) - np.log(or_low)) / (2 * norm.isf(0.025)))

        z_p = pd.Series(norm.isf(p.clip(lower=1e-300, upper=1) / 2), index=data.index)

        se_from_p = se.isna() & beta.notna() & (beta != 0) & (p > 0) & (p < 1)
        se = se.where(~se_from_p, beta.abs() / z_p)

    info["N_BETA_FROM_OR"] = int(beta_from_or.sum())
    info["N_SE_FROM_OR_CI"] = int(se_from_ci.sum())
    info["N_SE_FROM_BETA_P"] = int(se_from_p.sum())

    out["BETA"] = beta
    out["SE"] = se
    out["P"] = p

    eaf, _ = numeric(["EAF"])
    bad_eaf = (eaf < 0) | (eaf > 1)
    eaf = eaf.mask(bad_eaf)
    info["N_EAF_OUT_OF_RANGE_SET_MISSING"] = int(bad_eaf.sum())
    out["EAF"] = eaf
    out["MAF"] = np.minimum(eaf, 1 - eaf)

    # N: per-variant N when the GWAS has it, otherwise the study N_TOTAL.
    # N_EFF (binary) is kept in the study metadata, never substituted here.
    n, n_col = numeric(["N"])
    study_n = parse_number(trait.get("N_TOTAL"))
    n_from_study = n.isna() & bool(study_n)

    if study_n:
        n = n.fillna(float(study_n))

    out["N"] = n
    info["N_SOURCE"] = (
        ("variant" if n_col else "")
        + ("+" if n_col and n_from_study.any() else "")
        + ("study N_TOTAL" if n_from_study.any() else "")
    ) or "missing"
    info["N_FILLED_FROM_STUDY_TOTAL"] = int(n_from_study.sum())

    info_score, _ = numeric(["INFO"])
    bad_info = (info_score < 0) | (info_score > 2)
    info["N_INFO_OUT_OF_RANGE_SET_MISSING"] = int(bad_info.sum())
    out["INFO"] = info_score.mask(bad_info)

    # Z: BETA / SE, otherwise signed from P and the effect direction.
    z, _ = numeric(["Z"])

    with np.errstate(divide="ignore", invalid="ignore"):

        z_from_beta_se = z.isna() & beta.notna() & (se > 0)
        z = z.where(~z_from_beta_se, beta / se)

        z_from_p = z.isna() & beta.notna() & p.notna()
        z = z.where(~z_from_p, np.sign(beta) * z_p)

    out["Z"] = z.replace([np.inf, -np.inf], np.nan)
    info["N_Z_FROM_BETA_SE"] = int(z_from_beta_se.sum())
    info["N_Z_FROM_P"] = int(z_from_p.sum())

    for column in ("BETA", "SE", "P", "EAF", "MAF", "N", "INFO", "Z"):
        info[f"N_MISSING_{column}"] = int(out[column].isna().sum())

    # Effect scale of BETA: log(OR) for binary studies or when derived from OR.
    phenotype_type = trait.get("PHENOTYPE_TYPE", "Unknown")
    info["EFFECT_TYPE"] = (
        "log(OR)" if phenotype_type == "Binary" or (beta_col is None and beta_from_or.any())
        else effect_type(phenotype_type)
    )

    extra = [
        column for column in data.columns
        if str(column).upper() not in FINEMAP_COLUMNS
    ]

    return pd.concat([out, data[extra]], axis=1), info


# =============================================================================
# FIND STEP 02 DOWNLOAD MANIFEST
# =============================================================================

def find_download_manifest(
    phenotype_slug,
    ancestry_slug,
):

    manifest = (
        Path("01_gwas_catalog")
        / phenotype_slug
        / ancestry_slug
        / "GWAS_download_manifest.tsv"
    )

    if not manifest.exists():

        raise FileNotFoundError(
            "\nCannot find the phenotype/ancestry-specific "
            "Step 02 manifest:\n"
            f"{manifest.resolve()}\n\n"
            "Run Step02_Plan_QC_GWAS.py (without --qc) for this "
            "exact phenotype and ancestry, then submit the "
            "generated Step02 downloader."
        )

    return manifest.resolve()


# =============================================================================
# DETECT PHENOTYPE AND ANCESTRY
# =============================================================================

def determine_project(
    manifest_file,
    manifest,
):

    phenotype = None
    ancestry = None

    # -------------------------------------------------------------------------
    # Try manifest metadata first
    # -------------------------------------------------------------------------

    if (
        "PHENOTYPE"
        in manifest.columns
    ):

        values = (
            manifest[
                "PHENOTYPE"
            ]
            .dropna()
            .astype(str)
            .str.strip()
        )

        values = values[
            values != ""
        ]

        if len(values) > 0:

            phenotype = (
                values.iloc[0]
            )

    if (
        "ANCESTRY"
        in manifest.columns
    ):

        values = (
            manifest[
                "ANCESTRY"
            ]
            .dropna()
            .astype(str)
            .str.strip()
        )

        values = values[
            values != ""
        ]

        if len(values) > 0:

            ancestry = (
                values.iloc[0]
            )

    # -------------------------------------------------------------------------
    # Directory fallback
    #
    # Expected structure:
    #
    # 01_gwas_catalog/
    #     parkinson_s_disease/
    #         european/
    #             GWAS_download_manifest.tsv
    # -------------------------------------------------------------------------

    if phenotype is None:

        phenotype = (
            manifest_file
            .parent
            .parent
            .name
        )

    if ancestry is None:

        ancestry = (
            manifest_file
            .parent
            .name
        )

    return (
        phenotype,
        ancestry,
    )


# =============================================================================
# FIND DOWNLOADED GWAS FILE
# =============================================================================

def find_downloaded_file(
    raw_root,
    accession,
    expected_filename,
):

    # -------------------------------------------------------------------------
    # Expected Step 02 path
    # -------------------------------------------------------------------------

    expected = (
        raw_root
        / accession
        / expected_filename
    )

    if (
        expected.exists()
        and expected.is_file()
        and expected.stat().st_size > 0
    ):

        return (
            expected.resolve()
        )

    # -------------------------------------------------------------------------
    # Fallback:
    # Search inside accession directory.
    # -------------------------------------------------------------------------

    accession_directory = (
        raw_root
        / accession
    )

    if not accession_directory.exists():

        return None

    candidates = []

    valid_extensions = (
        ".tsv.gz",
        ".txt.gz",
        ".vcf.gz",
        ".tsv",
        ".txt",
        ".vcf",
    )

    for path in (
        accession_directory.iterdir()
    ):

        if not path.is_file():

            continue

        if (
            path.stat().st_size
            <= 0
        ):

            continue

        name = (
            path.name.lower()
        )

        # Ignore unfinished downloads.
        if name.endswith(
            ".part"
        ):

            continue

        if name.endswith(
            valid_extensions
        ):

            candidates.append(
                path
            )

    # If exactly one plausible GWAS file exists,
    # use it.
    if len(candidates) == 1:

        return (
            candidates[0]
            .resolve()
        )

    return None


# =============================================================================
# COLLECT EXISTING QC SUMMARIES
# =============================================================================

def collect_existing_qc_summaries(
    qc_root,
    destination,
):

    summary_files = list(
        qc_root.glob(
            "*/"
            "*_QC_summary.tsv"
        )
    )

    if not summary_files:

        return

    frames = []

    for path in summary_files:

        try:

            frame = pd.read_csv(
                path,
                sep="\t",
                low_memory=False,
            )

            frames.append(
                frame
            )

        except Exception:

            continue

    if not frames:

        return

    combined = pd.concat(
        frames,
        ignore_index=True,
    )

    if (
        "STUDY_ACCESSION"
        in combined.columns
    ):

        combined = (
            combined
            .drop_duplicates(
                subset=[
                    "STUDY_ACCESSION"
                ],
                keep="last",
            )
        )

    combined.to_csv(
        destination,
        sep="\t",
        index=False,
    )


# =============================================================================
# PLANNER MODE
# =============================================================================

def planner_mode(args):

    banner(
        "STEP 03 - GWASLAB QC PLANNER"
    )

    # -------------------------------------------------------------------------
    # Absolute paths
    # -------------------------------------------------------------------------

    working_directory = (
        Path.cwd()
        .resolve()
    )

    script_path = (
        Path(__file__)
        .resolve()
    )

    python_path = (
        working_directory
        / "envs"
        / "pipeline"
        / "bin"
        / "python"
    ).resolve()

    if not python_path.exists():

        raise FileNotFoundError(
            "\nPipeline Python is missing:\n"
            f"{python_path}\n\n"
            "Run the software setup before Step 03."
        )

    print(
        "Working directory:"
    )

    print(
        f"  {working_directory}"
    )

    print()

    print(
        "Python:"
    )

    print(
        f"  {python_path}"
    )

    print()

    print(
        "Step 03 script:"
    )

    print(
        f"  {script_path}"
    )

    # -------------------------------------------------------------------------
    # Requested project
    # -------------------------------------------------------------------------

    phenotype = args.phenotype.strip()

    if not phenotype:

        raise ValueError(
            "--phenotype cannot be empty."
        )

    ancestry_code, ancestry_label = (
        canonical_ancestry_pair(
            args.ancestry
        )
    )

    phenotype_slug = (
        slugify(
            phenotype
        )
    )

    ancestry_slug = (
        slugify(
            ancestry_label
        )
    )

    # -------------------------------------------------------------------------
    # Locate the EXACT Step 02 manifest
    # -------------------------------------------------------------------------

    download_manifest_file = (
        find_download_manifest(
            phenotype_slug=phenotype_slug,
            ancestry_slug=ancestry_slug,
        )
    )

    print()

    print(
        "Step 02 manifest:"
    )

    print(
        f"  {download_manifest_file}"
    )

    download_manifest = pd.read_csv(
        download_manifest_file,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    required_columns = [
        "STUDY_ACCESSION",
        "FILE_NAME",
    ]

    for column in required_columns:

        if (
            column
            not in download_manifest.columns
        ):

            raise RuntimeError(
                "\nStep 02 manifest is missing "
                f"required column:\n{column}"
            )

    # -------------------------------------------------------------------------
    # Validate Step 01 project metadata
    # -------------------------------------------------------------------------

    if (
        "PHENOTYPE"
        in download_manifest.columns
    ):

        manifest_phenotypes = {
            normalize_text(value)
            for value in (
                download_manifest[
                    "PHENOTYPE"
                ]
                .dropna()
                .astype(str)
            )
            if str(value).strip()
        }

        if (
            manifest_phenotypes
            and manifest_phenotypes
            != {
                normalize_text(
                    phenotype
                )
            }
        ):

            raise RuntimeError(
                "\nManifest phenotype does not match "
                "the requested phenotype.\n"
                f"Requested : {phenotype}\n"
                f"Manifest  : "
                f"{sorted(manifest_phenotypes)}"
            )

    manifest_ancestry_codes = set()

    if (
        "ANCESTRY_CODE"
        in download_manifest.columns
    ):

        manifest_ancestry_codes = {
            str(value).strip().upper()
            for value in (
                download_manifest[
                    "ANCESTRY_CODE"
                ]
                .dropna()
            )
            if str(value).strip()
        }

    elif (
        "ANCESTRY"
        in download_manifest.columns
    ):

        for value in (
            download_manifest[
                "ANCESTRY"
            ]
            .dropna()
            .astype(str)
        ):

            if not value.strip():

                continue

            try:

                code, _ = (
                    canonical_ancestry_pair(
                        value
                    )
                )

                manifest_ancestry_codes.add(
                    code
                )

            except ValueError:

                pass

    if (
        manifest_ancestry_codes
        and manifest_ancestry_codes
        != {
            ancestry_code
        }
    ):

        raise RuntimeError(
            "\nManifest ancestry does not match "
            "the requested ancestry.\n"
            f"Requested : {ancestry_code} "
            f"({ancestry_label})\n"
            f"Manifest  : "
            f"{sorted(manifest_ancestry_codes)}"
        )

    ancestry = ancestry_label

    banner(
        "PROJECT"
    )

    print(
        f"Phenotype : {phenotype}"
    )

    print(
        f"Ancestry  : {ancestry_label}"
    )

    print(
        f"Code      : {ancestry_code}"
    )

    print(
        f"Partition : {args.partition or 'scheduler default'}"
    )

    print(
        f"Walltime  : {args.time}"
    )

    print(
        f"Studies   : "
        f"{len(download_manifest)}"
    )

    # -------------------------------------------------------------------------
    # Paths
    # -------------------------------------------------------------------------

    metadata_directory = (
        download_manifest_file
        .parent
    )

    raw_root = (
        working_directory
        / "02_summary_stats"
        / phenotype_slug
        / ancestry_slug
        / "raw"
    )

    qc_root = (
        working_directory
        / "02_summary_stats"
        / phenotype_slug
        / ancestry_slug
        / "qc"
    )

    log_root = (
        working_directory
        / "logs"
        / "step03_qc"
        / phenotype_slug
        / ancestry_slug
    )

    qc_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    log_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    # -------------------------------------------------------------------------
    # Collect summaries from previous QC runs
    # -------------------------------------------------------------------------

    combined_summary = (
        metadata_directory
        / "GWAS_QC_summary.tsv"
    )

    collect_existing_qc_summaries(
        qc_root=qc_root,
        destination=combined_summary,
    )

    # -------------------------------------------------------------------------
    # Find currently downloaded GWAS
    # -------------------------------------------------------------------------

    banner(
        "CHECKING STEP 02 DOWNLOADS"
    )

    qc_rows = []

    missing_rows = []

    for _, row in (
        download_manifest.iterrows()
    ):

        accession = str(
            row[
                "STUDY_ACCESSION"
            ]
        ).strip()

        expected_filename = str(
            row[
                "FILE_NAME"
            ]
        ).strip()

        input_file = (
            find_downloaded_file(
                raw_root=raw_root,
                accession=accession,
                expected_filename=expected_filename,
            )
        )

        if input_file is None:

            print(
                f"MISSING  {accession}"
            )

            missing_rows.append(
                {
                    "STUDY_ACCESSION":
                        accession,

                    "EXPECTED_FILE":
                        str(
                            raw_root
                            / accession
                            / expected_filename
                        ),
                }
            )

            continue

        print(
            f"READY    {accession}"
        )

        output_directory = (
            qc_root
            / accession
        )

        qc_rows.append(
            build_qc_row(
                task_id=len(qc_rows) + 1,
                accession=accession,
                input_file=input_file,
                output_directory=output_directory,
                phenotype=phenotype,
                ancestry_code=ancestry_code,
                ancestry_label=ancestry_label,
                extra=normalise_trait({
                    field: row.get(field, "")
                    for field in TRAIT_FIELDS + ["PHENOTYPE_TYPE_SOURCE", "TRAIT_TYPE", "TRAIT_TYPE_SOURCE"]
                    if field in row
                }, accession),
            )
        )

    # -------------------------------------------------------------------------
    # Missing files report
    #
    # IMPORTANT:
    # Missing files DO NOT stop Step 03 anymore.
    # -------------------------------------------------------------------------

    if missing_rows:

        missing_report = (
            metadata_directory
            / "GWAS_QC_missing_downloads.tsv"
        )

        pd.DataFrame(
            missing_rows
        ).to_csv(
            missing_report,
            sep="\t",
            index=False,
        )

        banner(
            "SOME DOWNLOADS ARE STILL INCOMPLETE"
        )

        print(
            f"Total expected : "
            f"{len(download_manifest)}"
        )

        print(
            f"Ready for QC   : "
            f"{len(qc_rows)}"
        )

        print(
            f"Still missing  : "
            f"{len(missing_rows)}"
        )

        print()

        print(
            "Missing report:"
        )

        print(
            f"  {missing_report}"
        )

        print()

        print(
            "Continuing with all currently "
            "available GWAS files."
        )

        print()

        print(
            "You can run this planner again later "
            "when more downloads finish."
        )

    # -------------------------------------------------------------------------
    # Must have at least one file
    # -------------------------------------------------------------------------

    if not qc_rows:

        raise RuntimeError(
            "\nNo completed GWAS files are currently "
            "available for QC."
        )

    # -------------------------------------------------------------------------
    # Create Step 03 manifest
    # -------------------------------------------------------------------------

    qc_manifest = pd.DataFrame(
        qc_rows
    )

    qc_manifest_file = (
        metadata_directory
        / "GWAS_QC_manifest.tsv"
    ).resolve()

    qc_manifest.to_csv(
        qc_manifest_file,
        sep="\t",
        index=False,
    )

    number_of_jobs = len(
        qc_manifest
    )

    # -------------------------------------------------------------------------
    # Generate SLURM bash
    # -------------------------------------------------------------------------

    bash_file = (
        working_directory
        / (
            f"Step02_QC_"
            f"{phenotype_slug}_"
            f"{ancestry_slug}.sh"
        )
    )

    job_name = (
        f"QC_"
        f"{phenotype_slug}_"
        f"{ancestry_slug}"
    )

    # There is intentionally NO concurrency cap.
    #
    # Example:
    #
    # #SBATCH --array=1-12
    #
    # NOT:
    #
    # #SBATCH --array=1-12%4

    qc_resources = args.slurm_resources

    # SLURM settings from HPC_CPU.sh (command-line overrides already applied).
    sbatch_header = sbatch_header_block(
        args,
        job_name=job_name,
        output=f"{log_root}/qc.%A_%a.out",
        error=f"{log_root}/qc.%A_%a.err",
        n_tasks=number_of_jobs,
        max_parallel=qc_resources["max_parallel"],
    )

    bash_text = f"""#!/bin/bash
{sbatch_header}

set -euo pipefail

cd "{working_directory}"

export GWAS_QC_MANIFEST="{qc_manifest_file}"

echo "=============================================================="
echo "STEP 03 - GWASLAB QC"
echo "=============================================================="
echo "Job ID        : $SLURM_JOB_ID"
echo "Array task ID : $SLURM_ARRAY_TASK_ID"
echo "Phenotype     : {phenotype}"
echo "Ancestry      : {ancestry_label}"
echo "Ancestry code : {ancestry_code}"
echo "Partition     : {args.partition or 'scheduler default'}"
echo "Walltime      : {args.time}"
echo "Hostname      : $(hostname)"
echo "Python        : {python_path}"
echo "Script        : {script_path}"
echo "=============================================================="

"{python_path}" "{script_path}"

echo
echo "=============================================================="
echo "STEP 03 TASK COMPLETE"
echo "Task: $SLURM_ARRAY_TASK_ID"
echo "=============================================================="
"""

    bash_file.write_text(
        bash_text,
        encoding="utf-8",
    )

    bash_file.chmod(
        0o755
    )

    # -------------------------------------------------------------------------
    # Final report
    # -------------------------------------------------------------------------

    banner(
        "STEP 03 BASH GENERATED"
    )

    print(
        f"Phenotype         : {phenotype}"
    )

    print(
        f"Ancestry          : {ancestry_label}"
    )

    print(
        f"Ancestry code     : {ancestry_code}"
    )

    print(
        f"Partition         : {args.partition or 'scheduler default'}"
    )

    print(
        f"Walltime          : {args.time}"
    )

    print(
        f"Pipeline Python   : {python_path}"
    )

    print(
        f"Total Step 02 GWAS: "
        f"{len(download_manifest)}"
    )

    print(
        f"Available now     : "
        f"{number_of_jobs}"
    )

    print(
        f"Still missing     : "
        f"{len(missing_rows)}"
    )

    print(
        f"SLURM array       : "
        f"1-{number_of_jobs}"
    )

    print(
        f"Concurrency cap   : {qc_resources['max_parallel'] or 'NONE'}"
    )

    print(
        f"Memory per task   : {qc_resources['memory']}"
    )

    print(
        f"CPUs per task     : "
        f"{qc_resources['cpus']}"
    )

    print()

    print(
        "QC manifest:"
    )

    print(
        f"  {qc_manifest_file}"
    )

    print()

    print(
        "Generated bash:"
    )

    print(
        f"  {bash_file}"
    )

    print()

    print(
        "RUN:"
    )

    print()

    print(
        f"  sbatch {bash_file.name}"
    )

    print()

    banner(
        "QC ARRAY TASKS"
    )

    print(
        qc_manifest[
            [
                "QC_TASK_ID",
                "STUDY_ACCESSION",
                "INPUT_FILE",
            ]
        ].to_string(
            index=False
        )
    )


# =============================================================================
# WORKER MODE
# =============================================================================

def worker_mode():

    # -------------------------------------------------------------------------
    # SLURM automatically gives us the task ID.
    # -------------------------------------------------------------------------

    task_value = (
        os.environ.get(
            "SLURM_ARRAY_TASK_ID"
        )
    )

    if task_value is None:

        raise RuntimeError(
            "\nSLURM_ARRAY_TASK_ID is not defined."
        )

    task_id = int(
        task_value
    )

    # -------------------------------------------------------------------------
    # Generated bash automatically provides manifest.
    # -------------------------------------------------------------------------

    manifest_value = (
        os.environ.get(
            "GWAS_QC_MANIFEST"
        )
    )

    if not manifest_value:

        raise RuntimeError(
            "\nGWAS_QC_MANIFEST is not defined."
        )

    manifest_file = (
        Path(
            manifest_value
        )
        .resolve()
    )

    if not manifest_file.exists():

        raise FileNotFoundError(
            "\nQC manifest not found:\n"
            f"{manifest_file}"
        )

    # -------------------------------------------------------------------------
    # Load manifest
    # -------------------------------------------------------------------------

    manifest = pd.read_csv(
        manifest_file,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    required_columns = [
        "QC_TASK_ID",
        "STUDY_ACCESSION",
        "INPUT_FILE",
        "OUTPUT_DIR",
        "PHENOTYPE",
        "ANCESTRY",
        "ANCESTRY_CODE",
        "ANCESTRY_LABEL",
    ]

    for column in required_columns:

        if (
            column
            not in manifest.columns
        ):

            raise RuntimeError(
                "\nQC manifest missing column:\n"
                f"{column}"
            )

    task_numbers = pd.to_numeric(
        manifest[
            "QC_TASK_ID"
        ],
        errors="coerce",
    )

    selected = manifest[
        task_numbers
        == task_id
    ]

    if len(selected) != 1:

        raise RuntimeError(
            "\nExpected exactly one manifest row "
            f"for array task {task_id}. "
            f"Found {len(selected)}."
        )

    row = (
        selected.iloc[0]
    )

    # -------------------------------------------------------------------------
    # Current GWAS
    # -------------------------------------------------------------------------

    accession = str(
        row[
            "STUDY_ACCESSION"
        ]
    ).strip()

    phenotype = str(
        row[
            "PHENOTYPE"
        ]
    ).strip()

    ancestry = str(
        row[
            "ANCESTRY"
        ]
    ).strip()

    ancestry_code = str(
        row[
            "ANCESTRY_CODE"
        ]
    ).strip().upper()

    ancestry_label = str(
        row[
            "ANCESTRY_LABEL"
        ]
    ).strip()

    input_file = (
        Path(
            row[
                "INPUT_FILE"
            ]
        )
        .resolve()
    )

    output_directory = (
        Path(
            row[
                "OUTPUT_DIR"
            ]
        )
        .resolve()
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    # -------------------------------------------------------------------------
    # Output paths
    # -------------------------------------------------------------------------

    full_output = (
        output_directory
        / (
            f"{accession}"
            "_GRCh38_QC.tsv.gz"
        )
    )

    significant_output = (
        output_directory
        / (
            f"{accession}"
            "_GRCh38_significant.tsv.gz"
        )
    )

    summary_tsv = (
        output_directory
        / (
            f"{accession}"
            "_QC_summary.tsv"
        )
    )

    summary_json = (
        output_directory
        / (
            f"{accession}"
            "_QC_summary.json"
        )
    )

    failed_output = (
        output_directory
        / (
            f"{accession}"
            "_FAILED.txt"
        )
    )

    # -------------------------------------------------------------------------
    # Resume support
    #
    # If this GWAS was already processed successfully,
    # do nothing.
    # -------------------------------------------------------------------------

    if (
        full_output.exists()
        and full_output.stat().st_size > 0
        and significant_output.exists()
        and significant_output.stat().st_size > 0
        and summary_tsv.exists()
        and summary_tsv.stat().st_size > 0
    ):

        banner(
            "GWAS ALREADY COMPLETE"
        )

        print(
            f"Accession: {accession}"
        )

        print(
            "Skipping QC."
        )

        return

    # -------------------------------------------------------------------------
    # Verify raw GWAS
    # -------------------------------------------------------------------------

    if not input_file.exists():

        raise FileNotFoundError(
            "\nGWAS input file not found:\n"
            f"{input_file}"
        )

    if (
        input_file.stat().st_size
        <= 0
    ):

        raise RuntimeError(
            "\nGWAS input file is empty:\n"
            f"{input_file}"
        )

    # -------------------------------------------------------------------------
    # Load GWASLab only on compute worker
    # -------------------------------------------------------------------------

    ensure_packages({"scipy": "scipy", "gwaslab": "gwaslab==4.2.3"})

    import gwaslab as gl

    threads = int(
        os.environ.get(
            "SLURM_CPUS_PER_TASK",
            CPUS_PER_TASK,
        )
    )

    banner(
        "STEP 03 - GWASLAB QC WORKER"
    )

    print(
        f"Task       : {task_id}"
    )

    print(
        f"Accession  : {accession}"
    )

    print(
        f"Phenotype  : {phenotype}"
    )

    print(
        f"Ancestry   : {ancestry_label}"
    )

    print(
        f"Code       : {ancestry_code}"
    )

    print(
        f"Input      : {input_file}"
    )

    print(
        f"Output     : {output_directory}"
    )

    print(
        f"Threads    : {threads}"
    )

    try:

        # =====================================================================
        # LOAD SUMMARY STATISTICS
        # =====================================================================

        banner(
            "LOADING GWAS"
        )

        ss = gl.Sumstats(
            str(input_file),
            fmt="auto",
            build="38",
            verbose=True,
        )

        n_before = len(
            ss.data
        )

        # Missing per-variant N -> study N_TOTAL BEFORE basic_check, which
        # would otherwise drop those variants (fine-mapping needs them all).
        trait = normalise_trait({
            field: row.get(field, "")
            for field in TRAIT_FIELDS + ["PHENOTYPE_TYPE_SOURCE", "TRAIT_TYPE", "TRAIT_TYPE_SOURCE"]
            if field in row
        }, accession)

        n_prefilled = 0
        raw_n_col = find_column(ss.data, ["N"])
        study_n = parse_number(trait.get("N_TOTAL"))

        if raw_n_col is not None and study_n:
            values = pd.to_numeric(ss.data[raw_n_col], errors="coerce")
            n_prefilled = int(values.isna().sum())
            if n_prefilled:
                ss.data[raw_n_col] = values.fillna(study_n)
                print(f"N missing for {n_prefilled:,} variants: set to the study N_TOTAL ({study_n:,})")

        print()

        print(
            f"Variants before QC: "
            f"{n_before:,}"
        )

        # =====================================================================
        # GWASLAB STANDARDIZATION / QC
        # =====================================================================

        banner(
            "GWASLAB BASIC CHECK"
        )

        ss.basic_check(
            remove=True,
            remove_dup=True,
            normalize=True,
            threads=threads,
            sanity_check_stats_kwargs={
                "coltocheck": GWASLAB_RANGE_CHECK,
            },
            verbose=True,
        )

        data = (
            ss.data.copy()
        )

        n_after_basic_check = len(
            data
        )

        # =====================================================================
        # DETECT STANDARDIZED COLUMNS
        # =====================================================================

        chr_col = find_column(
            data,
            [
                "CHR",
            ],
        )

        pos_col = find_column(
            data,
            [
                "POS",
            ],
        )

        p_col = find_column(
            data,
            [
                "P",
            ],
        )

        ea_col = find_column(
            data,
            [
                "EA",
            ],
        )

        nea_col = find_column(
            data,
            [
                "NEA",
            ],
        )

        beta_col = find_column(
            data,
            [
                "BETA",
            ],
        )

        se_col = find_column(
            data,
            [
                "SE",
            ],
        )

        z_col = find_column(
            data,
            [
                "Z",
            ],
        )

        eaf_col = find_column(
            data,
            [
                "EAF",
            ],
        )

        n_col = find_column(
            data,
            [
                "N",
            ],
        )

        rsid_col = find_column(
            data,
            [
                "rsID",
                "RSID",
            ],
        )

        snpid_col = find_column(
            data,
            [
                "SNPID",
            ],
        )

        # =====================================================================
        # ESSENTIAL FIELDS
        # =====================================================================

        missing_essential = []

        if chr_col is None:

            missing_essential.append(
                "CHR"
            )

        if pos_col is None:

            missing_essential.append(
                "POS"
            )

        if p_col is None:

            missing_essential.append(
                "P"
            )

        if missing_essential:

            raise RuntimeError(
                "\nMissing essential fields after "
                "GWASLab standardization:\n"
                + ", ".join(
                    missing_essential
                )
            )

        # =====================================================================
        # CLEAN P VALUES
        # =====================================================================

        data[
            p_col
        ] = pd.to_numeric(
            data[
                p_col
            ],
            errors="coerce",
        )

        n_missing_p = int(
            data[
                p_col
            ].isna().sum()
        )

        data = data[
            data[
                p_col
            ].notna()
        ].copy()

        valid_p = (
            (
                data[
                    p_col
                ]
                >= 0
            )
            &
            (
                data[
                    p_col
                ]
                <= 1
            )
        )

        n_invalid_p = int(
            (~valid_p).sum()
        )

        data = data[
            valid_p
        ].copy()

        # =====================================================================
        # STANDARDISED FINE-MAPPING COLUMNS
        # CHR POS SNP EA NEA BETA SE P EAF MAF N INFO Z (all variants kept)
        # =====================================================================

        banner(
            "FINE-MAPPING COLUMNS"
        )

        data, finemap_info = standardize_finemap_columns(
            data,
            trait,
        )

        if n_prefilled:
            finemap_info["N_FILLED_FROM_STUDY_TOTAL"] += n_prefilled
            finemap_info["N_SOURCE"] = "variant+study N_TOTAL"

        for key, value in finemap_info.items():
            print(f"{key:<24}: {value}")

        # The detected columns are now the standardised ones.
        chr_col, pos_col, p_col = "CHR", "POS", "P"
        ea_col = "EA" if data["EA"].notna().any() else None
        nea_col = "NEA" if data["NEA"].notna().any() else None
        beta_col = "BETA" if data["BETA"].notna().any() else None
        se_col = "SE" if data["SE"].notna().any() else None
        z_col = "Z" if data["Z"].notna().any() else None
        eaf_col = "EAF" if data["EAF"].notna().any() else None
        n_col = "N" if data["N"].notna().any() else None

        # =====================================================================
        # ADD PIPELINE METADATA
        # =====================================================================

        data[
            "GWAS_ACCESSION"
        ] = accession

        data[
            "ANALYSIS_PHENOTYPE"
        ] = phenotype

        data[
            "ANALYSIS_ANCESTRY"
        ] = ancestry_label

        data[
            "ANALYSIS_ANCESTRY_CODE"
        ] = ancestry_code

        data[
            "PIPELINE_BUILD"
        ] = "GRCh38"

        # =====================================================================
        # SORT
        # =====================================================================

        data = data.sort_values(
            [
                chr_col,
                pos_col,
            ],
            kind="stable",
        )

        n_after = len(
            data
        )

        # =====================================================================
        # GENOME-WIDE SIGNIFICANT VARIANTS
        # =====================================================================

        significant = data[
            data[
                p_col
            ]
            < P_THRESHOLD
        ].copy()

        significant = (
            significant
            .sort_values(
                p_col,
                ascending=True,
                kind="stable",
            )
        )

        n_significant = len(
            significant
        )

        # =====================================================================
        # DOWNSTREAM READINESS
        # =====================================================================

        has_chr = (
            chr_col is not None
        )

        has_pos = (
            pos_col is not None
        )

        has_p = (
            p_col is not None
        )

        has_ea = (
            ea_col is not None
        )

        has_nea = (
            nea_col is not None
        )

        has_beta = (
            beta_col is not None
        )

        has_se = (
            se_col is not None
        )

        has_z = (
            z_col is not None
        )

        has_eaf = (
            eaf_col is not None
        )

        has_n = (
            n_col is not None
        )

        has_rsid = (
            rsid_col is not None
        )

        has_snpid = (
            snpid_col is not None
        )

        # ---------------------------------------------------------------------
        # Clumping readiness
        # ---------------------------------------------------------------------

        clumping_ready = (
            has_chr
            and has_pos
            and has_p
            and (
                has_rsid
                or has_snpid
                or (
                    has_ea
                    and has_nea
                )
            )
        )

        # ---------------------------------------------------------------------
        # SuSiE readiness
        #
        # We need effect information:
        #
        # BETA + SE
        #
        # OR
        #
        # Z
        # ---------------------------------------------------------------------

        susie_ready = (
            has_chr
            and has_pos
            and has_p
            and has_ea
            and has_nea
            and (
                (
                    has_beta
                    and has_se
                )
                or has_z
            )
        )

        if susie_ready:

            qc_status = "PASS"

        elif clumping_ready:

            qc_status = "WARNING"

        else:

            qc_status = "FAIL"

        # =====================================================================
        # SAVE FULL CLEAN GWAS
        # =====================================================================

        banner(
            "SAVING CLEAN GWAS"
        )

        data.to_csv(
            full_output,
            sep="\t",
            index=False,
            compression="gzip",
        )

        # =====================================================================
        # SAVE SIGNIFICANT VARIANTS
        # =====================================================================

        banner(
            "SAVING GENOME-WIDE SIGNIFICANT VARIANTS"
        )

        significant.to_csv(
            significant_output,
            sep="\t",
            index=False,
            compression="gzip",
        )

        # =====================================================================
        # SUMMARY
        # =====================================================================

        minimum_p = None

        if n_after > 0:

            minimum_p = float(
                data[
                    p_col
                ].min()
            )

        summary = {

            "QC_TASK_ID":
                task_id,

            "STUDY_ACCESSION":
                accession,

            "PHENOTYPE":
                phenotype,

            "ANCESTRY":
                ancestry_label,

            "ANCESTRY_CODE":
                ancestry_code,

            "ANCESTRY_LABEL":
                ancestry_label,

            # Binary vs continuous phenotype (read by fine-mapping / coloc).
            **trait,

            # Derivations behind the standardised fine-mapping columns.
            **finemap_info,

            "BUILD":
                "GRCh38",

            "INPUT_FILE":
                str(input_file),

            "FULL_QC_FILE":
                str(full_output),

            "SIGNIFICANT_FILE":
                str(
                    significant_output
                ),

            "N_BEFORE_QC":
                n_before,

            "N_AFTER_GWASLAB_BASIC_CHECK":
                n_after_basic_check,

            "N_AFTER_QC":
                n_after,

            "N_REMOVED":
                (
                    n_before
                    - n_after
                ),

            "N_MISSING_P_REMOVED":
                n_missing_p,

            "N_INVALID_P_REMOVED":
                n_invalid_p,

            "N_GENOME_WIDE_SIGNIFICANT":
                n_significant,

            "P_THRESHOLD":
                P_THRESHOLD,

            "MIN_P":
                minimum_p,

            "HAS_CHR":
                has_chr,

            "HAS_POS":
                has_pos,

            "HAS_P":
                has_p,

            "HAS_EA":
                has_ea,

            "HAS_NEA":
                has_nea,

            "HAS_BETA":
                has_beta,

            "HAS_SE":
                has_se,

            "HAS_Z":
                has_z,

            "HAS_EAF":
                has_eaf,

            "HAS_N":
                has_n,

            "HAS_RSID":
                has_rsid,

            "HAS_SNPID":
                has_snpid,

            "CLUMPING_READY":
                clumping_ready,

            "SUSIE_READY":
                susie_ready,

            "QC_STATUS":
                qc_status,

            "THREADS":
                threads,

            "QC_COMPLETED_UTC":
                datetime.now(
                    timezone.utc
                ).isoformat(),
        }

        # =====================================================================
        # SAVE SUMMARY TSV
        # =====================================================================

        pd.DataFrame(
            [summary]
        ).to_csv(
            summary_tsv,
            sep="\t",
            index=False,
        )

        # =====================================================================
        # SAVE SUMMARY JSON
        # =====================================================================

        with open(
            summary_json,
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                summary,
                handle,
                indent=2,
                default=str,
            )

        # =====================================================================
        # SAVE STUDY METADATA (fine-mapping)
        # =====================================================================

        metadata = study_metadata_row(
            trait,
            accession,
            ancestry_code,
            ancestry_label,
            effect=finemap_info["EFFECT_TYPE"],
            extra={
                "N_SOURCE": finemap_info["N_SOURCE"],
                "N_VARIANTS": n_after,
                "QC_STATUS": qc_status,
                "QC_FILE": str(full_output),
            },
        )

        pd.DataFrame([metadata]).to_csv(
            output_directory / f"{accession}_study_metadata.tsv",
            sep="\t",
            index=False,
        )

        with open(
            output_directory / f"{accession}_study_metadata.json",
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(metadata, handle, indent=2, default=str)

        # =====================================================================
        # Remove stale failure file
        # =====================================================================

        if failed_output.exists():

            failed_output.unlink()

        # =====================================================================
        # FINAL REPORT
        # =====================================================================

        banner(
            "GWAS QC COMPLETE"
        )

        print(
            f"Accession               : "
            f"{accession}"
        )

        print(
            f"Variants before QC      : "
            f"{n_before:,}"
        )

        print(
            f"Variants after QC       : "
            f"{n_after:,}"
        )

        print(
            f"Variants removed        : "
            f"{n_before - n_after:,}"
        )

        print(
            f"Genome-wide significant : "
            f"{n_significant:,}"
        )

        print(
            f"Minimum P               : "
            f"{minimum_p}"
        )

        print(
            f"Clumping ready          : "
            f"{clumping_ready}"
        )

        print(
            f"SuSiE ready             : "
            f"{susie_ready}"
        )

        print(
            f"QC status               : "
            f"{qc_status}"
        )

        print()

        print(
            "Full cleaned GWAS:"
        )

        print(
            f"  {full_output}"
        )

        print()

        print(
            "Significant variants:"
        )

        print(
            f"  {significant_output}"
        )

        print()

        print(
            "QC summary:"
        )

        print(
            f"  {summary_tsv}"
        )

    # =========================================================================
    # FAILURE HANDLING
    # =========================================================================

    except Exception as error:

        failure_text = (
            "STEP 03 GWAS QC FAILED\n\n"
            f"Task ID: {task_id}\n"
            f"Accession: {accession}\n"
            f"Phenotype: {phenotype}\n"
            f"Ancestry: {ancestry_label} ({ancestry_code})\n"
            f"Input file: {input_file}\n\n"
            f"Error:\n"
            f"{type(error).__name__}: "
            f"{error}\n\n"
            f"{traceback.format_exc()}"
        )

        failed_output.write_text(
            failure_text,
            encoding="utf-8",
        )

        print(
            failure_text,
            file=sys.stderr,
        )

        raise


# =============================================================================
# DIRECT SINGLE-INDEX MODE
# =============================================================================

def direct_index_mode(args):
    """--index N: the whole job for study N of the download manifest (ARRAY_ID N),
    on this machine: download (resumes a .part file; skipped when complete), then
    the same GWASLab QC worker the SLURM array runs. Only the planner (without
    --qc) has to have run before."""

    if args.index < 1:
        raise SystemExit("--index must be >= 1.")

    P = project_paths(args)
    manifest = read_download_manifest(P)
    ids = pd.to_numeric(manifest["ARRAY_ID"], errors="coerce")
    selected = manifest[ids == args.index]

    if selected.empty:
        listing = "\n".join(f"  {int(i):>3}  {a}" for i, a in zip(ids, manifest["STUDY_ACCESSION"]) if pd.notna(i))
        raise SystemExit(f"--index {args.index} is not in {P['download_manifest']}\nAvailable:\n{listing}")

    row = selected.iloc[0]
    accession = str(row["STUDY_ACCESSION"]).strip()
    file_name = str(row["FILE_NAME"]).strip()
    url = str(row.get("DOWNLOAD_URL", "")).strip()

    banner(f"STEP 02 - INDEX {args.index}: {accession} (download + QC on this machine)")
    print(f"Phenotype : {P['phenotype']}")
    print(f"Ancestry  : {P['label']} ({P['code']})")
    print(f"File      : {file_name}")
    print(f"CPUs      : {os.environ.get('SLURM_CPUS_PER_TASK', CPUS_PER_TASK)}")
    if not os.environ.get("SLURM_JOB_ID"):
        print("NOTE: running here, not as a SLURM job. On a login node use e.g.\n"
              f"  srun --cpus-per-task=8 --mem=64G --time=4:00:00 envs/pipeline/bin/python "
              f"{Path(__file__).name} --phenotype \"{P['phenotype']}\" --ancestry {P['code']} --index {args.index}")

    # 1. download -----------------------------------------------------------------
    raw = find_downloaded_file(P["raw"], accession, file_name)

    if raw is None:
        if not url:
            raise SystemExit(f"No DOWNLOAD_URL for {accession} in the download manifest.")
        target = P["raw"] / accession / file_name
        banner("DOWNLOADING")
        print(f"{url}\n-> {target}")
        download_study_file(url, target)
        raw = target
    else:
        print(f"Already downloaded: {raw}")

    # 2. QC (one-row manifest for the normal worker) --------------------------------
    qc_row = build_qc_row(
        task_id=1,
        accession=accession,
        input_file=raw,
        output_directory=P["qc"] / accession,
        phenotype=P["phenotype"],
        ancestry_code=P["code"],
        ancestry_label=P["label"],
        extra=normalise_trait(manifest_trait(row), accession),
    )
    single = P["analysis"] / f".index_{args.index}_qc_manifest.tsv"
    pd.DataFrame([qc_row]).to_csv(single, sep="\t", index=False)

    saved = {k: os.environ.get(k) for k in ("SLURM_ARRAY_TASK_ID", "GWAS_QC_MANIFEST")}
    os.environ["SLURM_ARRAY_TASK_ID"] = "1"
    os.environ["GWAS_QC_MANIFEST"] = str(single)

    try:
        worker_mode()
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        single.unlink(missing_ok=True)

    print()
    print(f"Index {args.index} ({accession}) finished. Progress of all studies:")
    print(f"  envs/pipeline/bin/python {Path(__file__).name} --phenotype \"{P['phenotype']}\" "
          f"--ancestry {P['code']} --status")


# =============================================================================
# --status  (read-only: downloads and QC per study)
# =============================================================================

def project_paths(args):

    ancestry_label = canonical_ancestry(args.ancestry)
    phenotype_slug = slugify(args.phenotype.strip())
    ancestry_slug = slugify(ancestry_label)
    stats_root = PROJECT_ROOT / "02_summary_stats" / phenotype_slug / ancestry_slug
    analysis = PROJECT_ROOT / "01_gwas_catalog" / phenotype_slug / ancestry_slug

    return {
        "phenotype": args.phenotype.strip(),
        "label": ancestry_label,
        "code": ancestry_code(ancestry_label),
        "p_slug": phenotype_slug,
        "a_slug": ancestry_slug,
        "analysis": analysis,
        "raw": stats_root / "raw",
        "qc": stats_root / "qc",
        "download_manifest": analysis / "GWAS_download_manifest.tsv",
        "stats_table": analysis / "GWAS_publication_stats.tsv",
        "stats_cache": analysis / ".GWAS_stats_cache.json",
    }


def qc_outputs(qc_root, accession):

    directory = qc_root / accession

    return {
        "full": directory / f"{accession}_GRCh38_QC.tsv.gz",
        "significant": directory / f"{accession}_GRCh38_significant.tsv.gz",
        "summary": directory / f"{accession}_QC_summary.tsv",
        "failed": directory / f"{accession}_FAILED.txt",
    }


def human_bytes(size):

    size = float(size or 0)

    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024

    return f"{size:.1f} PB"


def study_state(P, row):

    accession = str(row["STUDY_ACCESSION"]).strip()
    raw = find_downloaded_file(P["raw"], accession, str(row["FILE_NAME"]).strip())
    part = P["raw"] / accession / (str(row["FILE_NAME"]).strip() + ".part")

    if raw is not None:
        download, size = "DOWNLOADED", raw.stat().st_size
    elif part.is_file():
        download, size = "PARTIAL", part.stat().st_size
    else:
        download, size = "PENDING", 0

    out = qc_outputs(P["qc"], accession)
    detail = ""

    if all(out[k].is_file() and out[k].stat().st_size > 0 for k in ("full", "significant", "summary")):
        summary = pd.read_csv(out["summary"], sep="\t").iloc[0].to_dict()
        qc = f"DONE ({summary.get('QC_STATUS', '?')})"
    elif out["failed"].is_file():
        qc = "FAILED"
        lines = [x for x in out["failed"].read_text(encoding="utf-8", errors="replace").splitlines() if x.strip()]
        index = next((i for i, x in enumerate(lines) if x.startswith("Error:")), None)
        detail = lines[index + 1][:80] if index is not None and index + 1 < len(lines) else ""
    else:
        qc = "PENDING"

    return {"raw": raw, "download": download, "size": size, "qc": qc, "detail": detail}


def read_download_manifest(P):

    if not P["download_manifest"].is_file():
        raise SystemExit(
            f"No download manifest: {P['download_manifest']}\n"
            "Run the planner first (without --status / --stats)."
        )

    return pd.read_csv(P["download_manifest"], sep="\t", dtype=str, low_memory=False).fillna("")


# =============================================================================
# STUDY METADATA (fine-mapping): GWAS_study_metadata.tsv
# =============================================================================

def study_metadata_row(trait, accession, ancestry_code, ancestry_label, effect=None, extra=None):
    """GCST_ID PHENOTYPE_TYPE N_CASES N_CONTROLS N_TOTAL N_EFF ANCESTRY GENOME_BUILD
    EFFECT_TYPE, then provenance. Both N_TOTAL and N_EFF are kept; each
    fine-mapping method chooses which one it uses."""

    trait = normalise_trait(trait, accession)

    row = {
        "GCST_ID": accession,
        "PHENOTYPE_TYPE": trait["PHENOTYPE_TYPE"],
        "N_CASES": trait.get("N_CASES", ""),
        "N_CONTROLS": trait.get("N_CONTROLS", ""),
        "N_TOTAL": trait.get("N_TOTAL", ""),
        "N_EFF": trait["N_EFF"],
        "ANCESTRY": ancestry_code,
        "GENOME_BUILD": "GRCh38",
        "EFFECT_TYPE": effect or effect_type(trait["PHENOTYPE_TYPE"]),
        "ANCESTRY_LABEL": ancestry_label,
        "CASE_FRACTION": trait.get("CASE_FRACTION", ""),
        "EFFECT_REPORTED": trait.get("EFFECT_REPORTED", ""),
        "CLASSIFICATION_SOURCE": trait.get("CLASSIFICATION_SOURCE", ""),
        "CLASSIFICATION_CONFIDENCE": trait.get("CLASSIFICATION_CONFIDENCE", ""),
        "NEEDS_REVIEW": trait.get("NEEDS_REVIEW", False),
        "CLASSIFICATION_EVIDENCE": trait.get("CLASSIFICATION_EVIDENCE", ""),
        "LD_REFERENCE": f"resources/1000G/{ancestry_code}",
    }

    row.update(extra or {})

    return row


def manifest_trait(row):

    return {
        field: row.get(field, "")
        for field in TRAIT_FIELDS + ["PHENOTYPE_TYPE_SOURCE", "TRAIT_TYPE", "TRAIT_TYPE_SOURCE"]
        if field in row
    }


def review_file(analysis_dir):

    return Path(analysis_dir) / "phenotype_type_review.tsv"


def load_reviews(analysis_dir):
    """VERIFIED_PHENOTYPE_TYPE (Binary / Continuous) filled in by hand."""

    path = review_file(analysis_dir)

    if not path.is_file():
        return {}

    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    return {str(r["STUDY_ACCESSION"]).strip(): r.to_dict() for _, r in table.iterrows()}


def atomic_tsv(table, path):
    """Write via a unique temporary file + rename: parallel array tasks that
    refresh the same table never leave a half-written file."""

    path = Path(path)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    table.to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, path)


def write_review_file(analysis_dir, manifest):
    """Studies whose phenotype type needs manual verification. Verdicts already
    filled in are kept; studies verified by hand stay listed."""

    path = review_file(analysis_dir)
    old = load_reviews(analysis_dir)
    rows = []

    for _, item in manifest.iterrows():

        accession = str(item["STUDY_ACCESSION"]).strip()
        flagged = (
            str(item.get("NEEDS_REVIEW", "")).strip().lower() == "true"
            or phenotype_type_value(item.get("PHENOTYPE_TYPE") or item.get("TRAIT_TYPE")) == "Unknown"
        )

        if not flagged and accession not in old:
            continue

        row = {column: item.get(column, "") for column in REVIEW_COLUMNS}
        row["STUDY_ACCESSION"] = accession
        for column in ("VERIFIED_PHENOTYPE_TYPE", "VERIFIED_N_CASES", "VERIFIED_N_CONTROLS", "REVIEW_NOTE"):
            row[column] = old.get(accession, {}).get(column, "")
        rows.append(row)

    if not rows:
        if path.is_file() and not old:
            path.unlink()
        return []

    atomic_tsv(pd.DataFrame(rows, columns=REVIEW_COLUMNS), path)

    return [r for r in rows if phenotype_type_value(r["VERIFIED_PHENOTYPE_TYPE"]) == "Unknown"]


def write_study_metadata_table(analysis_dir, qc_root, manifest, ancestry_code, ancestry_label):
    """One row per study: the per-study QC metadata when QC has run, otherwise
    the download manifest (EFFECT_TYPE then follows from PHENOTYPE_TYPE)."""

    rows = []

    for _, item in manifest.iterrows():

        accession = str(item["STUDY_ACCESSION"]).strip()
        per_study = Path(qc_root) / accession / f"{accession}_study_metadata.tsv"

        # Phenotype type from the manifest + review verdicts (so a later
        # review also updates finished studies); QC fields from the QC run.
        qc = {"QC_STATUS": "PENDING"}
        effect = None

        if per_study.is_file() and per_study.stat().st_size > 0:
            done = pd.read_csv(per_study, sep="\t", dtype=str).fillna("").iloc[0].to_dict()
            qc = {k: done.get(k, "") for k in ("N_SOURCE", "N_VARIANTS", "QC_STATUS", "QC_FILE")}
            if done.get("EFFECT_TYPE") and done.get("PHENOTYPE_TYPE") == normalise_trait(
                    manifest_trait(item), accession)["PHENOTYPE_TYPE"]:
                effect = done["EFFECT_TYPE"]

        rows.append(study_metadata_row(
            manifest_trait(item),
            accession, ancestry_code, ancestry_label,
            effect=effect,
            extra={"PMID": item.get("PUBMEDID", ""), "DOI": item.get("DOI", ""), **qc},
        ))

    table = pd.DataFrame(rows)

    if table.empty:
        table = pd.DataFrame(columns=STUDY_METADATA_COLUMNS)

    order = STUDY_METADATA_COLUMNS + [c for c in table.columns if c not in STUDY_METADATA_COLUMNS]
    output = Path(analysis_dir) / "GWAS_study_metadata.tsv"
    atomic_tsv(table[order], output)

    return output


def refresh_study_metadata(args):
    """GWAS_study_metadata.tsv from the download manifest + finished QC."""

    P = project_paths(args)

    if not P["download_manifest"].is_file():
        return

    manifest = read_download_manifest(P)
    pending = write_review_file(P["analysis"], manifest)
    REVIEWS.clear()
    REVIEWS.update(load_reviews(P["analysis"]))

    output = write_study_metadata_table(
        P["analysis"], P["qc"], manifest, P["code"], P["label"]
    )
    print(f"Study metadata (fine-mapping): {output}")

    if pending:
        print()
        print(f"{len(pending)} studies need a manual phenotype-type check (Binary or Continuous):")
        for r in pending:
            print(f"  {r['STUDY_ACCESSION']:<14} {str(r['DISEASE_TRAIT'])[:45]:<45}  "
                  f"PMID {r['PUBMEDID'] or '-'}  DOI {r['DOI'] or '-'}")
        print(f"Fill VERIFIED_PHENOTYPE_TYPE (Binary / Continuous) in:\n  {review_file(P['analysis'])}")
        print("then re-run with --qc (unverified studies stay PHENOTYPE_TYPE=Unknown).")


def status_mode(args):

    P = project_paths(args)
    manifest = read_download_manifest(P)

    banner(f"STEP 02 STATUS - {P['phenotype']} / {P['label']} ({P['code']})")

    print(f"{'ID':>3}  {'STUDY':<14}{'TRAIT':<14}{'DOWNLOAD':<12}{'SIZE':>10}  {'QC':<15}DETAIL")

    downloads, qcs = {}, {}

    for _, row in manifest.iterrows():
        s = study_state(P, row)
        downloads[s["download"]] = downloads.get(s["download"], 0) + 1
        key = s["qc"].split(" ")[0]
        qcs[key] = qcs.get(key, 0) + 1
        trait = normalise_trait(manifest_trait(row), row["STUDY_ACCESSION"])
        phenotype_type = trait["PHENOTYPE_TYPE"] + (" *" if trait["NEEDS_REVIEW"] else "")
        print(f"{row.get('ARRAY_ID', ''):>3}  {row['STUDY_ACCESSION']:<14}{phenotype_type:<14}"
              f"{s['download']:<12}{human_bytes(s['size']) if s['size'] else '-':>10}  {s['qc']:<15}{s['detail']}")

    collect_existing_qc_summaries(P["qc"], P["analysis"] / "GWAS_QC_summary.tsv")

    print()
    print("Downloads : " + "  ".join(f"{k}={v}" for k, v in sorted(downloads.items())))
    print("QC        : " + "  ".join(f"{k}={v}" for k, v in sorted(qcs.items())))
    print()

    me = Path(__file__).name
    base = f"envs/pipeline/bin/python {me} --phenotype \"{P['phenotype']}\" --ancestry {P['code']}"

    if downloads.get("PENDING", 0) + downloads.get("PARTIAL", 0):
        print("Downloads still to finish (partial files resume):")
        print(f"  sbatch Step02_Download_GWAS_{P['p_slug']}_{P['a_slug']}.sh")
    if downloads.get("DOWNLOADED", 0) and (qcs.get("PENDING", 0) or qcs.get("FAILED", 0)):
        print("QC for the downloaded GWAS (finished studies are skipped):")
        print(f"  {base} --qc")
        print(f"  sbatch Step02_QC_{P['p_slug']}_{P['a_slug']}.sh")
    if not (downloads.get("PENDING", 0) or downloads.get("PARTIAL", 0) or qcs.get("PENDING", 0) or qcs.get("FAILED", 0)):
        print("Every study is downloaded and QC'd.")
        print(f"Publication table: {base} --stats")


# =============================================================================
# --stats  (publication table)
# =============================================================================

def join_unique(values):

    seen = []

    for value in values:
        for part in re.split(r"[|,]", str(value)):
            part = part.strip()
            if part and part.lower() not in {"nr", "nan", "na"} and part not in seen:
                seen.append(part)

    return ", ".join(seen)


def raw_file_stats(path, cache):
    """Columns + variant count of the downloaded file (counted once, then cached)."""

    key = str(path)
    signature = [path.stat().st_size, int(path.stat().st_mtime)]

    if key in cache and cache[key].get("signature") == signature:
        return cache[key]

    opener = gzip.open if path.name.endswith((".gz", ".bgz")) else open

    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = sum(1 for _ in handle)

    cache[key] = {"signature": signature, "columns": header, "n_variants": rows}

    return cache[key]


def qc_file_stats(path, cache):
    """Genomic inflation (lambda GC), N and chromosomes from the QC'd file (cached)."""

    key = str(path)
    signature = [path.stat().st_size, int(path.stat().st_mtime)]

    if key in cache and cache[key].get("signature") == signature:
        return cache[key]

    from scipy.stats import chi2

    header = pd.read_csv(path, sep="\t", nrows=0).columns
    data = pd.read_csv(path, sep="\t", usecols=[c for c in ("CHR", "P", "N", "EAF") if c in header],
                       low_memory=False)
    p = pd.to_numeric(data["P"], errors="coerce").dropna().clip(lower=1e-300, upper=1)

    result = {
        "signature": signature,
        "lambda_gc": round(float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)), 4) if len(p) else None,
        "n_chromosomes": int(data["CHR"].nunique()) if "CHR" in data else None,
    }

    if "N" in data:
        n = pd.to_numeric(data["N"], errors="coerce").dropna()
        result["n_median"] = int(n.median()) if len(n) else None

    if "EAF" in data:
        eaf = pd.to_numeric(data["EAF"], errors="coerce").dropna()
        result["n_maf_below_1pct"] = int((np.minimum(eaf, 1 - eaf) < 0.01).sum())

    cache[key] = result

    return result


def stats_mode(args):

    P = project_paths(args)
    manifest = read_download_manifest(P)

    session = make_session()
    cache_dir = PROJECT_ROOT / "reference_metadata" / "gwas_catalog"
    download_metadata(session, STUDIES_URL, cache_dir / "gwas_catalog_studies.tsv", args.refresh)
    download_metadata(session, ANCESTRY_URL, cache_dir / "gwas_catalog_ancestry.tsv", args.refresh)
    studies = pd.read_csv(cache_dir / "gwas_catalog_studies.tsv", sep="\t", dtype=str, low_memory=False)
    ancestry = pd.read_csv(cache_dir / "gwas_catalog_ancestry.tsv", sep="\t", dtype=str, low_memory=False)

    try:
        cache = json.loads(P["stats_cache"].read_text(encoding="utf-8")) if P["stats_cache"].is_file() else {}
    except (OSError, ValueError):
        cache = {}

    banner(f"STEP 02 PUBLICATION STATISTICS - {P['phenotype']} / {P['label']} ({P['code']})")
    print("Variant counts and lambda GC are computed once per file and cached.")

    rows = []

    for _, m in manifest.iterrows():

        accession = str(m["STUDY_ACCESSION"]).strip()
        r = {"ID": m.get("ARRAY_ID", ""), "STUDY_ACCESSION": accession}

        catalog = studies[studies["STUDY ACCESSION"].astype(str) == accession]
        if not catalog.empty:
            c = catalog.iloc[0]
            for key, column in (("FIRST_AUTHOR", "FIRST AUTHOR"), ("PUBLICATION_DATE", "DATE"),
                                ("JOURNAL", "JOURNAL"), ("PUBMEDID", "PUBMEDID"), ("TITLE", "STUDY"),
                                ("DISEASE_TRAIT", "DISEASE/TRAIT"), ("MAPPED_TRAIT", "MAPPED_TRAIT"),
                                ("MAPPED_TRAIT_URI", "MAPPED_TRAIT_URI"),
                                ("INITIAL_SAMPLE", "INITIAL SAMPLE SIZE"),
                                ("REPLICATION_SAMPLE", "REPLICATION SAMPLE SIZE"),
                                ("PLATFORM", "PLATFORM [SNPS PASSING QC]")):
                r[key] = "" if column not in c or pd.isna(c[column]) else str(c[column])

        # API v2: cohorts, imputation, SNP count, EFO terms (cached per study).
        api_key = f"api::{accession}"
        if api_key not in cache:
            cache[api_key] = fetch_json(session, f"{GWAS_API_V2}/{accession}") or {}
        api = cache[api_key]
        r.update({
            "COHORTS": ", ".join(api.get("cohort") or []),
            "IMPUTED": api.get("imputed", ""),
            "SNP_COUNT_REPORTED": api.get("snp_count", ""),
            "EFO_TERMS": ", ".join(f"{x.get('efo_trait')} ({x.get('efo_id')})" for x in api.get("efo_traits") or []),
        })

        anc = ancestry[ancestry["STUDY ACCESSION"].astype(str) == accession]
        stage = anc["STAGE"].fillna("").str.lower() if "STAGE" in anc.columns else pd.Series("", index=anc.index)
        initial, replication = anc[stage.str.contains("initial")], anc[stage.str.contains("replication")]
        r.update({
            "POPULATION": join_unique(initial.get("BROAD ANCESTRAL CATEGORY", [])),
            "COUNTRY_OF_RECRUITMENT": join_unique(initial.get("COUNTRY OF RECRUITMENT", [])),
            "N_REPLICATION": sum(parse_number(x) for x in replication.get("NUMBER OF INDIVIDUALS", [])) or "",
        })

        # Binary / continuous from the planner (API + sumstats metadata).
        r.update(normalise_trait(manifest_trait(m), accession))
        r["DOI"] = m.get("DOI", "")

        s = study_state(P, m)
        r.update({"DOWNLOAD_STATUS": s["download"], "QC_STATE": s["qc"], "FILE_NAME": m.get("FILE_NAME", ""),
                  "FILE_SIZE_BYTES": s["size"] or ""})

        out = qc_outputs(P["qc"], accession)
        summary = pd.read_csv(out["summary"], sep="\t").iloc[0].to_dict() if s["qc"].startswith("DONE") else {}

        if s["raw"] is not None:
            raw = raw_file_stats(s["raw"], cache)
            r.update({"N_COLUMNS": len(raw["columns"]), "COLUMNS": ", ".join(raw["columns"]),
                      "N_VARIANTS_RAW": summary.get("N_BEFORE_QC", raw["n_variants"])})

        if summary:
            r.update({
                "N_VARIANTS_QC": summary.get("N_AFTER_QC", ""),
                "N_GENOME_WIDE_SIGNIFICANT": summary.get("N_GENOME_WIDE_SIGNIFICANT", ""),
                "MIN_P": summary.get("MIN_P", ""),
                "QC_STATUS": summary.get("QC_STATUS", ""),
                "EFFECT_TYPE": summary.get("EFFECT_TYPE", r["EFFECT_TYPE"]),
                "N_SOURCE": summary.get("N_SOURCE", ""),
                "N_MISSING_Z": summary.get("N_MISSING_Z", ""),
                "AVAILABLE_FIELDS": ", ".join(k[4:] for k, v in summary.items()
                                              if k.startswith("HAS_") and str(v) == "True"),
            })
            q = qc_file_stats(out["full"], cache)
            r.update({"LAMBDA_GC": q.get("lambda_gc", ""), "N_CHROMOSOMES": q.get("n_chromosomes", ""),
                      "N_MEDIAN_IN_FILE": q.get("n_median", ""),
                      "N_VARIANTS_MAF_BELOW_1PCT": q.get("n_maf_below_1pct", "")})

        rows.append(r)
        print_study(r)

    try:
        P["stats_cache"].write_text(json.dumps(cache, default=str), encoding="utf-8")
    except OSError:
        pass

    table = pd.DataFrame(rows)
    table.to_csv(P["stats_table"], sep="\t", index=False)

    banner("SUMMARY TABLE")
    show = [c for c in ("STUDY_ACCESSION", "FIRST_AUTHOR", "PUBLICATION_DATE", "PUBMEDID", "DOI", "PHENOTYPE_TYPE", "N_TOTAL", "N_CASES",
                        "N_CONTROLS", "N_EFF", "N_VARIANTS_QC", "N_GENOME_WIDE_SIGNIFICANT", "LAMBDA_GC", "QC_STATUS")
            if c in table.columns]
    display = table[show].copy()
    for column in ("N_TOTAL", "N_CASES", "N_CONTROLS", "N_EFF", "N_VARIANTS_QC", "N_GENOME_WIDE_SIGNIFICANT"):
        if column in display.columns:
            display[column] = pd.to_numeric(display[column], errors="coerce").round().astype("Int64")
    print(display.astype(object).where(display.notna(), "").to_string(index=False))
    print(f"\nFull table (all columns): {P['stats_table']}")


def print_study(r):

    def fmt(value):
        if value in ("", None) or (isinstance(value, float) and np.isnan(value)):
            return "-"
        try:
            return f"{int(float(value)):,}"
        except (TypeError, ValueError):
            return str(value)

    print()
    print(f"{r['ID']}. {r['STUDY_ACCESSION']}  {r.get('FIRST_AUTHOR', '')} ({str(r.get('PUBLICATION_DATE', ''))[:4]}) "
          f"{r.get('JOURNAL', '')}  PMID {r.get('PUBMEDID', '')}  DOI {r.get('DOI', '') or '-'}")
    print(f"   Title        : {str(r.get('TITLE', ''))[:110]}")
    print(f"   Trait        : {r.get('DISEASE_TRAIT', '')}  [{r.get('EFO_TERMS', '') or r.get('MAPPED_TRAIT', '')}]")
    print(f"   Phenotype    : {r.get('PHENOTYPE_TYPE') or 'Unknown'} ({r.get('CLASSIFICATION_CONFIDENCE') or '-'})"
          f" | effect {r.get('EFFECT_TYPE') or '-'} (reported: {r.get('EFFECT_REPORTED') or '-'})"
          + ("  ** NEEDS REVIEW **" if r.get("NEEDS_REVIEW") else ""))
    print(f"   Evidence     : {r.get('CLASSIFICATION_EVIDENCE') or '-'} [{r.get('CLASSIFICATION_SOURCE') or '-'}]")
    print(f"   Samples      : N {fmt(r.get('N_TOTAL'))} | cases {fmt(r.get('N_CASES'))} | "
          f"controls {fmt(r.get('N_CONTROLS'))} | case fraction {r.get('CASE_FRACTION') or '-'} | "
          f"N_eff {fmt(r.get('N_EFF'))} | replication N {fmt(r.get('N_REPLICATION'))}")
    print(f"   Population   : {r.get('POPULATION', '') or '-'}; recruited in: {r.get('COUNTRY_OF_RECRUITMENT', '') or '-'}")
    if r.get("COHORTS"):
        print(f"   Cohorts      : {str(r['COHORTS'])[:110]}")
    print(f"   Platform     : {r.get('PLATFORM', '') or '-'} | imputed {r.get('IMPUTED', '-')} | "
          f"SNPs reported {fmt(r.get('SNP_COUNT_REPORTED'))}")
    print(f"   File         : {r.get('FILE_NAME', '')}  {r.get('DOWNLOAD_STATUS', '')}"
          f"{', ' + str(r.get('N_COLUMNS')) + ' columns' if r.get('N_COLUMNS') else ''}")
    if r.get("COLUMNS"):
        print(f"   Columns      : {str(r['COLUMNS'])[:110]}")
    if r.get("N_VARIANTS_RAW") not in (None, ""):
        print(f"   Variants     : raw {fmt(r.get('N_VARIANTS_RAW'))} -> after QC {fmt(r.get('N_VARIANTS_QC'))}"
              f" | genome-wide significant {fmt(r.get('N_GENOME_WIDE_SIGNIFICANT'))} | min P {r.get('MIN_P') or '-'}")
    if r.get("LAMBDA_GC") not in (None, ""):
        print(f"   QC           : {r.get('QC_STATUS', '')} | lambda GC {r.get('LAMBDA_GC')} | "
              f"chromosomes {fmt(r.get('N_CHROMOSOMES'))} | fields: {r.get('AVAILABLE_FIELDS', '')}")
        print(f"   Fine-mapping : N from {r.get('N_SOURCE') or '-'} | variants without Z {fmt(r.get('N_MISSING_Z'))}")
    else:
        print(f"   QC           : {r.get('QC_STATE', '')}")


# =============================================================================
# ARGUMENTS / MAIN
# =============================================================================

def pmid_to_doi_mode(args):

    dois = pmid_to_doi(make_session(), args.pmid_to_doi)

    print("PMID\tDOI")
    for pmid in args.pmid_to_doi:
        print(f"{clean_pmid(pmid)}\t{dois.get(clean_pmid(pmid), '') or 'NOT FOUND'}")


def arguments(argv=None):

    parser = argparse.ArgumentParser(
        description="GWAS discovery + download planner and GWASLab QC for one phenotype x ancestry."
    )

    parser.add_argument("--phenotype", default=None,
                        help="Any phenotype, searched as given in the GWAS Catalog.")
    parser.add_argument("--ancestry", default=None,
                        help="EUR, AFR, EAS, SAS, AMR, MID or European, East Asian, ...")
    parser.add_argument("--alias", action="append", default=[],
                        help="Extra search term (repeatable), e.g. --alias \"Parkinson disease\".")
    parser.add_argument("--max-studies", type=int, default=0,
                        help="Maximum eligible GWAS in the download manifest (0 = all).")
    parser.add_argument("--refresh", action="store_true",
                        help="Re-download the GWAS Catalog metadata.")
    parser.add_argument("--phenotype-type", choices=["auto", "binary", "continuous"], default="auto",
                        help="Binary (cases/controls) or Continuous. auto = decide per study from the "
                             "GWAS Catalog API / sumstats metadata (default).")

    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--qc", action="store_true",
                       help="QC planner: QC manifest + SLURM array for the GWAS downloaded so far.")
    modes.add_argument("--index", type=int, default=None,
                       help="Study N of the download manifest (ARRAY_ID): download it if needed, then QC it, "
                            "on this machine.")
    modes.add_argument("--status", action="store_true",
                       help="READ-ONLY: downloaded / partial / pending and QC done / failed / pending.")
    modes.add_argument("--stats", action="store_true",
                       help="Publication statistics per study (writes GWAS_publication_stats.tsv).")
    modes.add_argument("--pmid-to-doi", nargs="+", metavar="PMID", default=None,
                       help="Only convert PubMed IDs to DOIs (Europe PMC / NCBI) and exit.")

    slurm = parser.add_argument_group("SLURM (default: the #SBATCH header of HPC_CPU.sh)")
    slurm.add_argument("--template", default="HPC_CPU.sh")
    slurm.add_argument("--partition", default=None)
    slurm.add_argument("--qos", default=None)
    slurm.add_argument("--account", default=None)
    slurm.add_argument("--cpus", type=int, default=None)
    slurm.add_argument("--mem", default=None)
    slurm.add_argument("--time", type=validate_slurm_time, default=None)
    slurm.add_argument("--max-parallel", type=int, default=0,
                       help="Maximum simultaneously running array tasks (0 = no throttle).")

    args = parser.parse_args(argv)

    if not args.pmid_to_doi and not (args.phenotype and args.ancestry):
        parser.error("--phenotype and --ancestry are required (except with --pmid-to-doi)")

    return resolve_slurm(args)


def main():

    global ARGS

    # Generated QC SLURM array task (the download array is plain bash).
    if (
        os.environ.get("SLURM_ARRAY_TASK_ID")
        and os.environ.get("GWAS_QC_MANIFEST")
    ):
        worker_mode()
        return

    args = arguments()
    ARGS = args

    if args.pmid_to_doi:
        pmid_to_doi_mode(args)
        return

    REVIEWS.update(load_reviews(project_paths(args)["analysis"]))

    if args.status:
        status_mode(args)
    elif args.stats:
        stats_mode(args)
    elif args.index is not None:
        direct_index_mode(args)
    elif args.qc:
        planner_mode(args)
    else:
        plan(args)

    refresh_study_metadata(args)


if __name__ == "__main__":

    main()
