#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2m STEP 02 - FIND, DOWNLOAD AND QC GWAS SUMMARY STATISTICS
===============================================================================

One script for one phenotype x ancestry (formerly Step01_Plan_GWAS.py +
Step03_QC_One_GWAS.py):

  1. PLAN (default): search the GWAS Catalog for the phenotype, keep studies
     with full summary statistics, strict ancestry match and a harmonised file,
     rank them by sample size and write ONE SLURM array script in which each
     task DOWNLOADS one GWAS and then runs GWASLab QC on it.
  2. --status : what is downloaded / partial / pending, and QC done / failed /
     pending, per study (read-only).
  3. --stats  : publication table per study - authors, journal, year, PMID,
     trait (EFO), population and countries, total N, cases, controls,
     effective N, platform, number of variants (raw / after QC), file columns,
     genome-wide significant hits, minimum P, genomic inflation (lambda GC),
     available fields (read-only; writes GWAS_publication_stats.tsv).
  4. --index N: run task N (download + QC) directly on this machine.

Commands (from the project directory):

    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR
    sbatch Step02_GWAS_migraine_european.sh

    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --status
    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --stats
    envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype migraine --ancestry EUR --index 1

    # options: --alias "migraine disorder" --max-studies 5 --refresh
    #          --partition --qos --account --cpus --mem --time --parallel

SLURM settings come from the #SBATCH header of HPC_CPU.sh (or --template);
workers always run envs/pipeline/bin/python.

Outputs (unchanged, read by Step04 / Step05):

    01_gwas_catalog/<phenotype>/<ancestry>/
        candidate_studies.tsv, study_selection.tsv, discovery_summary.json
        GWAS_download_manifest.tsv, GWAS_QC_manifest.tsv, GWAS_QC_summary.tsv
        GWAS_publication_stats.tsv                       (--stats)
    02_summary_stats/<phenotype>/<ancestry>/
        raw/<GCST>/<harmonised file>
        qc/<GCST>/<GCST>_GRCh38_QC.tsv.gz
                  <GCST>_GRCh38_significant.tsv.gz
                  <GCST>_QC_summary.tsv / .json   (QC_STATUS, CLUMPING_READY, SUSIE_READY)
                  <GCST>_FAILED.txt                (only if QC failed)
===============================================================================
"""

from __future__ import annotations

import argparse
import gzip
import html
import json
import os
import re
import shutil
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd


# =============================================================================
# CONSTANTS
# =============================================================================

ROOT = Path.cwd().resolve()
SCRIPT = Path(__file__).resolve()
PIPELINE_PYTHON = ROOT / "envs" / "pipeline" / "bin" / "python"

STUDIES_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
               "gwas-catalog-download-studies-v1.0.3.1.txt")
ANCESTRY_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
                "gwas-catalog-download-ancestries-v1.0.3.1.txt")
SUMMARY_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
CATALOG_CACHE = ROOT / "reference_metadata" / "gwas_catalog"

P_THRESHOLD = 5e-8
DEFAULT_THREADS = 4
TASK_MANIFEST_ENV = "GWAS2M_STEP02_MANIFEST"

PHENOTYPE_ALIASES = {
    "cad": ["coronary artery disease", "coronary heart disease"],
    "chd": ["coronary heart disease", "coronary artery disease"],
    "t2d": ["type 2 diabetes", "type 2 diabetes mellitus"],
    "mdd": ["major depressive disorder", "major depression", "depression"],
}

ANCESTRY_ALIASES = {
    "eur": "European", "european": "European",
    "afr": "African", "african": "African",
    "eas": "East Asian", "east asian": "East Asian",
    "sas": "South Asian", "south asian": "South Asian",
    "amr": "Hispanic or Latin American", "hispanic": "Hispanic or Latin American",
    "latino": "Hispanic or Latin American", "hispanic or latin american": "Hispanic or Latin American",
    "mid": "Greater Middle Eastern", "middle eastern": "Greater Middle Eastern",
    "greater middle eastern": "Greater Middle Eastern",
}

ANCESTRY_CODES = {
    "European": "EUR", "African": "AFR", "East Asian": "EAS", "South Asian": "SAS",
    "Hispanic or Latin American": "AMR", "Greater Middle Eastern": "MID",
}


# =============================================================================
# SMALL HELPERS
# =============================================================================

def banner(text: str) -> None:
    print()
    print("=" * 80)
    print(text)
    print("=" * 80)


def normalize(value) -> str:
    if value is None:
        return ""
    value = re.sub(r"[^a-z0-9]+", " ", str(value).lower())
    return re.sub(r"\s+", " ", value).strip()


def slugify(value) -> str:
    return normalize(value).replace(" ", "_")


def canonical_ancestry(value) -> tuple[str, str]:
    """(code, label), e.g. ('EUR', 'European')."""
    label = ANCESTRY_ALIASES.get(normalize(value), str(value))
    return ANCESTRY_CODES.get(label, slugify(label).upper()), label


def parse_number(value) -> int:
    numbers = re.findall(r"\d[\d,]*", str(value or ""))
    try:
        return int(numbers[0].replace(",", "")) if numbers else 0
    except ValueError:
        return 0


def estimate_sample_size(text) -> int:
    values = []
    for number in re.findall(r"\d[\d,]*", str(text or "")):
        try:
            x = int(number.replace(",", ""))
        except ValueError:
            continue
        if x >= 100:
            values.append(x)
    return sum(values)


def has_full_summary_stats(value) -> bool:
    value = str(value or "").strip().lower()
    return bool(value) and value not in {
        "na", "n/a", "nan", "none", "no", "false", "0", "not available", "nr"}


def human_bytes(size) -> str:
    size = float(size or 0)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"


def find_column(frame: pd.DataFrame, candidates) -> str | None:
    lookup = {str(c).lower(): c for c in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


# =============================================================================
# PROJECT PATHS (same layout Step04 / Step05 read)
# =============================================================================

def project(phenotype: str, ancestry: str) -> dict:
    phenotype = phenotype.strip()
    if not phenotype:
        raise SystemExit("--phenotype cannot be empty.")
    code, label = canonical_ancestry(ancestry)
    p, a = slugify(phenotype), slugify(label)
    analysis = ROOT / "01_gwas_catalog" / p / a
    stats_root = ROOT / "02_summary_stats" / p / a
    return {
        "phenotype": phenotype, "code": code, "label": label, "p_slug": p, "a_slug": a,
        "analysis": analysis,
        "raw": stats_root / "raw",
        "qc": stats_root / "qc",
        "logs": ROOT / "logs" / "step02_gwas" / p / a,
        "download_manifest": analysis / "GWAS_download_manifest.tsv",
        "task_manifest": analysis / "GWAS_QC_manifest.tsv",
        "qc_summary": analysis / "GWAS_QC_summary.tsv",
        "stats_table": analysis / "GWAS_publication_stats.tsv",
        "stats_cache": analysis / ".GWAS_stats_cache.json",
        "script": ROOT / f"Step02_GWAS_{p}_{a}.sh",
    }


def qc_paths(output_dir: Path, accession: str) -> dict:
    return {
        "full": output_dir / f"{accession}_GRCh38_QC.tsv.gz",
        "significant": output_dir / f"{accession}_GRCh38_significant.tsv.gz",
        "summary_tsv": output_dir / f"{accession}_QC_summary.tsv",
        "summary_json": output_dir / f"{accession}_QC_summary.json",
        "failed": output_dir / f"{accession}_FAILED.txt",
    }


def qc_complete(paths: dict) -> bool:
    return all(paths[k].is_file() and paths[k].stat().st_size > 0
               for k in ("full", "significant", "summary_tsv"))


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


# =============================================================================
# GWAS CATALOG DISCOVERY (logic of the former Step01_Plan_GWAS.py)
# =============================================================================

def make_session():
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    session = requests.Session()
    retries = Retry(total=8, connect=8, read=8, status=8, backoff_factor=1.5,
                    status_forcelist=[429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.headers.update({"User-Agent": "GWAS2m-Step02/1.0"})
    return session


def download_metadata(session, url: str, path: Path, refresh: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0 and not refresh:
        print(f"Using cached metadata: {path.relative_to(ROOT)}")
        return
    print(f"Downloading:\n  {url}")
    response = session.get(url, timeout=(30, 300))
    response.raise_for_status()
    temporary = Path(str(path) + ".part")
    temporary.write_bytes(response.content)
    os.replace(temporary, path)


def load_catalog(refresh: bool = False, session=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    studies_file = CATALOG_CACHE / "gwas_catalog_studies.tsv"
    ancestry_file = CATALOG_CACHE / "gwas_catalog_ancestry.tsv"
    if refresh or not studies_file.is_file() or not ancestry_file.is_file():
        session = session or make_session()
        download_metadata(session, STUDIES_URL, studies_file, refresh)
        download_metadata(session, ANCESTRY_URL, ancestry_file, refresh)
    return read_tsv(studies_file), read_tsv(ancestry_file)


def build_terms(phenotype: str, user_aliases) -> list[str]:
    terms = [phenotype] + PHENOTYPE_ALIASES.get(normalize(phenotype), []) + list(user_aliases)
    final, seen = [], set()
    for term in terms:
        n = normalize(term)
        if n and n not in seen:
            seen.add(n)
            final.append(term)
    return final


def contains_term(series: pd.Series, term: str) -> pd.Series:
    term = normalize(term)
    if not term:
        return pd.Series(False, index=series.index)
    if " " not in term and len(term) <= 5:  # short terms: whole-word match only
        pattern = rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])"
    else:
        pattern = re.escape(term)
    return series.str.contains(pattern, regex=True, na=False)


def find_phenotype_studies(studies: pd.DataFrame, terms) -> pd.DataFrame:
    empty = pd.Series("", index=studies.index)
    disease = studies["DISEASE/TRAIT"].fillna("").map(normalize)
    mapped = studies["MAPPED_TRAIT"].fillna("").map(normalize) if "MAPPED_TRAIT" in studies else empty
    title = studies["STUDY"].fillna("").map(normalize) if "STUDY" in studies else empty
    scores = pd.Series(0, index=studies.index, dtype=int)
    for term in terms:
        t = normalize(term)
        for mask, score in ((disease == t, 100), (mapped == t, 95), (contains_term(disease, term), 80),
                            (contains_term(mapped, term), 70), (contains_term(title, term), 40)):
            scores.loc[mask] = np.maximum(scores.loc[mask], score)
    result = studies[scores > 0].copy()
    result["MATCH_SCORE"] = scores.loc[result.index]
    return result.drop_duplicates(subset=["STUDY ACCESSION"])


def get_ancestry_information(ancestry: pd.DataFrame, target_label: str) -> dict:
    """STRICT ancestry: every discovery (initial-stage) category must equal the request."""
    target = normalize(target_label)
    if "STAGE" in ancestry.columns:
        initial = ancestry[ancestry["STAGE"].fillna("").astype(str).str.lower().str.contains("initial", regex=False)]
        if not initial.empty:
            ancestry = initial
    records = {}
    for accession, group in ancestry.groupby("STUDY ACCESSION"):
        categories = {normalize(x) for x in group["BROAD ANCESTRAL CATEGORY"].dropna() if normalize(x)}
        total_n = sum(parse_number(x) for x in group["NUMBER OF INDIVIDUALS"]) \
            if "NUMBER OF INDIVIDUALS" in group.columns else 0
        records[str(accession)] = {
            "STRICT_MATCH": categories == {target},
            "ANCESTRY_CATEGORIES": " | ".join(sorted(categories)),
            "ANCESTRY_N": total_n,
        }
    return records


def gwas_bucket(accession: str) -> str:
    numeric = accession.replace("GCST", "")
    number, width = int(numeric), max(len(numeric), 8)
    start = ((number - 1) // 1000) * 1000 + 1
    return f"GCST{start:0{width}d}-GCST{start + 999:0{width}d}"


def find_harmonised_file(session, accession: str) -> dict | None:
    import requests

    directory = f"{SUMMARY_ROOT}/{gwas_bucket(accession)}/{accession}/harmonised/"
    try:
        response = session.get(directory, timeout=(20, 60))
    except requests.RequestException:
        return None
    if response.status_code != 200:
        return None
    files = []
    for href in re.findall(r'href=["\']([^"\']+)["\']', response.text, flags=re.I):
        name = html.unescape(href).split("/")[-1]
        lower = name.lower()
        if not name or lower.endswith((".tbi", ".log", ".yaml", ".yml", ".md5", ".json")):
            continue
        if ".h.tsv" in lower or lower.endswith((".tsv.gz", ".tsv")):
            files.append(name)
    if not files:
        return None

    def priority(name: str) -> int:
        lower = name.lower()
        for rank, suffix in enumerate((".h.tsv.gz", ".h.tsv", ".tsv.gz", ".tsv")):
            if lower.endswith(suffix):
                return rank
        return 10

    filename = sorted(files, key=priority)[0]
    return {"URL": urljoin(directory, filename), "FILENAME": filename, "DIRECTORY": directory}


def new_selection_table(candidates, phenotype, ancestry_label, ancestry_info) -> pd.DataFrame:
    """Audit table: every phenotype-matched study and why it was (not) selected."""
    rows = []
    for _, row in candidates.iterrows():
        accession = str(row["STUDY ACCESSION"])
        info = ancestry_info.get(accession, {})
        rows.append({
            "STUDY_ACCESSION": accession, "REQUESTED_PHENOTYPE": phenotype,
            "DISEASE_TRAIT": row.get("DISEASE/TRAIT", ""), "MAPPED_TRAIT": row.get("MAPPED_TRAIT", ""),
            "MATCH_SCORE": int(row["MATCH_SCORE"]), "REQUESTED_ANCESTRY": ancestry_label,
            "ANCESTRY_CATEGORIES": info.get("ANCESTRY_CATEGORIES", ""),
            "ANCESTRY_MATCH": bool(info.get("STRICT_MATCH", False)),
            "INITIAL_SAMPLE_SIZE": row.get("INITIAL SAMPLE SIZE", ""),
            "SUMMARY_STATS_AVAILABLE": has_full_summary_stats(row.get("FULL SUMMARY STATISTICS")),
            "HARMONISED_FILE_AVAILABLE": "NOT_CHECKED", "DOWNLOAD_URL": "",
            "ELIGIBLE": False, "EXCLUSION_REASON": "", "RANK": "",
        })
    selection = pd.DataFrame(rows, columns=[
        "STUDY_ACCESSION", "REQUESTED_PHENOTYPE", "DISEASE_TRAIT", "MAPPED_TRAIT", "MATCH_SCORE",
        "REQUESTED_ANCESTRY", "ANCESTRY_CATEGORIES", "ANCESTRY_MATCH", "INITIAL_SAMPLE_SIZE",
        "SUMMARY_STATS_AVAILABLE", "HARMONISED_FILE_AVAILABLE", "DOWNLOAD_URL", "ELIGIBLE",
        "EXCLUSION_REASON", "RANK"])
    selection.index = selection["STUDY_ACCESSION"].tolist()
    # Reasons follow the filter order: full summary statistics -> strict ancestry -> harmonised file.
    selection.loc[~selection["SUMMARY_STATS_AVAILABLE"].astype(bool), "EXCLUSION_REASON"] = \
        "NO_FULL_SUMMARY_STATISTICS"
    no_ancestry = (selection["EXCLUSION_REASON"] == "") & ~selection["ANCESTRY_MATCH"].astype(bool)
    selection.loc[no_ancestry & (selection["ANCESTRY_CATEGORIES"] == ""), "EXCLUSION_REASON"] = \
        "ANCESTRY_NOT_REPORTED"
    selection.loc[no_ancestry & (selection["ANCESTRY_CATEGORIES"] != ""), "EXCLUSION_REASON"] = \
        "ANCESTRY_MISMATCH"
    return selection


def write_study_selection(state: dict) -> None:
    analysis = state.get("analysis_dir")
    if analysis is None:
        return
    selection = state.get("selection")
    if selection is None:
        selection = pd.DataFrame(columns=["STUDY_ACCESSION", "REQUESTED_PHENOTYPE", "ELIGIBLE", "EXCLUSION_REASON"])
    selection = selection.copy()
    unchecked = (selection["EXCLUSION_REASON"] == "") & ~selection["ELIGIBLE"].astype(bool)
    selection.loc[unchecked, "EXCLUSION_REASON"] = "NOT_CHECKED_MAX_STUDIES_LIMIT"
    selection.to_csv(analysis / "study_selection.tsv", sep="\t", index=False)
    summary = {k: v for k, v in state.items() if k not in {"analysis_dir", "selection"}}
    (analysis / "discovery_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")


def discover(args, P: dict, state: dict) -> pd.DataFrame:
    """Eligible studies (ranked), as the former Step01 planner selected them."""
    terms = build_terms(P["phenotype"], args.alias)
    P["analysis"].mkdir(parents=True, exist_ok=True)
    state.update({"analysis_dir": P["analysis"], "REQUESTED_PHENOTYPE": P["phenotype"],
                  "REQUESTED_ANCESTRY": P["label"], "ANCESTRY_CODE": P["code"],
                  "SEARCH_TERMS": terms, "MAX_STUDIES": args.max_studies})
    print("Search terms: " + ", ".join(terms))

    banner("GWAS CATALOG METADATA")
    session = make_session()
    studies, ancestry = load_catalog(args.refresh, session)
    state.update({"GWAS_CATALOG_STUDIES_URL": STUDIES_URL, "GWAS_CATALOG_ANCESTRY_URL": ANCESTRY_URL,
                  "N_CATALOG_STUDIES": int(len(studies))})
    ancestry_info = get_ancestry_information(ancestry, P["label"])

    banner("PHENOTYPE FILTER")
    candidates = find_phenotype_studies(studies, terms)
    print(f"Phenotype matches              : {len(candidates):,}")
    state["selection"] = new_selection_table(candidates, P["phenotype"], P["label"], ancestry_info)
    state["N_PHENOTYPE_MATCHED"] = int(len(candidates))
    if candidates.empty:
        raise RuntimeError(f"No GWAS studies found for phenotype: {P['phenotype']}")
    if "FULL SUMMARY STATISTICS" not in candidates.columns:
        raise RuntimeError("GWAS Catalog studies file does not contain 'FULL SUMMARY STATISTICS'.")

    candidates = candidates[candidates["FULL SUMMARY STATISTICS"].map(has_full_summary_stats)].copy()
    print(f"With full summary statistics   : {len(candidates):,}")
    state["N_WITH_FULL_SUMMARY_STATISTICS"] = int(len(candidates))
    if candidates.empty:
        raise RuntimeError("No phenotype-matching studies have full summary statistics.")

    strict = {acc for acc, info in ancestry_info.items() if info["STRICT_MATCH"]}
    candidates = candidates[candidates["STUDY ACCESSION"].astype(str).isin(strict)].copy()
    print(f"Strict {P['label']} ancestry      : {len(candidates):,}")
    state["N_FULL_SUMMARY_STATISTICS_AND_ANCESTRY_MATCHED"] = int(len(candidates))
    if candidates.empty:
        raise RuntimeError(f"No strict {P['label']} GWAS studies found.")

    accessions = candidates["STUDY ACCESSION"].astype(str)
    candidates["ANCESTRY_N"] = [int(ancestry_info.get(a, {}).get("ANCESTRY_N", 0)) for a in accessions]
    candidates["ANCESTRY_CATEGORIES"] = [ancestry_info.get(a, {}).get("ANCESTRY_CATEGORIES", "") for a in accessions]
    candidates["ESTIMATED_INITIAL_N"] = candidates["INITIAL SAMPLE SIZE"].map(estimate_sample_size) \
        if "INITIAL SAMPLE SIZE" in candidates.columns else 0
    candidates["RANKING_N"] = candidates[["ANCESTRY_N", "ESTIMATED_INITIAL_N"]].max(axis=1)
    candidates = candidates.sort_values(by=["MATCH_SCORE", "RANKING_N"], ascending=[False, False])
    for rank, accession in enumerate(candidates["STUDY ACCESSION"].astype(str), start=1):
        state["selection"].at[accession, "RANK"] = rank
    candidate_file = P["analysis"] / "candidate_studies.tsv"
    candidates.to_csv(candidate_file, sep="\t", index=False)
    state["CANDIDATE_FILE"] = str(candidate_file)

    banner("CHECKING HARMONISED GWAS FILES")
    eligible = []
    for _, row in candidates.iterrows():
        accession = str(row["STUDY ACCESSION"])
        print(f"{accession}: ", end="", flush=True)
        result = find_harmonised_file(session, accession)
        if result is None:
            print("no harmonised file")
            state["selection"].at[accession, "HARMONISED_FILE_AVAILABLE"] = False
            state["selection"].at[accession, "EXCLUSION_REASON"] = "NO_HARMONISED_FILE"
            continue
        print(result["FILENAME"])
        state["selection"].at[accession, "HARMONISED_FILE_AVAILABLE"] = True
        state["selection"].at[accession, "DOWNLOAD_URL"] = result["URL"]
        state["selection"].at[accession, "ELIGIBLE"] = True
        eligible.append((row, result))
        if args.max_studies > 0 and len(eligible) >= args.max_studies:
            break
    if not eligible:
        raise RuntimeError("No matching harmonised GWAS files were found.")

    manifest = pd.DataFrame([{
        "ARRAY_ID": i, "STUDY_ACCESSION": str(row["STUDY ACCESSION"]),
        "DOWNLOAD_URL": download["URL"], "FILE_NAME": download["FILENAME"],
        "PHENOTYPE": P["phenotype"], "ANCESTRY": P["label"], "ANCESTRY_CODE": P["code"],
        "ANCESTRY_LABEL": P["label"], "RANKING_N": int(row["RANKING_N"]),
        "MATCH_SCORE": int(row["MATCH_SCORE"]), "DISEASE_TRAIT": row.get("DISEASE/TRAIT", ""),
        "MAPPED_TRAIT": row.get("MAPPED_TRAIT", ""), "FIRST_AUTHOR": row.get("FIRST AUTHOR", ""),
        "PUBMEDID": row.get("PUBMEDID", ""), "INITIAL_SAMPLE_SIZE": row.get("INITIAL SAMPLE SIZE", ""),
    } for i, (row, download) in enumerate(eligible, start=1)])
    manifest.to_csv(P["download_manifest"], sep="\t", index=False)
    state["N_ELIGIBLE"] = int(len(manifest))
    state["DOWNLOAD_MANIFEST"] = str(P["download_manifest"])
    return manifest


# =============================================================================
# SLURM SCRIPT (settings from HPC_CPU.sh; one task = download + QC of one GWAS)
# =============================================================================

PER_JOB_KEYS = {"job-name", "output", "error", "array"}


def read_template(path: Path) -> tuple[str, dict]:
    shebang, values = "#!/bin/bash --login", {}
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line.startswith("#!"):
                shebang = line
                continue
            match = re.match(r"^#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
            if match and (match.group(2) or "").strip():
                values[match.group(1)] = match.group(2).strip()
    return shebang, values


def write_slurm_script(args, P: dict, n_tasks: int) -> dict:
    shebang, template = read_template(ROOT / args.template)
    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": str(DEFAULT_THREADS), "mem": "64G", "time": "24:00:00"}
    values.update({k: v for k, v in template.items() if k not in PER_JOB_KEYS})
    for key, value in (("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time),
                       ("partition", args.partition), ("qos", args.qos), ("account", args.account)):
        if value not in (None, ""):
            values[key] = str(value)
    if n_tasks > 1000:
        raise SystemExit(f"{n_tasks} studies is above the 1000-job SLURM limit; use --max-studies.")
    throttle = f"%{args.parallel}" if args.parallel and args.parallel > 0 else ""
    P["logs"].mkdir(parents=True, exist_ok=True)
    directives = {**values, "job-name": f"GWAS_{P['p_slug']}_{P['a_slug']}"[:100],
                  "array": f"1-{n_tasks}{throttle}",
                  "output": f"{P['logs']}/step02.%A_%a.out", "error": f"{P['logs']}/step02.%A_%a.err"}
    order = ["job-name", "account", "partition", "qos", "nodes", "ntasks", "cpus-per-task", "mem", "time",
             "array", "output", "error"]
    lines = [f"#SBATCH --{k}={directives[k]}" for k in order if directives.get(k)]
    lines += [f"#SBATCH --{k}={v}" for k, v in directives.items() if k not in order and v]

    text = f"""{shebang}
{chr(10).join(lines)}

# GWAS2m STEP 02 - one array task = download + GWASLab QC of one GWAS
# Generated by Step02_Plan_QC_GWAS.py - re-run the planner instead of editing.
# Phenotype: {P['phenotype']}   Ancestry: {P['label']} ({P['code']})   Studies: {n_tasks}
#
#   sbatch {P['script'].name}
#   envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype "{P['phenotype']}" --ancestry {P['code']} --status
#   envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype "{P['phenotype']}" --ancestry {P['code']} --stats

set -euo pipefail
cd "{ROOT}"

export {TASK_MANIFEST_ENV}="{P['task_manifest']}"
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-{DEFAULT_THREADS}}}"
export OPENBLAS_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-{DEFAULT_THREADS}}}"

echo "Job ${{SLURM_ARRAY_JOB_ID:-}}_${{SLURM_ARRAY_TASK_ID:-}} on $(hostname) at $(date)"
"{PIPELINE_PYTHON}" "{SCRIPT}"
"""
    P["script"].write_text(text, encoding="utf-8", newline="\n")
    P["script"].chmod(0o755)
    return values


# =============================================================================
# PLAN
# =============================================================================

def plan(args, P: dict) -> None:
    banner("STEP 02 - FIND, DOWNLOAD AND QC GWAS")
    print(f"Phenotype : {P['phenotype']}")
    print(f"Ancestry  : {P['label']} ({P['code']})")
    if not PIPELINE_PYTHON.exists():
        raise SystemExit(f"Pipeline Python is missing: {PIPELINE_PYTHON}\nRun: python Step00_Installation.py")

    state: dict = {}
    try:
        manifest = discover(args, P, state)
    finally:
        write_study_selection(state)

    tasks = pd.DataFrame([{
        "QC_TASK_ID": int(row["ARRAY_ID"]), "STUDY_ACCESSION": row["STUDY_ACCESSION"],
        "DOWNLOAD_URL": row["DOWNLOAD_URL"], "FILE_NAME": row["FILE_NAME"],
        "INPUT_FILE": str(P["raw"] / row["STUDY_ACCESSION"] / row["FILE_NAME"]),
        "OUTPUT_DIR": str(P["qc"] / row["STUDY_ACCESSION"]),
        "PHENOTYPE": P["phenotype"], "ANCESTRY": P["label"], "ANCESTRY_CODE": P["code"],
        "ANCESTRY_LABEL": P["label"],
    } for _, row in manifest.iterrows()])
    tasks.to_csv(P["task_manifest"], sep="\t", index=False)
    values = write_slurm_script(args, P, len(tasks))

    banner("PLANNING COMPLETE")
    print(manifest[["ARRAY_ID", "STUDY_ACCESSION", "RANKING_N", "DISEASE_TRAIT", "FILE_NAME"]].to_string(index=False))
    print(f"\nStudies (array tasks) : {len(tasks)}")
    print(f"Task manifest         : {P['task_manifest'].relative_to(ROOT)}")
    print(f"SLURM per task        : " + "  ".join(f"--{k}={values[k]}" for k in
                                                   ("account", "partition", "qos", "cpus-per-task", "mem", "time")
                                                   if values.get(k)))
    print("\nNEXT:")
    print(f"  sbatch {P['script'].name}")
    print(f"  envs/pipeline/bin/python {SCRIPT.name} --phenotype \"{P['phenotype']}\" --ancestry {P['code']} --status")


# =============================================================================
# WORKER: DOWNLOAD + QC ONE GWAS (QC logic of the former Step03_QC_One_GWAS.py)
# =============================================================================

def tool(name: str) -> str | None:
    candidate = ROOT / "envs" / "pipeline" / "bin" / name
    return str(candidate) if candidate.exists() else shutil.which(name)


def download_gwas(url: str, target: Path) -> str:
    """Resumable wget/curl download into <target>.part, renamed when complete."""
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file() and target.stat().st_size > 0:
        return "SKIPPED_ALREADY_COMPLETE"
    temporary = Path(str(target) + ".part")
    wget, curl = tool("wget"), tool("curl")
    if wget:
        command = [wget, "--continue", "--tries=10", "--timeout=60", "-O", str(temporary), url]
    elif curl:
        command = [curl, "--location", "--fail", "--retry", "10", "--retry-delay", "5",
                   "--continue-at", "-", "--output", str(temporary), url]
    else:
        raise RuntimeError("Neither wget nor curl is available.")
    subprocess.run(command, check=True)
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f"Download produced no data: {url}")
    os.replace(temporary, target)
    return "COMPLETE"


def run_task(task_id: int, manifest_file: Path) -> None:
    manifest = read_tsv(manifest_file)
    selected = manifest[pd.to_numeric(manifest["QC_TASK_ID"], errors="coerce") == task_id]
    if len(selected) != 1:
        raise RuntimeError(f"Expected one manifest row for task {task_id}; found {len(selected)}.")
    row = selected.iloc[0]
    accession = str(row["STUDY_ACCESSION"]).strip()
    phenotype = str(row["PHENOTYPE"]).strip()
    ancestry_code = str(row["ANCESTRY_CODE"]).strip().upper()
    ancestry_label = str(row["ANCESTRY_LABEL"]).strip()
    input_file = Path(row["INPUT_FILE"])
    output_dir = Path(row["OUTPUT_DIR"])
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = qc_paths(output_dir, accession)

    banner(f"STEP 02 TASK {task_id}: {accession}")

    # ---------------- download ----------------
    print(f"Download   : {row['DOWNLOAD_URL']}")
    print(f"Raw file   : {input_file}")
    status = download_gwas(str(row["DOWNLOAD_URL"]), input_file)
    print(f"Download   : {status} ({human_bytes(input_file.stat().st_size)})")

    # ---------------- QC (resume: skip if already complete) ----------------
    if qc_complete(paths):
        print("QC already complete - skipping.")
        return
    qc_one_gwas(task_id, accession, phenotype, ancestry_code, ancestry_label, input_file, paths)


def qc_one_gwas(task_id, accession, phenotype, ancestry_code, ancestry_label, input_file: Path, paths: dict) -> None:
    import gwaslab as gl

    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", DEFAULT_THREADS))
    print(f"QC threads : {threads}")
    try:
        banner("LOADING GWAS (GWASLab)")
        ss = gl.Sumstats(str(input_file), fmt="auto", build="38", verbose=True)
        n_before = len(ss.data)
        print(f"\nVariants before QC: {n_before:,}")

        banner("GWASLAB BASIC CHECK")
        ss.basic_check(remove=True, remove_dup=True, normalize=True, threads=threads, verbose=True)
        data = ss.data.copy()
        n_after_basic_check = len(data)

        col = {name: find_column(data, options) for name, options in {
            "CHR": ["CHR"], "POS": ["POS"], "P": ["P"], "EA": ["EA"], "NEA": ["NEA"], "BETA": ["BETA"],
            "SE": ["SE"], "Z": ["Z"], "EAF": ["EAF"], "N": ["N"], "RSID": ["rsID", "RSID"], "SNPID": ["SNPID"],
        }.items()}
        missing = [k for k in ("CHR", "POS", "P") if col[k] is None]
        if missing:
            raise RuntimeError("Missing essential fields after GWASLab standardization: " + ", ".join(missing))

        p = col["P"]
        data[p] = pd.to_numeric(data[p], errors="coerce")
        n_missing_p = int(data[p].isna().sum())
        data = data[data[p].notna()].copy()
        valid_p = (data[p] >= 0) & (data[p] <= 1)
        n_invalid_p = int((~valid_p).sum())
        data = data[valid_p].copy()

        data["GWAS_ACCESSION"] = accession
        data["ANALYSIS_PHENOTYPE"] = phenotype
        data["ANALYSIS_ANCESTRY"] = ancestry_label
        data["ANALYSIS_ANCESTRY_CODE"] = ancestry_code
        data["PIPELINE_BUILD"] = "GRCh38"
        data = data.sort_values([col["CHR"], col["POS"]], kind="stable")
        n_after = len(data)

        significant = data[data[p] < P_THRESHOLD].sort_values(p, ascending=True, kind="stable")
        has = {k: v is not None for k, v in col.items()}
        clumping_ready = has["CHR"] and has["POS"] and has["P"] and (
            has["RSID"] or has["SNPID"] or (has["EA"] and has["NEA"]))
        susie_ready = has["CHR"] and has["POS"] and has["P"] and has["EA"] and has["NEA"] and (
            (has["BETA"] and has["SE"]) or has["Z"])
        qc_status = "PASS" if susie_ready else ("WARNING" if clumping_ready else "FAIL")

        banner("SAVING")
        data.to_csv(paths["full"], sep="\t", index=False, compression="gzip")
        significant.to_csv(paths["significant"], sep="\t", index=False, compression="gzip")

        summary = {
            "QC_TASK_ID": task_id, "STUDY_ACCESSION": accession, "PHENOTYPE": phenotype,
            "ANCESTRY": ancestry_label, "ANCESTRY_CODE": ancestry_code, "ANCESTRY_LABEL": ancestry_label,
            "BUILD": "GRCh38", "INPUT_FILE": str(input_file), "FULL_QC_FILE": str(paths["full"]),
            "SIGNIFICANT_FILE": str(paths["significant"]), "N_BEFORE_QC": n_before,
            "N_AFTER_GWASLAB_BASIC_CHECK": n_after_basic_check, "N_AFTER_QC": n_after,
            "N_REMOVED": n_before - n_after, "N_MISSING_P_REMOVED": n_missing_p,
            "N_INVALID_P_REMOVED": n_invalid_p, "N_GENOME_WIDE_SIGNIFICANT": len(significant),
            "P_THRESHOLD": P_THRESHOLD, "MIN_P": float(data[p].min()) if n_after else None,
            **{f"HAS_{k}": v for k, v in has.items()},
            "CLUMPING_READY": clumping_ready, "SUSIE_READY": susie_ready, "QC_STATUS": qc_status,
            "THREADS": threads, "QC_COMPLETED_UTC": datetime.now(timezone.utc).isoformat(),
        }
        pd.DataFrame([summary]).to_csv(paths["summary_tsv"], sep="\t", index=False)
        paths["summary_json"].write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
        if paths["failed"].exists():
            paths["failed"].unlink()

        banner("GWAS QC COMPLETE")
        print(f"Accession               : {accession}")
        print(f"Variants before / after : {n_before:,} / {n_after:,}")
        print(f"Genome-wide significant : {len(significant):,}")
        print(f"Clumping / SuSiE ready  : {clumping_ready} / {susie_ready}")
        print(f"QC status               : {qc_status}")
    except Exception as error:
        text = (f"STEP 02 GWAS QC FAILED\n\nTask ID: {task_id}\nAccession: {accession}\n"
                f"Phenotype: {phenotype}\nAncestry: {ancestry_label} ({ancestry_code})\n"
                f"Input file: {input_file}\n\nError:\n{type(error).__name__}: {error}\n\n"
                f"{traceback.format_exc()}")
        paths["failed"].write_text(text, encoding="utf-8")
        print(text, file=sys.stderr)
        raise


def direct_index(args, P: dict) -> None:
    if args.index < 1:
        raise SystemExit("--index must be >= 1.")
    if not P["task_manifest"].is_file():
        raise SystemExit(f"No task manifest: {P['task_manifest']}\nRun the planner first (without --index).")
    print("NOTE: --index runs on this machine; on an HPC login node prefer srun / sbatch.")
    run_task(args.index, P["task_manifest"])


# =============================================================================
# --status
# =============================================================================

def collect_qc_summaries(P: dict) -> pd.DataFrame:
    frames = []
    for path in sorted(P["qc"].glob("*/*_QC_summary.tsv")):
        try:
            frames.append(pd.read_csv(path, sep="\t", low_memory=False))
        except Exception:
            continue
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    if "STUDY_ACCESSION" in combined.columns:
        combined = combined.drop_duplicates(subset=["STUDY_ACCESSION"], keep="last")
    combined.to_csv(P["qc_summary"], sep="\t", index=False)
    return combined


def task_rows(P: dict) -> pd.DataFrame:
    if P["task_manifest"].is_file():
        return read_tsv(P["task_manifest"])
    if P["download_manifest"].is_file():  # older layout: download manifest only
        frame = read_tsv(P["download_manifest"])
        frame["INPUT_FILE"] = [str(P["raw"] / a / f) for a, f in zip(frame["STUDY_ACCESSION"], frame["FILE_NAME"])]
        frame["OUTPUT_DIR"] = [str(P["qc"] / a) for a in frame["STUDY_ACCESSION"]]
        frame["QC_TASK_ID"] = frame.get("ARRAY_ID", range(1, len(frame) + 1))
        return frame
    raise SystemExit(f"No plan for this phenotype/ancestry ({P['analysis'].relative_to(ROOT)}).\n"
                     "Run the planner first (without --status / --stats).")


def study_state(row) -> dict:
    accession = str(row["STUDY_ACCESSION"])
    raw = Path(row["INPUT_FILE"])
    part = Path(str(raw) + ".part")
    paths = qc_paths(Path(row["OUTPUT_DIR"]), accession)
    if raw.is_file() and raw.stat().st_size > 0:
        download, size = "DOWNLOADED", raw.stat().st_size
    elif part.is_file():
        download, size = "PARTIAL", part.stat().st_size
    else:
        download, size = "PENDING", 0
    detail = ""
    if qc_complete(paths):
        summary = pd.read_csv(paths["summary_tsv"], sep="\t").iloc[0].to_dict()
        qc = f"DONE ({summary.get('QC_STATUS', '?')})"
    elif paths["failed"].is_file():
        qc = "FAILED"
        lines = [x for x in paths["failed"].read_text(encoding="utf-8", errors="replace").splitlines() if x.strip()]
        error_index = next((i for i, x in enumerate(lines) if x.startswith("Error:")), None)
        detail = lines[error_index + 1][:90] if error_index is not None and error_index + 1 < len(lines) else ""
    else:
        qc = "PENDING"
    return {"download": download, "size": size, "qc": qc, "detail": detail}


def show_status(P: dict) -> None:
    tasks = task_rows(P)
    banner(f"STEP 02 STATUS - {P['phenotype']} / {P['label']} ({P['code']})")
    print(f"{'TASK':>4}  {'STUDY':<14}{'DOWNLOAD':<12}{'SIZE':>10}  {'QC':<16}DETAIL")
    counts = {"download": {}, "qc": {}}
    for _, row in tasks.iterrows():
        s = study_state(row)
        counts["download"][s["download"]] = counts["download"].get(s["download"], 0) + 1
        qc_key = s["qc"].split(" ")[0]
        counts["qc"][qc_key] = counts["qc"].get(qc_key, 0) + 1
        print(f"{row['QC_TASK_ID']:>4}  {row['STUDY_ACCESSION']:<14}{s['download']:<12}"
              f"{human_bytes(s['size']) if s['size'] else '-':>10}  {s['qc']:<16}{s['detail']}")
    collect_qc_summaries(P)
    print("\nDownloads : " + "  ".join(f"{k}={v}" for k, v in sorted(counts["download"].items())))
    print("QC        : " + "  ".join(f"{k}={v}" for k, v in sorted(counts["qc"].items())))
    pending = counts["qc"].get("PENDING", 0) + counts["qc"].get("FAILED", 0)
    if pending:
        print(f"\n{pending} task(s) still to run. Submit (finished tasks are skipped):")
        print(f"  sbatch {P['script'].name}")
    else:
        print("\nEvery study is downloaded and QC'd.")
        print(f"Publication table: envs/pipeline/bin/python {SCRIPT.name} "
              f"--phenotype \"{P['phenotype']}\" --ancestry {P['code']} --stats")


# =============================================================================
# --stats  (publication table)
# =============================================================================

def parse_cases_controls(text) -> tuple[int, int]:
    """Cases / controls from a GWAS Catalog sample description, e.g.
    '48,975 European ancestry cases, 540,381 European ancestry controls'."""
    text = str(text or "")
    cases = sum(int(x.replace(",", "")) for x in re.findall(r"(\d[\d,]*)\s+[^,;\d]*?\bcases?\b", text, flags=re.I))
    controls = sum(int(x.replace(",", "")) for x in
                   re.findall(r"(\d[\d,]*)\s+[^,;\d]*?\bcontrols?\b", text, flags=re.I))
    return cases, controls


def effective_n(cases: int, controls: int) -> int | None:
    return int(round(4 / (1 / cases + 1 / controls))) if cases and controls else None


def join_unique(values) -> str:
    seen = []
    for value in values:
        for part in re.split(r"[|,]", str(value)):
            part = part.strip()
            if part and part.lower() not in {"nr", "nan", "na"} and part not in seen:
                seen.append(part)
    return ", ".join(seen)


def catalog_record(accession: str, studies: pd.DataFrame, ancestry: pd.DataFrame) -> dict:
    out = {}
    row = studies[studies["STUDY ACCESSION"].astype(str) == accession]
    if not row.empty:
        r = row.iloc[0]
        for key, column in (("FIRST_AUTHOR", "FIRST AUTHOR"), ("PUBLICATION_DATE", "DATE"),
                            ("JOURNAL", "JOURNAL"), ("PUBMEDID", "PUBMEDID"), ("TITLE", "STUDY"),
                            ("DISEASE_TRAIT", "DISEASE/TRAIT"), ("MAPPED_TRAIT", "MAPPED_TRAIT"),
                            ("MAPPED_TRAIT_URI", "MAPPED_TRAIT_URI"), ("INITIAL_SAMPLE", "INITIAL SAMPLE SIZE"),
                            ("REPLICATION_SAMPLE", "REPLICATION SAMPLE SIZE"),
                            ("PLATFORM", "PLATFORM [SNPS PASSING QC]"), ("GENOTYPING", "GENOTYPING TECHNOLOGY"),
                            ("ASSOCIATION_COUNT", "ASSOCIATION COUNT"), ("COHORT", "COHORT")):
            out[key] = "" if column not in r or pd.isna(r[column]) else str(r[column])
    cases, controls = parse_cases_controls(out.get("INITIAL_SAMPLE", ""))
    anc = ancestry[ancestry["STUDY ACCESSION"].astype(str) == accession]
    stage = anc["STAGE"].fillna("").str.lower() if "STAGE" in anc.columns else pd.Series("", index=anc.index)
    initial, replication = anc[stage.str.contains("initial")], anc[stage.str.contains("replication")]

    def number(frame, column):
        return sum(parse_number(x) for x in frame[column]) if column in frame.columns else 0

    # Cases / controls: GWAS Catalog numeric columns when present, else parsed from the sample text.
    cases = number(initial, "NUMBER OF CASES") or cases
    controls = number(initial, "NUMBER OF CONTROLS") or controls
    out.update({
        "POPULATION": join_unique(initial.get("BROAD ANCESTRAL CATEGORY", [])),
        "COUNTRY_OF_RECRUITMENT": join_unique(initial.get("COUNTRY OF RECRUITMENT", [])),
        "COUNTRY_OF_ORIGIN": join_unique(initial.get("COUNTRY OF ORIGIN", [])),
        "ANCESTRY_DESCRIPTION": join_unique(initial.get("ADDITIONAL ANCESTRY DESCRIPTION", [])),
        "N_DISCOVERY": number(initial, "NUMBER OF INDIVIDUALS") or estimate_sample_size(out.get("INITIAL_SAMPLE")),
        "N_REPLICATION": number(replication, "NUMBER OF INDIVIDUALS") or "",
        "N_CASES": cases or "", "N_CONTROLS": controls or "",
        "EFFECTIVE_N": effective_n(cases, controls) or "",
        "TRAIT_TYPE": "case-control" if cases and controls else "quantitative / not reported",
    })
    return out


def raw_file_stats(path: Path, cache: dict) -> dict:
    """Header + variant count of the downloaded file (counted once, then cached)."""
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


def qc_file_stats(path: Path, cache: dict) -> dict:
    """Genomic inflation (lambda GC), N and chromosome count from the QC'd file (cached)."""
    key = str(path)
    signature = [path.stat().st_size, int(path.stat().st_mtime)]
    if key in cache and cache[key].get("signature") == signature:
        return cache[key]
    from scipy.stats import chi2

    header = pd.read_csv(path, sep="\t", nrows=0).columns
    use = [c for c in ("CHR", "P", "N", "EAF") if c in header]
    data = pd.read_csv(path, sep="\t", usecols=use, low_memory=False)
    p = pd.to_numeric(data["P"], errors="coerce").dropna().clip(lower=1e-300, upper=1)
    result = {"signature": signature,
              "lambda_gc": round(float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)), 4) if len(p) else None,
              "n_chromosomes": int(data["CHR"].nunique()) if "CHR" in data else None}
    if "N" in data:
        n = pd.to_numeric(data["N"], errors="coerce").dropna()
        result.update({"n_median": int(n.median()) if len(n) else None, "n_max": int(n.max()) if len(n) else None})
    if "EAF" in data:
        eaf = pd.to_numeric(data["EAF"], errors="coerce").dropna()
        maf = np.minimum(eaf, 1 - eaf)
        result["n_maf_below_1pct"] = int((maf < 0.01).sum())
    cache[key] = result
    return result


def show_stats(args, P: dict) -> None:
    tasks = task_rows(P)
    studies, ancestry = load_catalog(args.refresh)
    try:
        cache = json.loads(P["stats_cache"].read_text(encoding="utf-8")) if P["stats_cache"].is_file() else {}
    except (OSError, ValueError):
        cache = {}

    banner(f"STEP 02 PUBLICATION STATISTICS - {P['phenotype']} / {P['label']} ({P['code']})")
    print("Counting variants and computing lambda GC for downloaded / QC'd files "
          "(first time can take ~1 min per GWAS; cached afterwards).")
    rows = []
    for _, task in tasks.iterrows():
        accession = str(task["STUDY_ACCESSION"])
        record = {"TASK": task["QC_TASK_ID"], "STUDY_ACCESSION": accession,
                  **catalog_record(accession, studies, ancestry)}
        state = study_state(task)
        record.update({"DOWNLOAD_STATUS": state["download"], "QC_STATE": state["qc"],
                       "FILE_NAME": Path(task["INPUT_FILE"]).name, "FILE_SIZE_BYTES": state["size"] or ""})
        raw = Path(task["INPUT_FILE"])
        paths = qc_paths(Path(task["OUTPUT_DIR"]), accession)
        summary = pd.read_csv(paths["summary_tsv"], sep="\t").iloc[0].to_dict() if qc_complete(paths) else {}
        if state["download"] == "DOWNLOADED":
            raw_stats = raw_file_stats(raw, cache)
            record.update({"N_COLUMNS": len(raw_stats["columns"]), "COLUMNS": ", ".join(raw_stats["columns"]),
                           "N_VARIANTS_RAW": summary.get("N_BEFORE_QC", raw_stats["n_variants"])})
        if summary:
            fields = [k[4:] for k, v in summary.items() if k.startswith("HAS_") and str(v) == "True"]
            record.update({
                "N_VARIANTS_QC": summary.get("N_AFTER_QC", ""),
                "N_GENOME_WIDE_SIGNIFICANT": summary.get("N_GENOME_WIDE_SIGNIFICANT", ""),
                "MIN_P": summary.get("MIN_P", ""), "QC_STATUS": summary.get("QC_STATUS", ""),
                "AVAILABLE_FIELDS": ", ".join(fields), "BUILD": summary.get("BUILD", "GRCh38"),
            })
            qc = qc_file_stats(paths["full"], cache)
            record.update({"LAMBDA_GC": qc.get("lambda_gc", ""), "N_CHROMOSOMES": qc.get("n_chromosomes", ""),
                           "N_MEDIAN_SUMSTATS": qc.get("n_median", ""), "N_MAX_SUMSTATS": qc.get("n_max", ""),
                           "N_VARIANTS_MAF_BELOW_1PCT": qc.get("n_maf_below_1pct", "")})
        rows.append(record)
        print_study(record)

    try:
        P["stats_cache"].write_text(json.dumps(cache), encoding="utf-8")
    except OSError:
        pass
    table = pd.DataFrame(rows)
    table.to_csv(P["stats_table"], sep="\t", index=False)

    banner("SUMMARY TABLE")
    show = [c for c in ("STUDY_ACCESSION", "FIRST_AUTHOR", "PUBLICATION_DATE", "N_DISCOVERY", "N_CASES",
                        "N_CONTROLS", "N_VARIANTS_QC", "N_GENOME_WIDE_SIGNIFICANT", "LAMBDA_GC", "QC_STATUS")
            if c in table.columns]
    display = table[show].copy()
    for column in ("N_DISCOVERY", "N_CASES", "N_CONTROLS", "N_VARIANTS_QC", "N_GENOME_WIDE_SIGNIFICANT"):
        if column in display.columns:  # whole numbers, not 20000.0
            display[column] = pd.to_numeric(display[column], errors="coerce").round().astype("Int64")
    print(display.astype(object).where(display.notna(), "").to_string(index=False))
    print(f"\nFull table (all columns): {P['stats_table'].relative_to(ROOT)}")


def print_study(r: dict) -> None:
    def fmt(value):
        if value in ("", None) or (isinstance(value, float) and np.isnan(value)):
            return "-"
        try:
            return f"{int(float(value)):,}"
        except (TypeError, ValueError):
            return str(value)

    year = str(r.get("PUBLICATION_DATE", ""))[:4]
    print(f"\n{r['TASK']}. {r['STUDY_ACCESSION']}  {r.get('FIRST_AUTHOR', '')} ({year}) "
          f"{r.get('JOURNAL', '')}  PMID {r.get('PUBMEDID', '')}")
    print(f"   Title        : {str(r.get('TITLE', ''))[:110]}")
    print(f"   Trait        : {r.get('DISEASE_TRAIT', '')}  [mapped: {r.get('MAPPED_TRAIT', '')}"
          f"{' - ' + str(r.get('MAPPED_TRAIT_URI', '')).rsplit('/', 1)[-1] if r.get('MAPPED_TRAIT_URI') else ''}]")
    print(f"   Population   : {r.get('POPULATION', '') or '-'}; recruited in: {r.get('COUNTRY_OF_RECRUITMENT', '') or '-'}")
    print(f"   Samples      : discovery N {fmt(r.get('N_DISCOVERY'))} | cases {fmt(r.get('N_CASES'))} | "
          f"controls {fmt(r.get('N_CONTROLS'))} | effective N {fmt(r.get('EFFECTIVE_N'))} | "
          f"replication N {fmt(r.get('N_REPLICATION'))}  ({r.get('TRAIT_TYPE', '')})")
    print(f"   Description  : {str(r.get('INITIAL_SAMPLE', ''))[:110]}")
    print(f"   Platform     : {r.get('PLATFORM', '') or '-'}  ({r.get('GENOTYPING', '') or '-'})")
    print(f"   File         : {r.get('FILE_NAME', '')}  {r.get('DOWNLOAD_STATUS', '')}"
          f"{', ' + str(r.get('N_COLUMNS')) + ' columns' if r.get('N_COLUMNS') else ''}")
    if r.get("COLUMNS"):
        print(f"   Columns      : {str(r['COLUMNS'])[:110]}")
    if r.get("N_VARIANTS_RAW") not in (None, ""):
        print(f"   Variants     : raw {fmt(r.get('N_VARIANTS_RAW'))} -> after QC {fmt(r.get('N_VARIANTS_QC'))}"
              f" | genome-wide significant {fmt(r.get('N_GENOME_WIDE_SIGNIFICANT'))}"
              f" | min P {r.get('MIN_P', '-') or '-'}")
    if r.get("LAMBDA_GC") not in (None, ""):
        print(f"   QC           : {r.get('QC_STATUS', '')} | lambda GC {r.get('LAMBDA_GC')} | "
              f"chromosomes {fmt(r.get('N_CHROMOSOMES'))} | median N in file {fmt(r.get('N_MEDIAN_SUMSTATS'))} | "
              f"fields: {r.get('AVAILABLE_FIELDS', '')}")
    else:
        print(f"   QC           : {r.get('QC_STATE', '')}")


# =============================================================================
# ARGUMENTS / MAIN
# =============================================================================

def arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Find, download and QC GWAS summary statistics for one phenotype x ancestry.")
    parser.add_argument("--phenotype", required=True, help="e.g. migraine, CAD, \"parkinson's disease\"")
    parser.add_argument("--ancestry", required=True, help="EUR, AFR, EAS, SAS, AMR, MID or a label (European, ...)")
    parser.add_argument("--alias", action="append", default=[], help="Extra search term (repeatable).")
    parser.add_argument("--max-studies", type=int, default=0, help="Keep at most N eligible studies (0 = all).")
    parser.add_argument("--refresh", action="store_true", help="Re-download the GWAS Catalog metadata.")

    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--status", action="store_true", help="Downloaded / pending and QC state per study.")
    modes.add_argument("--stats", action="store_true", help="Publication statistics per study.")
    modes.add_argument("--index", type=int, default=None, help="Run task N (download + QC) on this machine.")

    slurm = parser.add_argument_group("SLURM (default: the #SBATCH header of HPC_CPU.sh)")
    slurm.add_argument("--template", default="HPC_CPU.sh")
    slurm.add_argument("--partition")
    slurm.add_argument("--qos")
    slurm.add_argument("--account")
    slurm.add_argument("--cpus", type=int)
    slurm.add_argument("--mem")
    slurm.add_argument("--time")
    slurm.add_argument("--parallel", type=int, default=0, help="Max simultaneously running tasks (0 = no cap).")
    args = parser.parse_args(argv)
    if args.max_studies < 0:
        parser.error("--max-studies must be >= 0")
    return args


def main() -> int:
    # Generated SLURM array task.
    if os.environ.get(TASK_MANIFEST_ENV) and os.environ.get("SLURM_ARRAY_TASK_ID"):
        run_task(int(os.environ["SLURM_ARRAY_TASK_ID"]), Path(os.environ[TASK_MANIFEST_ENV]))
        return 0

    args = arguments()
    P = project(args.phenotype, args.ancestry)
    if args.status:
        show_status(P)
    elif args.stats:
        show_stats(args, P)
    elif args.index is not None:
        direct_index(args, P)
    else:
        plan(args, P)
    return 0


if __name__ == "__main__":
    sys.exit(main())
