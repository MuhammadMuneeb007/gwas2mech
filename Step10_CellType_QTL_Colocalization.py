#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 10 - Single-cell / cell-type eQTL colocalisation (provider SuSiE, coloc.bf_bf), per locus
===============================================================================

    # what is needed and what is ready (tools, QTL resources, Step06 input); exit 1 if not ready
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources

    # plan: ONE job script, one SLURM array task per fine-mapped locus (max 1000 tasks;
    # more loci than that -> several loci per task), then submit it
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR
    sbatch Step10_CellTypeColoc_parkinson_s_disease_european.sh

    # study 2 (the Step02 number) on this machine: all its loci / one locus
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --index 2
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --index 2 --locus L001

    # progress and results
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --stats

    # redo finished loci (e.g. after changing Step10_Config.py); plan again without --force afterwards
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --force

    # ---- 1000 Genomes 30x LD panel (Steps 04-06 with --reference-panel 1000G30x) ----
    # reads 06_finemapping/<p>/<a>_1000G30x/, writes 10_celltype_coloc/<p>/<a>_1000G30x/
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --check-resources
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x
    sbatch Step10_CellTypeColoc_parkinson_s_disease_european_1000G30x.sh
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --index 2 --locus L001
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --status
    envs/pipeline/bin/python Step10_CellType_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

Method (per locus, nothing is downloaded and no QTL LD / QTL SuSiE is computed):
  1. GWAS side   the Step06 SuSiE-RSS fit of the locus (06_finemapping/<p>/<a>[_<PANEL>]/<GCST>/loci/<L>/);
                 only its credible-set components are used. No credible set -> NO_GWAS_CREDIBLE_SET.
  2. Discovery   dense QTL summary statistics in the locus window: every molecular trait with
                 P <= CANDIDATE_P (top MAX_CANDIDATES_PER_CONTEXT_LOCUS per context).
  3. CS filter   the provider (eQTL Catalogue) SuSiE credible sets of the trait; only components that
                 overlap the GWAS locus are kept (else EXCLUDED_*).
  4. QTL side    provider SuSiE log Bayes factors of those components, read for the locus window from
                 the position-sorted, tabix-indexed LBF file (prepared once by Step00).
  5. coloc::coloc.bf_bf(GWAS LBF, QTL LBF) for every signal pair; the best pair (max PP.H4) classifies
     the trait: STRONG (PP.H4 >= STRONG_PP4), SUGGESTIVE (>= SUGGESTIVE_PP4); PRIOR_ROBUST_STRONG when
     still strong at p12 = P12_LOW. Prior sensitivity at P12_LOW / P12 / P12_HIGH.

Resources (Step00_Download_Resources.py; the LBF files are position-sorted + tabix-indexed there):
    --groups eqtl_catalogue                  tabix_ftp_paths*.tsv (dataset list)
    --groups qtl_celltype --celltype-study-regex '^(OneK1K|Jerber_2021|Nathan_2022|Perez_2022|Randolph_2021|Cytoimmgen)$'
                                             dense + .tbi and provider credible sets + indexed LBF of the
                                             167 single-cell contexts (OneK1K, Jerber_2021, Nathan_2022,
                                             Perez_2022, Randolph_2021, Cytoimmgen)
  A context is used only when dense + .tbi, provider credible sets and the indexed LBF are complete;
  others are reported as not ready - never as "no colocalisation".

Parameters: Step10_Config.py (QTL types, studies, priors, thresholds, SLURM); --qtl-types / --study-regex /
--candidate-p / --strong-pp4 / --suggestive-pp4 override it. SLURM: HPC_CPU.sh header, then SLURM in
Step10_Config.py, then --partition --qos --account --cpus --mem --time --max-parallel.

Outputs (10_celltype_coloc/<phenotype>/<ancestry>[_<PANEL>]/)
    <GCST>/loci/<L>/locus_coloc_results.tsv   every candidate trait: status, PP.H0-H4, class, best pair
                    coloc_signal_pairs.tsv     every GWAS x QTL signal pair (all three p12 priors)
                    coloc_variant_pp_h4.tsv.gz per-variant PP.H4 (traits with PP.H4 >= SAVE_VARIANT_PP_MIN_PP4)
                    discovery_audit.tsv, coloc.log, provider_coloc_done.json
    coloc_manifest.tsv, qtl_contexts_snapshot.tsv (contexts frozen when the job script was made)
    --stats: coloc_stats.tsv, coloc_all_candidates.tsv, coloc_locus_qtl_summary.tsv,
             <GCST>/<GCST>_coloc_results.tsv
===============================================================================
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import runpy
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


STEP = {
    "NAME": "Step10",
    "TITLE": "single-cell / cell-type eQTL (eQTL Catalogue single-cell studies)",
    "CONFIG": "Step10_Config.py",
    "OUT_DIR": "10_celltype_coloc",
    "LOG_DIR": "step10_celltype_coloc",
    "SCRIPT_PREFIX": "Step10_CellTypeColoc",
    "JOB": "SCCOLOC",
    "DEFAULT_QTL_TYPES": ["eQTL"],
    "DEFAULT_STUDY_REGEX": r"^(OneK1K|Jerber_2021|Nathan_2022|Perez_2022|Randolph_2021|Cytoimmgen)$",
    "STEP00_GROUP": "qtl_celltype",
}

VERSION = "2.0.0-provider-susie-per-locus"
PROJECT_ROOT = Path.cwd().resolve()
MAX_ARRAY_TASKS = 1000
DEFAULT_PANEL = "1000G"
EQTL_ROOT = PROJECT_ROOT / "resources" / "coloc" / "eqtl_catalogue"
CATALOGUE_TABLES = ("tabix_ftp_paths.tsv", "tabix_ftp_paths_imported.tsv")
DONE_FILE = "provider_coloc_done.json"
FAILED_FILE = "LOCUS_FAILED.txt"
RESULTS_FILE = "locus_coloc_results.tsv"

# Parameters; the step's config file (STEP["CONFIG"]) overrides these, CLI flags override the config.
DEFAULT_CONFIG: dict[str, Any] = {
    "QTL_TYPES": STEP["DEFAULT_QTL_TYPES"],
    "STUDY_REGEX": STEP["DEFAULT_STUDY_REGEX"],
    "TISSUE_REGEX": "",
    "DATASET_IDS": [],
    "CANDIDATE_P": 1e-5,
    "MAX_CANDIDATES_PER_CONTEXT_LOCUS": 100,
    "P1": 1e-4,
    "P2": 1e-4,
    "P12": 5e-6,
    "P12_LOW": 1e-6,
    "P12_HIGH": 1e-5,
    "STRONG_PP4": 0.80,
    "SUGGESTIVE_PP4": 0.50,
    "OVERLAP_MIN": 0.5,
    "TRIM_BY_POSTERIOR": True,
    "SAVE_VARIANT_PP_MIN_PP4": 0.50,
    "SLURM": {},
}

ANCESTRY_ALIASES = {
    "eur": ("EUR", "European"), "european": ("EUR", "European"),
    "afr": ("AFR", "African"), "african": ("AFR", "African"),
    "eas": ("EAS", "East Asian"), "east asian": ("EAS", "East Asian"),
    "sas": ("SAS", "South Asian"), "south asian": ("SAS", "South Asian"),
    "amr": ("AMR", "Hispanic or Latin American"), "hispanic": ("AMR", "Hispanic or Latin American"),
    "latino": ("AMR", "Hispanic or Latin American"),
    "hispanic or latin american": ("AMR", "Hispanic or Latin American"),
}

NOT_TESTED_TERMINAL = {
    "EXCLUDED_NO_PROVIDER_CS_FOR_TRAIT", "EXCLUDED_PROVIDER_CS_OUTSIDE_LOCUS",
    "NOT_TESTED_NO_PROVIDER_LBF_ROWS", "NOT_TESTED_NO_PROVIDER_LBF_COLUMNS", "NOT_TESTED_NO_FINITE_LBF_SIGNAL",
    "NOT_TESTED_NO_SHARED_VARIANTS", "NOT_TESTED_INSUFFICIENT_POSTERIOR_OVERLAP",
}


# =============================================================================
# GENERAL HELPERS
# =============================================================================

def banner(text: str) -> None:
    print()
    print("=" * 112)
    print(text)
    print("=" * 112, flush=True)


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_text(value: Any) -> str:
    value = str(value).lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def slugify(value: Any) -> str:
    return normalize_text(value).replace(" ", "_")


def canonical_ancestry(value: str) -> tuple[str, str]:
    key = normalize_text(value)
    if key not in ANCESTRY_ALIASES:
        raise SystemExit(f"Unsupported ancestry {value!r}; use EUR, AFR, EAS, SAS or AMR.")
    return ANCESTRY_ALIASES[key]


def sha1_short(text: str, n: int = 16) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:n]


def nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def complete_file(path: Path) -> bool:
    """Non-empty and not still being downloaded (.aria2 / .part / .lock)."""
    return nonempty(path) and not any(Path(str(path) + s).exists() for s in (".aria2", ".part", ".lock"))


def safe_int(value: Any) -> int | None:
    try:
        if value is None or str(value).strip() == "":
            return None
        return int(float(value))
    except Exception:
        return None


def numeric(series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def first_existing(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    lookup = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if str(candidate).lower() in lookup:
            return lookup[str(candidate).lower()]
    return None


def read_json(path: Path) -> dict:
    if not nonempty(path):
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [clean_json(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(float(value)) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(clean_json(payload), indent=2, default=str, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    frame.to_csv(tmp, sep="\t", index=False, compression="gzip" if path.name.endswith(".gz") else None)
    os.replace(tmp, path)


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    if not nonempty(path):
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False, **kwargs)


def find_executable(name: str) -> str | None:
    local = PROJECT_ROOT / "envs" / "pipeline" / "bin" / name
    if local.exists() and os.access(local, os.X_OK):
        return str(local)
    return shutil.which(name)


def make_variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    c = re.sub(r"^chr", "", str(chrom).strip(), flags=re.I)
    try:
        c_int, p_int = int(float(c)), int(float(pos))
    except Exception:
        return ""
    ref, alt = str(ref).strip().upper(), str(alt).strip().upper()
    if not (1 <= c_int <= 22) or p_int <= 0 or ref in {"", "NAN", "NONE", "."} or alt in {"", "NAN", "NONE", "."}:
        return ""
    return f"{c_int}:{p_int}:{ref}:{alt}"


VARIANT_ID = re.compile(r"^(?:chr)?(\d{1,2})[:_](\d+)[:_]([A-Za-z]+)[:_]([A-Za-z]+)$")


def variant_keys(ids: pd.Series) -> pd.Series:
    """chr1_114643302_G_A / 1:114643302:G:A -> 1:114643302:G:A ('' if unparseable or not chr1-22)."""
    parts = ids.astype(str).str.extract(VARIANT_ID)
    chrom = pd.to_numeric(parts[0], errors="coerce")
    ok = parts.notna().all(axis=1) & chrom.between(1, 22)
    keys = (chrom.astype("Int64").astype(str) + ":" + parts[1] + ":" + parts[2].str.upper() + ":"
            + parts[3].str.upper())
    return keys.where(ok, "")


# =============================================================================
# CONFIGURATION AND SLURM
# =============================================================================

def load_config(path: Path) -> dict:
    cfg = {k: (list(v) if isinstance(v, list) else dict(v) if isinstance(v, dict) else v)
           for k, v in DEFAULT_CONFIG.items()}
    if path.is_file():
        values = runpy.run_path(str(path))
        unknown = [k for k in values if k.isupper() and k not in DEFAULT_CONFIG]
        if unknown:
            print(f"[CONFIG] ignoring unknown settings in {path.name}: {', '.join(unknown)}")
        cfg.update({k: v for k, v in values.items() if k in DEFAULT_CONFIG})
        cfg["_FILE"] = str(path)
    else:
        print(f"[CONFIG] {path} not found - using built-in defaults")
        cfg["_FILE"] = ""
    if isinstance(cfg["QTL_TYPES"], str):
        cfg["QTL_TYPES"] = [x.strip() for x in cfg["QTL_TYPES"].split(",") if x.strip()]
    return cfg


def apply_overrides(cfg: dict, args) -> dict:
    if args.qtl_types:
        cfg["QTL_TYPES"] = [x.strip() for x in args.qtl_types.split(",") if x.strip()]
    if args.study_regex is not None:
        cfg["STUDY_REGEX"] = args.study_regex
    for name in ("candidate_p", "strong_pp4", "suggestive_pp4"):
        if getattr(args, name) is not None:
            cfg[name.upper()] = getattr(args, name)
    return cfg


PER_JOB_KEYS = {"job-name", "output", "error", "array"}


def resolve_slurm(args, cfg: dict) -> dict:
    """Built-in < HPC_CPU.sh #SBATCH header < config SLURM {...} < command line."""
    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": "4", "mem": "32G", "time": "24:00:00"}
    template = PROJECT_ROOT / args.template
    if template.is_file():
        for line in template.read_text(encoding="utf-8").splitlines():
            match = re.match(r"^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
            if match and (match.group(2) or "").strip() and match.group(1) not in PER_JOB_KEYS:
                values[match.group(1)] = match.group(2).strip()
    for key, value in (cfg.get("SLURM") or {}).items():
        if value not in (None, "") and key not in PER_JOB_KEYS:
            values[str(key)] = str(value)
    for key, value in (("partition", args.partition), ("qos", args.qos), ("account", args.account),
                       ("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time)):
        if value not in (None, ""):
            values[key] = str(value)
    return values


# =============================================================================
# PATHS: Step02 study numbers, Step06 fine-mapping (per LD panel), this step's results
# =============================================================================

def results_slug(ancestry_slug: str, panel: str) -> str:
    return ancestry_slug if panel == DEFAULT_PANEL else f"{ancestry_slug}_{panel}"


def panel_option(panel: str) -> str:
    return "" if panel == DEFAULT_PANEL else f" --reference-panel {panel}"


def project_paths(args) -> dict:
    code, label = canonical_ancestry(args.ancestry)
    p, a = slugify(args.phenotype), slugify(label)
    r = results_slug(a, args.reference_panel)
    return {
        "phenotype": args.phenotype.strip(), "code": code, "label": label, "p": p, "a": a, "r": r,
        "panel": args.reference_panel,
        "S02": PROJECT_ROOT / "01_gwas_catalog" / p / a,          # study numbers, shared by every panel
        "S06": PROJECT_ROOT / "06_finemapping" / p / r,
        "OUT": PROJECT_ROOT / STEP["OUT_DIR"] / p / r,
        "LOG": PROJECT_ROOT / "logs" / STEP["LOG_DIR"] / p / r,
        "BASH": PROJECT_ROOT / f"{STEP['SCRIPT_PREFIX']}_{p}_{r}.sh",
    }


def base_command(P: dict) -> str:
    return (f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {shlex.quote(P['phenotype'])} "
            f"--ancestry {P['code']}{panel_option(P['panel'])}")


def step02_order(P: dict) -> dict[str, int]:
    path = P["S02"] / "GWAS_download_manifest.tsv"
    if not path.is_file():
        return {}
    m = pd.read_csv(path, sep="\t", dtype=str)
    ids = pd.to_numeric(m.get("ARRAY_ID"), errors="coerce")
    return {str(acc).strip(): int(i) for acc, i in zip(m["STUDY_ACCESSION"], ids) if pd.notna(i)}


def study_for_index(P: dict, index: int) -> str:
    order = step02_order(P)
    if not order:
        raise SystemExit(f"No Step02 download manifest under {P['S02']} - run the Step02 planner first.")
    for accession, i in order.items():
        if i == index:
            return accession
    listing = "\n".join(f"  {i:>3}  {a}" for a, i in sorted(order.items(), key=lambda x: x[1]))
    raise SystemExit(f"--index {index} is not in the Step02 manifest. Available:\n{listing}")


# =============================================================================
# STEP06 INPUT (GWAS side: the SuSiE-RSS fit of every fine-mapped locus)
# =============================================================================

def step06_studies(S06: Path) -> tuple[list[str], list[dict]]:
    studies, excluded = [], []
    if not S06.exists():
        return studies, excluded
    for d in sorted(p for p in S06.iterdir() if p.is_dir() and p.name.startswith("GCST")):
        summary = read_json(d / "finemapping_summary.json")
        status = str(summary.get("STATUS", "")).upper()
        if not summary:
            excluded.append({"STUDY_ACCESSION": d.name, "REASON": "Step06 not finished (no finemapping_summary.json)"})
        elif status == "NO_SIGNIFICANT_LOCI":
            excluded.append({"STUDY_ACCESSION": d.name, "REASON": "no genome-wide significant loci"})
        elif not list_loci(d):
            excluded.append({"STUDY_ACCESSION": d.name, "REASON": "no locus with a SuSiE fit"})
        else:
            studies.append(d.name)
    return studies, excluded


def list_loci(study_dir: Path) -> list[str]:
    loci_dir = study_dir / "loci"
    if not loci_dir.exists():
        return []
    return sorted(d.name for d in loci_dir.iterdir()
                  if d.is_dir() and re.fullmatch(r"L\d+", d.name) and (d / "susie" / "susie_fit.rds").exists())


def load_gwas_locus(study_dir: Path, locus_id: str) -> pd.DataFrame:
    locus_dir = study_dir / "loci" / locus_id
    source = next((p for p in (locus_dir / "reference_qc_matched_variants.tsv.gz", locus_dir / "matched_variants.tsv")
                   if nonempty(p)), None)
    if source is None:
        raise FileNotFoundError(f"No Step06 matched GWAS table for {locus_id} under {locus_dir}")
    x = read_tsv(source)
    cols = {name: first_existing(x.columns, aliases) for name, aliases in {
        "CHR": ["REFERENCE_CHR", "CHR"], "POS": ["REFERENCE_POS", "POS"], "REF": ["REFERENCE_REF"],
        "ALT": ["REFERENCE_ALT"], "REFERENCE_ID": ["REFERENCE_ID"], "P": ["P"]}.items()}
    missing = [k for k, v in cols.items() if v is None]
    if missing:
        raise RuntimeError(f"{source} lacks {missing}")
    out = pd.DataFrame({k: x[v] for k, v in cols.items()})
    out["VARIANT_KEY"] = [make_variant_key(c, p, r, a) for c, p, r, a in
                          zip(out["CHR"], out["POS"], out["REF"], out["ALT"])]
    out["POS"] = numeric(out["POS"])
    out["CHR"] = numeric(out["CHR"].astype(str).str.replace("chr", "", regex=False))
    out = out[out["VARIANT_KEY"].ne("") & out["REFERENCE_ID"].notna()].copy()
    return out.drop_duplicates("VARIANT_KEY").reset_index(drop=True)


def gwas_susie_inputs(study_dir: Path, locus_id: str, code: str, gwas: pd.DataFrame,
                      work: Path) -> tuple[Path, Path]:
    """Step06 SuSiE fit + its variant order as VARIANT_KEY (the LD .vars file order of the fit)."""
    locus = study_dir / "loci" / locus_id
    fit, var_file = locus / "susie" / "susie_fit.rds", locus / f"{code}_LD.unphased.vcor1.bin.vars"
    if not nonempty(fit) or not nonempty(var_file):
        raise FileNotFoundError(f"Step06 SuSiE fit or LD variant list missing for {locus_id}: {fit} / {var_file}")
    order = [v.strip() for v in var_file.read_text().splitlines() if v.strip()]
    mapping = gwas[["REFERENCE_ID", "VARIANT_KEY"]].drop_duplicates("REFERENCE_ID").set_index("REFERENCE_ID")
    missing = [v for v in order if v not in mapping.index]
    if missing:
        raise RuntimeError(f"{len(missing)} SuSiE variants of {locus_id} are not in the Step06 matched table "
                           f"(e.g. {missing[:3]}) - rerun Step06 for this study")
    path = work / "gwas_susie_variant_map.tsv"
    pd.DataFrame({"REFERENCE_ID": order, "VARIANT_KEY": mapping.loc[order, "VARIANT_KEY"].values}).to_csv(
        path, sep="\t", index=False)
    return fit, path


def gwas_cs_indices(rscript: str, fit: Path) -> list[int]:
    code = ("fit <- readRDS(commandArgs(trailingOnly=TRUE)[1]); "
            "idx <- tryCatch(fit$sets$cs_index, error=function(e) NULL); "
            "if (length(idx)) cat(paste(as.integer(idx), collapse=','))")
    proc = subprocess.run([rscript, "-e", code, str(fit)], text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Cannot read the Step06 SuSiE credible sets of {fit}: {proc.stderr.strip()[-500:]}")
    return sorted({int(v) for v in proc.stdout.strip().split(",") if v.strip().isdigit() and int(v) >= 1})


# =============================================================================
# QTL RESOURCES: eQTL Catalogue metadata, dense summary statistics, provider SuSiE (LBF + CS)
#   resources/coloc/eqtl_catalogue/tabix_ftp_paths*.tsv
#   resources/coloc/eqtl_catalogue/dense/<dataset_id>/<file>.tsv.gz (+.tbi)
#   resources/coloc/eqtl_catalogue/susie_provider/<dataset_id>/<id>.credible_sets.tsv.gz
#   resources/coloc/eqtl_catalogue/susie_provider/<dataset_id>/<id>.lbf_variable.sorted.tsv.gz (+.tbi)
#   (Step00_Download_Resources.py downloads them and position-sorts + tabix-indexes every LBF file)
# =============================================================================

def quant_method_to_qtl_type(value: Any) -> str | None:
    q = normalize_text(value)
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


def basename(url: Any) -> str:
    return str(url).strip().split("?", 1)[0].rstrip("/").rsplit("/", 1)[-1]


def usable(value: Any) -> bool:
    return str(value).strip().lower() not in {"", "na", "nan", "none", "."}


def lbf_index_path(path: Path) -> Path:
    """Position-sorted, bgzipped, tabix-indexed copy written by Step00 (same rule as Step00)."""
    name = path.name
    stem = name[:-7] if name.endswith(".txt.gz") else name[:-3] if name.endswith(".gz") else name
    return path.with_name(stem + ".sorted.tsv.gz")


def load_catalogue() -> pd.DataFrame:
    frames = []
    for name in CATALOGUE_TABLES:
        path = EQTL_ROOT / name
        if nonempty(path):
            frames.append(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False))
    if not frames:
        return pd.DataFrame()
    x = pd.concat(frames, ignore_index=True)
    for column in ("dataset_id", "study_label", "tissue_label", "sample_group", "condition_label", "sample_size",
                   "quant_method", "ftp_path", "ftp_cs_path", "ftp_lbf_path"):
        if column not in x.columns:
            x[column] = ""
    x["QTL_TYPE"] = x["quant_method"].map(quant_method_to_qtl_type)
    return x.drop_duplicates("dataset_id", keep="last").reset_index(drop=True)


def select_contexts(catalogue: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    if catalogue.empty:
        return catalogue
    types = {t.lower() for t in cfg["QTL_TYPES"]}
    x = catalogue[catalogue["QTL_TYPE"].fillna("").str.lower().isin(types)].copy()
    if cfg["DATASET_IDS"]:
        x = x[x["dataset_id"].isin({str(d) for d in cfg["DATASET_IDS"]})]
    elif cfg["STUDY_REGEX"]:
        pattern = re.compile(cfg["STUDY_REGEX"], flags=re.I)
        x = x[x["study_label"].map(lambda s: bool(pattern.search(str(s))))]
    if cfg["TISSUE_REGEX"]:
        pattern = re.compile(cfg["TISSUE_REGEX"], flags=re.I)
        x = x[(x["tissue_label"] + " " + x["sample_group"]).map(lambda s: bool(pattern.search(s)))]
    x = x[x["ftp_path"].map(usable) & x["ftp_lbf_path"].map(usable) & x["ftp_cs_path"].map(usable)]
    return x.reset_index(drop=True)


def resolve_resources(contexts: pd.DataFrame) -> pd.DataFrame:
    """Local files and readiness of every selected context (nothing is downloaded here)."""
    rows = []
    for _, c in contexts.iterrows():
        d = str(c["dataset_id"])
        dense = EQTL_ROOT / "dense" / d / basename(c["ftp_path"])
        if not nonempty(dense) and nonempty(EQTL_ROOT / "dense" / basename(c["ftp_path"])):
            dense = EQTL_ROOT / "dense" / basename(c["ftp_path"])
        cs = EQTL_ROOT / "susie_provider" / d / basename(c["ftp_cs_path"])
        lbf = EQTL_ROOT / "susie_provider" / d / basename(c["ftp_lbf_path"])
        lbf_idx = lbf_index_path(lbf)
        dense_ok = complete_file(dense) and nonempty(Path(str(dense) + ".tbi"))
        cs_ok = complete_file(cs)
        indexed = complete_file(lbf_idx) and nonempty(Path(str(lbf_idx) + ".tbi"))
        downloaded = complete_file(lbf)
        status = ("READY" if dense_ok and cs_ok and indexed else
                  "DENSE_MISSING" if not dense_ok else
                  "PROVIDER_CS_MISSING" if not cs_ok else
                  "LBF_NOT_INDEXED" if downloaded else "LBF_MISSING")
        rows.append({**c.to_dict(), "DENSE_PATH": str(dense), "CS_PATH": str(cs), "LBF_PATH": str(lbf),
                     "LBF_INDEX_PATH": str(lbf_idx), "DENSE_OK": dense_ok, "CS_OK": cs_ok,
                     "LBF_DOWNLOADED": downloaded or indexed, "LBF_INDEXED": indexed, "RESOURCE_STATUS": status})
    return pd.DataFrame(rows)


def context_hash(ready: pd.DataFrame) -> str:
    return sha1_short("|".join(sorted(ready["dataset_id"].astype(str)))) if len(ready) else "none"


def current_contexts(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(every selected context with its resource status, the READY ones)."""
    selected = select_contexts(load_catalogue(), cfg)
    if selected.empty:
        return selected, selected
    status = resolve_resources(selected)
    return status, status[status["RESOURCE_STATUS"] == "READY"].reset_index(drop=True)


def step00_command(cfg: dict) -> str:
    regex = cfg["STUDY_REGEX"]
    if STEP["STEP00_GROUP"] == "qtl_celltype":
        flags = f" --celltype-study-regex {shlex.quote(regex)}" if regex else ""
        types = ",".join(cfg["QTL_TYPES"])
        flags += f" --celltype-qtl-types {shlex.quote(types)}" if types.lower() != "eqtl" else ""
        groups = "eqtl_catalogue qtl_celltype"
    else:
        flags = "" if regex == STEP["DEFAULT_STUDY_REGEX"] else (
            f" --catalogue-study-regex {shlex.quote(regex)}" if regex else " --catalogue-all-studies")
        groups = "eqtl_catalogue qtl_catalogue qtl_catalogue_susie"
    if cfg["TISSUE_REGEX"] and STEP["STEP00_GROUP"] != "qtl_celltype":
        flags += f" --catalogue-tissue-regex {shlex.quote(cfg['TISSUE_REGEX'])}"
    return (f"envs/pipeline/bin/python Step00_Download_Resources.py --groups {groups}{flags} "
            f"--output Step00_Download_{STEP['NAME']}.sh\n  sbatch Step00_Download_{STEP['NAME']}.sh")


# =============================================================================
# TABIX READING (pysam; one open handle per file, reused across loci of a task)
# =============================================================================

_TABIX: OrderedDict[str, Any] = OrderedDict()
_HEADERS: dict[str, list[str]] = {}


def tabix_file(path: str):
    import pysam
    if path in _TABIX:
        _TABIX.move_to_end(path)
        return _TABIX[path]
    handle = pysam.TabixFile(path)
    _TABIX[path] = handle
    while len(_TABIX) > 64:
        _TABIX.popitem(last=False)[1].close()
    return handle


def file_header(path: str) -> list[str]:
    if path not in _HEADERS:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            _HEADERS[path] = handle.readline().rstrip("\n\r").lstrip("#").split("\t")
    return _HEADERS[path]


def tabix_fetch(path: str, chrom: int, start: int, end: int) -> list[list[str]]:
    handle = tabix_file(path)
    contigs = set(handle.contigs)
    contig = next((c for c in (str(chrom), f"chr{chrom}") if c in contigs), None)
    if contig is None:
        return []
    return [line.split("\t") for line in handle.fetch(contig, max(0, start - 1), end)]


# =============================================================================
# ONE LOCUS: discovery (dense QTL) -> provider credible-set filter -> provider LBF -> coloc.bf_bf
# =============================================================================

def discover_candidates(ready: pd.DataFrame, chrom: int, start: int, end: int,
                        cfg: dict) -> tuple[pd.DataFrame, list[dict]]:
    """Molecular traits with P <= CANDIDATE_P in the locus window, per ready context."""
    candidates, audit = [], []
    for _, c in ready.iterrows():
        status, n_rows, n_cand = "OK", 0, 0
        try:
            header = file_header(c["DENSE_PATH"])
            i_trait = header.index(first_existing(header, ["molecular_trait_id", "phenotype_id"]))
            i_p = header.index(first_existing(header, ["pvalue", "p_value", "pval"]))
            gene_col = first_existing(header, ["gene_id"])
            i_gene = header.index(gene_col) if gene_col else None
            best: dict[str, list] = {}
            for fields in tabix_fetch(c["DENSE_PATH"], chrom, start, end):
                n_rows += 1
                try:
                    p = float(fields[i_p])
                except (ValueError, IndexError):
                    continue
                if p <= cfg["CANDIDATE_P"]:
                    trait = fields[i_trait]
                    entry = best.setdefault(trait, [p, 0, fields[i_gene] if i_gene is not None else ""])
                    entry[0] = min(entry[0], p)
                    entry[1] += 1
            if not n_rows:
                status = "NO_REGION_ROWS"
            elif not best:
                status = "NO_TRAIT_P_LE_THRESHOLD"
            for trait, (p, n, gene) in sorted(best.items(), key=lambda kv: kv[1][0])[
                    :int(cfg["MAX_CANDIDATES_PER_CONTEXT_LOCUS"])]:
                n_cand += 1
                candidates.append({
                    "DATASET_ID": c["dataset_id"], "STUDY_LABEL": c["study_label"], "QTL_TYPE": c["QTL_TYPE"],
                    "TISSUE": c["tissue_label"], "SAMPLE_GROUP": c["sample_group"],
                    "CONDITION": c["condition_label"], "QTL_N": safe_int(c["sample_size"]),
                    "MOLECULAR_TRAIT_ID": trait, "GENE_ID": gene, "DISCOVERY_MIN_P": p, "N_DISCOVERY_ROWS": n,
                    "CS_PATH": c["CS_PATH"], "LBF_INDEX_PATH": c["LBF_INDEX_PATH"],
                })
        except Exception as exc:
            status = f"ERROR: {type(exc).__name__}: {exc}"
        audit.append({"DATASET_ID": c["dataset_id"], "STUDY_LABEL": c["study_label"], "QTL_TYPE": c["QTL_TYPE"],
                      "TISSUE": c["tissue_label"], "N_REGION_ROWS": n_rows, "N_CANDIDATES": n_cand,
                      "STATUS": status})
    return pd.DataFrame(candidates), audit


def provider_cs_components(cs_path: str, traits: set[str]) -> dict[str, dict[int, tuple[int, int, int]]]:
    """trait -> {SuSiE component L: (chrom, min_pos, max_pos)} from the provider credible-set file."""
    header = file_header(cs_path)
    trait_col = first_existing(header, ["molecular_trait_id"])
    cs_col = first_existing(header, ["cs_id", "cs_index"])
    var_col = first_existing(header, ["variant", "variant_id"])
    purity_col = first_existing(header, ["low_purity"])
    if not trait_col or not cs_col or not var_col:
        raise RuntimeError(f"Unexpected provider credible-set columns in {cs_path}: {header[:15]}")
    usecols = [c for c in (trait_col, cs_col, var_col, purity_col) if c]
    x = pd.read_csv(cs_path, sep="\t", usecols=usecols, dtype=str, keep_default_na=False)
    x = x[x[trait_col].isin(traits)]
    if purity_col:
        x = x[~x[purity_col].str.lower().isin({"true", "t", "1", "yes"})]
    if x.empty:
        return {}
    comp = x[cs_col].str.extract(r"(?:^|[_:.-])(?:L|CS)?(\d+)$", expand=False)
    parts = x[var_col].str.extract(VARIANT_ID)
    frame = pd.DataFrame({"TRAIT": x[trait_col].values, "COMP": pd.to_numeric(comp, errors="coerce").values,
                          "CHR": pd.to_numeric(parts[0], errors="coerce").values,
                          "POS": pd.to_numeric(parts[1], errors="coerce").values}).dropna()
    out: dict[str, dict[int, tuple[int, int, int]]] = {}
    for (trait, comp), g in frame.groupby(["TRAIT", "COMP"]):
        out.setdefault(str(trait), {})[int(comp)] = (int(g["CHR"].mode().iloc[0]), int(g["POS"].min()),
                                                     int(g["POS"].max()))
    return out


def provider_lbf_rows(lbf_index: str, chrom: int, start: int, end: int, traits: set[str]) -> pd.DataFrame:
    header = file_header(lbf_index)
    i_trait = header.index(first_existing(header, ["molecular_trait_id"]))
    rows = [f for f in tabix_fetch(lbf_index, chrom, start, end) if f[i_trait] in traits]
    if not rows:
        return pd.DataFrame(columns=header)
    return pd.DataFrame(rows, columns=header[:len(rows[0])])


R_WORKER = r'''
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 15) stop("Expected 15 arguments")
fit_file <- args[1]; map_file <- args[2]; qtl_file <- args[3]; cand_file <- args[4]
out_summary <- args[5]; out_status <- args[6]; out_varpp <- args[7]
p1 <- as.numeric(args[8]); p2 <- as.numeric(args[9]); p12 <- as.numeric(args[10])
p12_low <- as.numeric(args[11]); p12_high <- as.numeric(args[12])
overlap_min <- as.numeric(args[13]); trim <- as.logical(args[14]); save_min <- as.numeric(args[15])
suppressPackageStartupMessages(library(coloc))
if (!("coloc.bf_bf" %in% getNamespaceExports("coloc"))) stop("installed coloc has no coloc.bf_bf (need coloc >= 5)")

# GWAS side: Step06 susieR fit, signal x variant LBFs in the Step06 LD variant order.
fit <- readRDS(fit_file)
map <- read.delim(map_file, stringsAsFactors=FALSE, check.names=FALSE)
if (is.null(fit$lbf_variable)) stop("Step06 SuSiE object has no lbf_variable")
gbf <- as.matrix(fit$lbf_variable)
if (ncol(gbf) != nrow(map) && nrow(gbf) == nrow(map)) gbf <- t(gbf)
if (ncol(gbf) != nrow(map)) stop(sprintf("GWAS lbf_variable %dx%d does not match %d mapped variants",
                                         nrow(gbf), ncol(gbf), nrow(map)))
colnames(gbf) <- map$VARIANT_KEY
gidx <- tryCatch(as.integer(fit$sets$cs_index), error=function(e) integer())
gidx <- unique(gidx[is.finite(gidx) & gidx >= 1 & gidx <= nrow(gbf)])
if (length(gidx) == 0) stop("Step06 GWAS SuSiE has no credible-set component")
gbf <- gbf[gidx, , drop=FALSE]
rownames(gbf) <- paste0("L", gidx)
gbf <- gbf[, !is.na(colnames(gbf)) & nzchar(colnames(gbf)) & !duplicated(colnames(gbf)), drop=FALSE]
gbf <- gbf[apply(gbf, 1, function(z) any(is.finite(z))), , drop=FALSE]

# QTL side: provider SuSiE LBFs of every candidate, restricted to its provider credible-set components.
qtl <- read.delim(gzfile(qtl_file), stringsAsFactors=FALSE, check.names=FALSE)
cands <- read.delim(cand_file, stringsAsFactors=FALSE, check.names=FALSE, colClasses="character")
rows_of <- split(seq_len(nrow(qtl)), qtl$CANDIDATE_ID)
status <- list(); summaries <- list(); varpp <- list()
priors <- list(DEFAULT=p12, LOW_P12=p12_low, HIGH_P12=p12_high)

for (i in seq_len(nrow(cands))) {
  cid <- cands$CANDIDATE_ID[i]
  st <- list(CANDIDATE_ID=cid, R_STATUS="", N_GWAS_SIGNALS=nrow(gbf), N_QTL_SIGNALS=0L,
             N_QTL_LBF_VARIANTS=0L, N_COMMON_VARIANTS=0L, MESSAGE="")
  res <- tryCatch({
    rows <- rows_of[[cid]]
    comps <- suppressWarnings(as.integer(strsplit(cands$QTL_COMPONENTS[i], ",", fixed=TRUE)[[1]]))
    cols <- intersect(paste0("lbf_variable", comps[is.finite(comps)]), names(qtl))
    if (is.null(rows) || length(rows) == 0) { st$R_STATUS <- "NO_PROVIDER_LBF_ROWS"; NULL }
    else if (length(cols) == 0) { st$R_STATUS <- "NO_PROVIDER_LBF_COLUMNS"; NULL }
    else {
      sub <- qtl[rows, c("VARIANT_KEY", cols), drop=FALSE]
      sub <- sub[nzchar(sub$VARIANT_KEY) & !duplicated(sub$VARIANT_KEY), , drop=FALSE]
      qbf <- t(as.matrix(sub[, cols, drop=FALSE]))
      storage.mode(qbf) <- "double"
      colnames(qbf) <- sub$VARIANT_KEY
      rownames(qbf) <- sub("^lbf_variable", "L", cols)
      qbf <- qbf[apply(qbf, 1, function(z) any(is.finite(z))), , drop=FALSE]
      st$N_QTL_SIGNALS <- nrow(qbf)
      st$N_QTL_LBF_VARIANTS <- ncol(qbf)
      st$N_COMMON_VARIANTS <- length(intersect(colnames(gbf), colnames(qbf)))
      if (nrow(qbf) == 0) { st$R_STATUS <- "NO_FINITE_LBF_SIGNAL"; NULL }
      else if (st$N_COMMON_VARIANTS == 0) { st$R_STATUS <- "NO_SHARED_VARIANTS"; NULL }
      else {
        out <- list()
        for (label in names(priors)) {
          z <- suppressWarnings(coloc::coloc.bf_bf(gbf, qbf, p1=p1, p2=p2, p12=priors[[label]],
                                                   overlap.min=overlap_min, trim_by_posterior=trim))
          sm <- tryCatch(as.data.frame(z$summary), error=function(e) data.frame())
          if (nrow(sm) > 0) {
            sm$CANDIDATE_ID <- cid; sm$PRIOR_SET <- label; sm$P12 <- priors[[label]]
            out[[label]] <- list(summary=sm, obj=z)
          }
        }
        if (is.null(out[["DEFAULT"]])) { st$R_STATUS <- "INSUFFICIENT_POSTERIOR_OVERLAP"; NULL }
        else { st$R_STATUS <- "OK"; out }
      }
    }
  }, error=function(e) { st$R_STATUS <<- "ERROR"; st$MESSAGE <<- gsub("[\t\n]", " ", conditionMessage(e)); NULL })
  status[[length(status) + 1]] <- as.data.frame(st, stringsAsFactors=FALSE)
  if (!is.null(res)) {
    for (label in names(res)) summaries[[length(summaries) + 1]] <- res[[label]]$summary
    h4 <- suppressWarnings(max(as.numeric(res$DEFAULT$summary$PP.H4.abf), na.rm=TRUE))
    if (is.finite(h4) && h4 >= save_min) {
      rr <- tryCatch(as.data.frame(res$DEFAULT$obj$results), error=function(e) data.frame())
      if (nrow(rr) > 0 && "snp" %in% names(rr)) {
        for (pc in setdiff(names(rr), "snp")) {
          varpp[[length(varpp) + 1]] <- data.frame(CANDIDATE_ID=cid, PAIR=pc, VARIANT_KEY=rr$snp,
                                                   SNP_PP_H4=as.numeric(rr[[pc]]), stringsAsFactors=FALSE)
        }
      }
    }
  }
}
write.table(do.call(rbind, status), out_status, sep="\t", quote=FALSE, row.names=FALSE)
if (length(summaries)) {
  keep <- Reduce(intersect, lapply(summaries, names))
  write.table(do.call(rbind, lapply(summaries, function(s) s[, keep, drop=FALSE])), out_summary,
              sep="\t", quote=FALSE, row.names=FALSE)
} else {
  write.table(data.frame(CANDIDATE_ID=character()), out_summary, sep="\t", quote=FALSE, row.names=FALSE)
}
if (length(varpp)) {
  con <- gzfile(out_varpp, "w"); write.table(do.call(rbind, varpp), con, sep="\t", quote=FALSE, row.names=FALSE)
  close(con)
}
cat("CANDIDATES=", nrow(cands), " TESTED=", sum(sapply(status, function(s) s$R_STATUS == "OK")), "\n", sep="")
'''


def classify(summary: pd.DataFrame, strong: float, suggestive: float) -> dict:
    """Classify one candidate from its best (max PP.H4) signal pair - one consistent posterior vector."""
    h4 = numeric(summary["PP.H4.abf"])
    if not np.isfinite(h4.to_numpy(dtype=float)).any():
        return {"STATUS": "NOT_TESTED_NO_FINITE_H4"}
    best = summary.loc[h4.idxmax()]
    pp = [float(numeric(pd.Series([best.get(f"PP.H{i}.abf")])).iloc[0]) for i in range(5)]
    dominant = int(np.nanargmax(pp))
    pp4 = pp[4]
    if pp4 >= strong:
        label = "STRONG_SHARED_SIGNAL"
    elif pp4 >= suggestive:
        label = "SUGGESTIVE_SHARED_SIGNAL"
    else:
        label = {0: "NEITHER_ASSOCIATED_FAVORED", 1: "GWAS_ONLY_FAVORED", 2: "QTL_ONLY_FAVORED",
                 3: "BOTH_ASSOCIATED_DIFFERENT_SIGNALS_FAVORED", 4: "SHARED_SIGNAL_BELOW_THRESHOLD"}[dominant]
    out = {f"PP_H{i}": pp[i] for i in range(5)}
    out.update({"BEST_GWAS_VARIANT": str(best.get("hit1", "")), "BEST_QTL_VARIANT": str(best.get("hit2", "")),
                "BEST_GWAS_SIGNAL": str(best.get("idx1", "")), "BEST_QTL_SIGNAL": str(best.get("idx2", "")),
                "DOMINANT_HYPOTHESIS": f"H{dominant}", "N_SIGNAL_PAIRS": int(len(summary)),
                "COLOC_CLASS": label, "COLOCALIZED_STRONG": "YES" if pp4 >= strong else "NO",
                "SHARED_SIGNAL": "YES" if pp4 >= suggestive else "NO"})
    return out


def sensitivity_pp4(rows: pd.DataFrame, label: str, hit1: str, hit2: str) -> float:
    sub = rows[rows["PRIOR_SET"] == label]
    same = sub[(sub["hit1"].astype(str) == hit1) & (sub["hit2"].astype(str) == hit2)]
    sub = same if len(same) else sub
    vals = numeric(sub["PP.H4.abf"]) if len(sub) else pd.Series(dtype=float)
    return float(vals.max()) if vals.notna().any() else np.nan


def process_locus(row: pd.Series, args, cfg: dict, ready: pd.DataFrame, rscript: str) -> str:
    """Colocalise one Step06 locus with every ready QTL context. Returns the locus status."""
    P = project_paths(args)
    accession, locus_id = str(row["STUDY_ACCESSION"]), str(row["LOCUS_ID"])
    study_dir = P["S06"] / accession
    out = P["OUT"] / accession / "loci" / locus_id
    done_file = out / DONE_FILE
    c_hash = context_hash(ready)
    previous = read_json(done_file)
    if previous and not args.force and previous.get("CONTEXT_HASH") == c_hash:
        print(f"[RESUME] {accession} {locus_id}: {previous.get('STATUS')} (same {len(ready)} QTL contexts)")
        return str(previous.get("STATUS"))
    out.mkdir(parents=True, exist_ok=True)
    (out / FAILED_FILE).unlink(missing_ok=True)
    work = out / "work"
    work.mkdir(exist_ok=True)
    started = time.time()
    banner(f"{STEP['NAME']} - {accession} {locus_id} - {len(ready)} QTL contexts")

    def finish(status: str, results: pd.DataFrame, extra: dict) -> str:
        write_tsv(results, out / RESULTS_FILE)
        tested = results[results.get("COLOC_TESTED", pd.Series(dtype=str)).astype(str) == "YES"] \
            if len(results) else results
        write_json(done_file, {
            "VERSION": VERSION, "STEP": STEP["NAME"], "STATUS": status, "STUDY_ACCESSION": accession,
            "LOCUS_ID": locus_id, "CONTEXT_HASH": c_hash, "N_CONTEXTS_READY": int(len(ready)),
            "N_CANDIDATES": int(len(results)), "N_TESTED": int(len(tested)),
            "N_STRONG": int((tested.get("COLOCALIZED_STRONG", pd.Series(dtype=str)) == "YES").sum()) if len(tested) else 0,
            "SECONDS": round(time.time() - started, 1), "COMPLETED_UTC": utcnow(), **extra})
        if not args.keep_work:
            shutil.rmtree(work, ignore_errors=True)
        print(f"[OK] {accession} {locus_id}: {status}; {len(results)} candidates, {len(tested)} tested "
              f"({time.time() - started:.0f}s)")
        return status

    # 1. GWAS side from Step06.
    gwas = load_gwas_locus(study_dir, locus_id)
    chroms = sorted(set(gwas["CHR"].dropna().astype(int)))
    if len(chroms) != 1:
        raise RuntimeError(f"Locus spans chromosomes {chroms}")
    chrom, start, end = chroms[0], int(gwas["POS"].min()), int(gwas["POS"].max())
    fit, gmap = gwas_susie_inputs(study_dir, locus_id, P["code"], gwas, work)
    gidx = gwas_cs_indices(rscript, fit)
    region = {"CHR": chrom, "LOCUS_START": start, "LOCUS_END": end, "GWAS_CS_INDEXES": ",".join(f"L{i}" for i in gidx)}
    if not gidx:
        return finish("NO_GWAS_CREDIBLE_SET", pd.DataFrame(columns=RESULT_COLUMNS), region)

    # 2. Discovery in the dense QTL summary statistics of every ready context.
    candidates, audit = discover_candidates(ready, chrom, start, end, cfg)
    write_tsv(pd.DataFrame(audit), out / "discovery_audit.tsv")
    if candidates.empty:
        return finish("NO_QTL_CANDIDATES", pd.DataFrame(columns=RESULT_COLUMNS), region)
    candidates.insert(0, "CANDIDATE_ID", [sha1_short(f"{accession}|{locus_id}|{d}|{t}", 20) for d, t in
                                          zip(candidates["DATASET_ID"], candidates["MOLECULAR_TRAIT_ID"])])
    for key, value in (("STUDY_ACCESSION", accession), ("LOCUS_ID", locus_id), *region.items()):
        candidates[key] = value
    for column in RESULT_COLUMNS:
        if column not in candidates.columns:
            candidates[column] = np.nan if column.startswith(("PP_", "N_")) else ""
    candidates["COLOC_TESTED"] = "NO"

    # 3. Provider credible-set filter: only real provider SuSiE components overlapping the GWAS locus.
    # 4. Provider LBFs of the surviving traits, read from the position-indexed LBF file (locus window only).
    lbf_frames = []
    for dataset_id, group in candidates.groupby("DATASET_ID", sort=False):
        traits = set(group["MOLECULAR_TRAIT_ID"])
        comps = provider_cs_components(group["CS_PATH"].iloc[0], traits)
        keep_traits = set()
        for i, r in group.iterrows():
            trait_comps = comps.get(r["MOLECULAR_TRAIT_ID"], {})
            overlap = [k for k, (c, lo, hi) in sorted(trait_comps.items()) if c == chrom and hi >= start and lo <= end]
            if not trait_comps:
                candidates.at[i, "STATUS"] = "EXCLUDED_NO_PROVIDER_CS_FOR_TRAIT"
            elif not overlap:
                candidates.at[i, "STATUS"] = "EXCLUDED_PROVIDER_CS_OUTSIDE_LOCUS"
                candidates.at[i, "QTL_CS_INDEXES"] = ",".join(f"L{k}" for k in sorted(trait_comps))
            else:
                candidates.at[i, "QTL_CS_INDEXES"] = ",".join(f"L{k}" for k in overlap)
                keep_traits.add(r["MOLECULAR_TRAIT_ID"])
        if keep_traits:
            lbf = provider_lbf_rows(group["LBF_INDEX_PATH"].iloc[0], chrom, start, end, keep_traits)
            if len(lbf):
                trait_col = first_existing(lbf.columns, ["molecular_trait_id"])
                var_col = first_existing(lbf.columns, ["variant", "variant_id"])
                ids = group.set_index("MOLECULAR_TRAIT_ID")["CANDIDATE_ID"]
                frame = pd.DataFrame({"CANDIDATE_ID": lbf[trait_col].map(ids).values,
                                      "VARIANT_KEY": variant_keys(lbf[var_col]).values})
                for column in [c for c in lbf.columns if str(c).startswith("lbf_variable")]:
                    frame[column] = lbf[column].values
                lbf_frames.append(frame[frame["VARIANT_KEY"].ne("")])

    to_test = candidates[candidates["STATUS"] == ""].copy()
    if to_test.empty:
        return finish("COMPLETE", tidy_results(candidates), region)
    to_test["QTL_COMPONENTS"] = to_test["QTL_CS_INDEXES"].str.replace("L", "", regex=False)
    qtl_long = pd.concat(lbf_frames, ignore_index=True) if lbf_frames else pd.DataFrame(
        columns=["CANDIDATE_ID", "VARIANT_KEY", "lbf_variable1"])
    write_tsv(qtl_long, work / "qtl_provider_lbf.tsv.gz")
    write_tsv(to_test[["CANDIDATE_ID", "QTL_COMPONENTS"]], work / "candidates.tsv")

    # 5. coloc::coloc.bf_bf for every candidate of the locus in one R session.
    worker = work / "provider_coloc_bf_bf.R"
    worker.write_text(R_WORKER, encoding="utf-8")
    summary_file, status_file, varpp_file = out / "coloc_signal_pairs.tsv", work / "r_status.tsv", \
        out / "coloc_variant_pp_h4.tsv.gz"
    varpp_file.unlink(missing_ok=True)
    command = [rscript, str(worker), str(fit), str(gmap), str(work / "qtl_provider_lbf.tsv.gz"),
               str(work / "candidates.tsv"), str(summary_file), str(status_file), str(varpp_file),
               *(str(cfg[k]) for k in ("P1", "P2", "P12", "P12_LOW", "P12_HIGH", "OVERLAP_MIN")),
               "TRUE" if cfg["TRIM_BY_POSTERIOR"] else "FALSE", str(cfg["SAVE_VARIANT_PP_MIN_PP4"])]
    proc = subprocess.run(command, text=True, capture_output=True)
    (out / "coloc.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"R coloc worker failed (log {out / 'coloc.log'}):\n{proc.stderr[-3000:]}")

    r_status = read_tsv(status_file, dtype=str).set_index("CANDIDATE_ID")
    summary = read_tsv(summary_file)
    r_codes = {"NO_PROVIDER_LBF_ROWS": "NOT_TESTED_NO_PROVIDER_LBF_ROWS",
               "NO_PROVIDER_LBF_COLUMNS": "NOT_TESTED_NO_PROVIDER_LBF_COLUMNS",
               "NO_FINITE_LBF_SIGNAL": "NOT_TESTED_NO_FINITE_LBF_SIGNAL",
               "NO_SHARED_VARIANTS": "NOT_TESTED_NO_SHARED_VARIANTS",
               "INSUFFICIENT_POSTERIOR_OVERLAP": "NOT_TESTED_INSUFFICIENT_POSTERIOR_OVERLAP", "ERROR": "FAILED_R"}
    for i, r in to_test.iterrows():
        st = r_status.loc[r["CANDIDATE_ID"]] if r["CANDIDATE_ID"] in r_status.index else None
        if st is None:
            candidates.at[i, "STATUS"] = "FAILED_NO_R_RESULT"
            continue
        for key in ("N_GWAS_SIGNALS", "N_QTL_SIGNALS", "N_QTL_LBF_VARIANTS", "N_COMMON_VARIANTS"):
            candidates.at[i, key] = safe_int(st.get(key))
        if st["R_STATUS"] != "OK":
            candidates.at[i, "STATUS"] = r_codes.get(st["R_STATUS"], "FAILED_R")
            candidates.at[i, "ERROR"] = str(st.get("MESSAGE", "") or "")
            continue
        rows = summary[summary["CANDIDATE_ID"].astype(str) == r["CANDIDATE_ID"]]
        result = classify(rows[rows["PRIOR_SET"] == "DEFAULT"], cfg["STRONG_PP4"], cfg["SUGGESTIVE_PP4"])
        if result.get("STATUS"):
            candidates.at[i, "STATUS"] = result["STATUS"]
            continue
        for key, value in result.items():
            candidates.at[i, key] = value
        hit1, hit2 = result["BEST_GWAS_VARIANT"], result["BEST_QTL_VARIANT"]
        low, high = sensitivity_pp4(rows, "LOW_P12", hit1, hit2), sensitivity_pp4(rows, "HIGH_P12", hit1, hit2)
        candidates.at[i, "PP_H4_P12_LOW"], candidates.at[i, "PP_H4_P12_HIGH"] = low, high
        candidates.at[i, "PRIOR_ROBUST_STRONG"] = "YES" if np.isfinite(low) and low >= cfg["STRONG_PP4"] else "NO"
        candidates.at[i, "COLOC_TESTED"] = "YES"
        candidates.at[i, "STATUS"] = "TESTED"
    return finish("COMPLETE", tidy_results(candidates), region)


RESULT_COLUMNS = [
    "STUDY_ACCESSION", "LOCUS_ID", "CHR", "LOCUS_START", "LOCUS_END", "QTL_TYPE", "STUDY_LABEL", "TISSUE",
    "SAMPLE_GROUP", "CONDITION", "DATASET_ID", "MOLECULAR_TRAIT_ID", "GENE_ID", "QTL_N", "DISCOVERY_MIN_P",
    "N_DISCOVERY_ROWS", "GWAS_CS_INDEXES", "QTL_CS_INDEXES", "STATUS", "COLOC_TESTED", "PP_H4", "COLOC_CLASS",
    "COLOCALIZED_STRONG", "SHARED_SIGNAL", "PRIOR_ROBUST_STRONG", "PP_H4_P12_LOW", "PP_H4_P12_HIGH",
    "PP_H0", "PP_H1", "PP_H2", "PP_H3", "DOMINANT_HYPOTHESIS", "BEST_GWAS_VARIANT", "BEST_QTL_VARIANT",
    "BEST_GWAS_SIGNAL", "BEST_QTL_SIGNAL", "N_SIGNAL_PAIRS", "N_GWAS_SIGNALS", "N_QTL_SIGNALS",
    "N_QTL_LBF_VARIANTS", "N_COMMON_VARIANTS", "ERROR", "CANDIDATE_ID",
]


def tidy_results(candidates: pd.DataFrame) -> pd.DataFrame:
    x = candidates.copy()
    for column in RESULT_COLUMNS:
        if column not in x.columns:
            x[column] = np.nan if column.startswith(("PP_", "N_")) else ""
    x = x[RESULT_COLUMNS]
    order = numeric(x["PP_H4"]).fillna(-1)
    return x.assign(_O=order).sort_values(["_O", "DISCOVERY_MIN_P"], ascending=[False, True]).drop(columns="_O")


# =============================================================================
# CONTEXT SNAPSHOT (frozen when the job script is made, so every task tests the same contexts)
# =============================================================================

def snapshot_file(P: dict) -> Path:
    return P["OUT"] / "qtl_contexts_snapshot.tsv"


def load_snapshot(P: dict) -> pd.DataFrame:
    x = read_tsv(snapshot_file(P), dtype=str, keep_default_na=False)
    if x.empty:
        raise SystemExit(f"QTL context snapshot missing: {snapshot_file(P)} - run the planner again.")
    still = x["DENSE_PATH"].map(lambda p: nonempty(Path(p))) & x["LBF_INDEX_PATH"].map(lambda p: nonempty(Path(p)))
    if (~still).any():
        print(f"[WARN] {int((~still).sum())} snapshot contexts disappeared from disk - not tested")
    return x[still].reset_index(drop=True)


# =============================================================================
# MODES
# =============================================================================

def require_tools() -> str:
    rscript = find_executable("Rscript")
    if not rscript:
        raise SystemExit("Rscript not found - run Step00_Installation.py")
    import importlib.util
    if importlib.util.find_spec("pysam") is None:
        raise SystemExit("python pysam not installed - run Step00_Installation.py")
    return rscript


def run_locus(row: pd.Series, args, cfg, ready, rscript) -> bool:
    try:
        process_locus(row, args, cfg, ready, rscript)
        return True
    except Exception as exc:
        out = project_paths(args)["OUT"] / str(row["STUDY_ACCESSION"]) / "loci" / str(row["LOCUS_ID"])
        out.mkdir(parents=True, exist_ok=True)
        (out / FAILED_FILE).write_text(f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}", encoding="utf-8")
        print(f"[FAILED] {row['STUDY_ACCESSION']} {row['LOCUS_ID']}: {exc}", file=sys.stderr)
        return False


def config_path(args) -> Path:
    path = Path(args.config)
    return path if path.is_absolute() else PROJECT_ROOT / path


def passthrough(args) -> list[str]:
    options = ["--config", str(config_path(args))]
    if args.reference_panel != DEFAULT_PANEL:
        options += ["--reference-panel", args.reference_panel]
    for flag, value in (("--qtl-types", args.qtl_types), ("--study-regex", args.study_regex),
                        ("--candidate-p", args.candidate_p), ("--strong-pp4", args.strong_pp4),
                        ("--suggestive-pp4", args.suggestive_pp4)):
        if value is not None:
            options += [flag, str(value)]
    if args.force:
        options.append("--force")
    if args.keep_work:
        options.append("--keep-work")
    return options


def planner_mode(args, cfg) -> None:
    P = project_paths(args)
    studies, excluded = step06_studies(P["S06"])
    if args.study:
        studies = [s for s in studies if s == args.study]
    if not studies:
        raise SystemExit(f"No fine-mapped studies under {P['S06']} - run Step06{panel_option(P['panel'])} first.")
    status, ready = current_contexts(cfg)
    banner(f"{STEP['NAME']} PLANNER - {STEP['TITLE']} - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    print_context_summary(status, cfg)
    if ready.empty:
        raise SystemExit(f"\nNo QTL context is ready (dense + provider credible sets + indexed provider LBF).\n"
                         f"  {step00_command(cfg)}\nthen: {base_command(P)} --check-resources")
    order = step02_order(P)
    rows = [{"STUDY_INDEX": order.get(acc, ""), "STUDY_ACCESSION": acc, "LOCUS_ID": locus}
            for acc in studies for locus in list_loci(P["S06"] / acc)]
    n_tasks = min(len(rows), MAX_ARRAY_TASKS)
    manifest = pd.DataFrame(rows)
    manifest.insert(0, "TASK_ID", [i % n_tasks + 1 for i in range(len(rows))])   # >1000 loci: several per task
    P["OUT"].mkdir(parents=True, exist_ok=True)
    P["LOG"].mkdir(parents=True, exist_ok=True)
    manifest_file = P["OUT"] / "coloc_manifest.tsv"
    write_tsv(manifest, manifest_file)
    write_tsv(pd.DataFrame(excluded, columns=["STUDY_ACCESSION", "REASON"]), P["OUT"] / "coloc_excluded_studies.tsv")
    write_tsv(ready, snapshot_file(P))

    values = resolve_slurm(args, cfg)
    array = (f"1-{n_tasks}" if n_tasks > 1 else "1") + (f"%{args.max_parallel}" if args.max_parallel else "")
    directives = {**values, "job-name": f"{STEP['JOB']}_{P['p']}_{P['code'].lower()}"[:100], "array": array,
                  "output": f"{P['LOG']}/coloc.%A_%a.out", "error": f"{P['LOG']}/coloc.%A_%a.err"}
    order_keys = ["job-name", "account", "partition", "qos", "nodes", "ntasks", "cpus-per-task", "mem", "time",
                  "array", "output", "error"]
    header = "\n".join(f"#SBATCH --{k}={directives[k]}" for k in order_keys if directives.get(k))
    header += "".join(f"\n#SBATCH --{k}={v}" for k, v in directives.items() if k not in order_keys and v)
    python = PROJECT_ROOT / "envs" / "pipeline" / "bin" / "python"
    command = [str(python), str(Path(__file__).resolve()), "--phenotype", P["phenotype"], "--ancestry", P["code"],
               "--worker-manifest", str(manifest_file)] + passthrough(args)
    P["BASH"].write_text(f"""#!/bin/bash --login
{header}

set -euo pipefail
cd {shlex.quote(str(PROJECT_ROOT))}
export PATH={shlex.quote(str(PROJECT_ROOT / 'envs' / 'pipeline' / 'bin'))}:"$PATH"
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-1}}"

{' '.join(shlex.quote(x) for x in command)}
""", encoding="utf-8", newline="\n")
    P["BASH"].chmod(0o755)

    print(f"\nStudies         : {len(studies)}  ({len(excluded)} excluded, see coloc_excluded_studies.tsv)")
    print(f"Loci            : {len(manifest)} -> array {array}"
          + (f" ({len(manifest) / n_tasks:.1f} loci per task)" if len(manifest) > n_tasks else " (one locus per task)"))
    print(f"QTL contexts    : {len(ready)} ready, frozen in {snapshot_file(P).name}")
    print(f"Config          : {cfg['_FILE'] or 'built-in defaults'}")
    per = manifest.groupby("STUDY_ACCESSION").agg(IDX=("STUDY_INDEX", "first"), LOCI=("LOCUS_ID", "size"),
                                                   FIRST=("TASK_ID", "min"), LAST=("TASK_ID", "max"))
    for accession, r in per.iterrows():
        print(f"  {str(r['IDX']):>4}  {accession:<15} {int(r['LOCI']):>4} loci   tasks {int(r['FIRST'])}-{int(r['LAST'])}")
    print(f"\nManifest: {manifest_file}\n\nRUN:\n  sbatch {P['BASH'].name}")
    print(f"THEN:\n  {base_command(P)} --status\n  {base_command(P)} --stats")


def worker_mode(args, cfg) -> None:
    P = project_paths(args)
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    manifest = read_tsv(Path(args.worker_manifest), dtype=str)
    selected = manifest[numeric(manifest["TASK_ID"]) == task]
    if selected.empty:
        raise SystemExit(f"No manifest rows for array task {task}")
    rscript = require_tools()
    ready = load_snapshot(P)
    failed = [f"{r['STUDY_ACCESSION']} {r['LOCUS_ID']}" for _, r in selected.iterrows()
              if not run_locus(r, args, cfg, ready, rscript)]
    if failed:
        raise SystemExit("Colocalisation failed for: " + ", ".join(failed))


def index_mode(args, cfg) -> None:
    P = project_paths(args)
    accession = study_for_index(P, args.index)
    loci = list_loci(P["S06"] / accession)
    if not loci:
        raise SystemExit(f"Study {args.index} ({accession}) has no fine-mapped locus with a SuSiE fit. Run:\n"
                         f"  envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype "
                         f"{shlex.quote(P['phenotype'])} --ancestry {P['code']} --index {args.index}"
                         f"{panel_option(P['panel'])}")
    if args.locus:
        if args.locus not in loci:
            raise SystemExit(f"{args.locus} is not a fine-mapped locus of {accession}: {', '.join(loci)}")
        loci = [args.locus]
    rscript = require_tools()
    status, ready = current_contexts(cfg)
    print_context_summary(status, cfg)
    if ready.empty:
        raise SystemExit(f"No QTL context is ready.\n  {step00_command(cfg)}")
    failed = [locus for locus in loci if not run_locus(
        pd.Series({"STUDY_ACCESSION": accession, "LOCUS_ID": locus}), args, cfg, ready, rscript)]
    if failed:
        raise SystemExit(f"{accession}: colocalisation failed for {', '.join(failed)}")


def print_context_summary(status: pd.DataFrame, cfg: dict) -> None:
    if status.empty:
        print("QTL contexts    : none selected - check the eQTL Catalogue tables "
              "(Step00 --groups eqtl_catalogue) and the config")
        return
    table = status.groupby(["study_label", "QTL_TYPE"]).agg(
        SELECTED=("dataset_id", "size"), DENSE=("DENSE_OK", "sum"), PROVIDER_CS=("CS_OK", "sum"),
        LBF_DOWNLOADED=("LBF_DOWNLOADED", "sum"), LBF_INDEXED=("LBF_INDEXED", "sum"),
        READY=("RESOURCE_STATUS", lambda s: int((s == "READY").sum()))).reset_index()
    table = table.rename(columns={"study_label": "STUDY"})
    with pd.option_context("display.width", 200, "display.max_rows", None):
        print(table.to_string(index=False))
    n_ready = int((status["RESOURCE_STATUS"] == "READY").sum())
    print(f"QTL contexts ready: {n_ready}/{len(status)}  (QTL types {','.join(cfg['QTL_TYPES'])}; "
          f"studies /{cfg['STUDY_REGEX'] or 'all'}/{'; ' + str(len(cfg['DATASET_IDS'])) + ' dataset ids' if cfg['DATASET_IDS'] else ''})")
    not_ready = status["RESOURCE_STATUS"].value_counts().drop("READY", errors="ignore")
    if len(not_ready):
        print("Not ready (skipped, never counted as 'no colocalisation'): "
              + ", ".join(f"{k} {v}" for k, v in not_ready.items()))


def check_resources_mode(args, cfg) -> None:
    P = project_paths(args)
    problems = []
    banner(f"{STEP['NAME']} - RESOURCE CHECK - {STEP['TITLE']}")
    print(f"Config file     : {cfg['_FILE'] or 'missing - built-in defaults'}")
    print(f"Parameters      : CANDIDATE_P {cfg['CANDIDATE_P']}, priors p1 {cfg['P1']} p2 {cfg['P2']} p12 {cfg['P12']} "
          f"(sensitivity {cfg['P12_LOW']} / {cfg['P12_HIGH']}), strong PP.H4 >= {cfg['STRONG_PP4']}")

    banner("TOOLS")
    rscript = find_executable("Rscript")
    print(f"  {'OK     ' if rscript else 'MISSING'}  Rscript                {rscript or ''}")
    if not rscript:
        problems.append("Rscript not found (Step00_Installation.py)")
    else:
        proc = subprocess.run([rscript, "-e", 'cat(as.character(packageVersion("coloc")), '
                               '"coloc.bf_bf" %in% getNamespaceExports("coloc"))'], text=True, capture_output=True)
        ok = proc.returncode == 0 and proc.stdout.strip().endswith("TRUE")
        print(f"  {'OK     ' if ok else 'MISSING'}  R coloc (coloc.bf_bf)  {proc.stdout.strip() if ok else proc.stderr.strip()[-200:]}")
        if not ok:
            problems.append("R package coloc >= 5 with coloc.bf_bf (Step00_Installation.py)")
    try:
        import pysam
        print(f"  OK       python pysam           {pysam.__version__}")
    except ImportError:
        print("  MISSING  python pysam")
        problems.append("python pysam (Step00_Installation.py)")

    banner("QTL RESOURCES (eQTL Catalogue: dense sumstats + provider SuSiE credible sets + indexed LBF)")
    tables = [EQTL_ROOT / name for name in CATALOGUE_TABLES]
    for t in tables:
        print(f"  {'OK     ' if nonempty(t) else 'MISSING'}  {t.relative_to(PROJECT_ROOT)}")
    status, ready = current_contexts(cfg)
    print()
    print_context_summary(status, cfg)
    if ready.empty:
        problems.append("no QTL context ready")
    if len(status) and (status["RESOURCE_STATUS"] == "LBF_NOT_INDEXED").any():
        print("\nLBF_NOT_INDEXED = downloaded but not yet position-sorted + tabix-indexed. Run the Step00 command "
              "below again: finished downloads are only indexed (and the unsorted copy removed).")

    banner(f"STEP06 INPUT - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    studies, excluded = step06_studies(P["S06"])
    n_loci = sum(len(list_loci(P["S06"] / s)) for s in studies)
    print(f"  Fine-mapping folder : {P['S06']} ({'exists' if P['S06'].exists() else 'MISSING'})")
    print(f"  Studies / loci      : {len(studies)} studies, {n_loci} loci with a SuSiE fit ({len(excluded)} studies excluded)")
    print(f"  Results will go to  : {P['OUT']}")
    if not n_loci:
        problems.append(f"no Step06 loci - run Step06{panel_option(P['panel'])} first")

    banner("CHECK RESULT")
    if problems:
        print(f"NOT READY - {len(problems)} problem(s):")
        for p in problems:
            print(f"  - {p}")
        if any("QTL" in p for p in problems) or (len(status) and len(ready) < len(status)):
            print(f"\nDownload + index the QTL resources with:\n  {step00_command(cfg)}")
        raise SystemExit(1)
    print(f"READY - {len(ready)} QTL contexts, {n_loci} loci.")
    if len(ready) < len(status):
        print(f"({len(status) - len(ready)} contexts not ready yet; get them with:\n  {step00_command(cfg)})")
    print(f"NEXT:\n  {base_command(P)}\n  sbatch {P['BASH'].name}")


def locus_state(P: dict, accession: str, locus: str) -> tuple[str, dict]:
    d = P["OUT"] / accession / "loci" / locus
    done = read_json(d / DONE_FILE)
    if done:
        return str(done.get("STATUS", "COMPLETE")), done
    if (d / FAILED_FILE).exists():
        return "FAILED", {"ERROR": (d / FAILED_FILE).read_text(encoding="utf-8", errors="replace").splitlines()[0][:100]}
    return "PENDING", {}


def status_mode(args, cfg) -> None:
    P = project_paths(args)
    studies, _ = step06_studies(P["S06"])
    order = step02_order(P)
    banner(f"{STEP['NAME']} STATUS - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    snap = read_tsv(snapshot_file(P), dtype=str)
    print(f"QTL contexts in the submitted snapshot: {len(snap)}" if len(snap) else
          "No job script made yet (run the planner).")
    rows, failures = [], []
    for accession in studies:
        counts = {"LOCI": 0, "DONE": 0, "NO_GWAS_CS": 0, "NO_CANDIDATES": 0, "FAILED": 0, "PENDING": 0,
                  "CANDIDATES": 0, "TESTED": 0, "STRONG": 0}
        for locus in list_loci(P["S06"] / accession):
            state, info = locus_state(P, accession, locus)
            counts["LOCI"] += 1
            key = {"COMPLETE": "DONE", "NO_GWAS_CREDIBLE_SET": "NO_GWAS_CS", "NO_QTL_CANDIDATES": "NO_CANDIDATES",
                   "FAILED": "FAILED"}.get(state, "PENDING")
            counts[key] += 1
            if key == "FAILED":
                failures.append(f"{accession} {locus}: {info.get('ERROR', '')}")
            for k in ("CANDIDATES", "TESTED", "STRONG"):
                counts[k] += int(info.get(f"N_{k}", 0) or 0)
        rows.append({"IDX": order.get(accession, ""), "STUDY": accession, **counts})
    if not rows:
        print(f"No fine-mapped studies under {P['S06']}.")
        return
    table = pd.DataFrame(rows)
    with pd.option_context("display.width", 200, "display.max_rows", None):
        print(table.to_string(index=False))
    total = table.drop(columns=["IDX", "STUDY"]).sum()
    print("\nLoci: " + "  ".join(f"{k}={int(total[k])}" for k in ("LOCI", "DONE", "NO_GWAS_CS", "NO_CANDIDATES",
                                                                       "FAILED", "PENDING")))
    for line in failures[:20]:
        print(f"  FAILED {line}")
    if int(total["FAILED"]) or int(total["PENDING"]):
        print(f"\nRe-run failed / pending loci: sbatch {P['BASH'].name}   (finished loci are skipped)")


def collect_results(P: dict, studies: list[str]) -> pd.DataFrame:
    frames = []
    for accession in studies:
        for locus in list_loci(P["S06"] / accession):
            f = P["OUT"] / accession / "loci" / locus / RESULTS_FILE
            if nonempty(f):
                x = read_tsv(f)
                if len(x):
                    frames.append(x)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=RESULT_COLUMNS)


def gene_label(row) -> str:
    gene = str(row.get("GENE_ID", "") or "")
    return gene if gene and gene.lower() != "nan" else str(row.get("MOLECULAR_TRAIT_ID", ""))


def stats_mode(args, cfg) -> None:
    P = project_paths(args)
    studies, _ = step06_studies(P["S06"])
    order = step02_order(P)
    banner(f"{STEP['NAME']} RESULTS - {STEP['TITLE']} - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    x = collect_results(P, studies)
    if x.empty:
        print("No colocalisation results yet (see --status).")
        return
    tested = x[x["COLOC_TESTED"].astype(str) == "YES"].copy()
    tested["PP_H4"] = numeric(tested["PP_H4"])
    strong, suggestive = cfg["STRONG_PP4"], cfg["SUGGESTIVE_PP4"]

    rows = []
    for accession in studies:
        sx, st = x[x["STUDY_ACCESSION"] == accession], tested[tested["STUDY_ACCESSION"] == accession]
        states = [locus_state(P, accession, locus)[0] for locus in list_loci(P["S06"] / accession)]
        best = st.sort_values("PP_H4", ascending=False).head(1)
        rows.append({
            "IDX": order.get(accession, ""), "STUDY": accession, "LOCI": len(states),
            "LOCI_DONE": sum(s in ("COMPLETE", "NO_GWAS_CREDIBLE_SET", "NO_QTL_CANDIDATES") for s in states),
            "CANDIDATES": len(sx), "CS_EXCLUDED": int(sx["STATUS"].astype(str).str.startswith("EXCLUDED").sum()),
            "TESTED": len(st), "STRONG": int((st["PP_H4"] >= strong).sum()),
            "SUGGESTIVE": int(((st["PP_H4"] >= suggestive) & (st["PP_H4"] < strong)).sum()),
            "PRIOR_ROBUST": int((st["PRIOR_ROBUST_STRONG"].astype(str) == "YES").sum()),
            "LOCI_WITH_STRONG": int(st.loc[st["PP_H4"] >= strong, "LOCUS_ID"].nunique()),
            "TOP (PP.H4)": (f"{gene_label(best.iloc[0])} "
                            f"{str(best.iloc[0]['TISSUE'])[:22]} ({best.iloc[0]['PP_H4']:.2f})") if len(best) else ""})
    table = pd.DataFrame(rows)
    print(f"\nPER STUDY (STRONG = PP.H4 >= {strong}; PRIOR_ROBUST = still strong at p12 = {cfg['P12_LOW']})")
    with pd.option_context("display.width", 250, "display.max_rows", None, "display.max_colwidth", 45):
        print(table.to_string(index=False))

    top = tested.sort_values("PP_H4", ascending=False).head(25)
    if len(top):
        print("\nTOP COLOCALISATIONS (best signal pair per molecular trait)")
        show = top.assign(IDX=top["STUDY_ACCESSION"].map(order).fillna(""),
                          GENE=[gene_label(r) for _, r in top.iterrows()], PP_H4=top["PP_H4"].round(3))
        cols = {"IDX": "IDX", "STUDY_ACCESSION": "STUDY", "LOCUS_ID": "LOCUS", "QTL_TYPE": "QTL",
                "STUDY_LABEL": "QTL_STUDY", "TISSUE": "TISSUE", "GENE": "GENE / TRAIT", "PP_H4": "PP.H4",
                "PRIOR_ROBUST_STRONG": "ROBUST", "COLOC_CLASS": "CLASS"}
        with pd.option_context("display.width", 250, "display.max_rows", None, "display.max_colwidth", 30):
            print(show[list(cols)].rename(columns=cols).to_string(index=False))
    counts = x["STATUS"].astype(str).value_counts()
    print("\nCandidate status: " + "  ".join(f"{k}={v}" for k, v in counts.items()))

    # Tables on disk: every candidate, per study, and one row per locus x QTL type (best trait).
    write_tsv(x, P["OUT"] / "coloc_all_candidates.tsv")
    for accession, g in x.groupby("STUDY_ACCESSION"):
        write_tsv(g, P["OUT"] / accession / f"{accession}_coloc_results.tsv")
    locus_rows = []
    for (accession, locus, qtl), g in tested.groupby(["STUDY_ACCESSION", "LOCUS_ID", "QTL_TYPE"]):
        b = g.sort_values("PP_H4", ascending=False).iloc[0]
        locus_rows.append({"STUDY_ACCESSION": accession, "LOCUS_ID": locus, "QTL_TYPE": qtl, "N_TESTED": len(g),
                           "N_STRONG": int((g["PP_H4"] >= strong).sum()), "BEST_PP_H4": b["PP_H4"],
                           "BEST_TRAIT": b["MOLECULAR_TRAIT_ID"], "BEST_GENE": b["GENE_ID"], "BEST_TISSUE": b["TISSUE"],
                           "BEST_DATASET": b["DATASET_ID"], "BEST_GWAS_VARIANT": b["BEST_GWAS_VARIANT"],
                           "BEST_QTL_VARIANT": b["BEST_QTL_VARIANT"], "PRIOR_ROBUST_STRONG": b["PRIOR_ROBUST_STRONG"]})
    write_tsv(pd.DataFrame(locus_rows), P["OUT"] / "coloc_locus_qtl_summary.tsv")
    write_tsv(table, P["OUT"] / "coloc_stats.tsv")
    print(f"\nTables: {P['OUT'] / 'coloc_stats.tsv'}\n        {P['OUT'] / 'coloc_all_candidates.tsv'}\n"
          f"        {P['OUT'] / 'coloc_locus_qtl_summary.tsv'}")


def arguments():
    parser = argparse.ArgumentParser(description=f"{STEP['NAME']}: {STEP['TITLE']} - provider-SuSiE coloc.bf_bf, "
                                                 "one SLURM task per fine-mapped locus.")
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--index", type=int, help="Study N of the Step02 download manifest: run its loci here.")
    modes.add_argument("--status", action="store_true", help="READ-ONLY: loci done / failed / pending per study.")
    modes.add_argument("--stats", action="store_true", help="Colocalisation results (tables on screen + TSV).")
    modes.add_argument("--check-resources", action="store_true",
                       help="Verify tools, QTL resources (dense, provider CS, indexed LBF) and Step06 input.")
    modes.add_argument("--worker-manifest", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--locus", default=None, help="With --index: only this locus (e.g. L001).")
    parser.add_argument("--study", default=None, help="Planner: only this GCST accession.")
    parser.add_argument("--reference-panel", default=DEFAULT_PANEL,
                        help="LD panel of Steps 04-06: 1000G (default), 1000G30x or custom; reads "
                             f"06_finemapping/<p>/<a>_<PANEL>/, writes {STEP['OUT_DIR']}/<p>/<a>_<PANEL>/.")
    parser.add_argument("--config", default=STEP["CONFIG"], help=f"Parameter file (default {STEP['CONFIG']}).")
    parser.add_argument("--qtl-types", default=None, help="Override QTL_TYPES, e.g. eQTL,sQTL.")
    parser.add_argument("--study-regex", default=None, help="Override STUDY_REGEX (eQTL Catalogue study_label).")
    parser.add_argument("--candidate-p", type=float, default=None, help="Override CANDIDATE_P.")
    parser.add_argument("--strong-pp4", type=float, default=None, help="Override STRONG_PP4.")
    parser.add_argument("--suggestive-pp4", type=float, default=None, help="Override SUGGESTIVE_PP4.")
    parser.add_argument("--force", action="store_true", help="Redo loci that are already done.")
    parser.add_argument("--keep-work", action="store_true", help="Keep the per-locus work/ files (LBF tables, R input).")
    slurm = parser.add_argument_group("SLURM (default: #SBATCH header of HPC_CPU.sh, then SLURM in the config)")
    slurm.add_argument("--template", default="HPC_CPU.sh")
    slurm.add_argument("--partition", default=None)
    slurm.add_argument("--qos", default=None)
    slurm.add_argument("--account", default=None)
    slurm.add_argument("--cpus", type=int, default=None)
    slurm.add_argument("--mem", "--memory", dest="mem", default=None)
    slurm.add_argument("--time", default=None)
    slurm.add_argument("--max-parallel", type=int, default=0, help="Locus tasks running at once (0 = no limit).")
    args = parser.parse_args()
    if args.locus and args.index is None:
        parser.error("--locus needs --index")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.reference_panel):
        parser.error(f"--reference-panel {args.reference_panel!r}: use letters, digits, '_', '-' or '.'")
    return args


def main() -> None:
    args = arguments()
    cfg = apply_overrides(load_config(config_path(args)), args)
    if args.check_resources:
        check_resources_mode(args, cfg)
    elif args.status:
        status_mode(args, cfg)
    elif args.stats:
        stats_mode(args, cfg)
    elif args.worker_manifest:
        worker_mode(args, cfg)
    elif args.index is not None:
        index_mode(args, cfg)
    else:
        planner_mode(args, cfg)


if __name__ == "__main__":
    main()
