#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 11 - External molecular-QTL colocalisation (coloc.abf), per locus
===============================================================================

    # 1. inventory resources/molecular_qtl/ -> config/Step11_Custom_QTL_Registry.tsv (every row DISABLED)
    #    + config/Step11_Column_Audit.tsv (columns / first line of every file). Re-scanning keeps your edits.
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --scan

    # 2. edit the registry: verify the column mapping, fill the QTL sample size n (or n = AN), set enabled = 1
    # 3. what is ready (registry eligibility, tools, Step06 input); exit 1 if nothing is ready
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources

    # 4. plan: ONE job script, one SLURM array task per fine-mapped locus (max 1000; more loci -> several per task)
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR
    sbatch Step11_CustomQTL_parkinson_s_disease_european.sh

    # study 2 on this machine: all its loci / one locus
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --index 2
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --index 2 --locus L001

    # progress and results; redo finished loci
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --stats
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --force

    # ---- 1000 Genomes 30x LD panel: reads 06_finemapping/<p>/<a>_1000G30x/, writes 11_custom_qtl_coloc/<p>/<a>_1000G30x/
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --check-resources
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x
    sbatch Step11_CustomQTL_parkinson_s_disease_european_1000G30x.sh
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --index 2 --locus L001
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --status
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

Why a separate step: Steps 09 / 10 use provider SuSiE fine-mapping (eQTL Catalogue LBF + credible sets).
The local molecular-QTL resources (resources/molecular_qtl/: ADSP FunGen brain pQTL / sQTL / haQTL /
mQTL, Kumasaka caQTL, ...) only have summary statistics, so they are colocalised here with
coloc::coloc.abf - independent of Steps 09 / 10, reported as its own evidence channel.

Per locus:
  1. GWAS side   the Step06 matched regional table (beta / SE, or P), N and case fraction from Step06
                 finemapping_summary.json (PHENOTYPE_TYPE Binary -> cc, Continuous -> quant).
  2. Discovery   every READY registry dataset, read by tabix for the locus window (only the files of the
                 locus chromosome); molecular traits with P <= CANDIDATE_P (top MAX_TRAITS).
  3. Matching    exact GRCh38 CHR:POS:REF:ALT (REF/ALT-swapped pairs too: ABF uses beta^2/var only);
                 >= MIN_SHARED shared variants and >= MIN_OVERLAP of the smaller side.
  4. coloc.abf   beta + varbeta (QTL: + sdY, or MAF + N) or P + MAF; priors p1 / p2 / p12 with
                 sensitivity at P12_LOW / P12_HIGH. STRONG = PP.H4 >= STRONG_PP4.

Safety: a disabled / incomplete registry row, missing N / sdY / MAF, headerless BED tracks without a column
mapping, or too little overlap is NOT_TESTED - never a negative colocalisation. coloc.abf assumes at most
one causal variant per trait and locus: results are exploratory until multi-signal validation.

Registry columns (config/Step11_Custom_QTL_Registry.tsv, one row per file; files of one dataset share
the dataset name and carry their chromosome):
  enabled dataset qtl_type context source path chrom build header variant_col trait_col p_col beta_col
  se_col maf_col ac_col an_col n sdy notes
  n = QTL sample size (number), or AN to use median(AN)/2 of the region rows; header = 1 named columns,
  0 headerless (columns are c1, c2, ... in the column audit).

Outputs (11_custom_qtl_coloc/<phenotype>/<ancestry>[_<PANEL>]/)
    <GCST>/loci/<L>/locus_coloc_results.tsv    every tested trait (same columns as Step09, COLOC_METHOD=ABF)
                    discovery_audit.tsv          per dataset: region rows, candidates, status
                    custom_coloc_done.json
    coloc_manifest.tsv, registry_snapshot.tsv (READY rows frozen when the job script was made)
    --stats: coloc_stats.tsv, coloc_all_candidates.tsv

Parameters: Step11_Config.py; SLURM from the HPC_CPU.sh header, then SLURM in the config, then
--partition --qos --account --cpus --mem --time --max-parallel.
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
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


VERSION = "1.0.0-custom-qtl-abf-per-locus"
PROJECT_ROOT = Path.cwd().resolve()
MAX_ARRAY_TASKS = 1000
DEFAULT_PANEL = "1000G"
OUT_DIR = "11_custom_qtl_coloc"
LOG_DIR = "step11_custom_qtl_coloc"
SCRIPT_PREFIX = "Step11_CustomQTL"
DONE_FILE = "custom_coloc_done.json"
FAILED_FILE = "LOCUS_FAILED.txt"
RESULTS_FILE = "locus_coloc_results.tsv"

DEFAULT_CONFIG: dict[str, Any] = {
    "RESOURCE_ROOT": "resources/molecular_qtl",
    "REGISTRY": "config/Step11_Custom_QTL_Registry.tsv",
    "QTL_TYPES": ["sQTL", "pQTL", "haQTL", "caQTL", "mQTL"],
    "CANDIDATE_P": 1e-5,
    "MAX_TRAITS": 100,
    "MIN_SHARED": 50,
    "MIN_OVERLAP": 0.5,
    "P1": 1e-4,
    "P2": 1e-4,
    "P12": 5e-6,
    "P12_LOW": 1e-6,
    "P12_HIGH": 1e-5,
    "STRONG_PP4": 0.80,
    "SUGGESTIVE_PP4": 0.50,
    "MAX_REGION_ROWS": 300_000,
    "R_TIMEOUT": 3600,
    "SLURM": {},
}

REGISTRY = ["enabled", "dataset", "qtl_type", "context", "source", "path", "chrom", "build", "header",
            "variant_col", "trait_col", "p_col", "beta_col", "se_col", "maf_col", "ac_col", "an_col", "n", "sdy",
            "notes"]
ALIAS = {
    "variant_col": ("variant", "variantId", "variant_id", "snp", "SNPID", "variant_key"),
    "trait_col": ("molecular_trait_id", "target", "target_id", "phenotype_id", "gene_id", "feature_id"),
    "p_col": ("pvalue", "p_value", "p", "pval", "P", "pval_nominal"),
    "beta_col": ("beta", "BETA", "slope", "effect_size"),
    "se_col": ("se", "SE", "standard_error", "slope_se"),
    "maf_col": ("maf", "MAF", "minor_allele_frequency"),
    "ac_col": ("ac", "AC"),
    "an_col": ("an", "AN"),
}
QTL_TYPE_NAMES = {"caqtl": "caQTL", "haqtl": "haQTL", "mqtl": "mQTL", "meqtl": "mQTL", "pqtl": "pQTL",
                  "sqtl": "sQTL", "eqtl": "eQTL"}
VARRE = re.compile(r"^(?:chr)?([0-9]{1,2})[_:]([0-9]+)[_:]([ACGTN]+)[_:]([ACGTN]+)$", re.I)
ADSP_NAME = re.compile(r"^FunGenADv1\.(?P<study>[^.]+)\.(?P<qtl>[^.]+)\.(?P<context>[^.]+)\.(?P<id>[A-Za-z0-9]+)_v1_"
                       r"(?P<chrom>[0-9]{1,2})_(?P<kind>snp|indel)_all\.[0-9]+\.bed\.gz$")
ENSG = re.compile(r"ENSG\d{11}")

ANCESTRY_ALIASES = {
    "eur": ("EUR", "European"), "european": ("EUR", "European"),
    "afr": ("AFR", "African"), "african": ("AFR", "African"),
    "eas": ("EAS", "East Asian"), "east asian": ("EAS", "East Asian"),
    "sas": ("SAS", "South Asian"), "south asian": ("SAS", "South Asian"),
    "amr": ("AMR", "Hispanic or Latin American"), "hispanic": ("AMR", "Hispanic or Latin American"),
    "latino": ("AMR", "Hispanic or Latin American"),
    "hispanic or latin american": ("AMR", "Hispanic or Latin American"),
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
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", str(value).lower())).strip()


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


def safe_float(value: Any) -> float:
    try:
        x = float(value)
        return x if np.isfinite(x) else np.nan
    except Exception:
        return np.nan


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


def variant_key(value: Any) -> str:
    m = VARRE.match(str(value).strip())
    return f"{int(m[1])}:{int(m[2])}:{m[3].upper()}:{m[4].upper()}" if m and 1 <= int(m[1]) <= 22 else ""


def swapped(key: str) -> str:
    c, p, r, a = key.split(":")
    return f"{c}:{p}:{a}:{r}"


def project_path(value: str) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else PROJECT_ROOT / path


# =============================================================================
# CONFIGURATION, PATHS, SLURM
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
    cfg["QTL_TYPES"] = [QTL_TYPE_NAMES.get(t.lower(), t) for t in cfg["QTL_TYPES"]]
    return cfg


def results_slug(ancestry_slug: str, panel: str) -> str:
    return ancestry_slug if panel == DEFAULT_PANEL else f"{ancestry_slug}_{panel}"


def panel_option(panel: str) -> str:
    return "" if panel == DEFAULT_PANEL else f" --reference-panel {panel}"


def project_paths(args) -> dict:
    code, label = canonical_ancestry(args.ancestry)
    p, a = slugify(args.phenotype), slugify(label)
    r = results_slug(a, args.reference_panel)
    return {"phenotype": args.phenotype.strip(), "code": code, "label": label, "p": p, "a": a, "r": r,
            "panel": args.reference_panel, "S02": PROJECT_ROOT / "01_gwas_catalog" / p / a,
            "S06": PROJECT_ROOT / "06_finemapping" / p / r, "OUT": PROJECT_ROOT / OUT_DIR / p / r,
            "LOG": PROJECT_ROOT / "logs" / LOG_DIR / p / r, "BASH": PROJECT_ROOT / f"{SCRIPT_PREFIX}_{p}_{r}.sh"}


def base_command(P: dict) -> str:
    return (f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {shlex.quote(P['phenotype'])} "
            f"--ancestry {P['code']}{panel_option(P['panel'])}")


def resolve_slurm(args, cfg: dict) -> dict:
    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": "2", "mem": "16G", "time": "12:00:00"}
    template = PROJECT_ROOT / args.template
    if template.is_file():
        for line in template.read_text(encoding="utf-8").splitlines():
            match = re.match(r"^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
            if match and (match.group(2) or "").strip() and match.group(1) not in {"job-name", "output", "error", "array"}:
                values[match.group(1)] = match.group(2).strip()
    for key, value in (cfg.get("SLURM") or {}).items():
        if value not in (None, ""):
            values[str(key)] = str(value)
    for key, value in (("partition", args.partition), ("qos", args.qos), ("account", args.account),
                       ("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time)):
        if value not in (None, ""):
            values[key] = str(value)
    return values


def step02_order(P: dict) -> dict[str, int]:
    path = P["S02"] / "GWAS_download_manifest.tsv"
    if not path.is_file():
        return {}
    m = pd.read_csv(path, sep="\t", dtype=str)
    ids = pd.to_numeric(m.get("ARRAY_ID"), errors="coerce")
    return {str(acc).strip(): int(i) for acc, i in zip(m["STUDY_ACCESSION"], ids) if pd.notna(i)}


def study_for_index(P: dict, index: int) -> str:
    for accession, i in step02_order(P).items():
        if i == index:
            return accession
    raise SystemExit(f"--index {index} is not in the Step02 download manifest under {P['S02']}")


# =============================================================================
# STEP06 INPUT (GWAS side)
# =============================================================================

def list_loci(study_dir: Path) -> list[str]:
    loci_dir = study_dir / "loci"
    if not loci_dir.exists():
        return []
    return sorted(d.name for d in loci_dir.iterdir() if d.is_dir() and re.fullmatch(r"L\d+", d.name)
                  and any(nonempty(d / f) for f in ("reference_qc_matched_variants.tsv.gz", "matched_variants.tsv")))


def step06_studies(S06: Path) -> list[str]:
    if not S06.exists():
        return []
    return sorted(d.name for d in S06.iterdir() if d.is_dir() and d.name.startswith("GCST")
                  and nonempty(d / "finemapping_summary.json") and list_loci(d))


def gwas_meta(study_dir: Path) -> dict:
    """coloc dataset type / N / case fraction of the GWAS from Step06 finemapping_summary.json."""
    s = read_json(study_dir / "finemapping_summary.json")
    kind = str(s.get("PHENOTYPE_TYPE", "") or "")
    cases, controls = safe_float(s.get("N_CASES")), safe_float(s.get("N_CONTROLS"))
    n = safe_float(s.get("GWAS_N"))
    if not np.isfinite(n):
        n = safe_float(s.get("N_TOTAL"))
    if kind == "Binary" and np.isfinite(cases) and np.isfinite(controls) and cases + controls > 0:
        return {"TYPE": "cc", "N": cases + controls if not np.isfinite(n) else n,
                "S": cases / (cases + controls), "SOURCE": "Step06 N_CASES / N_CONTROLS"}
    if kind == "Continuous" and np.isfinite(n):
        return {"TYPE": "quant", "N": n, "S": np.nan, "SOURCE": "Step06 GWAS_N"}
    return {"TYPE": "", "N": n, "S": np.nan, "SOURCE": f"PHENOTYPE_TYPE {kind or 'unknown'} / N missing"}


def load_gwas(study_dir: Path, locus: str) -> pd.DataFrame:
    base = study_dir / "loci" / locus
    src = next((p for p in (base / "reference_qc_matched_variants.tsv.gz", base / "matched_variants.tsv")
                if nonempty(p)), None)
    if src is None:
        raise FileNotFoundError(f"Missing Step06 matched table under {base}")
    g = read_tsv(src)

    def get(*names):
        for c in names:
            if c in g.columns:
                return g[c]
        return pd.Series(np.nan, index=g.index)

    chrom = get("REFERENCE_CHR", "CHR").astype(str).str.replace("chr", "", case=False)
    pos = pd.to_numeric(get("REFERENCE_POS", "POS"), errors="coerce")
    keys = [variant_key(f"{c}:{int(x)}:{r}:{a}") if pd.notna(x) else ""
            for c, x, r, a in zip(chrom, pos, get("REFERENCE_REF"), get("REFERENCE_ALT"))]
    out = pd.DataFrame({"SNP": keys, "CHR": pd.to_numeric(chrom, errors="coerce"), "POS": pos,
                        "P": pd.to_numeric(get("P"), errors="coerce"),
                        "B": pd.to_numeric(get("BETA_REF", "BETA"), errors="coerce"),
                        "SE": pd.to_numeric(get("SE"), errors="coerce"),
                        "MAF": pd.to_numeric(get("MAF", "EAF"), errors="coerce")})
    out["MAF"] = np.minimum(out["MAF"], 1 - out["MAF"])
    out = out[(out.SNP != "") & out.P.gt(0) & out.P.le(1)].sort_values("P").drop_duplicates("SNP")
    if out.empty:
        raise RuntimeError(f"No valid GWAS variants in {src}")
    return out.reset_index(drop=True)


# =============================================================================
# REGISTRY: --scan inventory and eligibility
# =============================================================================

def text_header(path: Path) -> tuple[list[str], int, str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        for _ in range(15):
            line = fh.readline()
            if not line:
                break
            if not line.strip() or line.startswith("##") or line.startswith("track "):
                continue
            cells = line.rstrip("\n\r").lstrip("#").split("\t")
            named = any(z.lower() in {"variant", "variantid", "variant_id", "molecular_trait_id", "chromstart",
                                      "pvalue", "beta", "chrom", "p_value", "snp"} for z in cells)
            return (cells if named else [f"c{i + 1}" for i in range(len(cells))], int(named), line[:240].strip())
    return [], 0, ""


def first_col(cols: Iterable[str], names: Iterable[str]) -> str:
    bylow = {c.lower(): c for c in cols}
    return next((bylow[s.lower()] for s in names if s.lower() in bylow), "")


def describe_file(path: Path) -> dict:
    """dataset / qtl type / context / source / chromosome from the folder layout and file name."""
    parts = [p.lower() for p in path.parts]
    kind = next((QTL_TYPE_NAMES[x] for x in ("caqtl", "haqtl", "mqtl", "pqtl", "sqtl", "eqtl") if x in parts), "otherQTL")
    m = ADSP_NAME.match(path.name)
    if m:
        return {"dataset": f"ADSP_{m['study']}_{m['context']}_{kind}_{m['id']}", "qtl_type": kind,
                "context": m["context"], "source": f"ADSP FunGen {m['study']}", "chrom": m["chrom"]}
    stem = re.sub(r"\.(all\.)?(tsv|bed)\.gz$", "", path.name)
    source = path.parent.name
    return {"dataset": f"{source}_{stem}", "qtl_type": kind, "context": "unknown", "source": source, "chrom": ""}


def scan_mode(args, cfg: dict) -> None:
    root = project_path(cfg["RESOURCE_ROOT"])
    registry = project_path(cfg["REGISTRY"])
    files = sorted(p for p in root.rglob("*") if p.is_file()
                   and (p.name.endswith(".all.tsv.gz") or re.search(r"_all\.[0-9]+\.bed\.gz$", p.name)
                        or p.name.endswith(".all.bed.gz")))
    banner(f"STEP 11 - SCAN {root}")
    if not files:
        raise SystemExit(f"No association files (*.all.tsv.gz / *_all.<date>.bed.gz) under {root}\n"
                         "Download them with: envs/pipeline/bin/python Step00_Download_Resources.py "
                         "--groups molecular_qtl meqtl")
    previous = read_tsv(registry, dtype=str, keep_default_na=False) if registry.exists() and not args.overwrite \
        else pd.DataFrame(columns=REGISTRY)
    kept = {str(r["path"]): r for r in previous.to_dict(orient="records")}
    rows, audits = [], []
    for i, path in enumerate(files, 1):
        rel = str(path.relative_to(PROJECT_ROOT)) if path.is_relative_to(PROJECT_ROOT) else str(path)
        cols, header, sample = text_header(path)
        info = describe_file(path)
        if rel in kept:                     # keep your manual edits
            row = {k: kept[rel].get(k, "") for k in REGISTRY}
        else:
            row = {k: "" for k in REGISTRY}
            row.update(enabled="0", path=rel, build="GRCh38", header=str(header), **info)
            for k, names in ALIAS.items():
                row[k] = first_col(cols, names)
            notes = ["indexed" if Path(str(path) + ".tbi").exists() else "NO_TABIX_INDEX"]
            if not header:
                notes.append("HEADERLESS_NEEDS_SCHEMA (map c1, c2, ... from the column audit)")
            if row["an_col"]:
                notes.append("n can be AN (median AN/2 of the region)")
            if re.search(r"QTD0005[0-9]{2}", path.name) or "sun2018" in rel.lower():
                notes.append("ALSO IN STEP09 (eQTL Catalogue Sun_2018) - enabling it counts it twice")
            row["notes"] = "; ".join(notes)
        rows.append(row)
        audits.append({"dataset": row["dataset"], "path": rel, "qtl_type": row["qtl_type"], "chrom": row["chrom"],
                       "has_tbi": Path(str(path) + ".tbi").exists(), "header": bool(header),
                       "columns": ",".join(cols), "first_line": sample})
        if i % 200 == 0:
            print(f"  scanned {i}/{len(files)}", flush=True)
    frame = pd.DataFrame(rows, columns=REGISTRY)
    registry.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(registry, sep="\t", index=False)
    audit = registry.with_name("Step11_Column_Audit.tsv")
    pd.DataFrame(audits).to_csv(audit, sep="\t", index=False)
    new = sum(1 for r in rows if r["path"] not in kept)
    print(f"Files            : {len(rows)} ({new} new, {len(rows) - new} kept with your edits)")
    print(f"Datasets         : {frame['dataset'].nunique()}  by type: {dict(Counter(frame.drop_duplicates('dataset')['qtl_type']))}")
    print(f"Registry         : {registry}")
    print(f"Column audit     : {audit}")
    print("\nEvery new row is DISABLED. For each dataset: check the columns in the audit (headerless BED files:\n"
          "map c1, c2, ...), fill n (QTL sample size, or AN), then set enabled = 1. Then:\n"
          f"  {base_command(project_paths(args))} --check-resources")


def eligibility(row: dict) -> str:
    if str(row.get("enabled", "")).strip().lower() not in {"1", "true", "yes"}:
        return "DISABLED"
    path = project_path(row["path"])
    if not path.is_file():
        return "FILE_MISSING"
    if not Path(str(path) + ".tbi").exists():
        return "INDEX_MISSING"
    if str(row.get("build", "")).strip().upper() not in {"GRCH38", "HG38"}:
        return "NOT_GRCH38"
    if any(not str(row.get(c, "")).strip() for c in ("variant_col", "trait_col", "p_col")):
        return "NEEDS_COLUMN_MAPPING"
    n = str(row.get("n", "")).strip()
    if n.upper() == "AN":
        if not str(row.get("an_col", "")).strip():
            return "N_IS_AN_BUT_NO_AN_COLUMN"
    elif not np.isfinite(safe_float(n)) or safe_float(n) <= 0:
        return "NEEDS_QTL_N"
    cols, header, _ = text_header(path)
    if str(row.get("header", "")).strip() != str(header):
        return "HEADER_FLAG_MISMATCH"
    mapped = [row.get(k, "") for k in ("variant_col", "trait_col", "p_col", "beta_col", "se_col", "maf_col",
                                         "ac_col", "an_col") if str(row.get(k, "")).strip()]
    if any(c not in cols for c in mapped):
        return "COLUMN_NOT_IN_FILE"
    if not ((str(row.get("beta_col", "")).strip() and str(row.get("se_col", "")).strip())
            or str(row.get("maf_col", "")).strip()
            or (str(row.get("ac_col", "")).strip() and str(row.get("an_col", "")).strip())):
        return "NEEDS_BETA_SE_OR_MAF"
    return "READY"


def load_registry(cfg: dict) -> pd.DataFrame:
    path = project_path(cfg["REGISTRY"])
    if not path.exists():
        return pd.DataFrame(columns=REGISTRY + ["ELIGIBILITY"])
    reg = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = set(REGISTRY) - set(reg.columns)
    if missing:
        raise SystemExit(f"Registry {path} lacks columns {sorted(missing)} - run --scan again")
    reg["qtl_type"] = reg["qtl_type"].map(lambda t: QTL_TYPE_NAMES.get(str(t).lower(), str(t)))
    reg["ELIGIBILITY"] = [eligibility(r) for r in reg.to_dict(orient="records")]
    wanted = {t.lower() for t in cfg["QTL_TYPES"]}
    reg["TYPE_SELECTED"] = reg["qtl_type"].str.lower().isin(wanted)
    return reg


def registry_hash(ready: pd.DataFrame) -> str:
    cols = [c for c in REGISTRY if c != "notes"]
    return sha1_short("\n".join("|".join(str(r[c]) for c in cols) for r in
                                ready.sort_values(["dataset", "path"]).to_dict(orient="records"))) if len(ready) else "none"


def print_registry_summary(reg: pd.DataFrame) -> None:
    if reg.empty:
        print("Registry: none yet - run --scan")
        return
    per = reg.groupby(["dataset", "qtl_type"]).agg(
        FILES=("path", "size"), READY=("ELIGIBILITY", lambda s: int((s == "READY").sum())),
        STATE=("ELIGIBILITY", lambda s: ", ".join(f"{k} {v}" for k, v in s.value_counts().items()))).reset_index()
    with pd.option_context("display.width", 220, "display.max_rows", None, "display.max_colwidth", 70):
        print(per.to_string(index=False))
    print("\nFiles: " + "  ".join(f"{k}={v}" for k, v in reg["ELIGIBILITY"].value_counts().items()))


# =============================================================================
# ONE LOCUS
# =============================================================================

def read_region(path: Path, names: list[str], chrom: int, start: int, end: int, max_rows: int) -> pd.DataFrame:
    import pysam
    rows = []
    with pysam.TabixFile(str(path)) as tab:
        contigs = set(tab.contigs)
        contig = next((z for z in (str(chrom), f"chr{chrom}") if z in contigs), None)
        if contig is None:
            return pd.DataFrame(columns=names)
        for raw in tab.fetch(contig, max(start - 1, 0), end):
            if raw.startswith("#"):
                continue
            v = raw.split("\t")
            if len(v) != len(names):
                continue
            rows.append(v)
            if len(rows) > max_rows:
                raise RuntimeError(f"TOO_MANY_REGION_ROWS (> {max_rows}); raise MAX_REGION_ROWS in Step11_Config.py")
    return pd.DataFrame(rows, columns=names)


def normalize_qtl(df: pd.DataFrame, r: dict) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["SNP", "TRAIT", "P", "B", "SE", "MAF", "AN"])

    def field(k):
        c = str(r.get(k, "")).strip()
        return df[c] if c and c in df.columns else pd.Series(np.nan, index=df.index)

    out = pd.DataFrame({"SNP": field("variant_col").map(variant_key), "TRAIT": field("trait_col").astype(str),
                        "P": pd.to_numeric(field("p_col"), errors="coerce"),
                        "B": pd.to_numeric(field("beta_col"), errors="coerce"),
                        "SE": pd.to_numeric(field("se_col"), errors="coerce"),
                        "MAF": pd.to_numeric(field("maf_col"), errors="coerce"),
                        "AN": pd.to_numeric(field("an_col"), errors="coerce")})
    if not str(r.get("maf_col", "")).strip() and str(r.get("ac_col", "")).strip() and str(r.get("an_col", "")).strip():
        ac = pd.to_numeric(field("ac_col"), errors="coerce")
        out["MAF"] = np.minimum(ac / out["AN"], 1 - ac / out["AN"])
    out["MAF"] = np.minimum(out["MAF"], 1 - out["MAF"])
    out = out[(out.SNP != "") & out.TRAIT.ne("") & out.TRAIT.ne("nan") & out.P.gt(0) & out.P.le(1)]
    return out.sort_values("P").drop_duplicates(["SNP", "TRAIT"]).reset_index(drop=True)


R_ABF = r'''
suppressPackageStartupMessages(library(coloc))
a <- commandArgs(trailingOnly=TRUE)
x <- read.delim(gzfile(a[1]), stringsAsFactors=FALSE, check.names=FALSE)
cands <- read.delim(a[2], stringsAsFactors=FALSE, check.names=FALSE)
out_file <- a[3]
p1 <- as.numeric(a[4]); p2 <- as.numeric(a[5]); p12s <- as.numeric(c(a[6], a[7], a[8]))
gtype <- a[9]; gs <- as.numeric(a[10]); gn <- as.numeric(a[11])
ok_maf <- function(m) all(is.finite(m) & m > 0 & m <= 0.5)
rows <- list()
for (i in seq_len(nrow(cands))) {
  cid <- cands$CANDIDATE_ID[i]
  r <- list(CANDIDATE_ID=cid, R_STATUS="OK", MESSAGE="", PP_H0=NA, PP_H1=NA, PP_H2=NA, PP_H3=NA, PP_H4=NA,
            PP_H4_P12_LOW=NA, PP_H4_P12_HIGH=NA, BEST_SNP="", BEST_SNP_PP_H4=NA, N_SNPS=0, GWAS_MODE="", QTL_MODE="")
  res <- tryCatch({
    d <- x[x$CANDIDATE_ID == cid, , drop=FALSE]
    d <- d[!duplicated(d$SNP), , drop=FALSE]
    r$N_SNPS <- nrow(d)
    nq <- as.numeric(cands$QTL_N[i]); sy <- suppressWarnings(as.numeric(cands$QTL_SDY[i]))
    if (!is.finite(nq) || nq <= 0) stop("QTL sample size missing")
    d1 <- list(type=gtype, N=gn, snp=d$SNP)
    if (gtype == "cc") d1$s <- gs
    if (all(is.finite(d$B_GWAS) & is.finite(d$SE_GWAS) & d$SE_GWAS > 0)) {
      d1$beta <- d$B_GWAS; d1$varbeta <- d$SE_GWAS^2; r$GWAS_MODE <- "beta"
      if (gtype == "quant") { if (ok_maf(d$MAF)) d1$MAF <- d$MAF else stop("quantitative GWAS needs MAF") }
    } else if (ok_maf(d$MAF)) { d1$pvalues <- d$P_GWAS; d1$MAF <- d$MAF; r$GWAS_MODE <- "pvalue+MAF"
    } else stop("GWAS needs complete beta/SE or p-values + MAF")
    d2 <- list(type="quant", N=nq, snp=d$SNP)
    if (all(is.finite(d$B_QTL) & is.finite(d$SE_QTL) & d$SE_QTL > 0)) {
      d2$beta <- d$B_QTL; d2$varbeta <- d$SE_QTL^2; r$QTL_MODE <- "beta"
      if (is.finite(sy) && sy > 0) d2$sdY <- sy
      else if (ok_maf(d$MAF)) d2$MAF <- d$MAF
      else stop("QTL beta/SE needs MAF (with N) or a known sdY")
    } else if (ok_maf(d$MAF)) { d2$pvalues <- d$P_QTL; d2$MAF <- d$MAF; r$QTL_MODE <- "pvalue+MAF"
    } else stop("QTL needs beta/SE + MAF or sdY, or p-values + MAF")
    z <- lapply(p12s, function(p12) suppressWarnings(coloc::coloc.abf(d1, d2, p1=p1, p2=p2, p12=p12)))
    pp <- vapply(z, function(o) as.numeric(o$summary[paste0("PP.H", 0:4, ".abf")]), numeric(5))
    r$PP_H0 <- pp[1, 1]; r$PP_H1 <- pp[2, 1]; r$PP_H2 <- pp[3, 1]; r$PP_H3 <- pp[4, 1]; r$PP_H4 <- pp[5, 1]
    r$PP_H4_P12_LOW <- pp[5, 2]; r$PP_H4_P12_HIGH <- pp[5, 3]
    rr <- z[[1]]$results
    if (!is.null(rr) && "SNP.PP.H4" %in% names(rr)) {
      j <- which.max(rr$SNP.PP.H4); r$BEST_SNP <- as.character(rr$snp[j]); r$BEST_SNP_PP_H4 <- rr$SNP.PP.H4[j]
    }
    r
  }, error=function(e) { r$R_STATUS <- "ERROR"; r$MESSAGE <- gsub("[\t\n]", " ", conditionMessage(e)); r })
  rows[[length(rows) + 1]] <- as.data.frame(res, stringsAsFactors=FALSE)
}
write.table(do.call(rbind, rows), out_file, sep="\t", quote=FALSE, row.names=FALSE)
'''

RESULT_COLUMNS = [
    "STUDY_ACCESSION", "LOCUS_ID", "CHR", "LOCUS_START", "LOCUS_END", "QTL_TYPE", "STUDY_LABEL", "TISSUE",
    "SAMPLE_GROUP", "CONDITION", "DATASET_ID", "MOLECULAR_TRAIT_ID", "GENE_ID", "QTL_N", "DISCOVERY_MIN_P",
    "N_DISCOVERY_ROWS", "GWAS_CS_INDEXES", "QTL_CS_INDEXES", "STATUS", "COLOC_TESTED", "PP_H4", "COLOC_CLASS",
    "COLOCALIZED_STRONG", "SHARED_SIGNAL", "PRIOR_ROBUST_STRONG", "PP_H4_P12_LOW", "PP_H4_P12_HIGH",
    "PP_H0", "PP_H1", "PP_H2", "PP_H3", "DOMINANT_HYPOTHESIS", "BEST_GWAS_VARIANT", "BEST_QTL_VARIANT",
    "BEST_SNP_PP_H4", "N_SIGNAL_PAIRS", "N_GWAS_SIGNALS", "N_QTL_SIGNALS", "N_GWAS_VARIANTS", "N_QTL_VARIANTS",
    "N_COMMON_VARIANTS", "GWAS_MODE", "QTL_MODE", "COLOC_METHOD", "ERROR", "CANDIDATE_ID",
]


def classify(pp: list[float], strong: float, suggestive: float) -> dict:
    dominant = int(np.nanargmax(pp))
    pp4 = pp[4]
    label = ("STRONG_SHARED_SIGNAL" if pp4 >= strong else "SUGGESTIVE_SHARED_SIGNAL" if pp4 >= suggestive else
             {0: "NEITHER_ASSOCIATED_FAVORED", 1: "GWAS_ONLY_FAVORED", 2: "QTL_ONLY_FAVORED",
              3: "BOTH_ASSOCIATED_DIFFERENT_SIGNALS_FAVORED", 4: "SHARED_SIGNAL_BELOW_THRESHOLD"}[dominant])
    return {"COLOC_CLASS": label, "DOMINANT_HYPOTHESIS": f"H{dominant}",
            "COLOCALIZED_STRONG": "YES" if pp4 >= strong else "NO", "SHARED_SIGNAL": "YES" if pp4 >= suggestive else "NO"}


def process_locus(row: pd.Series, args, cfg: dict, ready: pd.DataFrame, rscript: str) -> str:
    P = project_paths(args)
    accession, locus = str(row["STUDY_ACCESSION"]), str(row["LOCUS_ID"])
    study_dir = P["S06"] / accession
    out = P["OUT"] / accession / "loci" / locus
    done_file = out / DONE_FILE
    r_hash = registry_hash(ready)
    previous = read_json(done_file)
    if previous and not args.force and previous.get("REGISTRY_HASH") == r_hash:
        print(f"[RESUME] {accession} {locus}: {previous.get('STATUS')} (same registry)")
        return str(previous.get("STATUS"))
    out.mkdir(parents=True, exist_ok=True)
    (out / FAILED_FILE).unlink(missing_ok=True)
    work = out / "work"
    work.mkdir(exist_ok=True)
    started = time.time()
    banner(f"STEP 11 - {accession} {locus} - {ready['dataset'].nunique()} QTL datasets (coloc.abf)")

    def finish(status: str, results: pd.DataFrame, audit: list[dict], extra: dict) -> str:
        write_tsv(results.reindex(columns=RESULT_COLUMNS), out / RESULTS_FILE)
        write_tsv(pd.DataFrame(audit), out / "discovery_audit.tsv")
        tested = results[results.get("COLOC_TESTED", pd.Series(dtype=str)).astype(str) == "YES"] if len(results) else results
        h4 = pd.to_numeric(tested.get("PP_H4", pd.Series(dtype=float)), errors="coerce")
        write_json(done_file, {"VERSION": VERSION, "STEP": "Step11", "STATUS": status, "STUDY_ACCESSION": accession,
                               "LOCUS_ID": locus, "REGISTRY_HASH": r_hash, "N_DATASETS": int(ready["dataset"].nunique()),
                               "N_CANDIDATES": int(len(results)), "N_TESTED": int(len(tested)),
                               "N_STRONG": int((h4 >= cfg["STRONG_PP4"]).sum()),
                               "SECONDS": round(time.time() - started, 1), "COMPLETED_UTC": utcnow(), **extra})
        if not args.keep_work:
            shutil.rmtree(work, ignore_errors=True)
        print(f"[OK] {accession} {locus}: {status}; {len(results)} candidates, {len(tested)} tested "
              f"({time.time() - started:.0f}s)")
        return status

    gwas = load_gwas(study_dir, locus)
    meta = gwas_meta(study_dir)
    chroms = gwas.CHR.dropna().unique()
    if len(chroms) != 1:
        raise RuntimeError(f"GWAS locus spans chromosomes {sorted(chroms)}")
    chrom, start, end = int(chroms[0]), int(gwas.POS.min()), int(gwas.POS.max())
    region = {"CHR": chrom, "LOCUS_START": start, "LOCUS_END": end, "GWAS_TYPE": meta["TYPE"],
              "GWAS_N": meta["N"], "GWAS_CASE_FRACTION": meta["S"], "GWAS_META_SOURCE": meta["SOURCE"]}
    if not meta["TYPE"] or not np.isfinite(meta["N"]):
        return finish("NOT_TESTED_GWAS_N_OR_TYPE_MISSING", pd.DataFrame(), [], region)

    gw = gwas.set_index("SNP")
    long_rows, cands, audit = [], [], []
    for dataset, files in ready.groupby("dataset", sort=True):
        files = files[files["chrom"].astype(str).isin(["", str(chrom)])]
        first = files.iloc[0].to_dict() if len(files) else ready[ready["dataset"] == dataset].iloc[0].to_dict()
        entry = {"DATASET_ID": dataset, "QTL_TYPE": first["qtl_type"], "CONTEXT": first["context"], "N_FILES": len(files),
                 "N_REGION_ROWS": 0, "N_CANDIDATES": 0, "STATUS": "OK"}
        try:
            frames = []
            for r in files.to_dict(orient="records"):
                path = project_path(r["path"])
                cols, _, _ = text_header(path)
                raw = read_region(path, cols, chrom, start, end, int(cfg["MAX_REGION_ROWS"]))
                entry["N_REGION_ROWS"] += len(raw)
                frames.append(normalize_qtl(raw, r))
            q = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
            if q.empty:
                entry["STATUS"] = "NO_REGION_ROWS" if len(files) else "NO_FILE_FOR_CHROMOSOME"
                audit.append(entry)
                continue
            best = q.groupby("TRAIT")["P"].min().sort_values()
            traits = best[best <= cfg["CANDIDATE_P"]].head(int(cfg["MAX_TRAITS"])).index.tolist()
            if not traits:
                entry["STATUS"] = "NO_TRAIT_P_LE_THRESHOLD"
                audit.append(entry)
                continue
            entry["N_CANDIDATES"] = len(traits)
            for trait in traits:
                qt = q[q.TRAIT == trait].copy()
                # exact ref/alt; REF/ALT-swapped pairs too (coloc.abf uses beta^2 / varbeta, the sign is irrelevant)
                direct = qt.SNP.isin(gw.index)
                qt.loc[~direct, "SNP"] = qt.loc[~direct, "SNP"].map(swapped)
                m = qt[qt.SNP.isin(gw.index)].drop_duplicates("SNP")
                cid = sha1_short(f"{accession}|{locus}|{dataset}|{trait}", 20)
                an = pd.to_numeric(qt["AN"], errors="coerce")
                n_field = str(first.get("n", "")).strip()
                qtl_n = float(np.nanmedian(an) / 2) if n_field.upper() == "AN" and an.notna().any() else safe_float(n_field)
                gene = ENSG.search(trait)
                base = {"CANDIDATE_ID": cid, "STUDY_ACCESSION": accession, "LOCUS_ID": locus, **region,
                        "QTL_TYPE": first["qtl_type"], "STUDY_LABEL": first["source"], "TISSUE": first["context"],
                        "SAMPLE_GROUP": "", "CONDITION": "", "DATASET_ID": dataset, "MOLECULAR_TRAIT_ID": trait,
                        "GENE_ID": gene.group(0) if gene else "", "QTL_N": qtl_n, "DISCOVERY_MIN_P": float(best[trait]),
                        "N_DISCOVERY_ROWS": int(len(qt)), "N_GWAS_VARIANTS": int(len(gwas)),
                        "N_QTL_VARIANTS": int(len(qt)), "N_COMMON_VARIANTS": int(len(m)),
                        "N_SIGNAL_PAIRS": 1, "N_GWAS_SIGNALS": 1, "N_QTL_SIGNALS": 1,
                        "COLOC_METHOD": "coloc.abf (single causal variant)", "COLOC_TESTED": "NO"}
                if len(m) < int(cfg["MIN_SHARED"]) or len(m) / max(1, min(len(gwas), len(qt))) < float(cfg["MIN_OVERLAP"]):
                    cands.append({**base, "STATUS": "NOT_TESTED_INSUFFICIENT_VARIANT_OVERLAP"})
                    continue
                g = gw.loc[m.SNP]
                maf = m["MAF"].to_numpy(dtype=float)
                maf = np.where(np.isfinite(maf), maf, g["MAF"].to_numpy(dtype=float))   # QTL MAF first, GWAS MAF fallback
                long_rows.append(pd.DataFrame({"CANDIDATE_ID": cid, "SNP": m.SNP.values, "P_GWAS": g["P"].values,
                                               "B_GWAS": g["B"].values, "SE_GWAS": g["SE"].values, "P_QTL": m.P.values,
                                               "B_QTL": m.B.values, "SE_QTL": m.SE.values, "MAF": maf}))
                cands.append({**base, "STATUS": "PENDING_R", "QTL_SDY": first.get("sdy", "") or "NA"})
        except Exception as exc:
            entry["STATUS"] = f"NOT_TESTED_RESOURCE_ERROR: {type(exc).__name__}: {str(exc)[:200]}"
        audit.append(entry)

    results = pd.DataFrame(cands)
    if results.empty:
        return finish("NO_QTL_CANDIDATES", pd.DataFrame(columns=RESULT_COLUMNS), audit, region)
    for column in RESULT_COLUMNS:
        if column not in results.columns:
            results[column] = np.nan if column.startswith(("PP_", "N_", "BEST_SNP")) else ""
    todo = results[results["STATUS"] == "PENDING_R"]
    if len(todo):
        write_tsv(pd.concat(long_rows, ignore_index=True), work / "abf_input.tsv.gz")
        write_tsv(todo[["CANDIDATE_ID", "QTL_N", "QTL_SDY"]].fillna("NA"), work / "candidates.tsv")
        worker = work / "coloc_abf.R"
        worker.write_text(R_ABF, encoding="utf-8")
        command = [rscript, str(worker), str(work / "abf_input.tsv.gz"), str(work / "candidates.tsv"),
                   str(work / "abf_results.tsv"), str(cfg["P1"]), str(cfg["P2"]), str(cfg["P12"]),
                   str(cfg["P12_LOW"]), str(cfg["P12_HIGH"]), meta["TYPE"],
                   str(meta["S"] if np.isfinite(meta["S"]) else "NA"), str(meta["N"])]
        proc = subprocess.run(command, text=True, capture_output=True, timeout=int(cfg["R_TIMEOUT"]))
        (out / "coloc.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
        if proc.returncode != 0:
            raise RuntimeError(f"R coloc.abf worker failed (log {out / 'coloc.log'}):\n{proc.stderr[-3000:]}")
        abf = read_tsv(work / "abf_results.tsv", dtype=str).set_index("CANDIDATE_ID")
        for i, r in todo.iterrows():
            a = abf.loc[r["CANDIDATE_ID"]] if r["CANDIDATE_ID"] in abf.index else None
            if a is None or a["R_STATUS"] != "OK":
                results.at[i, "STATUS"] = "NOT_TESTED_INVALID_STATS"
                results.at[i, "ERROR"] = "" if a is None else str(a.get("MESSAGE", ""))
                continue
            pp = [safe_float(a[f"PP_H{h}"]) for h in range(5)]
            for h in range(5):
                results.at[i, f"PP_H{h}"] = pp[h]
            low, high = safe_float(a["PP_H4_P12_LOW"]), safe_float(a["PP_H4_P12_HIGH"])
            results.at[i, "PP_H4_P12_LOW"], results.at[i, "PP_H4_P12_HIGH"] = low, high
            for key, value in classify(pp, cfg["STRONG_PP4"], cfg["SUGGESTIVE_PP4"]).items():
                results.at[i, key] = value
            results.at[i, "PRIOR_ROBUST_STRONG"] = "YES" if np.isfinite(low) and low >= cfg["STRONG_PP4"] else "NO"
            results.at[i, "BEST_GWAS_VARIANT"] = results.at[i, "BEST_QTL_VARIANT"] = str(a.get("BEST_SNP", ""))
            results.at[i, "BEST_SNP_PP_H4"] = safe_float(a.get("BEST_SNP_PP_H4"))
            results.at[i, "GWAS_MODE"], results.at[i, "QTL_MODE"] = a.get("GWAS_MODE", ""), a.get("QTL_MODE", "")
            results.at[i, "COLOC_TESTED"] = "YES"
            results.at[i, "STATUS"] = "TESTED"
    order = pd.to_numeric(results.get("PP_H4"), errors="coerce").fillna(-1) if "PP_H4" in results else 0
    results = results.assign(_O=order).sort_values(["_O", "DISCOVERY_MIN_P"], ascending=[False, True]).drop(columns="_O")
    return finish("COMPLETE", results, audit, region)


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


def ready_registry(cfg: dict) -> pd.DataFrame:
    reg = load_registry(cfg)
    return reg[(reg["ELIGIBILITY"] == "READY") & reg["TYPE_SELECTED"]].reset_index(drop=True)


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
    for flag, value in (("--qtl-types", args.qtl_types), ("--candidate-p", args.candidate_p)):
        if value is not None:
            options += [flag, str(value)]
    if args.force:
        options.append("--force")
    if args.keep_work:
        options.append("--keep-work")
    return options


def planner_mode(args, cfg) -> None:
    P = project_paths(args)
    studies = step06_studies(P["S06"])
    if args.study:
        studies = [s for s in studies if s == args.study]
    if not studies:
        raise SystemExit(f"No fine-mapped studies under {P['S06']} - run Step06{panel_option(P['panel'])} first.")
    reg = load_registry(cfg)
    banner(f"STEP 11 PLANNER - external molecular QTL (coloc.abf) - {P['phenotype']} / {P['label']} - "
           f"LD panel {P['panel']}")
    print_registry_summary(reg)
    ready = reg[(reg["ELIGIBILITY"] == "READY") & reg["TYPE_SELECTED"]]
    if ready.empty:
        raise SystemExit(f"\nNo READY dataset in {cfg['REGISTRY']}. Run --scan, edit the registry (columns, n, enabled=1), "
                         f"then --check-resources.")
    order = step02_order(P)
    rows = [{"STUDY_INDEX": order.get(acc, ""), "STUDY_ACCESSION": acc, "LOCUS_ID": locus}
            for acc in studies for locus in list_loci(P["S06"] / acc)]
    n_tasks = min(len(rows), MAX_ARRAY_TASKS)
    manifest = pd.DataFrame(rows)
    manifest.insert(0, "TASK_ID", [i % n_tasks + 1 for i in range(len(rows))])
    P["OUT"].mkdir(parents=True, exist_ok=True)
    P["LOG"].mkdir(parents=True, exist_ok=True)
    manifest_file = P["OUT"] / "coloc_manifest.tsv"
    write_tsv(manifest, manifest_file)
    write_tsv(ready.drop(columns=["ELIGIBILITY", "TYPE_SELECTED"]), P["OUT"] / "registry_snapshot.tsv")

    values = resolve_slurm(args, cfg)
    array = (f"1-{n_tasks}" if n_tasks > 1 else "1") + (f"%{args.max_parallel}" if args.max_parallel else "")
    directives = {**values, "job-name": f"CQTL_{P['p']}_{P['code'].lower()}"[:100], "array": array,
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

{' '.join(shlex.quote(x) for x in command)}
""", encoding="utf-8", newline="\n")
    P["BASH"].chmod(0o755)
    print(f"\nStudies   : {len(studies)}\nLoci      : {len(manifest)} -> array {array}"
          + (f" ({len(manifest) / n_tasks:.1f} loci per task)" if len(manifest) > n_tasks else " (one locus per task)"))
    print(f"Datasets  : {ready['dataset'].nunique()} READY ({len(ready)} files), frozen in registry_snapshot.tsv")
    print(f"Manifest  : {manifest_file}\n\nRUN:\n  sbatch {P['BASH'].name}\nTHEN:\n  {base_command(P)} --status\n"
          f"  {base_command(P)} --stats")


def load_snapshot(P: dict) -> pd.DataFrame:
    snap = read_tsv(P["OUT"] / "registry_snapshot.tsv", dtype=str, keep_default_na=False)
    if snap.empty:
        raise SystemExit("registry_snapshot.tsv missing - run the planner again.")
    return snap


def worker_mode(args, cfg) -> None:
    P = project_paths(args)
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    manifest = read_tsv(Path(args.worker_manifest), dtype=str)
    selected = manifest[pd.to_numeric(manifest["TASK_ID"]) == task]
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
        raise SystemExit(f"Study {args.index} ({accession}) has no Step06 locus{panel_option(P['panel'])}.")
    if args.locus:
        if args.locus not in loci:
            raise SystemExit(f"{args.locus} is not a locus of {accession}: {', '.join(loci)}")
        loci = [args.locus]
    rscript = require_tools()
    ready = ready_registry(cfg)
    if ready.empty:
        raise SystemExit("No READY dataset in the registry (see --check-resources).")
    failed = [locus for locus in loci if not run_locus(pd.Series({"STUDY_ACCESSION": accession, "LOCUS_ID": locus}),
                                                      args, cfg, ready, rscript)]
    if failed:
        raise SystemExit(f"{accession}: colocalisation failed for {', '.join(failed)}")


def check_resources_mode(args, cfg) -> None:
    P = project_paths(args)
    problems = []
    banner("STEP 11 - RESOURCE CHECK - external molecular QTL (coloc.abf)")
    print(f"Config     : {cfg['_FILE'] or 'built-in defaults'}")
    print(f"Registry   : {project_path(cfg['REGISTRY'])}")
    print(f"QTL types  : {', '.join(cfg['QTL_TYPES'])}   CANDIDATE_P {cfg['CANDIDATE_P']}   MIN_SHARED {cfg['MIN_SHARED']}")
    rscript = find_executable("Rscript")
    if rscript:
        proc = subprocess.run([rscript, "-e", 'cat(as.character(packageVersion("coloc")))'], text=True,
                              capture_output=True)
        ok = proc.returncode == 0
        print(f"  {'OK     ' if ok else 'MISSING'}  R coloc        {proc.stdout.strip() if ok else proc.stderr.strip()[-150:]}")
        if not ok:
            problems.append("R package coloc")
    else:
        print("  MISSING  Rscript")
        problems.append("Rscript (Step00_Installation.py)")
    import importlib.util
    if importlib.util.find_spec("pysam") is None:
        print("  MISSING  python pysam")
        problems.append("python pysam")
    else:
        print("  OK       python pysam")

    banner("REGISTRY ELIGIBILITY (per dataset)")
    reg = load_registry(cfg)
    print_registry_summary(reg)
    ready = reg[(reg["ELIGIBILITY"] == "READY") & reg["TYPE_SELECTED"]]
    if reg.empty:
        problems.append(f"no registry - run: {base_command(P)} --scan")
    elif ready.empty:
        problems.append("no READY dataset: edit the registry (columns, n, enabled = 1)")
    P["OUT"].mkdir(parents=True, exist_ok=True)
    if len(reg):
        write_tsv(reg, P["OUT"] / "resource_plan.tsv")
        print(f"\nPer-file eligibility: {P['OUT'] / 'resource_plan.tsv'}")

    banner(f"STEP06 INPUT - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    studies = step06_studies(P["S06"])
    n_loci = sum(len(list_loci(P["S06"] / s)) for s in studies)
    print(f"  {P['S06']}: {len(studies)} studies, {n_loci} loci")
    for s in studies:
        meta = gwas_meta(P["S06"] / s)
        flag = "OK" if meta["TYPE"] and np.isfinite(meta["N"]) else "NOT TESTABLE"
        print(f"  {s:<15} {flag:<13} type {meta['TYPE'] or '-':<5} N {meta['N']:<10.0f} "
              f"case fraction {meta['S']:.3f}  ({meta['SOURCE']})" if np.isfinite(meta["N"]) else
              f"  {s:<15} {flag:<13} ({meta['SOURCE']})")
    if not n_loci:
        problems.append(f"no Step06 loci - run Step06{panel_option(P['panel'])} first")
    banner("CHECK RESULT")
    if problems:
        print(f"NOT READY - {len(problems)} problem(s):")
        for p in problems:
            print(f"  - {p}")
        raise SystemExit(1)
    print(f"READY - {ready['dataset'].nunique()} datasets, {n_loci} loci. Registry readiness does not guarantee that a "
          "dataset is statistically compatible with the GWAS (check N, ancestry, build).")
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
    order = step02_order(P)
    banner(f"STEP 11 STATUS - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    snap = read_tsv(P["OUT"] / "registry_snapshot.tsv", dtype=str)
    print(f"Datasets in the submitted snapshot: {snap['dataset'].nunique()}" if len(snap)
          else "No job script made yet (run the planner).")
    rows, failures = [], []
    for accession in step06_studies(P["S06"]):
        counts = {"LOCI": 0, "DONE": 0, "NO_CANDIDATES": 0, "NOT_TESTABLE": 0, "FAILED": 0, "PENDING": 0,
                  "CANDIDATES": 0, "TESTED": 0, "STRONG": 0}
        for locus in list_loci(P["S06"] / accession):
            state, info = locus_state(P, accession, locus)
            counts["LOCI"] += 1
            key = {"COMPLETE": "DONE", "NO_QTL_CANDIDATES": "NO_CANDIDATES", "FAILED": "FAILED",
                   "NOT_TESTED_GWAS_N_OR_TYPE_MISSING": "NOT_TESTABLE"}.get(state, "PENDING")
            counts[key] += 1
            if key == "FAILED":
                failures.append(f"{accession} {locus}: {info.get('ERROR', '')}")
            for k in ("CANDIDATES", "TESTED", "STRONG"):
                counts[k] += int(info.get(f"N_{k}", 0) or 0)
        rows.append({"IDX": order.get(accession, ""), "STUDY": accession, **counts})
    if not rows:
        print(f"No Step06 studies under {P['S06']}.")
        return
    table = pd.DataFrame(rows)
    with pd.option_context("display.width", 200, "display.max_rows", None):
        print(table.to_string(index=False))
    total = table.drop(columns=["IDX", "STUDY"]).sum()
    print("\nLoci: " + "  ".join(f"{k}={int(total[k])}" for k in ("LOCI", "DONE", "NO_CANDIDATES", "NOT_TESTABLE",
                                                                       "FAILED", "PENDING")))
    for line in failures[:20]:
        print(f"  FAILED {line}")
    if int(total["FAILED"]) or int(total["PENDING"]):
        print(f"\nRe-run failed / pending loci: sbatch {P['BASH'].name}   (finished loci are skipped)")


def stats_mode(args, cfg) -> None:
    P = project_paths(args)
    order = step02_order(P)
    banner(f"STEP 11 RESULTS - external molecular QTL (coloc.abf) - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    frames = [read_tsv(f) for acc in step06_studies(P["S06"]) for f in
              sorted((P["OUT"] / acc / "loci").glob(f"L*/{RESULTS_FILE}")) if nonempty(f)]
    frames = [f for f in frames if len(f)]
    if not frames:
        print("No results yet (see --status).")
        return
    x = pd.concat(frames, ignore_index=True)
    tested = x[x["COLOC_TESTED"].astype(str) == "YES"].copy()
    tested["PP_H4"] = pd.to_numeric(tested["PP_H4"], errors="coerce")
    strong, suggestive = cfg["STRONG_PP4"], cfg["SUGGESTIVE_PP4"]
    rows = []
    for (acc, qtl), g in x.groupby(["STUDY_ACCESSION", "QTL_TYPE"]):
        t = tested[(tested.STUDY_ACCESSION == acc) & (tested.QTL_TYPE == qtl)]
        b = t.sort_values("PP_H4", ascending=False).head(1)
        rows.append({"IDX": order.get(acc, ""), "STUDY": acc, "QTL_TYPE": qtl, "CANDIDATES": len(g), "TESTED": len(t),
                     "STRONG": int((t.PP_H4 >= strong).sum()),
                     "SUGGESTIVE": int(((t.PP_H4 >= suggestive) & (t.PP_H4 < strong)).sum()),
                     "PRIOR_ROBUST": int((t.PRIOR_ROBUST_STRONG.astype(str) == "YES").sum()),
                     "TOP (PP.H4)": f"{b.iloc[0]['MOLECULAR_TRAIT_ID']} {b.iloc[0]['TISSUE']} ({b.iloc[0]['PP_H4']:.2f})"
                     if len(b) else ""})
    table = pd.DataFrame(rows)
    print(f"\nPER STUDY x QTL TYPE (STRONG = PP.H4 >= {strong}; single-causal-variant coloc.abf, exploratory)")
    with pd.option_context("display.width", 250, "display.max_rows", None, "display.max_colwidth", 50):
        print(table.to_string(index=False))
        top = tested.sort_values("PP_H4", ascending=False).head(25)
        if len(top):
            print("\nTOP COLOCALISATIONS")
            show = top.assign(IDX=top["STUDY_ACCESSION"].map(order).fillna(""), PP_H4=top["PP_H4"].round(3))
            cols = ["IDX", "STUDY_ACCESSION", "LOCUS_ID", "QTL_TYPE", "STUDY_LABEL", "TISSUE", "MOLECULAR_TRAIT_ID",
                    "N_COMMON_VARIANTS", "PP_H4", "PRIOR_ROBUST_STRONG", "COLOC_CLASS"]
            print(show[cols].to_string(index=False))
    print("\nCandidate status: " + "  ".join(f"{k}={v}" for k, v in x["STATUS"].astype(str).value_counts().items()))
    write_tsv(x, P["OUT"] / "coloc_all_candidates.tsv")
    write_tsv(table, P["OUT"] / "coloc_stats.tsv")
    print(f"\nTables: {P['OUT'] / 'coloc_stats.tsv'}\n        {P['OUT'] / 'coloc_all_candidates.tsv'}")


def arguments():
    parser = argparse.ArgumentParser(description="Step11: external molecular-QTL coloc.abf, one SLURM task per locus.")
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--scan", action="store_true", help="Inventory resources/molecular_qtl -> registry (disabled rows).")
    modes.add_argument("--index", type=int, help="Study N of the Step02 download manifest: run its loci here.")
    modes.add_argument("--status", action="store_true", help="READ-ONLY: loci done / failed / pending per study.")
    modes.add_argument("--stats", action="store_true", help="Colocalisation results (tables on screen + TSV).")
    modes.add_argument("--check-resources", action="store_true",
                       help="Registry eligibility, tools and Step06 input (writes resource_plan.tsv; exit 1 if not ready).")
    modes.add_argument("--worker-manifest", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--locus", default=None, help="With --index: only this locus (e.g. L001).")
    parser.add_argument("--study", default=None, help="Planner: only this GCST accession.")
    parser.add_argument("--reference-panel", default=DEFAULT_PANEL,
                        help="LD panel of Steps 04-06 (1000G default, 1000G30x, ...): reads 06_finemapping/<p>/<a>_<PANEL>/.")
    parser.add_argument("--config", default="Step11_Config.py", help="Parameter file (default Step11_Config.py).")
    parser.add_argument("--qtl-types", default=None, help="Override QTL_TYPES, e.g. pQTL,sQTL.")
    parser.add_argument("--candidate-p", type=float, default=None, help="Override CANDIDATE_P.")
    parser.add_argument("--overwrite", action="store_true", help="--scan: rebuild the registry (drops your edits).")
    parser.add_argument("--force", action="store_true", help="Redo loci that are already done.")
    parser.add_argument("--keep-work", action="store_true", help="Keep per-locus work/ files.")
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
    cfg = load_config(config_path(args))
    if args.qtl_types:
        cfg["QTL_TYPES"] = [QTL_TYPE_NAMES.get(t.strip().lower(), t.strip()) for t in args.qtl_types.split(",") if t.strip()]
    if args.candidate_p is not None:
        cfg["CANDIDATE_P"] = args.candidate_p
    if args.scan:
        scan_mode(args, cfg)
    elif args.check_resources:
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
