#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 11 - External molecular-QTL colocalisation, per locus
                    (ADSP FunGen xQTL brain sQTL / pQTL / haQTL / mQTL, Kumasaka caQTL, ...)
===============================================================================

    # 1. resource table: datasets from the Step00 manifests (no recursive scan), columns from the file
    #    header, coordinates / alleles checked against GRCh38, QTL N from metadata -> validation status
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --scan

    # 2. what is ready for this phenotype: datasets, duplicates of Step09 / Step10, GWAS N and type,
    #    every Step06 locus with its GWAS input (exit 1 if nothing can be tested)
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources

    # 3. plan: ONE job script, one SLURM array task per Step06 locus (max 1000; more loci -> several per task)
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR
    sbatch Step11_CustomQTL_parkinson_s_disease_european.sh

    # study 2 on this machine: all its loci / one locus; progress, results, redo
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --index 2
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --index 2 --locus L001
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --stats
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --force

    # ---- 1000 Genomes 30x LD panel: loci from 06_finemapping/<p>/<a>_1000G30x/, results in 11_custom_qtl_coloc/<p>/<a>_1000G30x/
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --check-resources
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x
    sbatch Step11_CustomQTL_parkinson_s_disease_european_1000G30x.sh
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --index 2 --locus L001
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --status
    envs/pipeline/bin/python Step11_Custom_QTL_Colocalization.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

Inputs (nothing is rediscovered that an earlier step already recorded)
  Resources  Step00 manifests: resources/molecular_qtl/<qtl>/adsp_fungen_brain/ADSP_FunGen_xQTL.v1.<QTL>.all.manifest.txt
             (+ the .metadata.json / .metadata.txt / provenance.json next to them), eQTL-Catalogue-format files listed in
             Step11_Config.py (Kumasaka caQTL QTD100018 + peak metadata, Sun_2018 QTD000584), the eQTL Catalogue
             tables (tabix_ftp_paths*.tsv: sample_size, tissue) and, if present, resources/molecular_qtl/resource_registry.tsv
             and resources/coloc/custom_qtl_registry.tsv. A modality directory without association files is reported
             as NOT_DOWNLOADED, never as processed.
  GWAS       Step02 01_gwas_catalog/<p>/<a>/GWAS_study_metadata.tsv: PHENOTYPE_TYPE, N_CASES, N_CONTROLS, N_TOTAL, N_EFF
             Binary -> coloc type cc, N = N_CASES + N_CONTROLS, s = N_CASES / N  (N_EFF, the SuSiE-RSS N, is never used)
             Continuous -> type quant, N = N_TOTAL.  Unknown / NEEDS_REVIEW -> NOT_TESTED.
  Loci       Step06 <GCST>/loci_definition.tsv: EVERY locus (also those without a GWAS credible set).
             GWAS statistics: Step06 regional_gwas/<L>.tsv.gz (extracted BEFORE LD-panel matching); if absent, the
             Step02 QC file within the same locus boundaries (recorded as STEP02_QC_FALLBACK). The LD-filtered
             Step06 matched table is never used.

Per locus and dataset
  1. QTL rows of the locus window by tabix (BED tracks: POS = chromStart + 1, verified per dataset).
  2. Candidates: molecular traits with P <= CANDIDATE_P (top MAX_TRAITS).
  3. Harmonisation on CHR + POS + both alleles: ALIGNED (QTL REF/ALT = GWAS NEA/EA) or SWAPPED (GWAS beta sign
     flipped); no strand-complement matching (both GRCh38 forward strand); indels only with identical alleles;
     palindromic SNPs per PALINDROMIC_POLICY. GWAS and QTL MAF are kept apart: a variant without its own QTL
     frequency is excluded, never filled from the GWAS. Every QTL row keeps its original id, harmonised id, match
     status and exclusion reason (harmonised_variants.tsv.gz).
  4. Overlap: >= MIN_SHARED shared variants and >= MIN_OVERLAP of the smaller side.
  5. Method: datasets with provider SuSiE (credible sets + LBF) are colocalised with coloc.bf_bf in Step09 / Step10
     (multi-signal) and are NOT re-run here as ABF. All others: coloc.abf (single causal variant), priors p1 / p2 /
     p12 with sensitivity at P12_LOW / P12_HIGH, PP.H0-H4 reported. Loci with several independent GWAS signals or
     SuSiE credible sets, or in a complex-LD region (MHC, 17q21.31 MAPT, 8p23.1), are flagged
     (ABF_ASSUMPTION_WARNING) - the single-causal-variant assumption may not hold there.
  6. Datasets already in the Step09 / Step10 QTL snapshot are skipped (DUPLICATE_OF_STEP09/10); same cohort +
     modality + tissue is kept but marked INDEPENDENT_EVIDENCE = NO.

NOT_TESTED is never evidence against colocalisation: missing N / MAF, too little overlap, no candidate trait or
a missing input give an explicit reason per candidate, dataset and locus.

Outputs
  config/Step11_Custom_QTL_Registry.tsv   resource table: one row per dataset (study, modality, context, build,
                                          files, columns, coordinates, QTL N + source, validation status)
  config/Step11_Custom_QTL_Files.tsv      one row per file (manifest md5 / size, present, indexed)
  11_custom_qtl_coloc/<phenotype>/<ancestry>[_<PANEL>]/
    resource_plan.tsv, locus_plan.tsv     (--check-resources / planner)
    registry_snapshot.tsv, coloc_manifest.tsv
    <GCST>/loci/<L>/locus_coloc_results.tsv, harmonised_variants.tsv.gz, dataset_audit.tsv, custom_coloc_done.json
    --stats: coloc_stats.tsv, coloc_all_candidates.tsv, locus_reasons.tsv

Parameters: Step11_Config.py; SLURM from the HPC_CPU.sh header, then SLURM in the config, then
--partition --qos --account --cpus --mem --time --max-parallel.
===============================================================================
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
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
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


VERSION = "2.0.0-custom-qtl-metadata-driven"
PROJECT_ROOT = Path.cwd().resolve()
MAX_ARRAY_TASKS = 1000
DEFAULT_PANEL = "1000G"
OUT_DIR = "11_custom_qtl_coloc"
LOG_DIR = "step11_custom_qtl_coloc"
SCRIPT_PREFIX = "Step11_CustomQTL"
DONE_FILE = "custom_coloc_done.json"
FAILED_FILE = "LOCUS_FAILED.txt"
RESULTS_FILE = "locus_coloc_results.tsv"
HARMONISED_FILE = "harmonised_variants.tsv.gz"
AUDIT_FILE = "dataset_audit.tsv"

DEFAULT_CONFIG: dict[str, Any] = {
    "RESOURCE_ROOT": "resources/molecular_qtl",
    "REGISTRY": "config/Step11_Custom_QTL_Registry.tsv",
    "FILES_TABLE": "config/Step11_Custom_QTL_Files.tsv",
    "ADSP_MODALITIES": {"sQTL": "sqtl", "pQTL": "pqtl", "haQTL": "haqtl", "mQTL": "mqtl"},
    "ADSP_PROVIDER_DIR": "adsp_fungen_brain",
    "EQTL_CATALOGUE_FORMAT_DATASETS": [],
    "EXTERNAL_REGISTRIES": ["resources/molecular_qtl/resource_registry.tsv", "resources/coloc/custom_qtl_registry.tsv"],
    "EQTL_CATALOGUE_TABLES": ["resources/coloc/eqtl_catalogue/tabix_ftp_paths.tsv",
                              "resources/coloc/eqtl_catalogue/tabix_ftp_paths_imported.tsv",
                              "resources/molecular_qtl/manifests/tabix_ftp_paths.tsv",
                              "resources/molecular_qtl/manifests/tabix_ftp_paths_imported.tsv"],
    "NIAGADS_MANIFEST": "resources/molecular_qtl/manifests/NG00184.v1.csv",
    "EQTL_CATALOGUE_SUSIE_DIR": "resources/coloc/eqtl_catalogue/susie_provider",
    "REFERENCE_FASTA": "resources/genome/GRCh38/GRCh38.primary_assembly.genome.fa",
    "QTL_TYPES": ["sQTL", "pQTL", "haQTL", "caQTL", "mQTL"],
    "DISABLED_DATASETS": [],
    "QTL_N": {},
    "QTL_SDY": {},
    "COLUMN_OVERRIDES": {},
    "ALLOW_VERIFIED_AN_HALF": True,
    "ALLOW_UNREVIEWED_GWAS_TYPE": False,
    "ALLOW_DUPLICATES_OF_STEP09_10": False,
    "INCLUDE_INDELS": True,
    "PALINDROMIC_POLICY": "drop_high_maf",
    "PALINDROMIC_MAX_MAF": 0.42,
    "CANDIDATE_P": 1e-5,
    "MAX_TRAITS": 100,
    "MIN_SHARED": 50,
    "MIN_OVERLAP": 0.5,
    "MAX_REGION_ROWS": 3_000_000,
    "VALIDATION_ROWS": 2000,
    "REF_CHECK_MIN_MATCH": 0.95,
    "STATS_MAX_INCONSISTENT": 0.10,
    "COMPLEX_LD_REGIONS": {"MHC": (6, 28_510_120, 33_480_577), "17q21.31_MAPT_inversion": (17, 45_300_000, 46_700_000),
                           "8p23.1_inversion": (8, 8_100_000, 12_000_000)},
    "P1": 1e-4,
    "P2": 1e-4,
    "P12": 5e-6,
    "P12_LOW": 1e-6,
    "P12_HIGH": 1e-5,
    "STRONG_PP4": 0.80,
    "SUGGESTIVE_PP4": 0.50,
    "R_TIMEOUT": 3600,
    "SLURM": {},
}

# parameters that change a locus result (resume hash)
RESULT_KEYS = ("CANDIDATE_P", "MAX_TRAITS", "MIN_SHARED", "MIN_OVERLAP", "INCLUDE_INDELS", "PALINDROMIC_POLICY",
               "PALINDROMIC_MAX_MAF", "P1", "P2", "P12", "P12_LOW", "P12_HIGH", "STRONG_PP4", "SUGGESTIVE_PP4",
               "ALLOW_UNREVIEWED_GWAS_TYPE", "ALLOW_DUPLICATES_OF_STEP09_10", "COMPLEX_LD_REGIONS", "QTL_SDY")

QTL_TYPE_NAMES = {"caqtl": "caQTL", "haqtl": "haQTL", "mqtl": "mQTL", "meqtl": "mQTL", "pqtl": "pQTL",
                  "sqtl": "sQTL", "eqtl": "eQTL"}

# column roles; matched case-insensitively against the file's own header (never against headerless c1, c2, ...)
ROLE_ALIASES = {
    "chrom": ("chrom", "chr", "chromosome", "seqnames", "chr_name"),
    "start": ("chromstart", "start", "bed_start", "chrom_start"),
    "end": ("chromend", "end", "stop", "bed_end", "chrom_end"),
    "pos": ("position", "pos", "bp", "base_pair_location", "variant_pos", "pos_hg38"),
    "variant": ("variant", "variant_id", "variantid", "snp", "snp_id", "snpid", "variant_key", "id", "marker",
                "markername"),
    "ref": ("ref", "ref_allele", "reference_allele", "other_allele", "non_effect_allele", "nea", "allele0"),
    "alt": ("alt", "alt_allele", "alternate_allele", "effect_allele", "ea", "allele1"),
    "trait": ("molecular_trait_id", "phenotype_id", "molecular_trait", "target", "target_id", "feature_id", "feature",
              "trait_id", "peak_id", "intron_id", "cpg", "cpg_id", "probe_id", "protein_id", "region_id", "event_id"),
    "gene": ("gene_id", "gene", "target_gene", "ensembl_gene_id", "gene_name", "target_gene_symbol", "gene_symbol"),
    "beta": ("beta", "slope", "effect_size", "bhat", "beta_hat", "b", "effect"),
    "se": ("se", "slope_se", "standard_error", "sebeta", "sbhat", "se_beta", "stderr"),
    "p": ("pvalue", "p_value", "pval", "p", "pval_nominal", "nominal_p", "p.value", "pvalue_nominal"),
    "freq": ("maf", "af", "alt_af", "alt_freq", "freq", "eaf", "effect_allele_frequency", "alt_allele_freq",
             "allele_frequency", "af_alt"),
    "n": ("n", "sample_size", "n_samples", "num_samples", "nsamples", "n_obs", "n_total"),
    "an": ("an",),
    "ac": ("ac",),
    "z": ("z", "zscore", "z_score"),
}
ROLES = tuple(ROLE_ALIASES)
EQTL_CATALOGUE_SCHEMA = {"molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "pvalue", "beta",
                         "se", "maf", "an"}
VARRE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y)[_:]([0-9]+)[_:]([ACGTN]+)[_:]([ACGTN]+)$", re.I)
ADSP_NAME = re.compile(r"^FunGenADv\d+\.(?P<cohort>[^.]+)\.(?P<qtl>[^.]+)\.(?P<context>.+?)\.(?P<id>[A-Za-z]+\d+)"
                       r"_v\d+_(?P<chrom>[0-9XYM]{1,2})_(?P<kind>snp|indel)_all\.(?P<date>\d+)\.bed\.gz$")
ENSG = re.compile(r"ENSG\d{11}")
ACGT = re.compile(r"^[ACGT]+$")
ACGT_PAT = r"^[ACGT]+$"          # pandas .str.match
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

REGISTRY_COLUMNS = [
    "DATASET_ID", "SOURCE_STUDY", "PROVIDER", "COHORT", "MODALITY", "CONTEXT", "BUILD", "BUILD_EVIDENCE",
    "FILE_FORMAT", "DISCOVERED_BY", "METADATA_FILES", "N_FILES_EXPECTED", "N_FILES_PRESENT", "N_FILES_INDEXED",
    "CHROMOSOMES", "EXAMPLE_FILE", "HEADER", "COLUMNS", "COLUMN_SOURCE",
    *[f"COL_{r.upper()}" for r in ROLES], "COORDINATES", "ALLELE_SOURCE", "EFFECT_ALLELE", "FREQ_KIND",
    "TRAIT_METADATA", "QTL_N", "QTL_N_SOURCE", "QTL_N_PER_TRAIT", "AN_HALF_MEDIAN", "QTL_SDY", "QTL_SDY_SOURCE",
    "PROVIDER_FINEMAPPING", "METHOD", "EQTL_CATALOGUE_ID", "VALIDATION_STATUS", "VALIDATION_DETAILS", "WARNINGS",
    "READY", "SCANNED_UTC", "STEP11_VERSION",
]
FILES_COLUMNS = ["DATASET_ID", "CHROM", "KIND", "PATH", "MANIFEST", "MANIFEST_MD5", "MANIFEST_SIZE", "PRESENT",
                 "SIZE_OK", "INDEXED"]

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


def usable(value: Any) -> bool:
    return str(value).strip().lower() not in {"", "na", "nan", "none", ".", "null"}


def read_json(path: Path) -> Any:
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


def project_path(value: Any) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else PROJECT_ROOT / path


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def clean_chrom(value: Any) -> str:
    c = str(value).strip()
    c = c[3:] if c.lower().startswith("chr") else c
    return c.upper() if c.upper() in {"X", "Y", "M", "MT"} else c.lstrip("0") or "0"


def log10_p_from_z(z: np.ndarray) -> np.ndarray:
    """log10 of the two-sided normal P of z (stable for large |z|)."""
    try:
        from scipy.stats import norm
        return (np.log(2.0) + norm.logsf(np.abs(z))) / np.log(10.0)
    except ImportError:
        return np.log10(np.maximum([math.erfc(abs(v) / math.sqrt(2.0)) for v in z], 1e-320))


def stats_inconsistency(beta: pd.Series, se: pd.Series, p: pd.Series) -> tuple[float, int]:
    """Fraction of rows whose reported P is not reproduced by beta / SE (wrong column mapping or scale)."""
    b, s, pv = (pd.to_numeric(x, errors="coerce").to_numpy(dtype=float) for x in (beta, se, p))
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0) & np.isfinite(pv) & (pv > 1e-300) & (pv <= 1)
    if ok.sum() == 0:
        return np.nan, 0
    got = -np.log10(pv[ok])
    exp = -log10_p_from_z(b[ok] / s[ok])
    bad = np.abs(got - exp) > 0.3 + 0.1 * np.maximum(got, exp)
    return float(bad.mean()), int(ok.sum())


def is_palindromic(a: str, b: str) -> bool:
    return len(a) == 1 and len(b) == 1 and COMPLEMENT.get(a) == b


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
    if cfg["PALINDROMIC_POLICY"] not in {"keep", "drop_high_maf", "drop"}:
        raise SystemExit("PALINDROMIC_POLICY must be keep, drop_high_maf or drop")
    return cfg


def results_hash(cfg: dict) -> str:
    return sha1_short(json.dumps({k: cfg[k] for k in RESULT_KEYS}, sort_keys=True, default=str))


def results_slug(ancestry_slug: str, panel: str) -> str:
    return ancestry_slug if panel == DEFAULT_PANEL else f"{ancestry_slug}_{panel}"


def panel_option(panel: str) -> str:
    return "" if panel == DEFAULT_PANEL else f" --reference-panel {panel}"


def project_paths(args) -> dict:
    code, label = canonical_ancestry(args.ancestry)
    p, a = slugify(args.phenotype), slugify(label)
    r = results_slug(a, args.reference_panel)
    return {"phenotype": args.phenotype.strip(), "code": code, "label": label, "p": p, "a": a, "r": r,
            "panel": args.reference_panel, "S01": PROJECT_ROOT / "01_gwas_catalog" / p / a,
            "S02": PROJECT_ROOT / "02_summary_stats" / p / a / "qc",
            "S06": PROJECT_ROOT / "06_finemapping" / p / r, "S09": PROJECT_ROOT / "09_coloc" / p / r,
            "S10": PROJECT_ROOT / "10_celltype_coloc" / p / r, "OUT": PROJECT_ROOT / OUT_DIR / p / r,
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
    path = P["S01"] / "GWAS_download_manifest.tsv"
    if not path.is_file():
        return {}
    m = pd.read_csv(path, sep="\t", dtype=str)
    ids = pd.to_numeric(m.get("ARRAY_ID"), errors="coerce")
    return {str(acc).strip(): int(i) for acc, i in zip(m["STUDY_ACCESSION"], ids) if pd.notna(i)}


def study_for_index(P: dict, index: int) -> str:
    for accession, i in step02_order(P).items():
        if i == index:
            return accession
    raise SystemExit(f"--index {index} is not in the Step02 download manifest under {P['S01']}")


# =============================================================================
# GWAS SIDE: Step02 metadata, Step06 loci, regional statistics
# =============================================================================

def step02_metadata(P: dict) -> dict[str, dict]:
    path = P["S01"] / "GWAS_study_metadata.tsv"
    if not nonempty(path):
        return {}
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    key = "GCST_ID" if "GCST_ID" in table.columns else "STUDY_ACCESSION"
    return {str(r[key]).strip(): r for r in table.to_dict(orient="records")}


def gwas_meta(P: dict, accession: str, cfg: dict, meta_table: dict | None = None) -> dict:
    """coloc type / N / case fraction of the GWAS, each with its source. NOT_TESTED reason if not determinable."""
    table = step02_metadata(P) if meta_table is None else meta_table
    row = table.get(accession)
    source_file = rel(P["S01"] / "GWAS_study_metadata.tsv")
    summary = read_json(P["S06"] / accession / "finemapping_summary.json")
    out = {"GWAS_TYPE": "", "GWAS_N": np.nan, "GWAS_CASE_FRACTION": np.nan, "GWAS_N_SOURCE": "",
           "GWAS_N_EFF_NOT_USED": np.nan, "GWAS_PHENOTYPE_TYPE": "", "GWAS_META_REASON": "", "GWAS_META_WARNING": ""}
    if row is None:
        if not summary:
            out["GWAS_META_REASON"] = f"NOT_TESTED_GWAS_METADATA_MISSING ({accession} not in {source_file})"
            return out
        row, source_file = summary, rel(P["S06"] / accession / "finemapping_summary.json") + " (Step02 copy)"
    kind = str(row.get("PHENOTYPE_TYPE", "") or "").strip()
    review = str(row.get("NEEDS_REVIEW", "")).strip().lower() == "true"
    build = str(row.get("GENOME_BUILD", "GRCh38") or "GRCh38").strip()
    out["GWAS_PHENOTYPE_TYPE"] = kind
    out["GWAS_N_EFF_NOT_USED"] = safe_float(row.get("N_EFF"))
    if build.upper() not in {"GRCH38", "HG38"}:
        out["GWAS_META_REASON"] = f"NOT_TESTED_GWAS_NOT_GRCH38 (GENOME_BUILD={build})"
        return out
    if (kind not in {"Binary", "Continuous"} or review) and not cfg["ALLOW_UNREVIEWED_GWAS_TYPE"]:
        out["GWAS_META_REASON"] = (f"NOT_TESTED_GWAS_TYPE_UNVERIFIED (PHENOTYPE_TYPE={kind or 'Unknown'}"
                                   f"{' NEEDS_REVIEW' if review else ''}; set VERIFIED_PHENOTYPE_TYPE in Step02)")
        return out
    cases, controls, total = (safe_float(row.get(k)) for k in ("N_CASES", "N_CONTROLS", "N_TOTAL"))
    if kind == "Binary":
        if not (np.isfinite(cases) and np.isfinite(controls) and cases > 0 and controls > 0):
            out["GWAS_META_REASON"] = ("NOT_TESTED_GWAS_CASE_CONTROL_COUNTS_MISSING (coloc cc needs N_CASES and "
                                       "N_CONTROLS for N and the case fraction; N_EFF is not a substitute)")
            return out
        n = cases + controls
        out.update(GWAS_TYPE="cc", GWAS_N=n, GWAS_CASE_FRACTION=cases / n,
                   GWAS_N_SOURCE=f"N_CASES + N_CONTROLS ({source_file})")
        if np.isfinite(total) and abs(total - n) / n > 0.01:
            out["GWAS_META_WARNING"] = f"N_TOTAL {total:.0f} differs from N_CASES + N_CONTROLS {n:.0f}; used the sum"
        return out
    if not (np.isfinite(total) and total > 0):
        out["GWAS_META_REASON"] = "NOT_TESTED_GWAS_N_MISSING (Continuous trait without N_TOTAL)"
        return out
    out.update(GWAS_TYPE="quant", GWAS_N=total, GWAS_N_SOURCE=f"N_TOTAL ({source_file})")
    return out


def loci_table(study_dir: Path) -> pd.DataFrame:
    path = study_dir / "loci_definition.tsv"
    if not nonempty(path):
        return pd.DataFrame()
    x = pd.read_csv(path, sep="\t", dtype={"LOCUS_ID": str})
    need = {"LOCUS_ID", "CHR", "LOCUS_START", "LOCUS_END"}
    if not need <= set(x.columns):
        raise RuntimeError(f"{path} lacks {sorted(need - set(x.columns))}")
    return x


def list_loci(study_dir: Path) -> list[str]:
    try:
        return sorted(loci_table(study_dir).get("LOCUS_ID", pd.Series(dtype=str)).astype(str))
    except Exception:
        return []


def step06_studies(S06: Path) -> list[str]:
    if not S06.exists():
        return []
    return sorted(d.name for d in S06.iterdir() if d.is_dir() and d.name.startswith("GCST") and list_loci(d))


def gwas_cs_counts(study_dir: Path) -> dict[str, int]:
    path = study_dir / f"{study_dir.name}_95pct_credible_sets.tsv"
    if not nonempty(path):
        return {}
    cs = pd.read_csv(path, sep="\t", dtype=str)
    if cs.empty or not {"LOCUS_ID", "CREDIBLE_SET"} <= set(cs.columns):
        return {}
    return cs.drop_duplicates(["LOCUS_ID", "CREDIBLE_SET"]).groupby("LOCUS_ID").size().astype(int).to_dict()


GWAS_COLUMNS = ["CHR", "POS", "SNPID", "EA", "NEA", "BETA", "SE", "P", "Z", "EAF"]


def canon_gwas(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    def get(name):
        return frame[name] if name in columns else pd.Series(np.nan, index=frame.index)
    out = pd.DataFrame({"CHR": get("CHR").map(clean_chrom), "POS": pd.to_numeric(get("POS"), errors="coerce"),
                        "SNPID": get("SNPID").astype(str) if "SNPID" in columns else get("SNP").astype(str),
                        "EA": get("EA").astype(str).str.upper(), "NEA": get("NEA").astype(str).str.upper(),
                        "BETA": pd.to_numeric(get("BETA"), errors="coerce"),
                        "SE": pd.to_numeric(get("SE"), errors="coerce"),
                        "P": pd.to_numeric(get("P"), errors="coerce"),
                        "EAF": pd.to_numeric(get("EAF"), errors="coerce")})
    return out


def read_qc_region(qc_file: Path, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Step02 QC rows of one locus (streamed; the QC file is not tabix-indexed)."""
    header = list(pd.read_csv(qc_file, sep="\t", nrows=0, compression="infer").columns)
    use = [c for c in ("CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P", "EAF", "MAF", "Z") if c in header]
    missing = {"CHR", "POS", "EA", "NEA", "P"} - set(use)
    if missing:
        raise RuntimeError(f"Step02 QC file {qc_file} lacks {sorted(missing)}")
    parts = []
    for chunk in pd.read_csv(qc_file, sep="\t", usecols=use, chunksize=500_000, low_memory=False, compression="infer"):
        c = chunk["CHR"].map(clean_chrom)
        pos = pd.to_numeric(chunk["POS"], errors="coerce")
        hit = chunk[(c == chrom) & pos.between(start, end)]
        if len(hit):
            parts.append(hit)
    x = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=use)
    return canon_gwas(x, set(x.columns))


def load_gwas_region(P: dict, accession: str, locus: dict) -> tuple[pd.DataFrame, str, list[str]]:
    """(GWAS rows, input source, notes). Step06 regional file first, Step02 QC within the locus otherwise."""
    chrom, start, end = clean_chrom(locus["CHR"]), int(locus["LOCUS_START"]), int(locus["LOCUS_END"])
    regional = P["S06"] / accession / "regional_gwas" / f"{locus['LOCUS_ID']}.tsv.gz"
    qc_file = P["S02"] / accession / f"{accession}_GRCh38_QC.tsv.gz"
    notes = []
    if regional.is_file():
        x = read_tsv(regional)
        missing = {"CHR", "POS", "P", "EA", "NEA"} - set(x.columns)
        if missing:
            raise RuntimeError(f"{regional} lacks {sorted(missing)}")
        g = canon_gwas(x, set(x.columns))
        source = f"STEP06_REGIONAL ({rel(regional)})"
    elif nonempty(qc_file):
        g = read_qc_region(qc_file, chrom, start, end)
        source = f"STEP02_QC_FALLBACK ({rel(qc_file)}; regional_gwas/{locus['LOCUS_ID']}.tsv.gz absent)"
    else:
        return pd.DataFrame(columns=GWAS_COLUMNS), "", [f"no regional file ({rel(regional)}) and no Step02 QC file"]
    n0 = len(g)
    g = g[(g.CHR == chrom) & g.POS.between(start, end) & g.P.gt(0) & g.P.le(1)
          & g.EA.str.match(ACGT_PAT) & g.NEA.str.match(ACGT_PAT) & g.EA.ne(g.NEA)]
    if len(g) < n0:
        notes.append(f"{n0 - len(g)} GWAS rows outside the locus or with invalid P / alleles removed")
    g = g.assign(POS=g.POS.astype(int)).sort_values("P").drop_duplicates(["POS", "EA", "NEA"])
    return g.reset_index(drop=True), source, notes


def add_gwas_maf(P: dict, accession: str, locus: dict, g: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """GWAS allele frequency (only when the coloc mode needs it): from the rows themselves, else Step02 QC EAF."""
    if g["EAF"].notna().any():
        return g.assign(MAF_GWAS=np.minimum(g.EAF, 1 - g.EAF)), "GWAS EAF (same rows)"
    qc_file = P["S02"] / accession / f"{accession}_GRCh38_QC.tsv.gz"
    if not nonempty(qc_file):
        return g.assign(MAF_GWAS=np.nan), "none (no Step02 QC file)"
    q = read_qc_region(qc_file, clean_chrom(locus["CHR"]), int(locus["LOCUS_START"]), int(locus["LOCUS_END"]))
    q = q.dropna(subset=["EAF"]).drop_duplicates(["POS", "EA", "NEA"])
    m = g.drop(columns=["EAF"]).merge(q[["POS", "EA", "NEA", "EAF"]], on=["POS", "EA", "NEA"], how="left")
    return m.assign(MAF_GWAS=np.minimum(m.EAF, 1 - m.EAF)), f"Step02 QC EAF ({rel(qc_file)})"


# =============================================================================
# RESOURCES: files, headers, column mapping, metadata
# =============================================================================

def open_text(path: Path):
    return gzip.open(path, "rt", errors="replace") if str(path).endswith(".gz") else open(path, errors="replace")


def file_header(path: Path) -> tuple[list[str], bool, int]:
    """(column names, named header?, number of header / comment lines). Headerless -> c1, c2, ..."""
    skipped = 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##") or line.startswith(("track ", "browser ")) or not line.strip():
                skipped += 1
                continue
            cells = line.rstrip("\n\r").split("\t")
            if line.startswith("#"):
                return [c.lstrip("#").strip() for c in cells], True, skipped + 1
            numeric = sum(bool(re.fullmatch(r"[-+]?[0-9.eE+-]+", c)) for c in cells)
            known = {a for v in ROLE_ALIASES.values() for a in v}
            if numeric == 0 and sum(c.strip().lower() in known for c in cells) >= 2:
                return [c.strip() for c in cells], True, skipped + 1
            return [f"c{i + 1}" for i in range(len(cells))], False, skipped
    return [], False, skipped


def sample_rows(path: Path, names: list[str], n_header: int, limit: int) -> pd.DataFrame:
    rows = []
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            if i < n_header or line.startswith("#"):
                continue
            v = line.rstrip("\n\r").split("\t")
            if len(v) == len(names):
                rows.append(v)
            if len(rows) >= limit:
                break
    return pd.DataFrame(rows, columns=names)


def map_columns(cols: list[str], header: bool, override: dict) -> tuple[dict, str, list[str]]:
    """role -> column. Only from a named header or an explicit COLUMN_OVERRIDES entry."""
    mapping = {r: "" for r in ROLES}
    notes = []
    if header:
        low = {c.lower(): c for c in cols}
        for role, names in ROLE_ALIASES.items():
            mapping[role] = next((low[n] for n in names if n in low), "")
        if not mapping["ref"] and not mapping["alt"] and {"a1", "a2"} <= set(low):
            notes.append("allele columns a1 / a2 have no documented REF / ALT meaning (resolved by the FASTA check)")
        if mapping["trait"] == "" and mapping["gene"]:
            mapping["trait"] = mapping["gene"]
            notes.append(f"molecular trait = gene column {mapping['gene']}")
    source = "file header" if header else "none (headerless)"
    if override:
        bad = [k for k in override if k not in ROLES]
        if bad:
            raise SystemExit(f"COLUMN_OVERRIDES: unknown roles {bad}; use {', '.join(ROLES)}")
        missing = [v for v in override.values() if v and v not in cols]
        if missing:
            notes.append(f"COLUMN_OVERRIDES columns not in the file: {missing}")
        mapping.update({k: v for k, v in override.items() if v in cols or v == ""})
        source = ("file header + " if header else "") + "Step11_Config.py COLUMN_OVERRIDES"
    return mapping, source, notes


def metadata_records(paths: list[Path]) -> list[dict]:
    """Flat records from provider .metadata.txt (tab-separated) / .metadata.json / provenance.json."""
    records = []
    for path in paths:
        if not nonempty(path):
            continue
        if path.suffix == ".json":
            stack = [read_json(path)]
            while stack:
                item = stack.pop()
                if isinstance(item, dict):
                    flat = {k: v for k, v in item.items() if not isinstance(v, (dict, list))}
                    if flat:
                        records.append({**flat, "_FILE": rel(path)})
                    stack.extend(v for v in item.values() if isinstance(v, (dict, list)))
                elif isinstance(item, list):
                    stack.extend(item)
        else:
            try:
                t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment=None)
                records += [{**r, "_FILE": rel(path)} for r in t.to_dict(orient="records")]
            except Exception:
                continue
    return records


def key_norm(key: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(key).lower()).strip("_")


N_KEY = re.compile(r"^(n|n_?samples?|sample_?size|samplesize|number_of_(samples|subjects|individuals|donors|"
                   r"participants)|num_?samples|n_?subjects|n_?individuals|n_?donors)$")
BUILD_KEY = re.compile(r"^(genome_?build|assembly|build|reference_genome)$")
CONTEXT_KEY = re.compile(r"^(tissue|tissue_label|biosample|biosample_term_name|cell_?type|brain_?region|"
                         r"tissue_category|sample_group)$")


def metadata_for(records: list[dict], tokens: list[str]) -> dict:
    """Fields of the metadata records that mention one of the tokens (dataset id, file stem)."""
    tokens = [t for t in tokens if t]
    found = {"N": np.nan, "N_SOURCE": "", "BUILD": "", "CONTEXT": "", "MATCHED": 0}
    for r in records:
        text = " ".join(str(v) for v in r.values())
        if not any(t in text for t in tokens):
            continue
        found["MATCHED"] += 1
        for k, v in r.items():
            nk = key_norm(k)
            if N_KEY.match(nk) and not np.isfinite(found["N"]):
                x = safe_float(str(v).replace(",", "").strip())
                if np.isfinite(x) and x > 0 and float(x).is_integer():
                    found["N"], found["N_SOURCE"] = x, f"{r['_FILE']} field '{k}'"
            elif BUILD_KEY.match(nk) and not found["BUILD"] and usable(v):
                found["BUILD"] = str(v)
            elif CONTEXT_KEY.match(nk) and not found["CONTEXT"] and usable(v):
                found["CONTEXT"] = str(v)
    return found


def external_registry_records(cfg: dict) -> list[dict]:
    out = []
    for name in cfg["EXTERNAL_REGISTRIES"]:
        path = project_path(name)
        if nonempty(path):
            try:
                t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
                out += [{**r, "_FILE": rel(path)} for r in t.to_dict(orient="records")]
            except Exception as exc:
                print(f"[WARN] cannot read {path}: {exc}")
    return out


def eqtl_catalogue_table(cfg: dict) -> pd.DataFrame:
    frames = [pd.read_csv(project_path(p), sep="\t", dtype=str, keep_default_na=False)
              for p in cfg["EQTL_CATALOGUE_TABLES"] if nonempty(project_path(p))]
    if not frames:
        return pd.DataFrame(columns=["dataset_id"])
    x = pd.concat(frames, ignore_index=True)
    return x[x.get("dataset_id", pd.Series(dtype=str)).astype(str).ne("")].drop_duplicates("dataset_id", keep="last")


@lru_cache(maxsize=4)
def fasta_handle(path: str):
    import pysam
    return pysam.FastaFile(path)


def fasta_contig(fa, chrom: str) -> str | None:
    for name in (f"chr{chrom}", chrom):
        if name in fa.references:
            return name
    return None


def ref_match_fraction(fasta: str, chrom: pd.Series, pos: pd.Series, allele: pd.Series) -> tuple[float, int]:
    if not fasta:
        return np.nan, 0
    fa = fasta_handle(fasta)
    hits, n = 0, 0
    for c, p, a in zip(chrom, pos, allele):
        contig = fasta_contig(fa, clean_chrom(c))
        if contig is None or not isinstance(a, str) or not ACGT.match(a) or not np.isfinite(p):
            continue
        n += 1
        hits += fa.fetch(contig, int(p) - 1, int(p) - 1 + len(a)).upper() == a
    return (hits / n if n else np.nan), n


# =============================================================================
# RESOURCES: discovery (manifests, no recursive scan) and validation
# =============================================================================

def adsp_datasets(cfg: dict, records_ext: list[dict]) -> tuple[list[dict], list[dict]]:
    """Datasets + files of the ADSP FunGen xQTL archives from the Step00 manifests."""
    root = project_path(cfg["RESOURCE_ROOT"])
    datasets, files = [], []
    for modality, folder in cfg["ADSP_MODALITIES"].items():
        base = root / folder / cfg["ADSP_PROVIDER_DIR"]
        stem = f"ADSP_FunGen_xQTL.v1.{modality}.all"
        manifests = [p for p in (base / f"{stem}.manifest.txt", base / "all" / "ADSP_FunGen_xQTL" / "v1" /
                                 f"{stem}.manifest.txt") if nonempty(p)]
        meta_files = [p for p in (base / "all" / "ADSP_FunGen_xQTL" / "v1" / f"{stem}.metadata.txt",
                                  base / "all" / "ADSP_FunGen_xQTL" / "v1" / f"{stem}.metadata.json",
                                  base / "provenance.json") if nonempty(p)]
        listing = []
        if manifests:
            m = pd.read_csv(manifests[0], sep="\t", dtype=str, keep_default_na=False)
            name_col = next((c for c in m.columns if "file" in c.lower() and "name" in c.lower()), m.columns[0])
            md5_col = next((c for c in m.columns if "md5" in c.lower()), None)
            size_col = next((c for c in m.columns if "size" in c.lower()), None)
            for r in m.to_dict(orient="records"):
                name = str(r[name_col]).strip()
                if name.endswith(".bed.gz"):
                    listing.append((name, r.get(md5_col, "") if md5_col else "",
                                    r.get(size_col, "") if size_col else "", rel(manifests[0])))
        else:
            niagads = project_path(cfg["NIAGADS_MANIFEST"])
            if nonempty(niagads):         # fallback: the NG00184 file listing
                t = pd.read_csv(niagads, dtype=str, keep_default_na=False)
                for r in t.to_dict(orient="records"):
                    for v in r.values():
                        v = str(v)
                        if v.endswith(".bed.gz") and f"/{modality}/" in v and "ADSP_FunGen_xQTL/" in v:
                            listing.append((v[v.index("ADSP_FunGen_xQTL/"):], "", "", rel(niagads)))
        if not listing:
            datasets.append({"DATASET_ID": f"ADSP_{modality}", "SOURCE_STUDY": "ADSP FunGen xQTL (NG00184)",
                             "PROVIDER": "ADSP FunGen-AD / NIAGADS", "MODALITY": modality, "FILE_FORMAT": "bed3plus17_qtl",
                             "DISCOVERED_BY": "Step00 layout", "METADATA_FILES": ";".join(rel(p) for p in meta_files),
                             "_STATUS": "NOT_DOWNLOADED",
                             "_DETAIL": (f"{rel(base)} has no association files / manifest"
                                         + (" (provenance.json only)" if meta_files else "")
                                         + f" - Step00 --groups {'meqtl' if modality == 'mQTL' else 'molecular_qtl'}")})
            continue
        groups: dict[str, list] = {}
        for name, md5, size, manifest in listing:
            mm = ADSP_NAME.match(Path(name).name)
            if not mm:
                continue
            key = f"ADSP_{mm['cohort']}_{modality}_{mm['context']}_{mm['id']}"
            groups.setdefault(key, []).append((mm, name, md5, size, manifest))
        records = metadata_records(meta_files)
        for key, items in sorted(groups.items()):
            mm = items[0][0]
            for m2, name, md5, size, manifest in items:
                path = base / "all" / name
                present = path.is_file()
                files.append({"DATASET_ID": key, "CHROM": m2["chrom"], "KIND": m2["kind"], "PATH": rel(path),
                              "MANIFEST": manifest, "MANIFEST_MD5": md5, "MANIFEST_SIZE": size, "PRESENT": present,
                              "SIZE_OK": (not usable(size)) or (present and str(path.stat().st_size) == str(size).strip()),
                              "INDEXED": Path(str(path) + ".tbi").is_file()})
            meta = metadata_for(records + records_ext, [mm["id"], f"{mm['cohort']}.{mm['qtl']}.{mm['context']}"])
            datasets.append({"DATASET_ID": key, "SOURCE_STUDY": f"ADSP FunGen xQTL (NG00184) {mm['cohort']}",
                             "PROVIDER": "ADSP FunGen-AD / NIAGADS", "COHORT": mm["cohort"], "MODALITY": modality,
                             "CONTEXT": mm["context"] + (f" ({meta['CONTEXT']})" if meta["CONTEXT"] else ""),
                             "FILE_FORMAT": "bed3plus17_qtl", "DISCOVERED_BY": items[0][4],
                             "METADATA_FILES": ";".join(rel(p) for p in meta_files),
                             "_DECLARED_BUILD": meta["BUILD"] or ("hg38 (bed3plus17_qtl/hg38 directory)"
                                                                  if "/hg38/" in items[0][1] else ""),
                             "_META_N": meta["N"], "_META_N_SOURCE": meta["N_SOURCE"],
                             "_ID_TOKENS": [mm["id"]], "_MANIFEST_PATHS": [x[1] for x in items]})
    return datasets, files


def eqtl_format_datasets(cfg: dict, catalogue: pd.DataFrame, records_ext: list[dict]) -> tuple[list[dict], list[dict]]:
    datasets, files = [], []
    cat = catalogue.set_index("dataset_id") if len(catalogue) else pd.DataFrame()
    for spec in cfg["EQTL_CATALOGUE_FORMAT_DATASETS"]:
        did = str(spec["dataset_id"])
        path = project_path(spec["path"])
        crow = cat.loc[did].to_dict() if did in getattr(cat, "index", []) else {}
        meta_files = [p for p in (path.parent / "provenance.json",) if nonempty(p)]
        meta = metadata_for(metadata_records(meta_files) + records_ext, [did, path.name])
        n_cat = safe_float(crow.get("sample_size"))
        context = spec.get("context") or crow.get("tissue_label") or meta["CONTEXT"] or "unknown"
        if crow.get("condition_label") and usable(crow.get("condition_label")) and crow["condition_label"] != "naive":
            context += f" / {crow['condition_label']}"
        tmeta = project_path(spec["trait_metadata"]) if spec.get("trait_metadata") else None
        if not path.is_file():
            status, detail = "NOT_DOWNLOADED", f"{rel(path)} missing (Step00 --groups molecular_qtl)"
        else:
            status, detail = "", ""
        files.append({"DATASET_ID": did, "CHROM": "all", "KIND": "all", "PATH": rel(path), "MANIFEST": "Step11_Config.py",
                      "MANIFEST_MD5": "", "MANIFEST_SIZE": "", "PRESENT": path.is_file(), "SIZE_OK": True,
                      "INDEXED": Path(str(path) + ".tbi").is_file()})
        datasets.append({"DATASET_ID": did, "SOURCE_STUDY": spec.get("source_study") or crow.get("study_label", did),
                         "PROVIDER": spec.get("provider", "eQTL Catalogue format"),
                         "COHORT": spec.get("source_study") or crow.get("study_label", ""), "MODALITY": spec["modality"],
                         "CONTEXT": context, "FILE_FORMAT": "eQTL Catalogue summary statistics",
                         "DISCOVERED_BY": "Step11_Config.py EQTL_CATALOGUE_FORMAT_DATASETS",
                         "METADATA_FILES": ";".join([rel(p) for p in meta_files] +
                                                    ([rel(tmeta)] if tmeta and tmeta.is_file() else [])),
                         "TRAIT_METADATA": rel(tmeta) if tmeta and tmeta.is_file() else "",
                         "EQTL_CATALOGUE_ID": did, "_DECLARED_BUILD": spec.get("build", "GRCh38"),
                         "_META_N": n_cat if np.isfinite(n_cat) else meta["N"],
                         "_META_N_SOURCE": ("eQTL Catalogue tabix_ftp_paths*.tsv sample_size" if np.isfinite(n_cat)
                                            else meta["N_SOURCE"]),
                         "_STATUS": status, "_DETAIL": detail, "_ID_TOKENS": [did]})
    return datasets, files


def registry_extra_datasets(cfg: dict, known_paths: set[str], records_ext: list[dict]) -> tuple[list[dict], list[dict]]:
    """Tabix files named in resources/*/resource_registry.tsv-style tables that the manifests did not cover."""
    datasets, files = [], []
    for r in records_ext:
        for k, v in r.items():
            v = str(v).strip()
            if not (v.endswith((".tsv.gz", ".bed.gz")) and ("path" in k.lower() or "file" in k.lower())):
                continue
            path = project_path(v)
            if rel(path) in known_paths or not path.is_file() or not Path(str(path) + ".tbi").is_file():
                continue
            if not str(path.resolve()).startswith(str(project_path(cfg["RESOURCE_ROOT"]).resolve())):
                continue         # e.g. eQTL Catalogue files of Step09 / Step10
            kind = next((QTL_TYPE_NAMES[x] for x in QTL_TYPE_NAMES if x in v.lower().split("/")), "")
            ds_key = next((str(r[c]) for c in r if key_norm(c) in {"dataset", "dataset_id", "id", "name"}
                           and usable(r[c])), path.name.split(".")[0])
            did = f"REG_{ds_key}"
            meta = metadata_for([r], [v])
            known_paths.add(rel(path))
            files.append({"DATASET_ID": did, "CHROM": "all", "KIND": "all", "PATH": rel(path), "MANIFEST": r["_FILE"],
                          "MANIFEST_MD5": "", "MANIFEST_SIZE": "", "PRESENT": True, "SIZE_OK": True, "INDEXED": True})
            datasets.append({"DATASET_ID": did, "SOURCE_STUDY": ds_key, "PROVIDER": r["_FILE"], "MODALITY": kind,
                             "CONTEXT": meta["CONTEXT"] or "unknown", "FILE_FORMAT": "from header",
                             "DISCOVERED_BY": r["_FILE"], "METADATA_FILES": r["_FILE"],
                             "_DECLARED_BUILD": meta["BUILD"], "_META_N": meta["N"], "_META_N_SOURCE": meta["N_SOURCE"],
                             "_ID_TOKENS": [ds_key]})
    return datasets, files


def provider_finemapping(d: dict, cfg: dict) -> str:
    """SuSiE / fSuSiE fine-mapping supplied by the provider (credible sets + variant LBF) for this dataset?"""
    did = d.get("EQTL_CATALOGUE_ID", "")
    if did:
        base = project_path(cfg["EQTL_CATALOGUE_SUSIE_DIR"]) / did
        cs = list(base.glob("*.credible_sets.tsv.gz")) if base.is_dir() else []
        lbf = list(base.glob("*.lbf_variable.sorted.tsv.gz")) if base.is_dir() else []
        if cs and lbf and Path(str(lbf[0]) + ".tbi").is_file():
            return f"SuSiE (eQTL Catalogue credible sets + indexed LBF in {rel(base)})"
        return "none downloaded (eQTL Catalogue SuSiE files not in " + rel(base) + ")"
    terms = ("susie", "finemap", "fine_map", "credible", "lbf")
    hits = [p for p in d.get("_MANIFEST_PATHS", []) if any(t in p.lower() for t in terms)]
    return (f"listed in the manifest ({len(hits)} files) - not used: no compatible LBF reader" if hits
            else "none (the archive holds association statistics only)")


def validate_dataset(d: dict, files: pd.DataFrame, cfg: dict) -> dict:
    """Columns, coordinates, alleles, build, statistics and QTL N of one dataset -> registry row."""
    row = {c: "" for c in REGISTRY_COLUMNS}
    row.update({k: v for k, v in d.items() if k in REGISTRY_COLUMNS})
    row.update(SCANNED_UTC=utcnow(), STEP11_VERSION=VERSION)
    did = d["DATASET_ID"]
    warnings, details = [], []
    f = files[files.DATASET_ID == did]
    row["N_FILES_EXPECTED"], row["N_FILES_PRESENT"] = len(f), int(f.PRESENT.sum()) if len(f) else 0
    row["N_FILES_INDEXED"] = int((f.PRESENT & f.INDEXED).sum()) if len(f) else 0
    row["CHROMOSOMES"] = ",".join(sorted({str(c) for c in f[f.PRESENT].CHROM}, key=lambda c: (len(c), c)))

    def done(status: str, detail: str = "") -> dict:
        if detail:
            details.append(detail)
        row["VALIDATION_STATUS"] = status
        row["VALIDATION_DETAILS"] = "; ".join(details)
        row["WARNINGS"] = "; ".join(warnings)
        row["READY"] = "YES" if status in {"READY", "READY_WITH_WARNINGS"} else "NO"
        return row

    if did in {str(x) for x in cfg["DISABLED_DATASETS"]}:
        return done("DISABLED", "listed in Step11_Config.py DISABLED_DATASETS")
    if d.get("_STATUS"):
        return done(d["_STATUS"], d.get("_DETAIL", ""))
    if row["N_FILES_PRESENT"] == 0:
        return done("NOT_DOWNLOADED", "no file of the manifest is on disk")
    if row["N_FILES_PRESENT"] < len(f) or not f[f.PRESENT].SIZE_OK.all():
        return done("INCOMPLETE_DOWNLOAD", f"{row['N_FILES_PRESENT']}/{len(f)} files present, "
                                          f"{int((~f.SIZE_OK).sum())} with a size different from the manifest "
                                          "- re-run Step00 for this group")
    if row["N_FILES_INDEXED"] < row["N_FILES_PRESENT"]:
        return done("INDEX_MISSING", f"{row['N_FILES_PRESENT'] - row['N_FILES_INDEXED']} files without .tbi")

    present = f[f.PRESENT]
    example = present[present.KIND.isin(["snp", "all"])]
    example_path = project_path((example if len(example) else present).iloc[0]["PATH"])
    row["EXAMPLE_FILE"] = rel(example_path)
    cols, header, n_header = file_header(example_path)
    row["HEADER"], row["COLUMNS"] = ("named" if header else "none (headerless)"), ",".join(cols)
    override = (cfg["COLUMN_OVERRIDES"].get(did) or cfg["COLUMN_OVERRIDES"].get(d.get("FILE_FORMAT", "")) or {})
    mapping, source, notes = map_columns(cols, header, override)
    warnings += notes
    row["COLUMN_SOURCE"] = source
    for role in ROLES:
        row[f"COL_{role.upper()}"] = mapping[role]
    # the other files of the dataset must share the header
    for p in present.PATH.iloc[1:]:
        other, other_header, _ = file_header(project_path(p))
        if other != cols:
            return done("SCHEMA_DIFFERS_BETWEEN_FILES", f"{p} has columns {other[:6]}... vs {cols[:6]}...")
    if not header and not override:
        first = sample_rows(example_path, cols, n_header, 1)
        shown = " | ".join(f"{c}={v}" for c, v in zip(cols, first.iloc[0])) if len(first) else "no data row"
        return done("NEEDS_COLUMN_MAPPING", "headerless file: the column meaning is not documented in the file; "
                                            "add COLUMN_OVERRIDES for this dataset (or 'bed3plus17_qtl' for all ADSP "
                                            f"tracks) in Step11_Config.py. First row: {shown[:400]}")
    if not mapping["trait"] or not mapping["p"]:
        return done("NEEDS_COLUMN_MAPPING", f"no molecular-trait / P column among {cols}")
    if not (mapping["pos"] or mapping["start"] or mapping["variant"]) or not (mapping["chrom"] or mapping["variant"]):
        return done("NEEDS_COLUMN_MAPPING", "no chromosome / position / variant column")

    x = sample_rows(example_path, cols, n_header, int(cfg["VALIDATION_ROWS"]))
    if x.empty:
        return done("NO_ROWS", f"{rel(example_path)} has no data rows")
    vid = x[mapping["variant"]].astype(str).str.extract(VARRE) if mapping["variant"] else pd.DataFrame()
    vid_ok = len(vid) and vid[1].notna().mean() >= 0.95
    chrom = x[mapping["chrom"]].map(clean_chrom) if mapping["chrom"] else vid[0].map(clean_chrom)

    fasta = str(project_path(cfg["REFERENCE_FASTA"])) if nonempty(project_path(cfg["REFERENCE_FASTA"])) and \
        Path(str(project_path(cfg["REFERENCE_FASTA"])) + ".fai").is_file() else ""
    # coordinates: candidate 1-based positions, verified against the variant id or the reference genome
    candidates = {}
    if mapping["pos"]:
        candidates["POS column (1-based)"] = ("pos", pd.to_numeric(x[mapping["pos"]], errors="coerce"))
    if mapping["start"]:
        s = pd.to_numeric(x[mapping["start"]], errors="coerce")
        candidates["BED chromStart + 1 (0-based start)"] = ("start+1", s + 1)
        candidates["chromStart holds the 1-based position (non-standard BED)"] = ("start", s)
    if mapping["end"]:
        candidates["chromEnd (1-based, SNPs)"] = ("end", pd.to_numeric(x[mapping["end"]], errors="coerce"))
    chosen, evidence = None, ""

    # alleles
    low = {c.lower(): c for c in cols}
    if mapping["ref"] and mapping["alt"]:
        ref, alt, allele_source = x[mapping["ref"]].str.upper(), x[mapping["alt"]].str.upper(), \
            f"columns {mapping['ref']} / {mapping['alt']}"
    elif vid_ok:
        ref, alt, allele_source = vid[2].str.upper(), vid[3].str.upper(), f"variant id {mapping['variant']} (chr:pos:ref:alt)"
    elif header and {"a1", "a2"} <= set(low):
        # a1 / a2 have no fixed meaning: the one that is the GRCh38 base at the verified position is REF
        if not fasta:
            return done("NEEDS_ALLELES", "allele columns a1 / a2 without documented REF / ALT meaning and no GRCh38 "
                                         "FASTA to resolve them (Step00 --groups genome, or COLUMN_OVERRIDES ref / alt)")
        a1, a2 = x[low["a1"]].str.upper(), x[low["a2"]].str.upper()
        snp12 = a1.str.len().eq(1) & a2.str.len().eq(1)
        best = (0.0, None, None)
        for label, (code, pos) in candidates.items():
            for ref_col, alt_col, allele in ((low["a1"], low["a2"], a1), (low["a2"], low["a1"], a2)):
                frac, _ = ref_match_fraction(fasta, chrom[snp12], pos[snp12], allele[snp12])
                if np.isfinite(frac) and frac > best[0]:
                    best = (frac, (label, code), (ref_col, alt_col))
        if best[0] < cfg["REF_CHECK_MIN_MATCH"]:
            return done("NEEDS_ALLELES", f"neither a1 nor a2 is the GRCh38 base at any position convention "
                                         f"(best {best[0]:.1%})")
        mapping["ref"], mapping["alt"] = best[2]
        row["COL_REF"], row["COL_ALT"] = best[2]
        ref, alt = x[mapping["ref"]].str.upper(), x[mapping["alt"]].str.upper()
        allele_source = f"columns {best[2][0]} (= GRCh38 REF in {best[0]:.1%} of SNPs) / {best[2][1]}"
        chosen, evidence = best[1], f"REF matches GRCh38 at this position in {best[0]:.1%} of SNP rows"
    else:
        return done("NEEDS_ALLELES", "no REF / ALT columns and no chr:pos:ref:alt variant id")

    snp = ref.str.len().eq(1) & alt.str.len().eq(1)
    if vid_ok and chosen is None:
        vpos = pd.to_numeric(vid[1], errors="coerce")
        for label, (code, pos) in candidates.items():
            rows = snp if code == "end" else pd.Series(True, index=x.index)
            frac = float((pos[rows] == vpos[rows]).mean()) if rows.any() else 0.0
            if frac >= 0.99:
                chosen, evidence = (label, code), f"matches the variant-id position in {frac:.1%} of {int(rows.sum())} rows"
                break
        if chosen is None and not candidates:
            chosen, evidence = ("variant id position (1-based)", "vid"), "position taken from the variant id"
    if chosen is None and fasta:
        best = (0.0, None)
        for label, (code, pos) in candidates.items():
            frac, n = ref_match_fraction(fasta, chrom[snp], pos[snp], ref[snp])
            if np.isfinite(frac) and frac > best[0]:
                best = (frac, (label, code))
        if best[0] >= cfg["REF_CHECK_MIN_MATCH"]:
            chosen, evidence = best[1], f"REF matches GRCh38 at this position in {best[0]:.1%} of SNP rows"
    if chosen is None:
        return done("COORDINATES_AMBIGUOUS", f"cannot verify the position convention among {list(candidates)} "
                                             f"(no variant id agreement{'' if fasta else ', no FASTA'})")
    if chosen[1] == "end" and not snp.all():
        warnings.append("chromEnd convention verified on SNPs only; indel positions use chromStart + 1")
    row["COORDINATES"] = f"{chosen[0]}: {evidence}"
    pos = vpos if chosen[1] == "vid" else candidates[chosen[0]][1]

    # build + allele orientation against GRCh38
    frac_ref, n_ref = ref_match_fraction(fasta, chrom[snp], pos[snp], ref[snp]) if fasta else (np.nan, 0)
    declared = str(d.get("_DECLARED_BUILD", ""))
    if fasta and n_ref:
        if frac_ref >= cfg["REF_CHECK_MIN_MATCH"]:
            row["BUILD"], row["BUILD_EVIDENCE"] = "GRCh38", f"REF = GRCh38 base in {frac_ref:.1%} of {n_ref} SNPs"
        else:
            frac_alt, _ = ref_match_fraction(fasta, chrom[snp], pos[snp], alt[snp])
            if frac_alt >= cfg["REF_CHECK_MIN_MATCH"] and not (mapping["ref"] and mapping["alt"] and header
                                                              and mapping["ref"].lower() == "ref"):
                ref, alt = alt, ref
                allele_source += " [ORDER SWAPPED: the second allele is the GRCh38 reference]"
                warnings.append(f"allele order swapped: the second allele is the GRCh38 reference "
                                f"({frac_alt:.1%} of SNPs); beta still refers to the provider's effect allele")
                row["BUILD"], row["BUILD_EVIDENCE"] = "GRCh38", f"second allele = GRCh38 base in {frac_alt:.1%} of SNPs"
            else:
                return done("BUILD_OR_COORDINATE_MISMATCH",
                            f"REF matches GRCh38 in only {frac_ref:.1%} of {n_ref} SNPs (declared: {declared or '-'})")
    elif re.search(r"(?i)grch38|hg38", declared):
        row["BUILD"], row["BUILD_EVIDENCE"] = "GRCh38", f"declared ({declared}); not sequence-verified (no FASTA)"
        warnings.append("GRCh38 declared but not verified against the reference genome (Step00 --groups genome)")
    else:
        return done("BUILD_UNKNOWN", f"no GRCh38 evidence (declared: {declared or '-'}) and no FASTA to verify")
    row["ALLELE_SOURCE"] = allele_source
    row["EFFECT_ALLELE"] = (f"{mapping['alt']} (named effect / alt allele)" if mapping["alt"] else
                            "ALT of the variant id (assumed; coloc.abf uses beta^2 / varbeta and MAF, so the "
                            "direction does not change PP)")
    bad_alleles = int((~ref.str.match(ACGT_PAT) | ~alt.str.match(ACGT_PAT)).sum())
    if bad_alleles:
        warnings.append(f"{bad_alleles} rows with non-ACGT alleles are excluded")

    # statistics
    p = pd.to_numeric(x[mapping["p"]], errors="coerce")
    if not p.between(0, 1).mean() >= 0.99:
        return done("INVALID_STATS", f"P column {mapping['p']}: {float((~p.between(0, 1)).mean()):.1%} outside 0-1")
    has_beta = bool(mapping["beta"] and mapping["se"])
    if has_beta:
        se = pd.to_numeric(x[mapping["se"]], errors="coerce")
        if (se <= 0).mean() > 0.01:
            return done("INVALID_STATS", f"SE column {mapping['se']} has non-positive values")
        bad, n_ok = stats_inconsistency(x[mapping["beta"]], x[mapping["se"]], p)
        if np.isfinite(bad) and bad > cfg["STATS_MAX_INCONSISTENT"]:
            return done("STATS_INCONSISTENT", f"P is not reproduced by beta / SE in {bad:.1%} of {n_ok} rows "
                                              f"({mapping['beta']}, {mapping['se']}, {mapping['p']}: wrong columns?)")
        details.append(f"beta / SE reproduce P in {1 - bad:.1%} of {n_ok} rows" if np.isfinite(bad) else "")
    if mapping["freq"]:
        fr = pd.to_numeric(x[mapping["freq"]], errors="coerce")
        if not fr.between(0, 1).mean() >= 0.99:
            return done("INVALID_STATS", f"frequency column {mapping['freq']} outside 0-1")
        is_maf = mapping["freq"].lower() == "maf"
        if is_maf and (fr > 0.5 + 1e-9).mean() > 0.01:
            warnings.append(f"{mapping['freq']} has values > 0.5 (treated as allele frequency, folded)")
            is_maf = False
        row["FREQ_KIND"] = "MAF" if is_maf else f"allele frequency ({mapping['freq']}; folded to MAF)"
    elif mapping["ac"] and mapping["an"]:
        row["FREQ_KIND"] = f"AC / AN ({mapping['ac']} / {mapping['an']}; folded to MAF)"
    else:
        row["FREQ_KIND"] = "none"
    sdy = cfg["QTL_SDY"].get(did)
    if sdy:
        row["QTL_SDY"], row["QTL_SDY_SOURCE"] = str(sdy[0]), str(sdy[1])
    if not has_beta and row["FREQ_KIND"] == "none":
        return done("NO_EFFECT_OR_FREQUENCY", "coloc.abf needs beta + SE (+ QTL MAF or sdY) or P + QTL MAF")
    if has_beta and row["FREQ_KIND"] == "none" and not sdy:
        return done("NEEDS_QTL_MAF_OR_SDY", "beta / SE without a QTL allele frequency: set QTL_SDY for this "
                                            "dataset (documented trait SD) - the GWAS MAF is never substituted")

    # QTL sample size
    an_half = np.nan
    if mapping["an"]:
        an = pd.to_numeric(x[mapping["an"]], errors="coerce").dropna()
        if len(an):
            an_half = float(np.median(an) / 2)
            row["AN_HALF_MEDIAN"] = f"{an_half:g}"
    if mapping["n"]:
        n = pd.to_numeric(x[mapping["n"]], errors="coerce")
        if n.gt(0).mean() >= 0.99:
            row["QTL_N"], row["QTL_N_SOURCE"], row["QTL_N_PER_TRAIT"] = (
                f"{float(n.median()):g} (median)", f"file column {mapping['n']} (per trait)", "YES")
    if not row["QTL_N"] and np.isfinite(safe_float(d.get("_META_N"))):
        row["QTL_N"], row["QTL_N_SOURCE"], row["QTL_N_PER_TRAIT"] = f"{safe_float(d['_META_N']):g}", d["_META_N_SOURCE"], "NO"
    if not row["QTL_N"] and did in cfg["QTL_N"]:
        n, src = cfg["QTL_N"][did]
        row["QTL_N"], row["QTL_N_SOURCE"], row["QTL_N_PER_TRAIT"] = f"{float(n):g}", f"Step11_Config.py QTL_N: {src}", "NO"
    if not row["QTL_N"] and np.isfinite(an_half):
        lowcols = {c.lower() for c in cols}
        an_vals = pd.to_numeric(x[mapping["an"]], errors="coerce").dropna()
        even = bool(len(an_vals) and (an_vals % 2 == 0).all())
        if cfg["ALLOW_VERIFIED_AN_HALF"] and EQTL_CATALOGUE_SCHEMA <= lowcols and even:
            row["QTL_N"], row["QTL_N_PER_TRAIT"] = "AN/2", "YES"
            row["QTL_N_SOURCE"] = ("AN / 2 per trait: eQTL Catalogue schema documents 'an' as the total allele number "
                                   f"(diploid: 2 x genotyped samples); AN even in all {len(an_vals)} checked rows, "
                                   f"median AN/2 = {an_half:g}")
        else:
            warnings.append(f"AN present (median AN/2 = {an_half:g}) but not accepted as N "
                            f"({'schema not eQTL Catalogue' if not EQTL_CATALOGUE_SCHEMA <= lowcols else 'odd AN values' if not even else 'ALLOW_VERIFIED_AN_HALF = False'})")
    if row["QTL_N"] and row["QTL_N"] != "AN/2" and np.isfinite(an_half):
        ratio = an_half / safe_float(str(row["QTL_N"]).split()[0])
        if np.isfinite(ratio) and not 0.8 <= ratio <= 1.05:
            warnings.append(f"metadata N {row['QTL_N']} vs median AN/2 {an_half:g} (ratio {ratio:.2f})")
    if not row["QTL_N"]:
        return done("NEEDS_QTL_N", "no sample size in the file, the provider metadata, the eQTL Catalogue tables or "
                                   "Step11_Config.py QTL_N (never guessed)")

    row["PROVIDER_FINEMAPPING"] = provider_finemapping(d, cfg)
    if row["PROVIDER_FINEMAPPING"].startswith("SuSiE"):
        row["METHOD"] = "coloc.bf_bf in Step09 / Step10 (provider SuSiE, multi-signal) - not re-run as ABF here"
        return done("PROVIDER_SUSIE_USE_STEP09",
                    f"add '{did}' to DATASET_IDS in Step09_Config.py (or Step10_Config.py) for multi-signal coloc")
    row["METHOD"] = "coloc.abf (single causal variant per trait and locus)"
    return done("READY_WITH_WARNINGS" if warnings else "READY")


def scan_resources(cfg: dict, verbose: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    records_ext = external_registry_records(cfg)
    catalogue = eqtl_catalogue_table(cfg)
    a_ds, a_files = adsp_datasets(cfg, records_ext)
    e_ds, e_files = eqtl_format_datasets(cfg, catalogue, records_ext)
    known = {f["PATH"] for f in a_files + e_files}
    r_ds, r_files = registry_extra_datasets(cfg, known, records_ext)
    files = pd.DataFrame(a_files + e_files + r_files, columns=FILES_COLUMNS)
    rows = []
    for d in a_ds + e_ds + r_ds:
        try:
            rows.append(validate_dataset(d, files, cfg))
        except Exception as exc:
            row = {c: "" for c in REGISTRY_COLUMNS}
            row.update({k: v for k, v in d.items() if k in REGISTRY_COLUMNS})
            row.update(VALIDATION_STATUS="VALIDATION_ERROR", VALIDATION_DETAILS=f"{type(exc).__name__}: {exc}",
                       READY="NO", SCANNED_UTC=utcnow(), STEP11_VERSION=VERSION)
            rows.append(row)
        if verbose:
            r = rows[-1]
            print(f"  {r['DATASET_ID']:<52} {r['MODALITY']:<6} {r['VALIDATION_STATUS']}", flush=True)
    return pd.DataFrame(rows, columns=REGISTRY_COLUMNS), files


def load_registry(cfg: dict) -> pd.DataFrame:
    path = project_path(cfg["REGISTRY"])
    if not nonempty(path):
        return pd.DataFrame(columns=REGISTRY_COLUMNS)
    reg = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "DATASET_ID" not in reg.columns or "VALIDATION_STATUS" not in reg.columns:
        print(f"[INFO] {path} is from an older Step11 version (manual enable / disable) - rebuilding it")
        return pd.DataFrame(columns=REGISTRY_COLUMNS)
    return reg


def load_files_table(cfg: dict) -> pd.DataFrame:
    return read_tsv(project_path(cfg["FILES_TABLE"]), dtype=str, keep_default_na=False)


# =============================================================================
# DUPLICATES OF STEP09 / STEP10
# =============================================================================

def step09_10_snapshots(P: dict) -> pd.DataFrame:
    frames = []
    for step, folder in (("Step09", P["S09"]), ("Step10", P["S10"])):
        x = read_tsv(folder / "qtl_contexts_snapshot.tsv", dtype=str, keep_default_na=False)
        if len(x):
            frames.append(x.assign(_STEP=step))
    snap = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["dataset_id", "_STEP"])
    for column in ("dataset_id", "study_label", "tissue_label", "sample_group", "QTL_TYPE"):
        if column not in snap.columns:
            snap[column] = ""
    return snap


def duplicate_check(reg: pd.DataFrame, snap: pd.DataFrame) -> pd.DataFrame:
    """DUPLICATE_OF (same dataset id) and INDEPENDENT_EVIDENCE (same cohort + modality + tissue) per dataset."""
    reg = reg.copy()
    reg["DUPLICATE_OF"], reg["INDEPENDENT_EVIDENCE"], reg["OVERLAP_NOTE"] = "", "YES", ""
    if snap.empty:
        reg["OVERLAP_NOTE"] = "no Step09 / Step10 snapshot for this phenotype / panel (not yet planned)"
        return reg
    ids = dict(zip(snap["dataset_id"].astype(str), snap["_STEP"]))
    for i, r in reg.iterrows():
        did = str(r.get("EQTL_CATALOGUE_ID", "") or "")
        if did and did in ids:
            reg.at[i, "DUPLICATE_OF"] = f"{ids[did]} {did}"
            reg.at[i, "INDEPENDENT_EVIDENCE"] = "NO"
            continue
        cohort, modality = normalize_text(r.get("COHORT", "")), str(r.get("MODALITY", ""))
        context = normalize_text(str(r.get("CONTEXT", "")).split(" (")[0])
        if not cohort:
            continue
        same = snap[snap["study_label"].map(lambda s: cohort.split()[0] in normalize_text(s))
                    & snap["QTL_TYPE"].astype(str).eq(modality)]
        if context:
            same = same[(same["tissue_label"] + " " + same["sample_group"]).map(lambda s: context in normalize_text(s))]
        if len(same):
            reg.at[i, "INDEPENDENT_EVIDENCE"] = "NO"
            reg.at[i, "OVERLAP_NOTE"] = (f"same cohort / modality / tissue as {same.iloc[0]['_STEP']} "
                                         f"{same.iloc[0]['dataset_id']} ({same.iloc[0].get('study_label', '')}) - "
                                         "tested, but not independent evidence")
    return reg


def plan_datasets(P: dict, cfg: dict, reg: pd.DataFrame) -> pd.DataFrame:
    """Registry + this phenotype's duplicate check + type selection -> PLAN_STATUS per dataset."""
    reg = duplicate_check(reg, step09_10_snapshots(P))
    wanted = {t.lower() for t in cfg["QTL_TYPES"]}
    status = []
    for _, r in reg.iterrows():
        if r["READY"] != "YES":
            status.append(r["VALIDATION_STATUS"])
        elif str(r["MODALITY"]).lower() not in wanted:
            status.append("TYPE_NOT_SELECTED")
        elif r["DUPLICATE_OF"] and not cfg["ALLOW_DUPLICATES_OF_STEP09_10"]:
            status.append(f"DUPLICATE_OF_{r['DUPLICATE_OF'].split()[0].upper()}")
        else:
            status.append("PLANNED")
    reg["PLAN_STATUS"] = status
    return reg


def planned_hash(planned: pd.DataFrame) -> str:
    volatile = {"SCANNED_UTC", "WARNINGS", "VALIDATION_DETAILS", "STEP11_VERSION"}
    cols = [c for c in REGISTRY_COLUMNS + ["DUPLICATE_OF", "INDEPENDENT_EVIDENCE"]
            if c in planned.columns and c not in volatile]
    return sha1_short("\n".join("|".join(str(r[c]) for c in cols)
                                for r in planned.sort_values("DATASET_ID").to_dict(orient="records"))) \
        if len(planned) else "none"


# =============================================================================
# ONE LOCUS
# =============================================================================

def tabix_rows(path: Path, chrom: str, start: int, end: int):
    import pysam
    with pysam.TabixFile(str(path)) as tab:
        contigs = set(tab.contigs)
        contig = next((z for z in (chrom, f"chr{chrom}") if z in contigs), None)
        if contig is None:
            return
        for raw in tab.fetch(contig, max(start - 2, 0), end + 1):
            if not raw.startswith("#"):
                yield raw.split("\t")


def dataset_files_for(files: pd.DataFrame, did: str, chrom: str) -> list[Path]:
    f = files[(files.DATASET_ID == did) & files.PRESENT.astype(str).isin(["True", "true", "1"])]
    f = f[f.CHROM.astype(str).isin([chrom, "all"])]
    return [project_path(p) for p in f.PATH]


def region_qtl(r: dict, paths: list[Path], chrom: str, start: int, end: int, cfg: dict) -> tuple[pd.DataFrame, pd.Series, int]:
    """(rows of candidate traits, min P per trait, region rows). Two tabix passes: P per trait, then the candidates."""
    cols = r["COLUMNS"].split(",")
    idx = {role: cols.index(r[f"COL_{role.upper()}"]) for role in ROLES if r.get(f"COL_{role.upper()}")}
    best: dict[str, float] = {}
    n_rows = 0
    for path in paths:
        for v in tabix_rows(path, chrom, start, end):
            if len(v) != len(cols):
                continue
            n_rows += 1
            p = safe_float(v[idx["p"]])
            t = v[idx["trait"]]
            if np.isfinite(p) and (t not in best or p < best[t]):
                best[t] = p
    min_p = pd.Series(best, dtype=float).sort_values()
    traits = set(min_p[min_p <= cfg["CANDIDATE_P"]].head(int(cfg["MAX_TRAITS"])).index)
    if not traits:
        return pd.DataFrame(), min_p, n_rows
    rows = []
    for path in paths:
        for v in tabix_rows(path, chrom, start, end):
            if len(v) == len(cols) and v[idx["trait"]] in traits:
                rows.append(v)
                if len(rows) > int(cfg["MAX_REGION_ROWS"]):
                    raise RuntimeError(f"TOO_MANY_REGION_ROWS (> {cfg['MAX_REGION_ROWS']}); raise MAX_REGION_ROWS")
    return normalize_qtl(pd.DataFrame(rows, columns=cols), r, chrom, start, end), min_p, n_rows


def normalize_qtl(df: pd.DataFrame, r: dict, chrom: str, start: int, end: int) -> pd.DataFrame:
    def col(role):
        c = r.get(f"COL_{role.upper()}", "")
        return df[c] if c and c in df.columns else pd.Series(np.nan, index=df.index)

    vid = col("variant").astype(str).str.extract(VARRE) if r.get("COL_VARIANT") else None
    coords = str(r["COORDINATES"])
    if coords.startswith("POS column"):
        pos = pd.to_numeric(col("pos"), errors="coerce")
    elif coords.startswith("BED chromStart + 1"):
        pos = pd.to_numeric(col("start"), errors="coerce") + 1
    elif coords.startswith("chromStart holds"):
        pos = pd.to_numeric(col("start"), errors="coerce")
    elif coords.startswith("chromEnd"):
        ref_len = col("ref").astype(str).str.len() if r.get("COL_REF") else 1
        pos = np.where(pd.Series(ref_len, index=df.index).eq(1), pd.to_numeric(col("end"), errors="coerce"),
                       pd.to_numeric(col("start"), errors="coerce") + 1)
        pos = pd.Series(pos, index=df.index)
    else:
        pos = pd.to_numeric(vid[1], errors="coerce")
    if r.get("COL_REF") and r.get("COL_ALT"):
        ref, alt = col("ref").astype(str).str.upper(), col("alt").astype(str).str.upper()
    else:
        ref, alt = vid[2].str.upper(), vid[3].str.upper()
    if "ORDER SWAPPED" in str(r.get("ALLELE_SOURCE", "")):
        ref, alt = alt, ref
    freq = pd.to_numeric(col("freq"), errors="coerce")
    if not r.get("COL_FREQ") and r.get("COL_AC") and r.get("COL_AN"):
        freq = pd.to_numeric(col("ac"), errors="coerce") / pd.to_numeric(col("an"), errors="coerce")
    gene = col("gene").astype(str) if r.get("COL_GENE") else pd.Series("", index=df.index)
    out = pd.DataFrame({
        "TRAIT": col("trait").astype(str), "GENE": gene,
        "CHR": chrom, "POS": pos, "REF": ref, "ALT": alt,
        "QTL_ORIGINAL_ID": col("variant").astype(str) if r.get("COL_VARIANT") else
        (chrom + ":" + pos.astype("Int64").astype(str) + ":" + ref + ":" + alt),
        "P_QTL": pd.to_numeric(col("p"), errors="coerce"), "BETA_QTL": pd.to_numeric(col("beta"), errors="coerce"),
        "SE_QTL": pd.to_numeric(col("se"), errors="coerce"),
        "MAF_QTL": np.minimum(freq, 1 - freq), "N_QTL_ROW": pd.to_numeric(col("n"), errors="coerce"),
        "AN_QTL": pd.to_numeric(col("an"), errors="coerce")})
    out = out[out.POS.between(start, end)]
    out["POS"] = out["POS"].astype(int)
    return out.reset_index(drop=True)


@lru_cache(maxsize=4)
def trait_metadata(path: str) -> pd.DataFrame:
    p = project_path(path)
    if not nonempty(p):
        return pd.DataFrame()
    x = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    key = next((c for c in ("phenotype_id", "molecular_trait_id", "peak_id") if c in x.columns), x.columns[0])
    return x.drop_duplicates(key).set_index(key)


def harmonise(q: pd.DataFrame, g: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Every QTL row of one trait with its GWAS partner (if any), match status and exclusion reason."""
    m = q.merge(g[["POS", "EA", "NEA", "SNPID", "P", "BETA", "SE", "MAF_GWAS"]].rename(
        columns={"P": "P_GWAS", "BETA": "BETA_GWAS", "SE": "SE_GWAS", "SNPID": "GWAS_ORIGINAL_ID"}),
        on="POS", how="left")
    aligned = (m.REF == m.NEA) & (m.ALT == m.EA)
    swapped = (m.REF == m.EA) & (m.ALT == m.NEA)
    m["MATCH_STATUS"] = np.select([aligned, swapped, m.EA.isna()], ["ALIGNED", "SWAPPED", "NO_GWAS_VARIANT"],
                                  "ALLELE_MISMATCH")
    # several GWAS rows at one position (multi-allelic): keep the allele-matched one
    m = m.sort_values("MATCH_STATUS", key=lambda s: s.map({"ALIGNED": 0, "SWAPPED": 1, "ALLELE_MISMATCH": 2,
                                                          "NO_GWAS_VARIANT": 3}))
    m = m.drop_duplicates(["POS", "REF", "ALT"], keep="first").copy()
    m.loc[m.MATCH_STATUS.isin(["ALLELE_MISMATCH", "NO_GWAS_VARIANT"]),
          ["GWAS_ORIGINAL_ID", "P_GWAS", "BETA_GWAS", "SE_GWAS", "MAF_GWAS"]] = np.nan
    # GWAS effect re-expressed for the QTL ALT allele
    m["BETA_GWAS_ALIGNED"] = np.where(m.MATCH_STATUS == "SWAPPED", -m.BETA_GWAS, m.BETA_GWAS)
    m["HARMONISED_ID"] = m.CHR.astype(str) + ":" + m.POS.astype(str) + ":" + m.REF + ":" + m.ALT
    reason = pd.Series("", index=m.index)
    reason[~m.REF.str.match(ACGT_PAT) | ~m.ALT.str.match(ACGT_PAT)] = "NON_ACGT_ALLELES"
    indel = (m.REF.str.len() != 1) | (m.ALT.str.len() != 1)
    if not cfg["INCLUDE_INDELS"]:
        reason[(reason == "") & indel] = "INDEL_EXCLUDED_BY_CONFIG"
    reason[(reason == "") & (m.MATCH_STATUS == "NO_GWAS_VARIANT")] = "NO_GWAS_VARIANT_AT_POSITION"
    reason[(reason == "") & (m.MATCH_STATUS == "ALLELE_MISMATCH")] = \
        np.where(indel[(reason == "") & (m.MATCH_STATUS == "ALLELE_MISMATCH")],
                 "INDEL_ALLELE_REPRESENTATION_DIFFERS", "ALLELE_MISMATCH")
    pal = np.array([is_palindromic(a, b) for a, b in zip(m.REF, m.ALT)], dtype=bool)
    m["PALINDROMIC"] = pal
    if cfg["PALINDROMIC_POLICY"] == "drop":
        reason[(reason == "") & pal] = "PALINDROMIC_DROPPED"
    elif cfg["PALINDROMIC_POLICY"] == "drop_high_maf":
        high = ~(m.MAF_QTL <= cfg["PALINDROMIC_MAX_MAF"])
        reason[(reason == "") & pal & high] = "PALINDROMIC_HIGH_OR_UNKNOWN_MAF"
    reason[(reason == "") & ~(m.P_QTL.gt(0) & m.P_QTL.le(1))] = "QTL_P_INVALID"
    m["EXCLUSION_REASON"] = reason
    return m


def choose_modes(m: pd.DataFrame, gmeta: dict, r: dict) -> tuple[str, str, str, pd.Series]:
    """(GWAS mode, QTL mode, NOT_TESTED reason, per-variant exclusion) for the matched rows of one trait."""
    usable_rows = m.EXCLUSION_REASON.eq("")
    reason = pd.Series("", index=m.index)
    g_beta = m.BETA_GWAS.notna() & m.SE_GWAS.gt(0)
    if g_beta[usable_rows].mean() >= 0.95:
        gmode = "beta"
        reason[usable_rows & ~g_beta] = "GWAS_BETA_SE_MISSING"
        if gmeta["GWAS_TYPE"] == "quant":
            reason[usable_rows & (reason == "") & ~m.MAF_GWAS.between(0, 0.5, inclusive="right")] = "GWAS_MAF_MISSING"
    else:
        gmode = "pvalue+MAF"
        reason[usable_rows & ~m.MAF_GWAS.between(0, 0.5, inclusive="right")] = "GWAS_MAF_MISSING"
    q_beta = m.BETA_QTL.notna() & m.SE_QTL.gt(0)
    q_maf = m.MAF_QTL.between(0, 0.5, inclusive="right")
    sdy = safe_float(r.get("QTL_SDY"))
    if r.get("COL_BETA") and r.get("COL_SE") and q_beta[usable_rows].mean() >= 0.95:
        if np.isfinite(sdy) and sdy > 0:
            qmode = "beta+sdY"
            reason[usable_rows & (reason == "") & ~q_beta] = "QTL_BETA_SE_MISSING"
        else:
            qmode = "beta+MAF"
            reason[usable_rows & (reason == "") & ~q_beta] = "QTL_BETA_SE_MISSING"
            reason[usable_rows & (reason == "") & ~q_maf] = "QTL_MAF_MISSING"
    elif q_maf[usable_rows].any():
        qmode = "pvalue+MAF"
        reason[usable_rows & (reason == "") & ~q_maf] = "QTL_MAF_MISSING"
    else:
        return gmode, "", "NOT_TESTED_QTL_MAF_MISSING (no QTL beta/SE and no QTL frequency; GWAS MAF not used)", reason
    return gmode, qmode, "", reason


R_ABF = r'''
suppressPackageStartupMessages(library(coloc))
a <- commandArgs(trailingOnly=TRUE)
x <- read.delim(gzfile(a[1]), stringsAsFactors=FALSE, check.names=FALSE,
                colClasses=c(CANDIDATE_ID="character", SNP="character"))
cands <- read.delim(a[2], stringsAsFactors=FALSE, check.names=FALSE, colClasses="character")
out_file <- a[3]
p1 <- as.numeric(a[4]); p2 <- as.numeric(a[5]); p12s <- as.numeric(c(a[6], a[7], a[8]))
gtype <- a[9]; gs <- suppressWarnings(as.numeric(a[10])); gn <- as.numeric(a[11])
rows <- list()
for (i in seq_len(nrow(cands))) {
  cid <- cands$CANDIDATE_ID[i]
  r <- list(CANDIDATE_ID=cid, R_STATUS="OK", MESSAGE="", PP_H0=NA, PP_H1=NA, PP_H2=NA, PP_H3=NA, PP_H4=NA,
            PP_H4_P12_LOW=NA, PP_H4_P12_HIGH=NA, BEST_SNP="", BEST_SNP_PP_H4=NA, N_SNPS=0)
  res <- tryCatch({
    d <- x[x$CANDIDATE_ID == cid, , drop=FALSE]
    r$N_SNPS <- nrow(d)
    gm <- cands$GWAS_MODE[i]; qm <- cands$QTL_MODE[i]
    d1 <- list(type=gtype, N=gn, snp=d$SNP)
    if (gtype == "cc") d1$s <- gs
    if (gm == "beta") { d1$beta <- d$BETA_GWAS; d1$varbeta <- d$SE_GWAS^2; if (gtype == "quant") d1$MAF <- d$MAF_GWAS
    } else { d1$pvalues <- d$P_GWAS; d1$MAF <- d$MAF_GWAS }
    d2 <- list(type="quant", N=as.numeric(cands$QTL_N[i]), snp=d$SNP)
    if (qm == "beta+sdY") { d2$beta <- d$BETA_QTL; d2$varbeta <- d$SE_QTL^2; d2$sdY <- as.numeric(cands$QTL_SDY[i])
    } else if (qm == "beta+MAF") { d2$beta <- d$BETA_QTL; d2$varbeta <- d$SE_QTL^2; d2$MAF <- d$MAF_QTL
    } else { d2$pvalues <- d$P_QTL; d2$MAF <- d$MAF_QTL }
    z <- lapply(p12s, function(p12) suppressWarnings(suppressMessages(coloc::coloc.abf(d1, d2, p1=p1, p2=p2, p12=p12))))
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
    # Step11 provenance / harmonisation / assumption flags
    "METHOD_FAMILY", "NOT_TESTED_REASON", "SOURCE_STUDY", "COHORT", "QTL_N_SOURCE", "QTL_MAF_SOURCE",
    "GWAS_INPUT_SOURCE", "GWAS_TYPE", "GWAS_N", "GWAS_N_SOURCE", "GWAS_CASE_FRACTION", "GWAS_N_EFF_NOT_USED",
    "GWAS_MAF_SOURCE", "N_ALIGNED", "N_SWAPPED", "N_EXCLUDED", "EXCLUSIONS", "GWAS_LEAD_IN_SHARED",
    "QTL_LEAD_IN_SHARED", "GWAS_N_CREDIBLE_SETS", "GWAS_N_INDEPENDENT_SIGNALS", "MULTI_SIGNAL_FLAG",
    "COMPLEX_LD_REGION", "ABF_ASSUMPTION_WARNING", "DUPLICATE_OF", "INDEPENDENT_EVIDENCE", "TRAIT_ANNOTATION",
]


def classify(pp: list[float], strong: float, suggestive: float) -> dict:
    dominant = int(np.nanargmax(pp))
    pp4 = pp[4]
    label = ("STRONG_SHARED_SIGNAL" if pp4 >= strong else "SUGGESTIVE_SHARED_SIGNAL" if pp4 >= suggestive else
             {0: "NEITHER_ASSOCIATED_FAVORED", 1: "GWAS_ONLY_FAVORED", 2: "QTL_ONLY_FAVORED",
              3: "BOTH_ASSOCIATED_DIFFERENT_SIGNALS_FAVORED", 4: "SHARED_SIGNAL_BELOW_THRESHOLD"}[dominant])
    return {"COLOC_CLASS": label, "DOMINANT_HYPOTHESIS": f"H{dominant}",
            "COLOCALIZED_STRONG": "YES" if pp4 >= strong else "NO", "SHARED_SIGNAL": "YES" if pp4 >= suggestive else "NO"}


def locus_flags(locus: dict, n_cs: int | None, cfg: dict) -> dict:
    chrom, start, end = clean_chrom(locus["CHR"]), int(locus["LOCUS_START"]), int(locus["LOCUS_END"])
    n_ind = safe_float(locus.get("N_INDEPENDENT_SIGNALS"))
    n_ind = int(n_ind) if np.isfinite(n_ind) else 0
    complex_ld = [name for name, (c, s, e) in cfg["COMPLEX_LD_REGIONS"].items()
                  if str(c) == chrom and start <= e and end >= s]
    multi = []
    if n_ind > 1:
        multi.append(f"{n_ind} independent GWAS lead variants (Step05)")
    if n_cs is not None and n_cs > 1:
        multi.append(f"{n_cs} GWAS SuSiE credible sets (Step06)")
    warning = "; ".join(multi + [f"complex LD region {x}" for x in complex_ld])
    return {"GWAS_N_CREDIBLE_SETS": n_cs if n_cs is not None else "NO_GWAS_CREDIBLE_SET (ABF still applicable)",
            "GWAS_N_INDEPENDENT_SIGNALS": n_ind, "MULTI_SIGNAL_FLAG": "YES" if multi else "NO",
            "COMPLEX_LD_REGION": ",".join(complex_ld),
            "ABF_ASSUMPTION_WARNING": (warning + " - coloc.abf assumes one causal variant per trait; a multi-signal "
                                       "method is needed to confirm") if warning else ""}


def process_locus(row: pd.Series, args, cfg: dict, planned: pd.DataFrame, files: pd.DataFrame, rscript: str) -> str:
    P = project_paths(args)
    accession, locus_id = str(row["STUDY_ACCESSION"]), str(row["LOCUS_ID"])
    study_dir = P["S06"] / accession
    out = P["OUT"] / accession / "loci" / locus_id
    done_file = out / DONE_FILE
    run_hash = sha1_short(planned_hash(planned) + results_hash(cfg) + VERSION)
    previous = read_json(done_file)
    if previous and not args.force and previous.get("RUN_HASH") == run_hash:
        print(f"[RESUME] {accession} {locus_id}: {previous.get('STATUS')} (same datasets and parameters)")
        return str(previous.get("STATUS"))
    out.mkdir(parents=True, exist_ok=True)
    (out / FAILED_FILE).unlink(missing_ok=True)
    work = out / "work"
    work.mkdir(exist_ok=True)
    started = time.time()
    lt = loci_table(study_dir)
    hit = lt[lt.LOCUS_ID.astype(str) == locus_id]
    if hit.empty:
        raise RuntimeError(f"{locus_id} not in {study_dir / 'loci_definition.tsv'}")
    locus = hit.iloc[0].to_dict()
    chrom, start, end = clean_chrom(locus["CHR"]), int(locus["LOCUS_START"]), int(locus["LOCUS_END"])
    gmeta = gwas_meta(P, accession, cfg)
    flags = locus_flags(locus, gwas_cs_counts(study_dir).get(locus_id), cfg)
    banner(f"STEP 11 - {accession} {locus_id} chr{chrom}:{start:,}-{end:,} - {len(planned)} QTL datasets")
    region = {"CHR": chrom, "LOCUS_START": start, "LOCUS_END": end,
              **{k: gmeta[k] for k in ("GWAS_TYPE", "GWAS_N", "GWAS_N_SOURCE", "GWAS_CASE_FRACTION",
                                       "GWAS_N_EFF_NOT_USED")}, **flags}

    def finish(status: str, results: pd.DataFrame, audit: list[dict], harm: list[pd.DataFrame], extra: dict) -> str:
        write_tsv(results.reindex(columns=RESULT_COLUMNS), out / RESULTS_FILE)
        write_tsv(pd.DataFrame(audit), out / AUDIT_FILE)
        if harm:
            write_tsv(pd.concat(harm, ignore_index=True), out / HARMONISED_FILE)
        tested = results[results["COLOC_TESTED"].astype(str) == "YES"] if len(results) else results
        h4 = pd.to_numeric(tested.get("PP_H4", pd.Series(dtype=float)), errors="coerce")
        write_json(done_file, {"VERSION": VERSION, "STEP": "Step11", "STATUS": status, "STUDY_ACCESSION": accession,
                               "LOCUS_ID": locus_id, "RUN_HASH": run_hash, "N_DATASETS": int(len(planned)),
                               "N_CANDIDATES": int(len(results)), "N_TESTED": int(len(tested)),
                               "N_STRONG": int((h4 >= cfg["STRONG_PP4"]).sum()), **region, **extra,
                               "SECONDS": round(time.time() - started, 1), "COMPLETED_UTC": utcnow()})
        if not args.keep_work:
            shutil.rmtree(work, ignore_errors=True)
        print(f"[OK] {accession} {locus_id}: {status}; {len(results)} candidates, {len(tested)} tested "
              f"({time.time() - started:.0f}s)")
        return status

    if gmeta["GWAS_META_REASON"]:
        return finish(gmeta["GWAS_META_REASON"].split(" (")[0], pd.DataFrame(), [], [],
                      {"REASON": gmeta["GWAS_META_REASON"]})
    g, gsource, gnotes = load_gwas_region(P, accession, locus)
    region["GWAS_INPUT_SOURCE"] = gsource
    if g.empty:
        return finish("NOT_TESTED_NO_GWAS_INPUT" if not gsource else "NOT_TESTED_NO_GWAS_VARIANTS", pd.DataFrame(), [],
                      [], {"REASON": "; ".join(gnotes) or "no GWAS variants in the locus"})
    g_has_beta = g.BETA.notna().mean() >= 0.95 and g.SE.gt(0).mean() >= 0.95
    bad, n_ok = stats_inconsistency(g.BETA, g.SE, g.P) if g_has_beta else (np.nan, 0)
    if np.isfinite(bad) and bad > cfg["STATS_MAX_INCONSISTENT"]:
        gnotes.append(f"GWAS beta / SE do not reproduce P in {bad:.1%} of {n_ok} rows - P + MAF mode used")
        g_has_beta = False
        g = g.assign(BETA=np.nan, SE=np.nan)
    gmaf_source = "not needed (case-control GWAS with beta / SE)"
    if gmeta["GWAS_TYPE"] == "quant" or not g_has_beta:
        g, gmaf_source = add_gwas_maf(P, accession, locus, g)
    else:
        g = g.assign(MAF_GWAS=np.nan)
    region["GWAS_MAF_SOURCE"] = gmaf_source
    gwas_lead = g.sort_values("P").iloc[0]

    long_rows, cands, audit, harm = [], [], [], []
    for r in planned.to_dict(orient="records"):
        did = r["DATASET_ID"]
        entry = {"DATASET_ID": did, "QTL_TYPE": r["MODALITY"], "CONTEXT": r["CONTEXT"], "N_FILES": 0,
                 "N_REGION_ROWS": 0, "N_TRAITS_IN_REGION": 0, "N_CANDIDATES": 0, "MIN_P": np.nan, "STATUS": "OK"}
        try:
            paths = dataset_files_for(files, did, chrom)
            entry["N_FILES"] = len(paths)
            if not paths:
                entry["STATUS"] = "NOT_TESTED_NO_FILE_FOR_CHROMOSOME"
                audit.append(entry)
                continue
            q, min_p, n_rows = region_qtl(r, paths, chrom, start, end, cfg)
            entry.update(N_REGION_ROWS=n_rows, N_TRAITS_IN_REGION=len(min_p),
                         MIN_P=float(min_p.iloc[0]) if len(min_p) else np.nan)
            if n_rows == 0:
                entry["STATUS"] = "NO_QTL_ROWS_IN_LOCUS"
                audit.append(entry)
                continue
            if q.empty:
                entry["STATUS"] = f"NO_TRAIT_WITH_P_LE_{cfg['CANDIDATE_P']:g}"
                audit.append(entry)
                continue
            traits = list(dict.fromkeys(min_p[min_p.index.isin(set(q.TRAIT))].index))
            entry["N_CANDIDATES"] = len(traits)
            tmeta = trait_metadata(r["TRAIT_METADATA"]) if r.get("TRAIT_METADATA") else pd.DataFrame()
            for trait in traits:
                qt = q[q.TRAIT == trait].sort_values("P_QTL").drop_duplicates(["POS", "REF", "ALT"])
                cid = sha1_short(f"{accession}|{locus_id}|{did}|{trait}", 20)
                m = harmonise(qt, g, cfg)
                n_matched = int(m.EXCLUSION_REASON.eq("").sum())
                gmode, qmode, not_tested, var_reason = choose_modes(m, gmeta, r)
                m.loc[(m.EXCLUSION_REASON == "") & (var_reason != ""), "EXCLUSION_REASON"] = var_reason
                use = m[m.EXCLUSION_REASON.eq("")]
                # QTL sample size of this trait
                if r["QTL_N"] == "AN/2":
                    an = qt.AN_QTL.dropna()
                    qtl_n = float(np.median(an) / 2) if len(an) else np.nan
                elif r["QTL_N_PER_TRAIT"] == "YES" and qt.N_QTL_ROW.notna().any():
                    qtl_n = float(qt.N_QTL_ROW.median())
                else:
                    qtl_n = safe_float(str(r["QTL_N"]).split()[0])
                ann = ""
                if len(tmeta) and trait in tmeta.index:
                    t = tmeta.loc[trait]
                    ann = ";".join(f"{k}={t[k]}" for k in ("gene_id", "gene_name", "chromosome", "phenotype_pos",
                                                          "start", "end") if k in t.index and usable(t[k]))
                gene_text = " ".join([trait, str(qt.GENE.iloc[0]), ann])
                gene = ENSG.search(gene_text)
                excluded = m[m.EXCLUSION_REASON.ne("")]
                base = {"CANDIDATE_ID": cid, "STUDY_ACCESSION": accession, "LOCUS_ID": locus_id, **region,
                        "QTL_TYPE": r["MODALITY"], "STUDY_LABEL": r["SOURCE_STUDY"], "SOURCE_STUDY": r["SOURCE_STUDY"],
                        "COHORT": r["COHORT"], "TISSUE": r["CONTEXT"], "SAMPLE_GROUP": "", "CONDITION": "",
                        "DATASET_ID": did, "MOLECULAR_TRAIT_ID": trait, "GENE_ID": gene.group(0) if gene else "",
                        "TRAIT_ANNOTATION": ann, "QTL_N": qtl_n, "QTL_N_SOURCE": r["QTL_N_SOURCE"],
                        "QTL_MAF_SOURCE": r["FREQ_KIND"], "DISCOVERY_MIN_P": float(min_p[trait]),
                        "N_DISCOVERY_ROWS": int(len(qt)), "N_GWAS_VARIANTS": int(len(g)),
                        "N_QTL_VARIANTS": int(len(qt)), "N_COMMON_VARIANTS": int(len(use)),
                        "N_ALIGNED": int((use.MATCH_STATUS == "ALIGNED").sum()),
                        "N_SWAPPED": int((use.MATCH_STATUS == "SWAPPED").sum()), "N_EXCLUDED": int(len(excluded)),
                        "EXCLUSIONS": ";".join(f"{k}={v}" for k, v in excluded.EXCLUSION_REASON.value_counts().items()),
                        "GWAS_LEAD_IN_SHARED": "YES" if (use.POS == int(gwas_lead.POS)).any() else "NO",
                        "QTL_LEAD_IN_SHARED": "YES" if len(use) and use.P_QTL.min() <= qt.P_QTL.min() else "NO",
                        "N_SIGNAL_PAIRS": 1, "N_GWAS_SIGNALS": 1, "N_QTL_SIGNALS": 1, "GWAS_MODE": gmode,
                        "QTL_MODE": qmode, "COLOC_METHOD": "coloc.abf", "METHOD_FAMILY": "ABF_SINGLE_CAUSAL_VARIANT",
                        "DUPLICATE_OF": r.get("DUPLICATE_OF", ""), "INDEPENDENT_EVIDENCE": r.get("INDEPENDENT_EVIDENCE", ""),
                        "COLOC_TESTED": "NO", "GWAS_CS_INDEXES": "", "QTL_CS_INDEXES": ""}
                harm.append(m.assign(CANDIDATE_ID=cid, DATASET_ID=did, USED_IN_COLOC=m.EXCLUSION_REASON.eq(""))[
                    ["CANDIDATE_ID", "DATASET_ID", "TRAIT", "HARMONISED_ID", "QTL_ORIGINAL_ID", "GWAS_ORIGINAL_ID",
                     "MATCH_STATUS", "PALINDROMIC", "EXCLUSION_REASON", "USED_IN_COLOC", "P_GWAS", "BETA_GWAS",
                     "BETA_GWAS_ALIGNED", "SE_GWAS", "MAF_GWAS", "P_QTL", "BETA_QTL", "SE_QTL", "MAF_QTL"]])
                why = ""
                if n_matched < int(cfg["MIN_SHARED"]):
                    why = (f"NOT_TESTED_INSUFFICIENT_VARIANT_OVERLAP ({n_matched} allele-matched variants; "
                           f"need >= {cfg['MIN_SHARED']})")
                why = why or not_tested
                if not why and not (np.isfinite(qtl_n) and qtl_n > 0):
                    why = "NOT_TESTED_QTL_N_MISSING (trait has no usable N)"
                if not why and (len(use) < int(cfg["MIN_SHARED"]) or
                                len(use) / max(1, min(len(g), len(qt))) < float(cfg["MIN_OVERLAP"])):
                    why = (f"NOT_TESTED_INSUFFICIENT_VARIANT_OVERLAP ({len(use)} shared; need >= {cfg['MIN_SHARED']} "
                           f"and >= {cfg['MIN_OVERLAP']:.0%} of the smaller side)")
                if why:
                    cands.append({**base, "STATUS": why.split(" (")[0], "NOT_TESTED_REASON": why})
                    continue
                long_rows.append(pd.DataFrame({
                    "CANDIDATE_ID": cid, "SNP": use.HARMONISED_ID.values, "P_GWAS": use.P_GWAS.values,
                    "BETA_GWAS": use.BETA_GWAS_ALIGNED.values, "SE_GWAS": use.SE_GWAS.values,
                    "MAF_GWAS": use.MAF_GWAS.values, "P_QTL": use.P_QTL.values, "BETA_QTL": use.BETA_QTL.values,
                    "SE_QTL": use.SE_QTL.values, "MAF_QTL": use.MAF_QTL.values}))
                cands.append({**base, "STATUS": "PENDING_R", "QTL_SDY": r.get("QTL_SDY") or "NA"})
        except Exception as exc:
            entry["STATUS"] = f"NOT_TESTED_RESOURCE_ERROR: {type(exc).__name__}: {str(exc)[:200]}"
        audit.append(entry)

    extra = {"GWAS_INPUT_SOURCE": gsource, "GWAS_NOTES": "; ".join(gnotes)}
    results = pd.DataFrame(cands)
    if results.empty:
        return finish("NO_QTL_CANDIDATES", pd.DataFrame(columns=RESULT_COLUMNS), audit, harm, extra)
    for column in RESULT_COLUMNS:
        if column not in results.columns:
            results[column] = np.nan if column.startswith(("PP_", "BEST_SNP_PP")) else ""
    results = results.astype({c: object for c in ("COLOC_CLASS", "DOMINANT_HYPOTHESIS", "COLOCALIZED_STRONG",
                                                  "SHARED_SIGNAL", "PRIOR_ROBUST_STRONG", "BEST_GWAS_VARIANT",
                                                  "BEST_QTL_VARIANT", "ERROR", "NOT_TESTED_REASON")})
    todo = results[results["STATUS"] == "PENDING_R"]
    if len(todo):
        write_tsv(pd.concat(long_rows, ignore_index=True), work / "abf_input.tsv.gz")
        write_tsv(todo[["CANDIDATE_ID", "QTL_N", "QTL_SDY", "GWAS_MODE", "QTL_MODE"]].fillna("NA"),
                  work / "candidates.tsv")
        worker = work / "coloc_abf.R"
        worker.write_text(R_ABF, encoding="utf-8")
        s = gmeta["GWAS_CASE_FRACTION"]
        command = [rscript, str(worker), str(work / "abf_input.tsv.gz"), str(work / "candidates.tsv"),
                   str(work / "abf_results.tsv"), str(cfg["P1"]), str(cfg["P2"]), str(cfg["P12"]),
                   str(cfg["P12_LOW"]), str(cfg["P12_HIGH"]), gmeta["GWAS_TYPE"],
                   str(s if np.isfinite(s) else "NA"), str(gmeta["GWAS_N"])]
        proc = subprocess.run(command, text=True, capture_output=True, timeout=int(cfg["R_TIMEOUT"]))
        (out / "coloc.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
        if proc.returncode != 0:
            raise RuntimeError(f"R coloc.abf worker failed (log {out / 'coloc.log'}):\n{proc.stderr[-3000:]}")
        abf = read_tsv(work / "abf_results.tsv", dtype=str).set_index("CANDIDATE_ID")
        for i, r in todo.iterrows():
            a = abf.loc[r["CANDIDATE_ID"]] if r["CANDIDATE_ID"] in abf.index else None
            if a is None or a["R_STATUS"] != "OK":
                msg = "" if a is None else str(a.get("MESSAGE", ""))
                results.at[i, "STATUS"] = "NOT_TESTED_COLOC_ERROR"
                results.at[i, "NOT_TESTED_REASON"] = f"NOT_TESTED_COLOC_ERROR ({msg})"
                results.at[i, "ERROR"] = msg
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
            results.at[i, "COLOC_TESTED"] = "YES"
            results.at[i, "STATUS"] = "TESTED"
    order = pd.to_numeric(results["PP_H4"], errors="coerce").fillna(-1)
    results = results.assign(_O=order).sort_values(["_O", "DISCOVERY_MIN_P"], ascending=[False, True]).drop(columns="_O")
    status = "COMPLETE" if (results["COLOC_TESTED"] == "YES").any() else "COMPLETE_NO_TESTABLE_CANDIDATE"
    return finish(status, results, audit, harm, extra)


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


def run_locus(row: pd.Series, args, cfg, planned, files, rscript) -> bool:
    try:
        process_locus(row, args, cfg, planned, files, rscript)
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


def print_registry(reg: pd.DataFrame, plan: bool = False) -> None:
    if reg.empty:
        print("Resource table: empty - run --scan")
        return
    cols = ["DATASET_ID", "MODALITY", "CONTEXT", "QTL_N", "N_FILES_PRESENT", "VALIDATION_STATUS"]
    if plan:
        cols += ["PLAN_STATUS", "INDEPENDENT_EVIDENCE"]
    show = reg[cols].copy()
    show["CONTEXT"] = show["CONTEXT"].astype(str).str.slice(0, 28)
    show["QTL_N"] = show["QTL_N"].astype(str).str.slice(0, 14)
    with pd.option_context("display.width", 250, "display.max_rows", None, "display.max_colwidth", 60):
        print(show.to_string(index=False))
    for modality in ("sQTL", "pQTL", "haQTL", "caQTL", "mQTL"):
        sub = reg[reg.MODALITY == modality]
        ready = int((sub.READY == "YES").sum())
        state = "usable" if ready else "NOT usable"
        print(f"  {modality:<6} {state:<11} {ready} of {len(sub)} datasets ready"
              + (f"  ({', '.join(sorted(set(sub.VALIDATION_STATUS)))})" if len(sub) and not ready else ""))
    problems = reg[reg.READY != "YES"]
    if len(problems):
        print("\nNot ready (why):")
        for _, r in problems.iterrows():
            print(f"  {r['DATASET_ID']}: {r['VALIDATION_STATUS']} - {str(r['VALIDATION_DETAILS'])[:220]}")


def scan_mode(args, cfg: dict) -> None:
    banner(f"STEP 11 - RESOURCE TABLE (manifest-driven; root {project_path(cfg['RESOURCE_ROOT'])})")
    fasta = project_path(cfg["REFERENCE_FASTA"])
    print(f"GRCh38 FASTA : {fasta if nonempty(fasta) else 'MISSING - builds are declared, not verified (Step00 --groups genome)'}")
    reg, files = scan_resources(cfg)
    path = project_path(cfg["REGISTRY"])
    write_tsv(reg, path)
    write_tsv(files, project_path(cfg["FILES_TABLE"]))
    print()
    print_registry(reg)
    print(f"\nResource table : {path}\nFiles table    : {project_path(cfg['FILES_TABLE'])}")
    print(f"NEXT:\n  {base_command(project_paths(args))} --check-resources")


def current_registry(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    reg = load_registry(cfg)
    files = load_files_table(cfg)
    if reg.empty or files.empty:
        print("[INFO] no resource table yet - scanning the resources now")
        reg, files = scan_resources(cfg, verbose=False)
        write_tsv(reg, project_path(cfg["REGISTRY"]))
        write_tsv(files, project_path(cfg["FILES_TABLE"]))
    return reg, files


def locus_plan(P: dict, cfg: dict) -> pd.DataFrame:
    rows = []
    meta_table = step02_metadata(P)
    for accession in step06_studies(P["S06"]):
        study_dir = P["S06"] / accession
        gm = gwas_meta(P, accession, cfg, meta_table)
        cs = gwas_cs_counts(study_dir)
        qc = P["S02"] / accession / f"{accession}_GRCh38_QC.tsv.gz"
        for locus in loci_table(study_dir).to_dict(orient="records"):
            lid = str(locus["LOCUS_ID"])
            regional = study_dir / "regional_gwas" / f"{lid}.tsv.gz"
            source = "STEP06_REGIONAL" if regional.is_file() else "STEP02_QC_FALLBACK" if nonempty(qc) else "MISSING"
            flags = locus_flags(locus, cs.get(lid), cfg)
            reason = gm["GWAS_META_REASON"] or ("NOT_TESTED_NO_GWAS_INPUT" if source == "MISSING" else "")
            rows.append({"STUDY_ACCESSION": accession, "LOCUS_ID": lid, "CHR": locus["CHR"],
                         "LOCUS_START": locus["LOCUS_START"], "LOCUS_END": locus["LOCUS_END"],
                         "GWAS_INPUT": source, "GWAS_TYPE": gm["GWAS_TYPE"], "GWAS_N": gm["GWAS_N"],
                         "GWAS_CASE_FRACTION": gm["GWAS_CASE_FRACTION"], "GWAS_N_SOURCE": gm["GWAS_N_SOURCE"],
                         **{k: flags[k] for k in ("GWAS_N_CREDIBLE_SETS", "GWAS_N_INDEPENDENT_SIGNALS",
                                                  "MULTI_SIGNAL_FLAG", "COMPLEX_LD_REGION")},
                         "TESTABLE": "NO" if reason else "YES", "REASON": reason})
    return pd.DataFrame(rows)


def check_resources_mode(args, cfg) -> None:
    P = project_paths(args)
    problems = []
    banner("STEP 11 - RESOURCE CHECK - external molecular QTL")
    print(f"Config     : {cfg['_FILE'] or 'built-in defaults'}")
    print(f"QTL types  : {', '.join(cfg['QTL_TYPES'])}   CANDIDATE_P {cfg['CANDIDATE_P']}   MIN_SHARED "
          f"{cfg['MIN_SHARED']}   MIN_OVERLAP {cfg['MIN_OVERLAP']}   palindromic: {cfg['PALINDROMIC_POLICY']}")
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
    for module in ("pysam", "scipy"):
        found = importlib.util.find_spec(module) is not None
        print(f"  {'OK     ' if found else 'MISSING'}  python {module}")
        if not found and module == "pysam":
            problems.append("python pysam")
    fasta = project_path(cfg["REFERENCE_FASTA"])
    print(f"  {'OK     ' if nonempty(fasta) else 'MISSING'}  GRCh38 FASTA   {rel(fasta)}"
          + ("" if nonempty(fasta) else "  (builds / alleles not sequence-verified)"))

    banner(f"DATASETS - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    reg, files = current_registry(cfg)
    plan = plan_datasets(P, cfg, reg)
    print_registry(plan, plan=True)
    for _, r in plan[plan.OVERLAP_NOTE.astype(str).ne("")].iterrows():
        print(f"  NOTE {r['DATASET_ID']}: {r['OVERLAP_NOTE']}")
    if not (plan.PLAN_STATUS == "PLANNED").any():
        problems.append("no dataset is ready and planned (see VALIDATION_STATUS / PLAN_STATUS above)")
    P["OUT"].mkdir(parents=True, exist_ok=True)
    write_tsv(plan, P["OUT"] / "resource_plan.tsv")

    banner(f"GWAS AND LOCI - Step02 metadata + Step06 loci_definition ({P['S06']})")
    lp = locus_plan(P, cfg)
    if lp.empty:
        problems.append(f"no Step06 loci - run Step06{panel_option(P['panel'])} first")
    else:
        write_tsv(lp, P["OUT"] / "locus_plan.tsv")
        per = lp.groupby("STUDY_ACCESSION").agg(
            LOCI=("LOCUS_ID", "size"), TESTABLE=("TESTABLE", lambda s: int((s == "YES").sum())),
            REGIONAL=("GWAS_INPUT", lambda s: int((s == "STEP06_REGIONAL").sum())),
            QC_FALLBACK=("GWAS_INPUT", lambda s: int((s == "STEP02_QC_FALLBACK").sum())),
            NO_GWAS_CS=("GWAS_N_CREDIBLE_SETS", lambda s: int(s.astype(str).str.startswith("NO_").sum())),
            MULTI_SIGNAL=("MULTI_SIGNAL_FLAG", lambda s: int((s == "YES").sum())),
            TYPE=("GWAS_TYPE", "first"), N=("GWAS_N", "first"), CASE_FRACTION=("GWAS_CASE_FRACTION", "first"),
            REASON=("REASON", "first")).reset_index()
        per["CASE_FRACTION"] = pd.to_numeric(per["CASE_FRACTION"], errors="coerce").round(3)
        with pd.option_context("display.width", 250, "display.max_colwidth", 70):
            print(per.to_string(index=False))
        print(f"\nLoci: {len(lp)} from loci_definition.tsv ({int((lp.TESTABLE == 'YES').sum())} with GWAS input + N); "
              f"loci without a GWAS credible set stay in (coloc.abf does not need one).")
        print(f"GWAS N source: {lp.GWAS_N_SOURCE.iloc[0] or '-'}  (N_EFF is never used for coloc.abf)")
        if not (lp.TESTABLE == "YES").any():
            problems.append("no testable locus (GWAS metadata / input)")
    banner("CHECK RESULT")
    print(f"Tables: {P['OUT'] / 'resource_plan.tsv'}\n        {P['OUT'] / 'locus_plan.tsv'}")
    if problems:
        print(f"NOT READY - {len(problems)} problem(s):")
        for p in problems:
            print(f"  - {p}")
        raise SystemExit(1)
    print(f"READY - {int((plan.PLAN_STATUS == 'PLANNED').sum())} datasets x {len(lp)} loci.")
    print(f"NEXT:\n  {base_command(P)}\n  sbatch {P['BASH'].name}")


def planner_mode(args, cfg) -> None:
    P = project_paths(args)
    studies = step06_studies(P["S06"])
    if args.study:
        studies = [s for s in studies if s == args.study]
    if not studies:
        raise SystemExit(f"No Step06 loci_definition.tsv under {P['S06']} - run Step06{panel_option(P['panel'])} first.")
    banner(f"STEP 11 PLANNER - external molecular QTL - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    reg, files = current_registry(cfg)
    plan = plan_datasets(P, cfg, reg)
    print_registry(plan, plan=True)
    planned = plan[plan.PLAN_STATUS == "PLANNED"]
    if planned.empty:
        raise SystemExit(f"\nNo dataset ready and planned. See {base_command(P)} --check-resources")
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
    write_tsv(planned, P["OUT"] / "registry_snapshot.tsv")
    write_tsv(files[files.DATASET_ID.isin(set(planned.DATASET_ID))], P["OUT"] / "files_snapshot.tsv")
    write_tsv(plan, P["OUT"] / "resource_plan.tsv")
    write_tsv(locus_plan(P, cfg), P["OUT"] / "locus_plan.tsv")

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
    print(f"\nStudies   : {len(studies)}\nLoci      : {len(manifest)} (all of loci_definition.tsv) -> array {array}"
          + (f" ({len(manifest) / n_tasks:.1f} loci per task)" if len(manifest) > n_tasks else " (one locus per task)"))
    print(f"Datasets  : {len(planned)} planned, frozen in registry_snapshot.tsv")
    print(f"Manifest  : {manifest_file}\n\nRUN:\n  sbatch {P['BASH'].name}\nTHEN:\n  {base_command(P)} --status\n"
          f"  {base_command(P)} --stats")


def load_snapshot(P: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    snap = read_tsv(P["OUT"] / "registry_snapshot.tsv", dtype=str, keep_default_na=False)
    files = read_tsv(P["OUT"] / "files_snapshot.tsv", dtype=str, keep_default_na=False)
    if snap.empty or files.empty:
        raise SystemExit("registry_snapshot.tsv / files_snapshot.tsv missing - run the planner again.")
    return snap, files


def worker_mode(args, cfg) -> None:
    P = project_paths(args)
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    manifest = read_tsv(Path(args.worker_manifest), dtype=str)
    selected = manifest[pd.to_numeric(manifest["TASK_ID"]) == task]
    if selected.empty:
        raise SystemExit(f"No manifest rows for array task {task}")
    rscript = require_tools()
    planned, files = load_snapshot(P)
    failed = [f"{r['STUDY_ACCESSION']} {r['LOCUS_ID']}" for _, r in selected.iterrows()
              if not run_locus(r, args, cfg, planned, files, rscript)]
    if failed:
        raise SystemExit("Colocalisation failed for: " + ", ".join(failed))


def index_mode(args, cfg) -> None:
    P = project_paths(args)
    accession = study_for_index(P, args.index)
    loci = list_loci(P["S06"] / accession)
    if not loci:
        raise SystemExit(f"Study {args.index} ({accession}) has no Step06 loci_definition.tsv{panel_option(P['panel'])}.")
    if args.locus:
        if args.locus not in loci:
            raise SystemExit(f"{args.locus} is not a locus of {accession}: {', '.join(loci)}")
        loci = [args.locus]
    rscript = require_tools()
    reg, files = current_registry(cfg)
    plan = plan_datasets(P, cfg, reg)
    planned = plan[plan.PLAN_STATUS == "PLANNED"]
    if planned.empty:
        raise SystemExit("No dataset ready and planned (see --check-resources).")
    failed = [locus for locus in loci if not run_locus(pd.Series({"STUDY_ACCESSION": accession, "LOCUS_ID": locus}),
                                                      args, cfg, planned, files, rscript)]
    if failed:
        raise SystemExit(f"{accession}: colocalisation failed for {', '.join(failed)}")


def locus_state(P: dict, accession: str, locus: str) -> tuple[str, dict]:
    d = P["OUT"] / accession / "loci" / locus
    done = read_json(d / DONE_FILE)
    if done:
        return str(done.get("STATUS", "COMPLETE")), done
    if (d / FAILED_FILE).exists():
        return "FAILED", {"ERROR": (d / FAILED_FILE).read_text(encoding="utf-8", errors="replace").splitlines()[0][:100]}
    return "PENDING", {}


def status_key(state: str) -> str:
    if state == "COMPLETE":
        return "TESTED"
    if state in {"COMPLETE_NO_TESTABLE_CANDIDATE", "NO_QTL_CANDIDATES"}:
        return "NO_TESTABLE_QTL"
    if state.startswith("NOT_TESTED"):
        return "NOT_TESTABLE"
    return state if state in {"FAILED", "PENDING"} else "PENDING"


def status_mode(args, cfg) -> None:
    P = project_paths(args)
    order = step02_order(P)
    banner(f"STEP 11 STATUS - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    snap = read_tsv(P["OUT"] / "registry_snapshot.tsv", dtype=str)
    print(f"Datasets in the submitted snapshot: {len(snap)} ({', '.join(sorted(set(snap.MODALITY)))})" if len(snap)
          else "No job script made yet (run the planner).")
    rows, failures, reasons = [], [], Counter()
    for accession in step06_studies(P["S06"]):
        counts = {"LOCI": 0, "TESTED": 0, "NO_TESTABLE_QTL": 0, "NOT_TESTABLE": 0, "FAILED": 0, "PENDING": 0,
                  "CANDIDATES": 0, "COLOC_RUN": 0, "STRONG": 0}
        for locus in list_loci(P["S06"] / accession):
            state, info = locus_state(P, accession, locus)
            key = status_key(state)
            counts["LOCI"] += 1
            counts[key] += 1
            if key in {"NOT_TESTABLE", "NO_TESTABLE_QTL"}:
                reasons[state] += 1
            if key == "FAILED":
                failures.append(f"{accession} {locus}: {info.get('ERROR', '')}")
            counts["CANDIDATES"] += int(info.get("N_CANDIDATES", 0) or 0)
            counts["COLOC_RUN"] += int(info.get("N_TESTED", 0) or 0)
            counts["STRONG"] += int(info.get("N_STRONG", 0) or 0)
        rows.append({"IDX": order.get(accession, ""), "STUDY": accession, **counts})
    if not rows:
        print(f"No Step06 studies under {P['S06']}.")
        return
    table = pd.DataFrame(rows)
    with pd.option_context("display.width", 200, "display.max_rows", None):
        print(table.to_string(index=False))
    total = table.drop(columns=["IDX", "STUDY"]).sum()
    print("\nLoci: " + "  ".join(f"{k}={int(total[k])}" for k in ("LOCI", "TESTED", "NO_TESTABLE_QTL", "NOT_TESTABLE",
                                                                       "FAILED", "PENDING")))
    for state, n in reasons.most_common():
        print(f"  {n:>4} loci  {state}")
    print("(NO_TESTABLE_QTL / NOT_TESTABLE are not evidence against colocalisation)")
    for line in failures[:20]:
        print(f"  FAILED {line}")
    if int(total["FAILED"]) or int(total["PENDING"]):
        print(f"\nRe-run failed / pending loci: sbatch {P['BASH'].name}   (finished loci are skipped)")


def stats_mode(args, cfg) -> None:
    P = project_paths(args)
    order = step02_order(P)
    banner(f"STEP 11 RESULTS - external molecular QTL - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    frames, reasons = [], []
    for acc in step06_studies(P["S06"]):
        for locus in list_loci(P["S06"] / acc):
            state, info = locus_state(P, acc, locus)
            reasons.append({"STUDY_ACCESSION": acc, "LOCUS_ID": locus, "STATE": state,
                            "N_CANDIDATES": info.get("N_CANDIDATES", 0), "N_TESTED": info.get("N_TESTED", 0),
                            "GWAS_INPUT_SOURCE": str(info.get("GWAS_INPUT_SOURCE", "")).split(" (")[0],
                            "GWAS_N_CREDIBLE_SETS": info.get("GWAS_N_CREDIBLE_SETS", ""),
                            "MULTI_SIGNAL_FLAG": info.get("MULTI_SIGNAL_FLAG", ""),
                            "REASON": info.get("REASON", info.get("ERROR", ""))})
            f = P["OUT"] / acc / "loci" / locus / RESULTS_FILE
            if nonempty(f):
                x = read_tsv(f)
                if len(x):
                    frames.append(x)
    rt = pd.DataFrame(reasons)
    if len(rt):
        print("\nLOCI (every Step06 locus; NOT_TESTED / no candidate is not evidence against colocalisation)")
        print("  " + "  ".join(f"{k}={v}" for k, v in rt.STATE.value_counts().items()))
        write_tsv(rt, P["OUT"] / "locus_reasons.tsv")
    if not frames:
        print("No candidate results yet (see --status).")
        return
    x = pd.concat(frames, ignore_index=True)
    tested = x[x["COLOC_TESTED"].astype(str) == "YES"].copy()
    tested["PP_H4"] = pd.to_numeric(tested["PP_H4"], errors="coerce")
    strong, suggestive = cfg["STRONG_PP4"], cfg["SUGGESTIVE_PP4"]
    rows = []
    for (acc, qtl), grp in x.groupby(["STUDY_ACCESSION", "QTL_TYPE"]):
        t = tested[(tested.STUDY_ACCESSION == acc) & (tested.QTL_TYPE == qtl)]
        b = t.sort_values("PP_H4", ascending=False).head(1)
        rows.append({"IDX": order.get(acc, ""), "STUDY": acc, "QTL_TYPE": qtl, "METHOD": "coloc.abf",
                     "CANDIDATES": len(grp), "TESTED": len(t), "STRONG": int((t.PP_H4 >= strong).sum()),
                     "SUGGESTIVE": int(((t.PP_H4 >= suggestive) & (t.PP_H4 < strong)).sum()),
                     "PRIOR_ROBUST": int((t.PRIOR_ROBUST_STRONG.astype(str) == "YES").sum()),
                     "STRONG_MULTI_SIGNAL_LOCI": int(((t.PP_H4 >= strong) & (t.MULTI_SIGNAL_FLAG.astype(str) == "YES")).sum()),
                     "TOP (PP.H4)": f"{b.iloc[0]['MOLECULAR_TRAIT_ID']} {b.iloc[0]['TISSUE']} ({b.iloc[0]['PP_H4']:.2f})"
                     if len(b) else ""})
    table = pd.DataFrame(rows)
    print(f"\nPER STUDY x QTL TYPE (coloc.abf = single causal variant; STRONG = PP.H4 >= {strong})")
    with pd.option_context("display.width", 250, "display.max_rows", None, "display.max_colwidth", 50):
        print(table.to_string(index=False))
        top = tested.sort_values("PP_H4", ascending=False).head(25)
        if len(top):
            print("\nTOP COLOCALISATIONS (ABF_ASSUMPTION_WARNING = several GWAS signals / complex LD: confirm with a "
                  "multi-signal method)")
            show = top.assign(IDX=top["STUDY_ACCESSION"].map(order).fillna(""), PP_H4=top["PP_H4"].round(3),
                              WARN=top["ABF_ASSUMPTION_WARNING"].astype(str).ne("").map({True: "YES", False: ""}))
            cols = ["IDX", "STUDY_ACCESSION", "LOCUS_ID", "QTL_TYPE", "DATASET_ID", "MOLECULAR_TRAIT_ID",
                    "N_COMMON_VARIANTS", "PP_H4", "PRIOR_ROBUST_STRONG", "INDEPENDENT_EVIDENCE", "WARN"]
            print(show[cols].to_string(index=False))
    print("\nCandidate status: " + "  ".join(f"{k}={v}" for k, v in x["STATUS"].astype(str).value_counts().items()))
    write_tsv(x, P["OUT"] / "coloc_all_candidates.tsv")
    write_tsv(table, P["OUT"] / "coloc_stats.tsv")
    print(f"\nTables: {P['OUT'] / 'coloc_stats.tsv'}\n        {P['OUT'] / 'coloc_all_candidates.tsv'}\n"
          f"        {P['OUT'] / 'locus_reasons.tsv'}")


def arguments():
    parser = argparse.ArgumentParser(description="Step11: external molecular-QTL colocalisation, one SLURM task per locus.")
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--scan", action="store_true",
                       help="Resource table from the Step00 manifests + validation (config/Step11_Custom_QTL_Registry.tsv).")
    modes.add_argument("--index", type=int, help="Study N of the Step02 download manifest: run its loci here.")
    modes.add_argument("--status", action="store_true", help="READ-ONLY: loci done / not testable / failed / pending.")
    modes.add_argument("--stats", action="store_true", help="Colocalisation results (tables on screen + TSV).")
    modes.add_argument("--check-resources", action="store_true",
                       help="Datasets, duplicates of Step09/10, GWAS N and loci (resource_plan.tsv, locus_plan.tsv; "
                            "exit 1 if not ready).")
    modes.add_argument("--worker-manifest", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--locus", default=None, help="With --index: only this locus (e.g. L001).")
    parser.add_argument("--study", default=None, help="Planner: only this GCST accession.")
    parser.add_argument("--reference-panel", default=DEFAULT_PANEL,
                        help="LD panel of Steps 04-06 (1000G default, 1000G30x, ...): loci from 06_finemapping/<p>/<a>_<PANEL>/.")
    parser.add_argument("--config", default="Step11_Config.py", help="Parameter file (default Step11_Config.py).")
    parser.add_argument("--qtl-types", default=None, help="Override QTL_TYPES, e.g. pQTL,sQTL.")
    parser.add_argument("--candidate-p", type=float, default=None, help="Override CANDIDATE_P.")
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
