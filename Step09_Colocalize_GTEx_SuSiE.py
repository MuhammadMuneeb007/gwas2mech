#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 09 - Formal molecular-QTL colocalisation (coloc.susie), per locus
===============================================================================

    # plan: ONE SLURM ARRAY TASK PER FINE-MAPPED LOCUS (all studies; >1000 loci share tasks)
    envs/pipeline/bin/python Step09_Colocalize_GTEx_SuSiE.py --phenotype "Parkinson's disease" --ancestry EUR
    sbatch Step09_Coloc_parkinson_s_disease_european.sh

    # study 2 (the number the Step02 planner prints) on this machine, all its loci / one locus
    envs/pipeline/bin/python Step09_Colocalize_GTEx_SuSiE.py --phenotype "Parkinson's disease" --ancestry EUR --index 2
    envs/pipeline/bin/python Step09_Colocalize_GTEx_SuSiE.py --phenotype "Parkinson's disease" --ancestry EUR --index 2 --locus L001

    # resources, progress, results
    envs/pipeline/bin/python Step09_Colocalize_GTEx_SuSiE.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources
    envs/pipeline/bin/python Step09_Colocalize_GTEx_SuSiE.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step09_Colocalize_GTEx_SuSiE.py --phenotype "Parkinson's disease" --ancestry EUR --stats

For every fine-mapped GWAS locus, test whether the GWAS signal and a molecular
trait (eQTL, sQTL, pQTL, isoQTL, exonQTL, ... ) share a causal variant:

  1. GWAS side   = the Step06 SuSiE-RSS fit of that locus (reused, not refitted)
  2. QTL side    = dense regional summary statistics (eQTL Catalogue), every
                   molecular trait with P <= --candidate-p in the locus window
  3. LD          = signed r from the Step04 ancestry-matched 1000G panel
                   (ref-based, same REF orientation as Step06 Z and the QTL beta)
  4. coloc::runsusie (QTL) + coloc::coloc.susie (GWAS x QTL), PP.H4 per signal
     pair, with prior sensitivity p12 = 1e-6 / 5e-6 / 1e-5
  5. STRONG_SHARED_SIGNAL = PP.H4 >= 0.8 and H4 > H3; PRIOR_ROBUST_STRONG when
     it also holds at the low p12

Inputs (nothing is downloaded here)
    Step06   06_finemapping/<p>/<a>/<GCST>/loci/<L>/ (matched variants, SuSiE fit, LD ids)
             <GCST>/finemapping_summary.json (PHENOTYPE_TYPE, N, N_EFF)
    Step04   resources/1000G/<ANC>/chr<N>_<ANC>_GRCh38.pgen/.pvar/.psam
    Step00   Step00_Download_Resources.py --groups eqtl_catalogue qtl_catalogue
             resources/coloc/eqtl_catalogue/tabix_ftp_paths*.tsv  (dataset list)
             resources/coloc/eqtl_catalogue/dense/<dataset_id>/*.tsv.gz + .tbi
             default contexts: GTEx v8 (eQTL, sQTL, isoQTL, exonQTL) + Sun_2018 pQTL;
             --catalog-study-regex / --catalog-all-studies for others
    optional resources/coloc/custom_qtl_registry.tsv (any other dense QTL file,
             e.g. meQTL / caQTL / hQTL; LOCAL_PATH, column names, ENABLED=True)

A molecular-QTL class with no local resource is reported as
RESOURCE_NOT_CONFIGURED - never as a biological "no colocalisation".

Outputs (09_coloc/<p>/<a>/)
    coloc_manifest.tsv, coloc_excluded_studies.tsv, coloc_stats.tsv (--stats)
    <GCST>/<GCST>_formal_coloc_summary.tsv   every tested molecular trait
           <GCST>_locus_qtl_matrix.tsv       locus x QTL class evidence
           <GCST>_locus_evidence_table.tsv   one row per locus
           <GCST>_molecular_qtl_candidates.tsv, resource_status.tsv, coloc_summary.json
           loci/<L>/candidates/<id>/          coloc summary / variant results / priors / LD diagnostics

SLURM settings come from the #SBATCH header of HPC_CPU.sh; --partition --qos
--account --cpus --mem --time --max-parallel override them.
===============================================================================
"""

from __future__ import annotations

import argparse
import gzip
import io
from collections import OrderedDict
from contextlib import contextmanager
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

from typing import Any, Iterable

import numpy as np
import pandas as pd


STEP09_VERSION = "5.0.0-per-locus"
PROJECT_ROOT = Path.cwd().resolve()
MAX_ARRAY_TASKS = 1000

DEFAULT_CPUS = 4
DEFAULT_REMOTE_DELAY = 0.0
DEFAULT_LD_MAF = 0.01
DEFAULT_LD_MISMATCH_S = 0.10

DEFAULT_CANDIDATE_P = 1e-5
DEFAULT_MIN_COMMON_SNPS = 100
DEFAULT_MAX_COMMON_SNPS = 12000
DEFAULT_STRONG_PP4 = 0.80
DEFAULT_SUGGESTIVE_PP4 = 0.50
DEFAULT_MAX_CANDIDATES_PER_CONTEXT_LOCUS = 100

COLOC_PRIOR_P1 = 1e-4
COLOC_PRIOR_P2 = 1e-4
COLOC_PRIOR_P12 = 5e-6
COLOC_PRIOR_P12_LOW = 1e-6
COLOC_PRIOR_P12_HIGH = 1e-5

SUPPORTED_QTL_TYPES = [
    "eQTL",
    "sQTL",
    "pQTL",
    "meQTL",
    "caQTL",
    "hQTL",
    "isoQTL",
    "exonQTL",
    "apaQTL",
    "metaboliteQTL",
    "otherQTL",
]

ANCESTRY_ALIASES = {
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
}

EQTL_CATALOG_PUBLIC = {
    "TABIX_PATHS": {
        "url": (
            "https://raw.githubusercontent.com/eQTL-Catalogue/"
            "eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv"
        ),
        "path": "resources/coloc/eqtl_catalogue/tabix_ftp_paths.tsv",
    },
    "TABIX_IMPORTED": {
        "url": (
            "https://raw.githubusercontent.com/eQTL-Catalogue/"
            "eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv"
        ),
        "path": "resources/coloc/eqtl_catalogue/tabix_ftp_paths_imported.tsv",
    },
    "DATASET_METADATA_R8": {
        "url": (
            "https://raw.githubusercontent.com/eQTL-Catalogue/"
            "eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r8_beta.tsv"
        ),
        "path": "resources/coloc/eqtl_catalogue/dataset_metadata_r8_beta.tsv",
    },
    "POPULATION_ASSIGNMENTS": {
        "url": (
            "https://raw.githubusercontent.com/eQTL-Catalogue/"
            "eQTL-Catalogue-resources/master/data_tables/population_assignments.tsv"
        ),
        "path": "resources/coloc/eqtl_catalogue/population_assignments.tsv",
    },
}

CUSTOM_REGISTRY_COLUMNS = [
    "ENABLED",
    "DATASET",
    "QTL_TYPE",
    "TISSUE",
    "CONDITION",
    "SOURCE_URL",
    "LOCAL_PATH",
    "ACCESS_MODE",
    "SAMPLE_SIZE",
    "TRAIT_COLUMN",
    "GENE_COLUMN",
    "VARIANT_COLUMN",
    "CHR_COLUMN",
    "POS_COLUMN",
    "REF_COLUMN",
    "ALT_COLUMN",
    "BETA_COLUMN",
    "SE_COLUMN",
    "P_COLUMN",
    "MAF_COLUMN",
    "EFFECT_ALLELE_MODE",
    "EFFECT_ALLELE_COLUMN",
    "GENOME_BUILD",
    "NOTES",
]


def banner(text: str) -> None:
    print()
    print("=" * 124)
    print(text)
    print("=" * 124)


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
        raise ValueError(
            f"Unsupported ancestry {value!r}; use EUR, AFR, EAS, SAS, AMR."
        )
    return ANCESTRY_ALIASES[key]


def clean_chr(value: Any) -> int | None:
    if value is None:
        return None
    text = re.sub(r"^chr", "", str(value).strip(), flags=re.I)
    try:
        chrom = int(float(text))
    except Exception:
        return None
    return chrom if 1 <= chrom <= 22 else None


def clean_allele(value: Any) -> str | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    text = str(value).strip().upper()
    if text in {"", ".", "NA", "NAN", "NONE", "<NA>"}:
        return None
    return text


def make_variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    chrom = clean_chr(chrom)
    ref = clean_allele(ref)
    alt = clean_allele(alt)
    try:
        pos = int(float(pos))
    except Exception:
        return ""
    if chrom is None or pos <= 0 or ref is None or alt is None:
        return ""
    return f"{chrom}:{pos}:{ref}:{alt}"


def parse_variant_id(value: Any) -> tuple[int | None, int | None, str | None, str | None]:
    text = str(value).strip()
    match = re.match(
        r"^(?:chr)?(\d+)[:_](\d+)[:_]([^:_]+)[:_]([^:_]+)(?:[:_]b\d+)?$",
        text,
        flags=re.I,
    )
    if not match:
        return None, None, None, None
    chrom, pos, ref, alt = match.groups()
    return clean_chr(chrom), int(pos), clean_allele(ref), clean_allele(alt)


def safe_int(value: Any) -> int | None:
    try:
        if value is None or str(value).strip() == "":
            return None
        return int(float(value))
    except Exception:
        return None


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def boolish(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except Exception:
        pass
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def first_existing(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    lookup = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if str(candidate).lower() in lookup:
            return lookup[str(candidate).lower()]
    return None


def safe_json(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    lower = path.name.lower()
    if lower.endswith(".parquet") or lower.endswith(".pq"):
        return pd.read_parquet(path, **kwargs)
    if lower.endswith(".csv") or lower.endswith(".csv.gz"):
        return pd.read_csv(path, low_memory=False, **kwargs)
    return pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        low_memory=False,
        **kwargs,
    )


def validate_walltime(value: str) -> str:
    text = str(value).strip()
    if not re.fullmatch(r"(\d+-)?\d{1,3}(:\d{2}){0,2}", text):
        raise argparse.ArgumentTypeError(f"walltime {value!r} must look like 12:00:00 or 2-00:00:00")
    return text


def run_command(
    command: list[str | Path],
    *,
    log_file: Path | None = None,
    capture: bool = False,
    check: bool = True,
) -> str:
    command = [str(x) for x in command]
    print("$ " + " ".join(shlex.quote(x) for x in command), flush=True)
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    output = []
    handle = None
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handle = log_file.open("w", encoding="utf-8")
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            output.append(line)
            if not capture:
                print(line, end="", flush=True)
            if handle:
                handle.write(line)
                handle.flush()
    finally:
        if handle:
            handle.close()
    rc = proc.wait()
    text = "".join(output)
    if check and rc != 0:
        raise RuntimeError(
            f"Command failed with exit code {rc}:\n"
            + " ".join(shlex.quote(x) for x in command)
            + "\n"
            + text[-8000:]
        )
    return text


@contextmanager
def exclusive_file_lock(path: Path):
    """Cross-process advisory lock for shared downloads/aggregate writes on Linux HPC."""
    import fcntl   # Linux / HPC only

    lock_path = Path(str(path) + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("w") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


_REGION_CACHE: OrderedDict[tuple, pd.DataFrame] = OrderedDict()
_REGION_CACHE_LIMIT = 32
_CUSTOM_FILE_CACHE: OrderedDict[str, pd.DataFrame] = OrderedDict()
_CUSTOM_FILE_CACHE_LIMIT = 2
_LBF_TABLE_CACHE: OrderedDict[str, pd.DataFrame] = OrderedDict()
_LBF_TABLE_CACHE_LIMIT = 2


def _cache_put(cache: OrderedDict, key, value, limit: int) -> None:
    cache[key] = value
    cache.move_to_end(key)
    while len(cache) > limit:
        cache.popitem(last=False)


def cached_custom_file(path: Path) -> pd.DataFrame:
    key = str(path.resolve())
    if key in _CUSTOM_FILE_CACHE:
        _CUSTOM_FILE_CACHE.move_to_end(key)
        return _CUSTOM_FILE_CACHE[key]
    x = read_table(path)
    _cache_put(_CUSTOM_FILE_CACHE, key, x, _CUSTOM_FILE_CACHE_LIMIT)
    return x


def read_resource_env(root: Path) -> dict[str, str]:
    path = root / "resource_paths.env"
    result = {}
    if not path.exists():
        return result
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.replace("export ", "").strip()] = value.strip().strip('"').strip("'")
    return result


def find_executable(
    root: Path,
    explicit: str | None,
    env_name: str,
    command: str,
) -> Path:
    env = read_resource_env(root)
    candidates = [
        explicit,
        os.environ.get(env_name),
        env.get(env_name),
        shutil.which(command),
        root / "envs" / "pipeline" / "bin" / command,
    ]
    for item in candidates:
        if not item:
            continue
        path = Path(item).expanduser()
        if path.exists() and os.access(path, os.X_OK):
            return path.resolve()
    raise FileNotFoundError(
        f"Could not locate executable {command!r}; provide it explicitly."
    )


def write_custom_registry_template(root: Path, overwrite: bool = False) -> Path:
    path = root / "resources" / "coloc" / "custom_qtl_registry.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        return path

    base = {
        "ENABLED": False,
        "TISSUE": "",
        "CONDITION": "",
        "SOURCE_URL": "",
        "LOCAL_PATH": "",
        "ACCESS_MODE": "file",
        "SAMPLE_SIZE": "",
        "GENE_COLUMN": "",
        "VARIANT_COLUMN": "variant",
        "CHR_COLUMN": "",
        "POS_COLUMN": "",
        "REF_COLUMN": "",
        "ALT_COLUMN": "",
        "BETA_COLUMN": "beta",
        "SE_COLUMN": "se",
        "P_COLUMN": "pvalue",
        "MAF_COLUMN": "maf",
        "EFFECT_ALLELE_MODE": "ALT",
        "EFFECT_ALLELE_COLUMN": "",
        "GENOME_BUILD": "GRCh38",
        "NOTES": "",
    }

    specs = [
        ("YOUR_MEQTL_DATASET", "meQTL", "cpg_id", "DNA methylation QTL"),
        ("YOUR_CAQTL_DATASET", "caQTL", "peak_id", "chromatin accessibility QTL"),
        ("YOUR_HQTL_DATASET", "hQTL", "peak_id", "histone-mark QTL"),
        ("YOUR_METABOLITE_QTL_DATASET", "metaboliteQTL", "metabolite_id", "metabolite QTL"),
        ("YOUR_OTHER_QTL_DATASET", "otherQTL", "trait_id", "future/custom molecular QTL"),
    ]

    rows = []
    for dataset, qtl_type, trait_col, note in specs:
        row = dict(base)
        row.update({
            "DATASET": dataset,
            "QTL_TYPE": qtl_type,
            "TRAIT_COLUMN": trait_col,
            "NOTES": (
                f"Configure a licensed/public dense {note} summary-statistics source: "
                "set LOCAL_PATH (bgzip + tabix), the column names and ENABLED=True. Nothing is downloaded."
            ),
        })
        rows.append(row)

    pd.DataFrame(rows, columns=CUSTOM_REGISTRY_COLUMNS).to_csv(
        path, sep="\t", index=False
    )
    return path


def arguments():
    parser = argparse.ArgumentParser(
        description="Formal molecular-QTL colocalisation (coloc.susie), one SLURM task per fine-mapped locus."
    )
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)

    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--index", type=int,
                       help="Study N of the Step02 download manifest: run its loci on this machine.")
    modes.add_argument("--status", action="store_true", help="READ-ONLY: loci done / failed / pending.")
    modes.add_argument("--stats", action="store_true", help="Colocalisation results per study (coloc_stats.tsv).")
    modes.add_argument("--check-resources", action="store_true",
                       help="Which eQTL Catalogue / custom QTL datasets are available locally.")
    modes.add_argument("--worker-manifest", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--locus", default=None, help="With --index: only this locus (e.g. L001).")
    parser.add_argument("--study", default=None, help="Planner: only this GCST accession.")

    parser.add_argument("--registry", default="resources/coloc/custom_qtl_registry.tsv")
    parser.add_argument("--qtl-types", default=",".join(SUPPORTED_QTL_TYPES))
    parser.add_argument("--catalog-all-studies", action="store_true",
                        help="Every downloaded eQTL Catalogue study (default: GTEx v8 + Sun_2018).")
    parser.add_argument("--catalog-study-regex", default="")

    parser.add_argument("--candidate-p", type=float, default=DEFAULT_CANDIDATE_P)
    parser.add_argument("--min-common-snps", type=int, default=DEFAULT_MIN_COMMON_SNPS)
    parser.add_argument("--max-common-snps", type=int, default=DEFAULT_MAX_COMMON_SNPS)
    parser.add_argument("--strong-pp4", type=float, default=DEFAULT_STRONG_PP4)
    parser.add_argument("--suggestive-pp4", type=float, default=DEFAULT_SUGGESTIVE_PP4)
    parser.add_argument("--max-candidates-per-context-locus", type=int,
                        default=DEFAULT_MAX_CANDIDATES_PER_CONTEXT_LOCUS)
    parser.add_argument("--ld-maf", type=float, default=DEFAULT_LD_MAF)
    parser.add_argument("--ld-mismatch-s", type=float, default=DEFAULT_LD_MISMATCH_S)
    parser.add_argument("--plink2")
    parser.add_argument("--tabix")
    parser.add_argument("--rscript")
    parser.add_argument("--force", action="store_true", help="Re-run loci / candidates that are complete.")

    slurm = parser.add_argument_group("SLURM (default: the #SBATCH header of HPC_CPU.sh)")
    slurm.add_argument("--template", default="HPC_CPU.sh")
    slurm.add_argument("--partition", default=None)
    slurm.add_argument("--qos", default=None)
    slurm.add_argument("--account", default=None)
    slurm.add_argument("--cpus", type=int, default=None)
    slurm.add_argument("--mem", "--memory", dest="mem", default=None)
    slurm.add_argument("--time", type=validate_walltime, default=None)
    slurm.add_argument("--max-parallel", type=int, default=0,
                       help="Maximum simultaneously running locus tasks (0 = no throttle).")
    args = parser.parse_args()

    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": str(DEFAULT_CPUS), "mem": "32G", "time": "12:00:00"}
    template = PROJECT_ROOT / args.template
    if template.is_file():
        for line in template.read_text(encoding="utf-8").splitlines():
            match = re.match(r"^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
            if match and (match.group(2) or "").strip() and match.group(1) not in {"job-name", "output", "error", "array"}:
                values[match.group(1)] = match.group(2).strip()
    for key, value in (("partition", args.partition), ("qos", args.qos), ("account", args.account),
                       ("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time)):
        if value not in (None, ""):
            values[key] = str(value)
    args.slurm_values = values
    args.cpus = int(values["cpus-per-task"])
    if args.max_common_snps < args.min_common_snps:
        parser.error("--max-common-snps must be >= --min-common-snps")
    if args.locus and args.index is None:
        parser.error("--locus needs --index")
    return args


def project_paths(root: Path, phenotype: str, ancestry_label: str) -> dict[str, Path]:
    p = slugify(phenotype)
    a = slugify(ancestry_label)
    return {
        "S02": root / "01_gwas_catalog" / p / a,
        "S06": root / "06_finemapping" / p / a,
        "OUT": root / "09_coloc" / p / a,
        "LOG": root / "logs" / "step09_coloc" / p / a,
    }


def discover_step06_studies(s06: Path) -> tuple[list[str], list[dict]]:
    studies = []
    excluded = []
    if not s06.exists():
        return studies, excluded

    accessions = set()
    manifest = s06 / "finemapping_manifest.tsv"
    if manifest.exists():
        x = pd.read_csv(manifest, sep="\t", dtype=str, low_memory=False)
        if "STUDY_ACCESSION" in x.columns:
            accessions.update(
                v.strip()
                for v in x["STUDY_ACCESSION"].dropna().astype(str)
                if v.strip()
            )

    for directory in s06.glob("GCST*"):
        if directory.is_dir():
            accessions.add(directory.name)

    for accession in sorted(accessions):
        study_dir = s06 / accession
        summary = safe_json(study_dir / "finemapping_summary.json")
        fine_file = study_dir / f"{accession}_finemapped_variants.tsv.gz"
        status = str(summary.get("STATUS", "")).upper()
        n_success = safe_int(summary.get("N_LOCI_SUCCESS")) or 0

        reason = ""
        if not summary:
            reason = "Step06 summary missing/unreadable"
        elif status == "NO_SIGNIFICANT_LOCI":
            reason = "No significant loci"
        elif n_success <= 0 and status not in {"COMPLETE", "PARTIAL"}:
            reason = f"No successful loci; status={status}"
        elif not fine_file.exists() or fine_file.stat().st_size == 0:
            reason = "Combined Step06 fine-mapped file missing/empty"

        if reason:
            excluded.append({
                "STUDY_ACCESSION": accession,
                "STEP06_STATUS": status,
                "REASON": reason,
            })
        else:
            studies.append(accession)

    return studies, excluded


def list_loci(study_dir: Path) -> list[str]:
    """Step06 loci that have a SuSiE fit (failed loci cannot be colocalised)."""
    loci_dir = study_dir / "loci"
    if not loci_dir.exists():
        return []
    return sorted(
        d.name for d in loci_dir.iterdir()
        if d.is_dir() and re.fullmatch(r"L\d+", d.name) and (d / "susie" / "susie_fit.rds").exists()
    )


def load_gwas_metadata(s06: Path, accession: str) -> dict:
    """Phenotype type and N from Step06 (finemapping_summary.json, which carries
    Step02's PHENOTYPE_TYPE / N_CASES / N_CONTROLS / N_TOTAL / N_EFF), else
    from Step02 GWAS_study_metadata.tsv. The GWAS side of coloc.susie is the
    Step06 SuSiE fit, so N only documents the run."""
    meta = safe_json(s06 / accession / "finemapping_summary.json")
    step02 = s06.parents[2] / "01_gwas_catalog" / s06.parent.name / s06.name / "GWAS_study_metadata.tsv"
    if step02.exists() and not meta.get("PHENOTYPE_TYPE"):
        table = read_table(step02, dtype=str)
        row = table[table.get("GCST_ID", pd.Series(dtype=str)).astype(str) == accession]
        if not row.empty:
            meta = {**row.iloc[0].to_dict(), **meta}

    n_cases, n_controls = safe_int(meta.get("N_CASES")), safe_int(meta.get("N_CONTROLS"))
    phenotype_type = str(meta.get("PHENOTYPE_TYPE", "") or "")
    binary = phenotype_type == "Binary"
    return {
        "METADATA_FILE": str(s06 / accession / "finemapping_summary.json"),
        "PHENOTYPE_TYPE": phenotype_type or "Unknown",
        "COLOC_GWAS_N": safe_int(meta.get("GWAS_N")) or safe_int(meta.get("N_TOTAL")),
        "GWAS_N_SOURCE": meta.get("GWAS_N_SOURCE", ""),
        "N_CASES": n_cases if binary else None,
        "N_CONTROLS": n_controls if binary else None,
        "N_EFF": safe_int(meta.get("N_EFF")),
        "GWAS_TYPE": "cc" if binary else "quant",
        "CASE_PROPORTION": n_cases / (n_cases + n_controls) if binary and n_cases and n_controls else None,
        "METADATA_STATUS": "OK" if phenotype_type in ("Binary", "Continuous") else "PHENOTYPE_TYPE_UNKNOWN",
    }


def load_gwas_locus(study_dir: Path, locus_id: str) -> pd.DataFrame:
    locus_dir = study_dir / "loci" / locus_id
    paths = [
        locus_dir / "reference_qc_matched_variants.tsv.gz",
        locus_dir / "matched_variants.tsv",
    ]
    source = next(
        (p for p in paths if p.exists() and p.stat().st_size > 0),
        None,
    )
    if source is None:
        raise FileNotFoundError(
            f"No Step06 dense matched GWAS table for {locus_id} under {locus_dir}"
        )

    x = read_table(source)

    chr_col = first_existing(x.columns, ["REFERENCE_CHR", "CHR"])
    pos_col = first_existing(x.columns, ["REFERENCE_POS", "POS"])
    ref_col = first_existing(x.columns, ["REFERENCE_REF"])
    alt_col = first_existing(x.columns, ["REFERENCE_ALT"])
    rid_col = first_existing(x.columns, ["REFERENCE_ID"])
    beta_col = first_existing(x.columns, ["BETA_REF"])
    se_col = first_existing(x.columns, ["SE"])
    p_col = first_existing(x.columns, ["P"])

    required = {
        "CHR": chr_col, "POS": pos_col, "REF": ref_col, "ALT": alt_col,
        "REFERENCE_ID": rid_col, "BETA_REF": beta_col, "SE": se_col, "P": p_col,
    }
    missing = [name for name, col in required.items() if col is None]
    if missing:
        raise RuntimeError(
            f"{source} is missing fields needed for formal coloc: {missing}"
        )

    out = pd.DataFrame({
        "CHR": x[chr_col].map(clean_chr),
        "POS": numeric(x[pos_col]),
        "REF": x[ref_col].map(clean_allele),
        "ALT": x[alt_col].map(clean_allele),
        "REFERENCE_ID": x[rid_col].astype(str),
        "GWAS_BETA": numeric(x[beta_col]),
        "GWAS_SE": numeric(x[se_col]),
        "GWAS_P": numeric(x[p_col]),
    })

    out["VARIANT_KEY"] = [
        make_variant_key(c, p, r, a)
        for c, p, r, a in zip(out["CHR"], out["POS"], out["REF"], out["ALT"])
    ]

    out = out[
        out["VARIANT_KEY"].ne("")
        & out["GWAS_BETA"].notna()
        & out["GWAS_SE"].notna()
        & (out["GWAS_SE"] > 0)
        & out["REFERENCE_ID"].notna()
    ].copy()

    out = (
        out.sort_values(["CHR", "POS", "GWAS_P"], kind="stable")
        .drop_duplicates("VARIANT_KEY", keep="first")
        .reset_index(drop=True)
    )
    return out


def locus_region(gwas: pd.DataFrame) -> tuple[int, int, int]:
    chroms = sorted(set(gwas["CHR"].dropna().astype(int)))
    if len(chroms) != 1:
        raise RuntimeError(f"Locus spans multiple chromosomes: {chroms}")
    return chroms[0], int(gwas["POS"].min()), int(gwas["POS"].max())


def ensure_catalogue_metadata(root: Path, auto_download: bool = False) -> None:
    """The dataset lists come from Step00_Download_Resources.py --groups eqtl_catalogue."""
    missing = [spec["path"] for spec in EQTL_CATALOG_PUBLIC.values()
               if spec["path"].endswith("tabix_ftp_paths.tsv") and not (root / spec["path"]).exists()]
    if missing:
        print("[INFO] eQTL Catalogue dataset list missing: " + ", ".join(missing)
              + "\n       envs/pipeline/bin/python Step00_Download_Resources.py --groups eqtl_catalogue")


def normalize_catalogue_manifest(path: Path) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    rename = {}
    mapping = {
        "study_id": ["study_id"],
        "dataset_id": ["dataset_id"],
        "study_label": ["study_label", "study"],
        "sample_group": ["sample_group", "qtl_group"],
        "tissue_id": ["tissue_id"],
        "tissue_label": ["tissue_label", "tissue"],
        "condition_label": ["condition_label", "condition"],
        "sample_size": ["sample_size", "n"],
        "quant_method": ["quant_method", "quantification_method"],
        "ftp_path": ["ftp_path", "sumstats_path"],
        "ftp_cs_path": ["ftp_cs_path"],
        "ftp_lbf_path": ["ftp_lbf_path"],
    }
    for target, aliases in mapping.items():
        col = first_existing(x.columns, aliases)
        if col:
            rename[col] = target
    x = x.rename(columns=rename)

    for column in mapping:
        if column not in x.columns:
            x[column] = ""
    return x


def load_eqtl_catalogue(root: Path, auto_download: bool) -> pd.DataFrame:
    ensure_catalogue_metadata(root, auto_download)
    paths = [
        root / EQTL_CATALOG_PUBLIC["TABIX_PATHS"]["path"],
        root / EQTL_CATALOG_PUBLIC["TABIX_IMPORTED"]["path"],
    ]
    frames = [
        normalize_catalogue_manifest(path)
        for path in paths
        if path.exists() and path.stat().st_size > 0
    ]
    if not frames:
        return pd.DataFrame()
    return (
        pd.concat(frames, ignore_index=True)
        .drop_duplicates(["dataset_id", "ftp_path"], keep="last")
        .reset_index(drop=True)
    )


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


def select_catalogue_contexts(
    catalogue: pd.DataFrame,
    args,
    allowed_types: set[str],
) -> pd.DataFrame:
    if catalogue.empty:
        return catalogue

    x = catalogue.copy()
    x["QTL_TYPE"] = x["quant_method"].map(quant_method_to_qtl_type)

    unmapped = (
        x.loc[x["QTL_TYPE"].isna(), "quant_method"]
        .fillna("")
        .astype(str)
        .value_counts()
    )
    if not unmapped.empty:
        print("[INFO] eQTL Catalogue quant_method values not mapped to a QTL class:")
        for method, count in unmapped.items():
            print(f"       {method or '<EMPTY>'}: {count}")

    x = x[x["QTL_TYPE"].notna()].copy()
    x = x[
        x["QTL_TYPE"].astype(str).str.lower().isin(allowed_types)
    ]

    if args.catalog_study_regex:
        pattern = re.compile(args.catalog_study_regex, flags=re.I)
        x = x[
            x["study_label"].astype(str).map(
                lambda value: bool(pattern.search(value))
            )
        ]
    elif not args.catalog_all_studies:
        labels = x["study_label"].astype(str)
        x = x[
            labels.str.contains("GTEx", case=False, regex=False)
            | labels.str.fullmatch("Sun_2018", case=False)
        ].copy()

    x = x[x["ftp_path"].astype(str).str.strip().ne("")].copy()

    if not x.empty:
        print("[INFO] Formal-coloc source counts by study/QTL type:")
        counts = (
            x.groupby(["study_label", "QTL_TYPE"], dropna=False)
            .size()
            .reset_index(name="N_CONTEXTS")
        )
        print(counts.to_string(index=False))
        if x["study_label"].astype(str).str.contains("GTEx_v8", case=False, regex=False).any():
            print("[INFO] GTEx data in the eQTL Catalogue are GTEx v8 (labelled GTEx_v8).")

    return x.reset_index(drop=True)


_HEADER_CACHE: dict[str, list[str]] = {}


def _is_remote_source(value: str) -> bool:
    return bool(re.match(r"^(?:https?|ftp)://", str(value).strip(), flags=re.I))


def _catalogue_basename(value: str) -> str:
    text = str(value).strip().split("?", 1)[0].rstrip("/")
    return text.rsplit("/", 1)[-1]


def resolve_catalogue_local_sources(root: Path, contexts: pd.DataFrame) -> pd.DataFrame:
    """Resolve every selected eQTL-Catalogue context to local BGZF + TBI files.

    Normal Step09 analysis is deliberately network-free.  The manifest FTP URL
    is retained only as provenance; ACCESS_PATH is always a local filesystem path.
    """
    if contexts.empty:
        return contexts.copy()

    dense_root = root / "resources" / "coloc" / "eqtl_catalogue" / "dense"
    if not dense_root.exists():
        print(f"[LOCAL QTL] no dense eQTL Catalogue files under {dense_root}")
        return contexts.iloc[0:0].assign(ACCESS_PATH="", ACCESS_INDEX="", ACCESS_MODE="")

    by_name: dict[str, list[Path]] = {}
    for path in dense_root.rglob("*.tsv.gz"):
        if path.is_file() and path.stat().st_size > 0:
            by_name.setdefault(path.name, []).append(path.resolve())

    rows = []
    for _, row in contexts.iterrows():
        out = row.to_dict()
        dataset_id = str(row.get("dataset_id", "")).strip()
        remote = str(row.get("ftp_path", "")).strip()
        basename = _catalogue_basename(remote)

        candidates: list[Path] = []
        if dataset_id and basename:
            candidates.append(dense_root / dataset_id / basename)
        if basename:
            candidates.append(dense_root / basename)
            candidates.extend(by_name.get(basename, []))

        # Prefer a path whose parent is the dataset_id when duplicates exist.
        candidates = [Path(x) for x in dict.fromkeys(map(str, candidates))]
        candidates.sort(
            key=lambda x: (
                0 if dataset_id and x.parent.name == dataset_id else 1,
                len(str(x)),
            )
        )

        data = None
        index = None
        for candidate in candidates:
            if not candidate.exists() or candidate.stat().st_size == 0:
                continue
            tbi = Path(str(candidate) + ".tbi")
            if tbi.exists() and tbi.stat().st_size > 0:
                data = candidate.resolve()
                index = tbi.resolve()
                break

        out["REMOTE_SOURCE_URL"] = remote
        out["ACCESS_PATH"] = str(data) if data else ""
        out["ACCESS_INDEX"] = str(index) if index else ""
        out["ACCESS_MODE"] = "LOCAL_TABIX" if data and index else "MISSING_LOCAL"
        rows.append(out)

    result = pd.DataFrame(rows)
    missing = result[result["ACCESS_MODE"].ne("LOCAL_TABIX")]
    print(
        f"[LOCAL QTL] {len(result) - len(missing)}/{len(result)} selected eQTL Catalogue datasets "
        "are downloaded (.tsv.gz + .tbi)"
    )
    if not missing.empty:
        print(f"[LOCAL QTL] {len(missing)} not downloaded - not tested (RESOURCE_NOT_CONFIGURED); get them with\n"
              "            envs/pipeline/bin/python Step00_Download_Resources.py --groups qtl_catalogue")
    return result[result["ACCESS_MODE"].eq("LOCAL_TABIX")].reset_index(drop=True)


def tabix_header(tabix: Path, source: str) -> list[str]:
    source = str(source).strip()
    if _is_remote_source(source):
        raise RuntimeError(
            f"Remote tabix is disabled in this Step09 build: {source}"
        )
    path = Path(source)
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Local tabix source missing/empty: {path}")

    key = str(path.resolve())
    if key in _HEADER_CACHE:
        return _HEADER_CACHE[key]

    text = run_command([tabix, "-H", path], capture=True, check=False)
    header_lines = [line for line in text.splitlines() if line.startswith("#")]

    # Some BGZF files have a plain first header line rather than a tabix '#'-header.
    if not header_lines:
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                first = handle.readline().rstrip("\n\r")
            if first:
                header_lines = [first]
        except Exception:
            pass

    if not header_lines:
        raise RuntimeError(f"Could not retrieve header for local file: {path}")

    header = header_lines[-1].lstrip("#").split("\t")
    _HEADER_CACHE[key] = header
    return header


def check_local_tabix_connectivity(tabix: Path, contexts: pd.DataFrame) -> None:
    if contexts.empty:
        return
    source = str(contexts.iloc[0]["ACCESS_PATH"])
    header = tabix_header(tabix, source)
    if not header:
        raise RuntimeError(f"Local eQTL-Catalogue header is empty: {source}")
    print(f"[OK] Local tabix/header: {source}")


def tabix_region(
    tabix: Path,
    source: str,
    chrom: int,
    start: int,
    end: int,
) -> pd.DataFrame:
    source = str(source).strip()
    if _is_remote_source(source):
        raise RuntimeError(f"Remote tabix is disabled: {source}")
    path = Path(source)
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Local tabix source missing/empty: {path}")
    tbi = Path(str(path) + ".tbi")
    if not tbi.exists() or tbi.stat().st_size == 0:
        raise FileNotFoundError(f"Local tabix index missing/empty: {tbi}")

    key = (str(path.resolve()), int(chrom), int(start), int(end))
    if key in _REGION_CACHE:
        _REGION_CACHE.move_to_end(key)
        return _REGION_CACHE[key].copy()

    header = tabix_header(tabix, str(path))
    errors = []
    successful_query = False

    for region in (f"{chrom}:{start}-{end}", f"chr{chrom}:{start}-{end}"):
        proc = subprocess.run(
            [str(tabix), str(path), region],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=120,
        )
        if proc.returncode != 0:
            errors.append(proc.stderr.strip() or f"returncode={proc.returncode}")
            continue

        successful_query = True
        lines = [
            line for line in proc.stdout.splitlines()
            if line and not line.startswith("#")
        ]
        if lines:
            result = pd.read_csv(
                io.StringIO("\n".join(lines)),
                sep="\t",
                names=header,
                dtype=str,
                low_memory=False,
            )
            _cache_put(_REGION_CACHE, key, result, _REGION_CACHE_LIMIT)
            return result.copy()

    if not successful_query and errors:
        raise RuntimeError(
            f"tabix failed for local file {path}: {errors[-1][:500]}"
        )

    result = pd.DataFrame(columns=header)
    _cache_put(_REGION_CACHE, key, result, _REGION_CACHE_LIMIT)
    return result.copy()

def standardize_eqtl_catalogue_region(
    x: pd.DataFrame,
    context: pd.Series,
) -> pd.DataFrame:
    if x.empty:
        return pd.DataFrame()

    trait_col = first_existing(x.columns, ["molecular_trait_id", "phenotype_id", "trait_id"])
    gene_col = first_existing(x.columns, ["gene_id"])
    chr_col = first_existing(x.columns, ["chromosome", "chr", "chrom"])
    pos_col = first_existing(x.columns, ["position", "pos", "bp"])
    ref_col = first_existing(x.columns, ["ref", "reference_allele"])
    alt_col = first_existing(x.columns, ["alt", "alternate_allele"])
    variant_col = first_existing(x.columns, ["variant", "variant_id"])
    p_col = first_existing(x.columns, ["pvalue", "p_value", "pval", "pval_nominal"])
    beta_col = first_existing(x.columns, ["beta", "slope", "effect_size"])
    se_col = first_existing(x.columns, ["se", "slope_se", "standard_error"])
    maf_col = first_existing(x.columns, ["maf", "af", "allele_frequency"])
    an_col = first_existing(x.columns, ["an", "allele_number"])

    if trait_col is None or beta_col is None or se_col is None:
        return pd.DataFrame()

    out = pd.DataFrame(index=x.index)

    if all(c is not None for c in [chr_col, pos_col, ref_col, alt_col]):
        out["CHR"] = x[chr_col].map(clean_chr)
        out["POS"] = numeric(x[pos_col])
        out["REF"] = x[ref_col].map(clean_allele)
        out["ALT"] = x[alt_col].map(clean_allele)
    elif variant_col:
        parsed = x[variant_col].map(parse_variant_id)
        out["CHR"] = parsed.map(lambda z: z[0])
        out["POS"] = parsed.map(lambda z: z[1])
        out["REF"] = parsed.map(lambda z: z[2])
        out["ALT"] = parsed.map(lambda z: z[3])
    else:
        return pd.DataFrame()

    out["VARIANT_KEY"] = [
        make_variant_key(c, p, r, a)
        for c, p, r, a in zip(out["CHR"], out["POS"], out["REF"], out["ALT"])
    ]
    out["MOLECULAR_TRAIT_ID"] = x[trait_col].fillna("").astype(str)
    out["GENE_ID"] = x[gene_col].fillna("").astype(str) if gene_col else ""
    out["P"] = numeric(x[p_col]) if p_col else np.nan
    out["BETA_ALT"] = numeric(x[beta_col])
    out["SE"] = numeric(x[se_col])
    out["MAF"] = numeric(x[maf_col]) if maf_col else np.nan
    out["AN"] = numeric(x[an_col]) if an_col else np.nan

    # eQTL Catalogue states ALT is the effect allele.
    out["BETA_REF"] = -out["BETA_ALT"]

    out["DATASET_ID"] = str(context.get("dataset_id", ""))
    out["STUDY_LABEL"] = str(context.get("study_label", ""))
    out["TISSUE"] = str(context.get("tissue_label", ""))
    out["CONDITION"] = str(context.get("condition_label", ""))
    out["QTL_TYPE"] = str(context.get("QTL_TYPE", ""))
    metadata_n = safe_int(context.get("sample_size"))
    if out["AN"].notna().any():
        out["QTL_N"] = np.floor(out["AN"] / 2.0)
        out["QTL_N"] = out["QTL_N"].where(out["QTL_N"] > 0, metadata_n)
    else:
        out["QTL_N"] = metadata_n
    out["SOURCE_URL"] = str(context.get("ftp_path", ""))

    out = out[
        out["VARIANT_KEY"].ne("")
        & out["MOLECULAR_TRAIT_ID"].str.strip().ne("")
        & out["BETA_REF"].notna()
        & out["SE"].notna()
        & (out["SE"] > 0)
    ].copy()

    return (
        out.drop_duplicates(
            ["MOLECULAR_TRAIT_ID", "VARIANT_KEY"],
            keep="first",
        )
        .reset_index(drop=True)
    )


def discover_catalogue_candidates_for_locus(
    tabix: Path,
    contexts: pd.DataFrame,
    locus_id: str,
    gwas: pd.DataFrame,
    args,
) -> tuple[pd.DataFrame, list[dict]]:
    chrom, start, end = locus_region(gwas)
    candidates = []
    audit = []

    for i, context in contexts.iterrows():
        url = str(context["ACCESS_PATH"])
        qtl_type = str(context["QTL_TYPE"])
        status = "OK"
        error = ""
        n_region = 0
        n_candidates = 0

        print(
            f"[CATALOG {i+1}/{len(contexts)}] "
            f"{qtl_type:<8} | "
            f"{context.get('study_label',''):<18} | "
            f"{context.get('tissue_label',''):<35}",
            flush=True,
        )

        try:
            raw = tabix_region(tabix, url, chrom, start, end)
            standardized = standardize_eqtl_catalogue_region(raw, context)
            n_region = len(standardized)

            if standardized.empty:
                status = "NO_REGION_ROWS"
            else:
                selected = standardized[
                    standardized["P"].notna()
                    & (standardized["P"] <= args.candidate_p)
                ].copy()

                if selected.empty:
                    status = "NO_CANDIDATE_P_LE_THRESHOLD"
                else:
                    grouped = (
                        selected.groupby(
                            ["MOLECULAR_TRAIT_ID", "GENE_ID"],
                            dropna=False,
                            as_index=False,
                        )
                        .agg(
                            DISCOVERY_MIN_P=("P", "min"),
                            N_DISCOVERY_ROWS=("P", "size"),
                        )
                        .sort_values("DISCOVERY_MIN_P", kind="stable")
                        .head(args.max_candidates_per_context_locus)
                    )
                    n_candidates = len(grouped)

                    for _, item in grouped.iterrows():
                        candidates.append({
                            "SOURCE": "EQTL_CATALOGUE_DENSE_REGION",
                            "LOCUS_ID": locus_id,
                            "DATASET": str(context.get("FORMAL_DATASET", context.get("study_label", ""))),
                            "DATASET_ID": str(context.get("dataset_id", "")),
                            "QTL_TYPE": qtl_type,
                            "TISSUE": str(context.get("tissue_label", "")),
                            "CONDITION": str(context.get("condition_label", "")),
                            "MOLECULAR_TRAIT_ID": str(item["MOLECULAR_TRAIT_ID"]),
                            "GENE_ID": str(item["GENE_ID"]),
                            "DISCOVERY_MIN_P": item["DISCOVERY_MIN_P"],
                            "N_DISCOVERY_ROWS": item["N_DISCOVERY_ROWS"],
                            "QTL_N": safe_int(context.get("sample_size")),
                            "SOURCE_URL": str(context.get("REMOTE_SOURCE_URL", context.get("ftp_path", ""))),
                            "ACCESS_PATH": url,
                            "ACCESS_MODE": "local_tabix",
                        })

        except Exception as exc:
            status = "ERROR"
            error = f"{type(exc).__name__}: {exc}"

        audit.append({
            "LOCUS_ID": locus_id,
            "DATASET_ID": context.get("dataset_id", ""),
            "STUDY_LABEL": context.get("study_label", ""),
            "QTL_TYPE": qtl_type,
            "TISSUE": context.get("tissue_label", ""),
            "N_REGION_ROWS": n_region,
            "N_CANDIDATES": n_candidates,
            "STATUS": status,
            "ERROR": error,
        })

    return pd.DataFrame(candidates), audit


def load_custom_registry(root: Path, path_value: str) -> pd.DataFrame:
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = root / path

    if not path.exists():
        write_custom_registry_template(root)

    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=CUSTOM_REGISTRY_COLUMNS)

    x = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in CUSTOM_REGISTRY_COLUMNS:
        if col not in x.columns:
            x[col] = ""
    x["ENABLED_BOOL"] = x["ENABLED"].map(boolish)
    return x


def ensure_custom_registry_files(
    root: Path,
    registry: pd.DataFrame,
    auto_download: bool = False,
) -> pd.DataFrame:
    """Resolve custom registry LOCAL_PATH values without downloading anything."""
    if registry.empty:
        return registry

    x = registry.copy()
    for idx, row in x[x["ENABLED_BOOL"]].iterrows():
        local_text = str(row["LOCAL_PATH"]).strip()
        if not local_text:
            continue
        local = Path(local_text).expanduser()
        if not local.is_absolute():
            local = root / local
        x.at[idx, "LOCAL_PATH"] = str(local.resolve())
    return x

def custom_coordinates(x: pd.DataFrame, row: pd.Series) -> pd.DataFrame:
    chr_col = str(row["CHR_COLUMN"]).strip()
    pos_col = str(row["POS_COLUMN"]).strip()
    ref_col = str(row["REF_COLUMN"]).strip()
    alt_col = str(row["ALT_COLUMN"]).strip()
    variant_col = str(row["VARIANT_COLUMN"]).strip()

    out = pd.DataFrame(index=x.index)

    if all(
        c and c in x.columns
        for c in [chr_col, pos_col, ref_col, alt_col]
    ):
        out["CHR"] = x[chr_col].map(clean_chr)
        out["POS"] = numeric(x[pos_col])
        out["REF"] = x[ref_col].map(clean_allele)
        out["ALT"] = x[alt_col].map(clean_allele)
    elif variant_col and variant_col in x.columns:
        parsed = x[variant_col].map(parse_variant_id)
        out["CHR"] = parsed.map(lambda z: z[0])
        out["POS"] = parsed.map(lambda z: z[1])
        out["REF"] = parsed.map(lambda z: z[2])
        out["ALT"] = parsed.map(lambda z: z[3])
    else:
        raise RuntimeError(
            f"Custom registry row {row.get('DATASET')} has no usable coordinates."
        )

    out["VARIANT_KEY"] = [
        make_variant_key(c, p, r, a)
        for c, p, r, a in zip(out["CHR"], out["POS"], out["REF"], out["ALT"])
    ]
    return out


def standardize_custom_region(raw: pd.DataFrame, row: pd.Series) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame()

    trait_col = str(row["TRAIT_COLUMN"]).strip()
    gene_col = str(row["GENE_COLUMN"]).strip()
    beta_col = str(row["BETA_COLUMN"]).strip()
    se_col = str(row["SE_COLUMN"]).strip()
    p_col = str(row["P_COLUMN"]).strip()
    maf_col = str(row["MAF_COLUMN"]).strip()

    for needed in [trait_col, beta_col, se_col]:
        if not needed or needed not in raw.columns:
            raise RuntimeError(
                f"Custom dataset {row.get('DATASET')} missing column {needed!r}."
            )

    out = custom_coordinates(raw, row)
    out["MOLECULAR_TRAIT_ID"] = raw[trait_col].fillna("").astype(str)
    out["GENE_ID"] = (
        raw[gene_col].fillna("").astype(str)
        if gene_col and gene_col in raw.columns
        else ""
    )
    out["BETA_RAW"] = numeric(raw[beta_col])
    out["SE"] = numeric(raw[se_col])
    out["P"] = numeric(raw[p_col]) if p_col and p_col in raw.columns else np.nan
    out["MAF"] = numeric(raw[maf_col]) if maf_col and maf_col in raw.columns else np.nan

    mode = str(row["EFFECT_ALLELE_MODE"]).strip().upper() or "ALT"

    if mode == "ALT":
        out["BETA_REF"] = -out["BETA_RAW"]
    elif mode == "REF":
        out["BETA_REF"] = out["BETA_RAW"]
    elif mode == "COLUMN":
        ea_col = str(row["EFFECT_ALLELE_COLUMN"]).strip()
        if not ea_col or ea_col not in raw.columns:
            raise RuntimeError(
                "EFFECT_ALLELE_MODE=COLUMN requires EFFECT_ALLELE_COLUMN."
            )
        values = []
        for ea, ref, alt, beta in zip(
            raw[ea_col].map(clean_allele),
            out["REF"], out["ALT"], out["BETA_RAW"],
        ):
            if ea == ref:
                values.append(beta)
            elif ea == alt:
                values.append(-beta)
            else:
                values.append(np.nan)
        out["BETA_REF"] = values
    else:
        raise RuntimeError(f"Unsupported EFFECT_ALLELE_MODE={mode!r}")

    out["DATASET"] = str(row["DATASET"])
    out["QTL_TYPE"] = str(row["QTL_TYPE"])
    out["TISSUE"] = str(row["TISSUE"])
    out["CONDITION"] = str(row["CONDITION"])
    out["QTL_N"] = safe_int(row["SAMPLE_SIZE"])

    out = out[
        out["VARIANT_KEY"].ne("")
        & out["MOLECULAR_TRAIT_ID"].str.strip().ne("")
        & out["BETA_REF"].notna()
        & out["SE"].notna()
        & (out["SE"] > 0)
    ].copy()

    return (
        out.drop_duplicates(
            ["MOLECULAR_TRAIT_ID", "VARIANT_KEY"],
            keep="first",
        )
        .reset_index(drop=True)
    )


def load_custom_region(
    tabix: Path,
    row: pd.Series,
    chrom: int,
    start: int,
    end: int,
) -> pd.DataFrame:
    access = str(row["ACCESS_MODE"]).strip().lower() or "file"

    if access in {"tabix", "remote_tabix"}:
        source = str(row["LOCAL_PATH"]).strip()
        if not source:
            raise RuntimeError(
                f"{row.get('DATASET')} has no LOCAL_PATH; remote custom QTL access is disabled."
            )
        return standardize_custom_region(
            tabix_region(tabix, source, chrom, start, end),
            row,
        )

    path_text = str(row["LOCAL_PATH"]).strip()
    if not path_text:
        raise RuntimeError(f"{row.get('DATASET')} has no LOCAL_PATH.")

    path = Path(path_text)
    if not path.exists():
        raise FileNotFoundError(path)

    if Path(str(path) + ".tbi").exists():
        raw = tabix_region(tabix, str(path), chrom, start, end)
    else:
        raw = cached_custom_file(path)

    standardized = standardize_custom_region(raw, row)
    return standardized[
        (standardized["CHR"] == chrom)
        & (standardized["POS"] >= start)
        & (standardized["POS"] <= end)
    ].copy()


def discover_custom_candidates_for_locus(
    tabix: Path,
    registry: pd.DataFrame,
    locus_id: str,
    gwas: pd.DataFrame,
    args,
    allowed_types: set[str],
) -> tuple[pd.DataFrame, list[dict]]:
    chrom, start, end = locus_region(gwas)
    candidates = []
    audit = []

    for _, row in registry[registry["ENABLED_BOOL"]].iterrows():
        qtl_type = str(row["QTL_TYPE"])
        if qtl_type.lower() not in allowed_types:
            continue

        status = "OK"
        error = ""
        n_region = 0
        n_candidates = 0

        try:
            region = load_custom_region(tabix, row, chrom, start, end)
            n_region = len(region)

            selected = region[
                region["P"].notna()
                & (region["P"] <= args.candidate_p)
            ].copy()

            if region.empty:
                status = "NO_REGION_ROWS"
            elif selected.empty:
                status = "NO_CANDIDATE_P_LE_THRESHOLD"
            else:
                grouped = (
                    selected.groupby(
                        ["MOLECULAR_TRAIT_ID", "GENE_ID"],
                        dropna=False,
                        as_index=False,
                    )
                    .agg(
                        DISCOVERY_MIN_P=("P", "min"),
                        N_DISCOVERY_ROWS=("P", "size"),
                    )
                    .sort_values("DISCOVERY_MIN_P", kind="stable")
                    .head(args.max_candidates_per_context_locus)
                )
                n_candidates = len(grouped)

                for _, item in grouped.iterrows():
                    candidates.append({
                        "SOURCE": "CUSTOM_QTL_REGISTRY",
                        "LOCUS_ID": locus_id,
                        "DATASET": str(row["DATASET"]),
                        "DATASET_ID": str(row["DATASET"]),
                        "QTL_TYPE": qtl_type,
                        "TISSUE": str(row["TISSUE"]),
                        "CONDITION": str(row["CONDITION"]),
                        "MOLECULAR_TRAIT_ID": str(item["MOLECULAR_TRAIT_ID"]),
                        "GENE_ID": str(item["GENE_ID"]),
                        "DISCOVERY_MIN_P": item["DISCOVERY_MIN_P"],
                        "N_DISCOVERY_ROWS": item["N_DISCOVERY_ROWS"],
                        "QTL_N": safe_int(row["SAMPLE_SIZE"]),
                        "SOURCE_URL": str(row["SOURCE_URL"]),
                        "ACCESS_MODE": "custom_registry",
                    })

        except Exception as exc:
            status = "ERROR"
            error = f"{type(exc).__name__}: {exc}"

        audit.append({
            "LOCUS_ID": locus_id,
            "DATASET": row["DATASET"],
            "QTL_TYPE": qtl_type,
            "TISSUE": row["TISSUE"],
            "N_REGION_ROWS": n_region,
            "N_CANDIDATES": n_candidates,
            "STATUS": status,
            "ERROR": error,
        })

    return pd.DataFrame(candidates), audit


def context_from_catalogue_candidate(
    contexts: pd.DataFrame,
    candidate: pd.Series,
) -> pd.Series | None:
    if contexts.empty:
        return None
    match = contexts[
        contexts["dataset_id"].astype(str)
        == str(candidate.get("DATASET_ID", ""))
    ]
    return None if match.empty else match.iloc[0]


def dense_qtl_for_candidate(
    tabix: Path,
    candidate: pd.Series,
    gwas: pd.DataFrame,
    catalogue_contexts: pd.DataFrame,
    custom_registry: pd.DataFrame,
) -> pd.DataFrame:
    chrom, start, end = locus_region(gwas)
    source = str(candidate["SOURCE"])

    if source == "EQTL_CATALOGUE_DENSE_REGION":
        context = context_from_catalogue_candidate(catalogue_contexts, candidate)
        if context is None:
            return pd.DataFrame()
        region = standardize_eqtl_catalogue_region(
            tabix_region(tabix, str(context["ACCESS_PATH"]), chrom, start, end),
            context,
        )

    elif source == "CUSTOM_QTL_REGISTRY":
        enabled = custom_registry[custom_registry["ENABLED_BOOL"]].copy()
        match = enabled[
            (enabled["DATASET"].astype(str) == str(candidate["DATASET"]))
            & (
                enabled["QTL_TYPE"].astype(str).str.lower()
                == str(candidate["QTL_TYPE"]).lower()
            )
            & (enabled["TISSUE"].astype(str) == str(candidate["TISSUE"]))
        ]
        if match.empty:
            return pd.DataFrame()
        region = load_custom_region(tabix, match.iloc[0], chrom, start, end)
    else:
        return pd.DataFrame()

    if region.empty:
        return region

    return region[
        region["MOLECULAR_TRAIT_ID"].astype(str)
        == str(candidate["MOLECULAR_TRAIT_ID"])
    ].reset_index(drop=True)


def reference_prefix(root: Path, ancestry_code: str, chrom: int) -> Path:
    return (
        root / "resources" / "1000G" / ancestry_code
        / f"chr{chrom}_{ancestry_code}_GRCh38"
    )


def build_external_ld(
    root: Path,
    plink2: Path,
    ancestry_code: str,
    common: pd.DataFrame,
    out_dir: Path,
    threads: int,
    maf: float = DEFAULT_LD_MAF,
) -> tuple[Path, Path, list[str], Path]:
    """Build/cache ONE signed-LD matrix + frequency table per GWAS locus.

    Candidate-specific QTL tests subset this locus-wide matrix in R, so PLINK2
    is not rerun for every gene/protein/splicing trait.
    """
    chroms = sorted(set(common["CHR"].dropna().astype(int)))
    if len(chroms) != 1:
        raise RuntimeError(f"Cannot construct one LD matrix across chromosomes: {chroms}")
    chrom = chroms[0]

    prefix = reference_prefix(root, ancestry_code, chrom)
    for ext in [".pgen", ".pvar", ".psam"]:
        path = Path(str(prefix) + ext)
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing 1000G reference file: {path}")

    out_dir.mkdir(parents=True, exist_ok=True)
    ids_file = out_dir / "common_reference_ids.txt"
    ld_prefix = out_dir / "external_LD"
    ld_file = Path(str(ld_prefix) + ".unphased.vcor1.bin")
    vars_file = Path(str(ld_prefix) + ".unphased.vcor1.bin.vars")
    freq_file = out_dir / "external_FREQ.afreq"

    # Fast resume: reuse a validated locus-wide LD/frequency cache.
    if ld_file.exists() and vars_file.exists() and freq_file.exists():
        order = [x.strip() for x in vars_file.read_text().splitlines() if x.strip()]
        expected = len(order) * len(order) * 4
        if order and ld_file.stat().st_size == expected and freq_file.stat().st_size > 0:
            print(f"[CACHE] locus LD/frequency: {out_dir}")
            return ld_file, vars_file, order, freq_file

    ids = common["REFERENCE_ID"].dropna().astype(str).drop_duplicates()
    ids_file.write_text("\n".join(ids) + "\n", encoding="utf-8")

    run_command(
        [
            plink2,
            "--pfile", prefix,
            "--extract", ids_file,
            "--maf", str(maf),
            "--geno", "0.05",
            "--min-alleles", "2",
            "--max-alleles", "2",
            "--r-unphased", "square", "bin4", "ref-based", "yes-really",
            "--threads", str(threads),
            "--out", ld_prefix,
        ],
        log_file=out_dir / "plink_ld.log",
    )

    if not ld_file.exists() or not vars_file.exists():
        raise RuntimeError("PLINK2 did not generate signed-r LD output.")

    order = [line.strip() for line in vars_file.read_text().splitlines() if line.strip()]
    expected = len(order) * len(order) * 4
    actual = ld_file.stat().st_size
    if actual != expected:
        raise RuntimeError(
            f"LD validation failed: p={len(order)}, expected={expected} bytes, actual={actual}"
        )

    freq_prefix = out_dir / "external_FREQ"
    run_command(
        [
            plink2,
            "--pfile", prefix,
            "--extract", vars_file,
            "--freq",
            "--threads", str(threads),
            "--out", freq_prefix,
        ],
        log_file=out_dir / "plink_freq.log",
    )
    if not freq_file.exists() or freq_file.stat().st_size == 0:
        raise RuntimeError("PLINK2 did not generate reference allele frequencies.")

    return ld_file, vars_file, order, freq_file

def load_reference_maf(freq_file: Path) -> dict[str, float]:
    if not freq_file.exists() or freq_file.stat().st_size == 0:
        raise FileNotFoundError(f"Missing PLINK frequency file: {freq_file}")

    x = pd.read_csv(freq_file, sep=r"\s+", dtype=str, comment=None)
    if "ID" not in x.columns:
        raise RuntimeError(
            f"PLINK frequency file lacks ID column: {freq_file}; columns={list(x.columns)}"
        )
    freq_col = "ALT_FREQS" if "ALT_FREQS" in x.columns else (
        "ALT_FREQ" if "ALT_FREQ" in x.columns else None
    )
    if freq_col is None:
        raise RuntimeError(
            f"PLINK frequency file lacks ALT frequency column: {freq_file}; columns={list(x.columns)}"
        )

    af = pd.to_numeric(
        x[freq_col].astype(str).str.split(",").str[0],
        errors="coerce",
    )
    maf = np.minimum(af, 1.0 - af)
    return dict(zip(x["ID"].astype(str), maf))


def step06_susie_paths(study_dir: Path, locus_id: str, ancestry_code: str) -> dict[str, Path]:
    locus = study_dir / "loci" / locus_id
    return {
        "fit": locus / "susie" / "susie_fit.rds",
        "diag": locus / "susie" / "susie_diagnostics.tsv",
        "vars": locus / f"{ancestry_code}_LD.unphased.vcor1.bin.vars",
    }


def prepare_gwas_susie_map(
    study_dir: Path,
    locus_id: str,
    ancestry_code: str,
    gwas: pd.DataFrame,
    out_dir: Path,
) -> tuple[Path | None, Path | None]:
    paths = step06_susie_paths(study_dir, locus_id, ancestry_code)
    if not paths["fit"].exists() or not paths["vars"].exists():
        return None, None

    order = [x.strip() for x in paths["vars"].read_text().splitlines() if x.strip()]
    mapping = gwas[["REFERENCE_ID", "VARIANT_KEY"]].drop_duplicates("REFERENCE_ID")
    order_map = {v: i for i, v in enumerate(order)}
    mapping = mapping[mapping["REFERENCE_ID"].isin(order_map)].copy()
    mapping["_ORDER"] = mapping["REFERENCE_ID"].map(order_map)
    mapping = mapping.sort_values("_ORDER").drop(columns="_ORDER")

    if len(mapping) != len(order):
        return None, None

    path = out_dir / "gwas_susie_variant_map.tsv"
    mapping.to_csv(path, sep="\t", index=False)
    return paths["fit"], path


def read_step06_s_rss(study_dir: Path, locus_id: str) -> float:
    diag = study_dir / "loci" / locus_id / "susie" / "susie_diagnostics.tsv"
    x = read_table(diag)
    if x.empty or "ESTIMATE_S_RSS" not in x.columns:
        return np.nan
    vals = numeric(x["ESTIMATE_S_RSS"]).dropna()
    return float(vals.iloc[0]) if not vals.empty else np.nan


def catalogue_context_for_candidate(
    contexts: pd.DataFrame,
    candidate: pd.Series,
) -> pd.Series | None:
    if contexts.empty:
        return None
    x = contexts[
        contexts["dataset_id"].astype(str)
        == str(candidate.get("DATASET_ID", ""))
    ]
    return None if x.empty else x.iloc[0]


R_WORKER = r"""
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 18) stop("Expected 18 arguments")

stats_file <- args[1]
ld_file <- args[2]
vars_file <- args[3]
out_summary <- args[4]
out_results <- args[5]
out_sensitivity <- args[6]
out_diag <- args[7]
gwas_n <- as.numeric(args[8])
qtl_n <- as.numeric(args[9])
gwas_type <- args[10]
case_prop <- suppressWarnings(as.numeric(args[11]))
p1 <- as.numeric(args[12])
p2 <- as.numeric(args[13])
p12 <- as.numeric(args[14])
p12_low <- as.numeric(args[15])
p12_high <- as.numeric(args[16])
gwas_fit_file <- args[17]
gwas_map_file <- args[18]

suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(susieR))

stats <- read.delim(stats_file, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
full_vars <- scan(vars_file, what="character", quiet=TRUE)

# Load the locus-wide LD ONCE, then subset it to this molecular-trait candidate.
p <- length(full_vars)
con <- file(ld_file, "rb")
values <- readBin(con, what="numeric", n=p*p, size=4, endian="little")
close(con)
if (length(values) != p*p) stop("LD binary length mismatch")
Rfull <- matrix(values, nrow=p, ncol=p, byrow=TRUE)
rm(values); gc()
Rfull <- (Rfull + t(Rfull))/2
diag(Rfull) <- 1
dimnames(Rfull) <- list(full_vars, full_vars)

idx <- match(as.character(stats$REFERENCE_ID), full_vars)
keep <- which(!is.na(idx))
if (length(keep) < 20) stop("Too few candidate SNPs represented in locus LD")
stats <- stats[keep,,drop=FALSE]
idx <- idx[keep]

# Remove duplicated LD IDs while preserving common_stats order.
dup <- duplicated(idx)
if (any(dup)) {
    stats <- stats[!dup,,drop=FALSE]
    idx <- idx[!dup]
}
R <- Rfull[idx, idx, drop=FALSE]
rm(Rfull); gc()
vars <- as.character(stats$REFERENCE_ID)
dimnames(R) <- list(vars, vars)

bad <- which(!apply(is.finite(R), 1, all) | !apply(is.finite(R), 2, all))
if (length(bad) > 0) {
    cat("Filtering", length(bad), "variants with non-finite LD rows/columns\n")
    keep <- setdiff(seq_len(nrow(R)), bad)
    R <- R[keep, keep, drop=FALSE]
    stats <- stats[keep,,drop=FALSE]
    vars <- vars[keep]
}
if (nrow(R) < 20) stop("Too few SNPs after non-finite LD filtering")

ok <- is.finite(stats$QTL_BETA) &
      is.finite(stats$QTL_SE) & stats$QTL_SE > 0 &
      is.finite(stats$QTL_MAF) &
      stats$QTL_MAF > 0 & stats$QTL_MAF < 0.5
if (!all(ok)) {
    R <- R[ok, ok, drop=FALSE]
    stats <- stats[ok,,drop=FALSE]
    vars <- vars[ok]
}
if (nrow(R) < 20) stop("Too few finite QTL SNPs after MAF/statistic filtering")
diag(R) <- 1
dimnames(R) <- list(vars, vars)

# Reuse Step06 GWAS SuSiE fit instead of rerunning GWAS SuSiE for every QTL.
s1 <- readRDS(gwas_fit_file)
gmap <- read.delim(gwas_map_file, stringsAsFactors=FALSE, check.names=FALSE)
if (!"REFERENCE_ID" %in% names(gmap)) stop("GWAS SuSiE map lacks REFERENCE_ID")
gwas_snps <- as.character(gmap$REFERENCE_ID)
if (is.null(s1$lbf_variable)) stop("Step06 SuSiE fit has no lbf_variable")
if (ncol(s1$lbf_variable) > length(gwas_snps) + 1) {
    stop("Step06 SuSiE lbf_variable has more columns than the variant map")
}
colnames(s1$lbf_variable) <- c(gwas_snps, "null")[seq_len(ncol(s1$lbf_variable))]
if (!is.null(s1$alpha)) {
    if (ncol(s1$alpha) > length(gwas_snps) + 1) stop("Step06 alpha/map dimension mismatch")
    colnames(s1$alpha) <- c(gwas_snps, "null")[seq_len(ncol(s1$alpha))]
}
if (!is.null(s1$pip)) {
    if (length(s1$pip) != length(gwas_snps)) stop("Step06 PIP/map dimension mismatch")
    names(s1$pip) <- gwas_snps
}
if (!is.null(s1$sets$cs) && length(s1$sets$cs)) {
    s1$sets$cs <- lapply(s1$sets$cs, function(x) {
        good <- x[x >= 1 & x <= length(gwas_snps)]
        names(good) <- gwas_snps[good]
        good
    })
}
class(s1) <- unique(c("susie", class(s1)))

# Fine-map only the molecular QTL side for this candidate.
d2 <- list(
    beta=as.numeric(stats$QTL_BETA),
    varbeta=as.numeric(stats$QTL_SE)^2,
    snp=as.character(stats$REFERENCE_ID),
    position=as.integer(stats$POS),
    N=as.numeric(qtl_n),
    type="quant",
    MAF=as.numeric(stats$QTL_MAF),
    LD=R
)
coloc::check_dataset(d2, req=c("type","snp","LD","N"))
cat("QTL CHECK = PASS\n")

z2 <- d2$beta / sqrt(d2$varbeta)
s_qtl <- tryCatch(susieR::estimate_s_rss(z=z2, R=R, n=qtl_n), error=function(e) NA_real_)
diag_df <- data.frame(
    TRAIT="QTL",
    ESTIMATE_S_RSS=s_qtl,
    MAX_ABS_Z=max(abs(z2), na.rm=TRUE),
    N_SNPS=length(z2),
    stringsAsFactors=FALSE
)
write.table(diag_df, out_diag, sep="\t", quote=FALSE, row.names=FALSE)

cat("Running coloc::runsusie for molecular QTL only\n")
s2 <- coloc::runsusie(
    d2, suffix="QTL", maxit=200, repeat_until_convergence=TRUE
)
cat("QTL SuSiE COMPLETE; reusing Step06 GWAS SuSiE\n")

run_susie_coloc <- function(prior12, label) {
    res <- suppressWarnings(coloc::coloc.susie(s1, s2, p1=p1, p2=p2, p12=prior12))
    if (!is.list(res) || is.null(res$summary)) {
        sm <- data.frame(
            NO_SIGNAL_PAIR="YES", METHOD="coloc.susie_step06_reuse",
            PRIOR_SET=label, P12=prior12, stringsAsFactors=FALSE
        )
        return(list(result=list(results=NULL), summary=sm))
    }
    sm <- tryCatch(as.data.frame(res$summary), error=function(e) data.frame())
    if (nrow(sm) == 0) {
        sm <- data.frame(
            NO_SIGNAL_PAIR="YES", METHOD="coloc.susie_step06_reuse",
            PRIOR_SET=label, P12=prior12, stringsAsFactors=FALSE
        )
    } else {
        sm$NO_SIGNAL_PAIR <- "NO"
        sm$METHOD <- "coloc.susie_step06_reuse"
        sm$PRIOR_SET <- label
        sm$P12 <- prior12
    }
    list(result=res, summary=sm)
}

base <- run_susie_coloc(p12, "DEFAULT")
low <- run_susie_coloc(p12_low, "LOW_P12")
high <- run_susie_coloc(p12_high, "HIGH_P12")

write.table(base$summary, out_summary, sep="\t", quote=FALSE, row.names=FALSE)
results_df <- tryCatch(as.data.frame(base$result$results), error=function(e) data.frame())
write.table(results_df, out_results, sep="\t", quote=FALSE, row.names=FALSE)
sens <- rbind(base$summary, low$summary, high$summary)
write.table(sens, out_sensitivity, sep="\t", quote=FALSE, row.names=FALSE)

cat("\n===== COLOC SUMMARY =====\n")
print(base$summary)
cat("\n===== QTL LD MISMATCH DIAGNOSTICS =====\n")
print(diag_df)
"""


def write_r_worker(output_root: Path) -> Path:
    path = output_root / "_step09_coloc_susie_worker.R"
    path.write_text(R_WORKER, encoding="utf-8")
    return path


def verify_r_packages(rscript: Path) -> dict:
    text = run_command(
        [
            rscript,
            "-e",
            (
                'stopifnot(requireNamespace("coloc", quietly=TRUE)); '
                'stopifnot(requireNamespace("susieR", quietly=TRUE)); '
                'cat("coloc=", as.character(packageVersion("coloc")), "\\n", sep=""); '
                'cat("susieR=", as.character(packageVersion("susieR")), "\\n", sep="")'
            ),
        ],
        capture=True,
    )
    result = {}
    for line in text.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def posterior_column(columns: Iterable[str], hypothesis: str) -> str | None:
    target = hypothesis.upper()
    patterns = [f"PP.{target}", f"PP_{target}", f"PP.{target}.ABF", target]
    for column in columns:
        upper = str(column).upper()
        for pattern in patterns:
            if pattern.upper() in upper:
                return column
    return None


def _best_coloc_row(summary: pd.DataFrame) -> pd.Series | None:
    if summary.empty:
        return None
    h4_col = posterior_column(summary.columns, "H4")
    if h4_col is None:
        return summary.iloc[0]
    vals = numeric(summary[h4_col])
    if not vals.notna().any():
        return summary.iloc[0]
    return summary.loc[vals.idxmax()]


def _row_posterior(row: pd.Series | None, summary: pd.DataFrame, hypothesis: str) -> float:
    if row is None:
        return np.nan
    col = posterior_column(summary.columns, hypothesis)
    if col is None:
        return np.nan
    try:
        return float(row[col])
    except Exception:
        return np.nan


def _hit_value(row: pd.Series | None, names: list[str]) -> str:
    if row is None:
        return ""
    lookup = {str(c).lower(): c for c in row.index}
    for name in names:
        if name.lower() in lookup:
            value = row[lookup[name.lower()]]
            if pd.notna(value):
                return str(value)
    return ""


def classify_coloc(
    summary: pd.DataFrame,
    strong_pp4: float,
    suggestive_pp4: float,
) -> dict:
    row = _best_coloc_row(summary)
    values = {
        h: _row_posterior(row, summary, h)
        for h in ["H0", "H1", "H2", "H3", "H4"]
    }
    h3 = values["H3"]
    h4 = values["H4"]

    hit1 = _hit_value(row, ["hit1", "signal1", "cs1"])
    hit2 = _hit_value(row, ["hit2", "signal2", "cs2"])
    method = _hit_value(row, ["METHOD"]) or "coloc.susie"

    if pd.isna(h4):
        classification = "NO_PP_H4_RETURNED"
        strong = "NO"
        shared = "UNKNOWN"
    elif h4 >= strong_pp4 and (pd.isna(h3) or h4 > h3):
        classification = "STRONG_SHARED_SIGNAL"
        strong = "YES"
        shared = "YES"
    elif h4 >= suggestive_pp4:
        classification = "SUGGESTIVE_SHARED_SIGNAL"
        strong = "NO"
        shared = "SUGGESTIVE"
    elif pd.notna(h3) and h3 > h4:
        classification = "BOTH_ASSOCIATED_DIFFERENT_SIGNALS_FAVORED"
        strong = "NO"
        shared = "NO"
    else:
        classification = "NO_STRONG_SHARED_SIGNAL"
        strong = "NO"
        shared = "NO"

    return {
        "MAX_PP_H0": values["H0"],
        "MAX_PP_H1": values["H1"],
        "MAX_PP_H2": values["H2"],
        "MAX_PP_H3": values["H3"],
        "MAX_PP_H4": values["H4"],
        "BEST_HIT1": hit1,
        "BEST_HIT2": hit2,
        "BEST_SIGNAL_PAIR": (f"{hit1}|{hit2}" if hit1 or hit2 else ""),
        "COLOC_METHOD": method,
        "COLOCALIZED_STRONG": strong,
        "SHARED_SIGNAL": shared,
        "COLOC_CLASS": classification,
    }


def sensitivity_pp4_for_pair(
    sensitivity: pd.DataFrame,
    prior_label: str,
    hit1: str,
    hit2: str,
) -> float:
    if sensitivity.empty or "PRIOR_SET" not in sensitivity.columns:
        return np.nan
    sub = sensitivity[sensitivity["PRIOR_SET"].astype(str) == prior_label].copy()
    if sub.empty:
        return np.nan

    # Keep the same signal pair used for the default classification whenever possible.
    if hit1 and hit2:
        h1col = first_existing(sub.columns, ["hit1", "signal1", "cs1"])
        h2col = first_existing(sub.columns, ["hit2", "signal2", "cs2"])
        if h1col and h2col:
            exact = sub[
                (sub[h1col].astype(str) == hit1)
                & (sub[h2col].astype(str) == hit2)
            ]
            if not exact.empty:
                sub = exact

    row = _best_coloc_row(sub)
    return _row_posterior(row, sub, "H4")


def process_candidate(
    root: Path,
    accession: str,
    study_dir: Path,
    ancestry_code: str,
    gwas_meta: dict,
    candidate: pd.Series,
    catalogue_contexts: pd.DataFrame,
    custom_registry: pd.DataFrame,
    plink2: Path,
    tabix: Path,
    rscript: Path,
    r_worker: Path,
    locus_output_dir: Path,
    locus_resources: dict,
    args,
) -> dict:
    locus_id = str(candidate["LOCUS_ID"])
    qtl_type = str(candidate["QTL_TYPE"])
    dataset = str(candidate["DATASET"])
    tissue = str(candidate["TISSUE"])
    condition = str(candidate.get("CONDITION", ""))
    trait_id = str(candidate["MOLECULAR_TRAIT_ID"])
    gene_id = str(candidate.get("GENE_ID", ""))

    candidate_id = slugify(
        f"{locus_id}_{dataset}_{qtl_type}_{tissue}_{condition}_{trait_id}"
    )[:220]
    out = locus_output_dir / "candidates" / candidate_id
    out.mkdir(parents=True, exist_ok=True)
    result_json = out / "candidate_result.json"

    if not args.force and result_json.exists():
        previous = safe_json(result_json)
        if previous.get("STATUS") == "COMPLETE":
            print(f"[RESUME] {candidate_id}")
            return previous

    result = {
        "STUDY_ACCESSION": accession,
        "LOCUS_ID": locus_id,
        "QTL_TYPE": qtl_type,
        "DATASET": dataset,
        "DATASET_ID": candidate.get("DATASET_ID", ""),
        "TISSUE": tissue,
        "CONDITION": condition,
        "MOLECULAR_TRAIT_ID": trait_id,
        "GENE_ID": gene_id,
        "DISCOVERY_SOURCE": candidate.get("SOURCE", ""),
        "DISCOVERY_MIN_P": candidate.get("DISCOVERY_MIN_P", np.nan),
        "QTL_N": candidate.get("QTL_N", np.nan),
        "RESOURCE_AVAILABLE": "YES",
        "CANDIDATE_PRESENT": "YES",
        "COLOC_TESTED": "NO",
        "N_GWAS_VARIANTS": 0,
        "N_QTL_VARIANTS": 0,
        "N_COMMON_SNPS": 0,
        "LD_SOURCE": "",
        "QTL_MAF_SOURCE": "",
        "QTL_MAF_REFERENCE_FALLBACK_N": 0,
        "GWAS_STEP06_S_RSS": read_step06_s_rss(study_dir, locus_id),
        "GWAS_S_RSS": read_step06_s_rss(study_dir, locus_id),
        "QTL_S_RSS": np.nan,
        "LD_MISMATCH_WARNING": "NO",
        "MAX_PP_H0": np.nan,
        "MAX_PP_H1": np.nan,
        "MAX_PP_H2": np.nan,
        "MAX_PP_H3": np.nan,
        "MAX_PP_H4": np.nan,
        "PP_H4_P12_LOW": np.nan,
        "PP_H4_P12_DEFAULT": np.nan,
        "PP_H4_P12_HIGH": np.nan,
        "PRIOR_ROBUST_STRONG": "NO",
        "BEST_HIT1": "",
        "BEST_HIT2": "",
        "BEST_SIGNAL_PAIR": "",
        "COLOC_METHOD": "",
        "COLOCALIZED_STRONG": "NO",
        "SHARED_SIGNAL": "UNKNOWN",
        "COLOC_CLASS": "",
        "STATUS": "",
        "ERROR": "",
        "OUTPUT_DIR": str(out),
    }

    def finish(status: str | None = None):
        if status is not None:
            result["STATUS"] = status
        write_json(result_json, result)
        return result

    try:
        gwas = locus_resources["gwas"]
        result["N_GWAS_VARIANTS"] = len(gwas)

        # LOCAL dense path only: no eQTL-Catalogue LBF download/use.
        qtl = dense_qtl_for_candidate(
            tabix, candidate, gwas, catalogue_contexts, custom_registry
        )
        result["N_QTL_VARIANTS"] = len(qtl)
        if qtl.empty:
            return finish("NOT_TESTED_NO_DENSE_QTL_ROWS")

        qtl_n = None
        if "QTL_N" in qtl.columns:
            qn = numeric(qtl["QTL_N"]).dropna()
            qn = qn[qn > 0]
            if not qn.empty:
                qtl_n = int(round(float(qn.median())))
        if not qtl_n:
            qtl_n = safe_int(candidate.get("QTL_N"))
        if not qtl_n:
            return finish("NOT_TESTED_QTL_SAMPLE_SIZE_MISSING")
        result["QTL_N"] = qtl_n

        common = gwas.merge(
            qtl[["VARIANT_KEY", "BETA_REF", "SE", "P", "MAF"]].rename(
                columns={
                    "BETA_REF": "QTL_BETA",
                    "SE": "QTL_SE",
                    "P": "QTL_P",
                    "MAF": "QTL_MAF",
                }
            ),
            on="VARIANT_KEY",
            how="inner",
            validate="one_to_one",
        )
        common = common[
            common["GWAS_BETA"].notna()
            & common["GWAS_SE"].notna()
            & (common["GWAS_SE"] > 0)
            & common["QTL_BETA"].notna()
            & common["QTL_SE"].notna()
            & (common["QTL_SE"] > 0)
            & common["REFERENCE_ID"].notna()
        ].drop_duplicates("REFERENCE_ID", keep="first")
        result["N_COMMON_SNPS"] = len(common)

        if len(common) < args.min_common_snps:
            return finish("NOT_TESTED_TOO_FEW_COMMON_SNPS")
        if len(common) > args.max_common_snps:
            result["ERROR"] = (
                f"{len(common)} common SNPs > {args.max_common_snps}; no significance-based subsampling."
            )
            return finish("NOT_TESTED_TOO_MANY_COMMON_SNPS_NO_SUBSAMPLING")

        # The GWAS side is the Step06 SuSiE fit (fitted with Step06's N), so
        # N / case fraction are recorded here but not needed by coloc.susie.
        gwas_n = safe_int(gwas_meta.get("COLOC_GWAS_N")) or 0

        if "ld" not in locus_resources:
            ld_dir = locus_output_dir / "locus_ld"
            locus_resources["ld"] = build_external_ld(
                root, plink2, ancestry_code, gwas, ld_dir, args.cpus, maf=args.ld_maf
            )
        ld_file, vars_file, ld_order, freq_file = locus_resources["ld"]
        result["LD_SOURCE"] = f"1000G_{ancestry_code}_SIGNED_R_MAF_GE_{args.ld_maf}"

        available = set(ld_order)
        common = common[common["REFERENCE_ID"].isin(available)].copy()

        reference_maf = load_reference_maf(freq_file)
        common["GWAS_MAF"] = common["REFERENCE_ID"].astype(str).map(reference_maf)

        # Prefer cohort QTL MAF when supplied.  If the dense QTL file does not
        # carry MAF, use the same ancestry-matched 1000G reference MAF used for
        # LD.  coloc needs MAF+N for quantitative traits when sdY is unavailable.
        qtl_raw = numeric(common["QTL_MAF"])
        qtl_raw = np.minimum(qtl_raw, 1.0 - qtl_raw)
        valid_qtl = qtl_raw.notna() & (qtl_raw > 0) & (qtl_raw < 0.5)
        fallback_mask = ~valid_qtl
        common["QTL_MAF"] = qtl_raw.where(valid_qtl, common["GWAS_MAF"])
        result["QTL_MAF_REFERENCE_FALLBACK_N"] = int(fallback_mask.sum())
        result["QTL_MAF_SOURCE"] = (
            "QTL_FILE+1000G_FALLBACK" if valid_qtl.any() and fallback_mask.any()
            else "QTL_FILE" if valid_qtl.all()
            else f"1000G_{ancestry_code}_REFERENCE"
        )

        if gwas_meta.get("GWAS_TYPE") == "quant":
            common = common[common["GWAS_MAF"].notna()].copy()
        common = common[common["QTL_MAF"].notna()].copy()

        result["N_COMMON_SNPS"] = len(common)
        if len(common) < args.min_common_snps:
            return finish("NOT_TESTED_TOO_FEW_SNPS_AFTER_LD_OR_MAF")

        common_file = out / "common_stats.tsv.gz"
        common[
            [
                "REFERENCE_ID", "CHR", "POS", "REF", "ALT", "VARIANT_KEY",
                "GWAS_BETA", "GWAS_SE", "GWAS_P", "GWAS_MAF",
                "QTL_BETA", "QTL_SE", "QTL_P", "QTL_MAF",
            ]
        ].to_csv(common_file, sep="\t", index=False, compression="gzip")

        gwas_fit_file, gwas_map_file = prepare_gwas_susie_map(
            study_dir, locus_id, ancestry_code, gwas, locus_output_dir
        )
        if gwas_fit_file is None or gwas_map_file is None:
            return finish("NOT_TESTED_STEP06_SUSIE_FIT_OR_MAP_MISSING")

        summary_file = out / "coloc_summary.tsv"
        results_file = out / "coloc_results.tsv"
        sensitivity_file = out / "prior_sensitivity.tsv"
        diagnostics_file = out / "diagnostics.tsv"
        case_prop = gwas_meta.get("CASE_PROPORTION")

        run_command(
            [
                rscript, r_worker, common_file, ld_file, vars_file,
                summary_file, results_file, sensitivity_file, diagnostics_file,
                str(gwas_n), str(qtl_n), str(gwas_meta.get("GWAS_TYPE", "quant")),
                str(case_prop) if case_prop is not None else "NA",
                str(COLOC_PRIOR_P1), str(COLOC_PRIOR_P2), str(COLOC_PRIOR_P12),
                str(COLOC_PRIOR_P12_LOW), str(COLOC_PRIOR_P12_HIGH),
                gwas_fit_file, gwas_map_file,
            ],
            log_file=out / "coloc.log",
        )

        summary = read_table(summary_file)
        no_signal_pair = (
            not summary.empty
            and "NO_SIGNAL_PAIR" in summary.columns
            and summary["NO_SIGNAL_PAIR"].astype(str).eq("YES").all()
        )
        classification = classify_coloc(summary, args.strong_pp4, args.suggestive_pp4)
        result.update(classification)
        sensitivity = read_table(sensitivity_file)
        result["PP_H4_P12_LOW"] = sensitivity_pp4_for_pair(
            sensitivity, "LOW_P12",
            classification.get("BEST_HIT1", ""),
            classification.get("BEST_HIT2", ""),
        )
        result["PP_H4_P12_DEFAULT"] = sensitivity_pp4_for_pair(
            sensitivity, "DEFAULT",
            classification.get("BEST_HIT1", ""),
            classification.get("BEST_HIT2", ""),
        )
        result["PP_H4_P12_HIGH"] = sensitivity_pp4_for_pair(
            sensitivity, "HIGH_P12",
            classification.get("BEST_HIT1", ""),
            classification.get("BEST_HIT2", ""),
        )
        result["PRIOR_ROBUST_STRONG"] = (
            "YES"
            if pd.notna(result["PP_H4_P12_LOW"])
            and result["PP_H4_P12_LOW"] >= args.strong_pp4
            else "NO"
        )

        diagnostics = read_table(diagnostics_file)
        if not diagnostics.empty and "TRAIT" in diagnostics.columns:
            for trait, field in [("GWAS", "GWAS_S_RSS"), ("QTL", "QTL_S_RSS")]:
                sub = diagnostics[diagnostics["TRAIT"].astype(str) == trait]
                if not sub.empty and "ESTIMATE_S_RSS" in sub.columns:
                    vals = numeric(sub["ESTIMATE_S_RSS"]).dropna()
                    if not vals.empty:
                        result[field] = float(vals.iloc[0])
        svals = [result.get("GWAS_S_RSS"), result.get("QTL_S_RSS")]
        result["LD_MISMATCH_WARNING"] = (
            "YES"
            if any(pd.notna(v) and float(v) > args.ld_mismatch_s for v in svals)
            else "NO"
        )

        if results_file.exists() and results_file.stat().st_size > 0:
            rr = read_table(results_file)
            rr.to_csv(
                out / "coloc_variant_results.tsv.gz",
                sep="\t", index=False, compression="gzip"
            )
            results_file.unlink(missing_ok=True)

        result["COLOC_TESTED"] = "YES"
        result["STATUS"] = "COMPLETE_NO_SIGNAL_PAIR" if no_signal_pair else "COMPLETE"
        return finish()

    except Exception as exc:
        result["STATUS"] = "FAILED"
        result["ERROR"] = f"{type(exc).__name__}: {exc}"
        (out / "FAILED.txt").write_text(
            result["ERROR"] + "\n\n" + traceback.format_exc(),
            encoding="utf-8",
        )
        return finish()


def resource_status_table(
    allowed_types: set[str],
    catalogue_contexts: pd.DataFrame,
    custom_registry: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for qtl_type in SUPPORTED_QTL_TYPES:
        if qtl_type.lower() not in allowed_types:
            continue
        catalog_n = int((catalogue_contexts["QTL_TYPE"].astype(str).str.lower() == qtl_type.lower()).sum()) \
            if not catalogue_contexts.empty else 0
        custom_n = int((custom_registry["ENABLED_BOOL"]
                        & (custom_registry["QTL_TYPE"].astype(str).str.lower() == qtl_type.lower())).sum()) \
            if not custom_registry.empty else 0
        formal = catalog_n > 0 or custom_n > 0
        rows.append({
            "QTL_TYPE": qtl_type,
            "FORMAL_RESOURCE_AVAILABLE": formal,
            "N_EQTL_CATALOG_CONTEXTS": catalog_n,
            "N_CUSTOM_CONTEXTS": custom_n,
            "INTERPRETATION": "FORMAL_RESOURCE_AVAILABLE" if formal else "RESOURCE_NOT_CONFIGURED",
        })
    return pd.DataFrame(rows)


def print_candidates(candidates: pd.DataFrame) -> None:
    banner("MOLECULAR-QTL CANDIDATES — ALL CONFIGURED TISSUES / CONTEXTS")
    if candidates.empty:
        print("No candidates found.")
        return

    cols = [
        "LOCUS_ID", "QTL_TYPE", "DATASET", "TISSUE", "CONDITION",
        "MOLECULAR_TRAIT_ID", "GENE_ID", "DISCOVERY_MIN_P", "QTL_N", "SOURCE",
    ]
    cols = [c for c in cols if c in candidates.columns]

    with pd.option_context(
        "display.max_rows", 1000,
        "display.max_columns", None,
        "display.width", 300,
        "display.max_colwidth", 55,
        "display.float_format", lambda value: f"{value:.4g}",
    ):
        print(candidates[cols].to_string(index=False))


def print_formal_summary(summary: pd.DataFrame) -> None:
    banner("FORMAL COLOCALIZATION RESULTS")
    if summary.empty:
        print("No formal candidates/results.")
        return

    cols = [
        "LOCUS_ID", "QTL_TYPE", "DATASET", "TISSUE", "CONDITION",
        "MOLECULAR_TRAIT_ID", "GENE_ID", "COLOC_TESTED", "STATUS",
        "N_COMMON_SNPS", "MAX_PP_H3", "MAX_PP_H4",
        "PP_H4_P12_LOW", "PP_H4_P12_DEFAULT", "PP_H4_P12_HIGH",
        "PRIOR_ROBUST_STRONG", "BEST_HIT1", "BEST_HIT2",
        "GWAS_S_RSS", "QTL_S_RSS", "LD_MISMATCH_WARNING",
        "SHARED_SIGNAL", "COLOC_METHOD", "COLOC_CLASS",
    ]
    cols = [c for c in cols if c in summary.columns]

    with pd.option_context(
        "display.max_rows", 2000,
        "display.max_columns", None,
        "display.width", 340,
        "display.max_colwidth", 55,
        "display.float_format", lambda value: f"{value:.4g}",
    ):
        print(summary[cols].to_string(index=False))


def build_locus_qtl_matrix(
    summary: pd.DataFrame,
    loci: list[str],
    resource_status: pd.DataFrame,
    audit: pd.DataFrame,
    suggestive_pp4: float,
) -> pd.DataFrame:
    rows = []
    status_map = {str(row["QTL_TYPE"]): row for _, row in resource_status.iterrows()}

    for locus in loci:
        for qtl_type, rs in status_map.items():
            local = (
                summary[
                    (summary["LOCUS_ID"].astype(str) == locus)
                    & (summary["QTL_TYPE"].astype(str).str.lower() == qtl_type.lower())
                ] if not summary.empty else pd.DataFrame()
            )
            errors = (
                audit[
                    (audit["LOCUS_ID"].astype(str) == locus)
                    & (audit["QTL_TYPE"].astype(str).str.lower() == qtl_type.lower())
                    & (audit["STATUS"].astype(str) == "ERROR")
                ] if not audit.empty and "QTL_TYPE" in audit.columns else pd.DataFrame()
            )

            tested = int((local["COLOC_TESTED"] == "YES").sum()) if not local.empty else 0
            strong = int((local["COLOCALIZED_STRONG"] == "YES").sum()) if not local.empty else 0
            robust = int((local["PRIOR_ROBUST_STRONG"] == "YES").sum()) if not local.empty and "PRIOR_ROBUST_STRONG" in local.columns else 0
            max_h4 = (
                float(numeric(local["MAX_PP_H4"]).max())
                if not local.empty and numeric(local["MAX_PP_H4"]).notna().any()
                else np.nan
            )

            best_tissue = best_trait = best_gene = best_hit1 = best_hit2 = ""
            if not local.empty and numeric(local["MAX_PP_H4"]).notna().any():
                best = local.loc[numeric(local["MAX_PP_H4"]).idxmax()]
                best_tissue = str(best.get("TISSUE", ""))
                best_trait = str(best.get("MOLECULAR_TRAIT_ID", ""))
                best_gene = str(best.get("GENE_ID", ""))
                best_hit1 = str(best.get("BEST_HIT1", ""))
                best_hit2 = str(best.get("BEST_HIT2", ""))

            if strong > 0:
                evidence = "YES_STRONG_COLOC"
            elif tested > 0 and pd.notna(max_h4) and max_h4 >= suggestive_pp4:
                evidence = "SUGGESTIVE_COLOC"
            elif tested > 0:
                evidence = "NO_STRONG_COLOC"
            elif len(local) > 0:
                evidence = "CANDIDATE_NOT_FORMALLY_TESTED"
            elif len(errors) > 0:
                evidence = "RESOURCE_ERROR"
            elif bool(rs["FORMAL_RESOURCE_AVAILABLE"]):
                evidence = "NO_CANDIDATE_IN_CONFIGURED_RESOURCE"
            else:
                evidence = "RESOURCE_NOT_CONFIGURED"

            rows.append({
                "LOCUS_ID": locus,
                "QTL_TYPE": qtl_type,
                "FORMAL_RESOURCE_AVAILABLE": bool(rs["FORMAL_RESOURCE_AVAILABLE"]),
                "N_RESOURCE_ERRORS": len(errors),
                "N_CANDIDATES": len(local),
                "N_FORMALLY_TESTED": tested,
                "N_STRONG_COLOC": strong,
                "N_PRIOR_ROBUST_STRONG": robust,
                "MAX_PP_H4": max_h4,
                "BEST_TISSUE": best_tissue,
                "BEST_TRAIT": best_trait,
                "BEST_GENE": best_gene,
                "BEST_HIT1": best_hit1,
                "BEST_HIT2": best_hit2,
                "EVIDENCE_STATUS": evidence,
            })

    return pd.DataFrame(rows)


def build_wide_evidence_table(matrix: pd.DataFrame) -> pd.DataFrame:
    if matrix.empty:
        return matrix

    records = []
    for locus, group in matrix.groupby("LOCUS_ID", sort=False):
        record = {"LOCUS_ID": locus}
        for _, row in group.iterrows():
            prefix = re.sub(r"[^A-Za-z0-9]+", "_", str(row["QTL_TYPE"])).upper()
            record[f"{prefix}_STATUS"] = row["EVIDENCE_STATUS"]
            record[f"{prefix}_MAX_PP_H4"] = row["MAX_PP_H4"]
            record[f"{prefix}_BEST_TISSUE"] = row["BEST_TISSUE"]
            record[f"{prefix}_BEST_TRAIT"] = row["BEST_TRAIT"]
            record[f"{prefix}_BEST_GENE"] = row["BEST_GENE"]
        records.append(record)

    return pd.DataFrame(records)


def _combine_tables(paths: list[Path]) -> pd.DataFrame:
    frames = [read_table(p) for p in paths if p.exists() and p.stat().st_size > 0]
    frames = [x for x in frames if not x.empty]
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def aggregate_study_outputs(
    root: Path,
    phenotype: str,
    ancestry_label: str,
    accession: str,
    args,
) -> None:
    paths = project_paths(root, phenotype, ancestry_label)
    study_out = paths["OUT"] / accession
    locus_root = study_out / "loci"
    study_out.mkdir(parents=True, exist_ok=True)

    with exclusive_file_lock(study_out / ".aggregate"):
        candidate_files = sorted(locus_root.glob("L*/locus_candidates.tsv"))
        result_files = sorted(locus_root.glob("L*/locus_formal_coloc.tsv"))
        audit_files = sorted(locus_root.glob("L*/locus_resource_audit.tsv"))
        resource_files = sorted(locus_root.glob("L*/resource_status.tsv"))

        candidates = _combine_tables(candidate_files)
        summary = _combine_tables(result_files)
        audit = _combine_tables(audit_files)

        if resource_files:
            resource_status = read_table(resource_files[0])
        else:
            resource_status = pd.DataFrame()

        candidates.to_csv(
            study_out / f"{accession}_molecular_qtl_candidates.tsv",
            sep="\t", index=False
        )
        summary.to_csv(
            study_out / f"{accession}_formal_coloc_summary.tsv",
            sep="\t", index=False
        )
        audit.to_csv(study_out / "resource_scan_audit.tsv", sep="\t", index=False)
        resource_status.to_csv(study_out / "resource_status.tsv", sep="\t", index=False)

        manifest = read_table(paths["OUT"] / "coloc_manifest.tsv")
        study_manifest = (
            manifest[manifest["STUDY_ACCESSION"].astype(str) == accession]
            if not manifest.empty else pd.DataFrame()
        )
        loci = sorted(study_manifest["LOCUS_ID"].astype(str).unique()) if not study_manifest.empty else []

        matrix = build_locus_qtl_matrix(
            summary, loci, resource_status, audit, args.suggestive_pp4
        ) if not resource_status.empty else pd.DataFrame()
        matrix.to_csv(
            study_out / f"{accession}_locus_qtl_matrix.tsv",
            sep="\t", index=False
        )
        wide = build_wide_evidence_table(matrix)
        wide.to_csv(
            study_out / f"{accession}_locus_evidence_table.tsv",
            sep="\t", index=False
        )

        done_loci = 0
        for locus in loci:
            meta = safe_json(locus_root / locus / "locus_summary.json")
            if meta.get("STATUS") == "COMPLETE":
                done_loci += 1

        status = "COMPLETE" if loci and done_loci == len(loci) else "PARTIAL"
        write_json(
            study_out / "coloc_summary.json",
            {
                "STEP09_VERSION": STEP09_VERSION,
                "METHOD": "local dense QTL + Step06 GWAS SuSiE fit + QTL runsusie + coloc.susie",
                "STUDY_ACCESSION": accession,
                "PHENOTYPE": phenotype,
                "N_LOCI_EXPECTED": len(loci),
                "N_LOCI_COMPLETE": done_loci,
                "N_TESTS": int((summary["COLOC_TESTED"] == "YES").sum()) if not summary.empty else 0,
                "N_STRONG": int((summary["COLOCALIZED_STRONG"] == "YES").sum()) if not summary.empty else 0,
                "N_PRIOR_ROBUST_STRONG": int((summary["PRIOR_ROBUST_STRONG"] == "YES").sum()) if not summary.empty and "PRIOR_ROBUST_STRONG" in summary.columns else 0,
                "STATUS": status,
                "UPDATED_UTC": utcnow(),
            },
        )


def process_locus_task(root: Path, manifest_row: pd.Series, args) -> None:
    phenotype = str(manifest_row["PHENOTYPE"])
    ancestry_code = str(manifest_row["ANCESTRY_CODE"])
    ancestry_label = str(manifest_row["ANCESTRY_LABEL"])
    accession = str(manifest_row["STUDY_ACCESSION"])
    locus_id = str(manifest_row["LOCUS_ID"])

    paths = project_paths(root, phenotype, ancestry_label)
    study_dir = paths["S06"] / accession
    study_out = paths["OUT"] / accession
    locus_out = study_out / "loci" / locus_id
    locus_out.mkdir(parents=True, exist_ok=True)

    locus_meta = locus_out / "locus_summary.json"
    if not args.force and safe_json(locus_meta).get("STATUS") == "COMPLETE":
        banner(f"STEP09 LOCUS ALREADY COMPLETE — {accession} {locus_id}")
        aggregate_study_outputs(root, phenotype, ancestry_label, accession, args)
        return

    allowed_types = {v.strip().lower() for v in args.qtl_types.split(",") if v.strip()}
    plink2 = find_executable(root, args.plink2, "PLINK2", "plink2")
    tabix = find_executable(root, args.tabix, "TABIX", "tabix")
    rscript = find_executable(root, args.rscript, "RSCRIPT", "Rscript")
    verify_r_packages(rscript)

    # Local files only (Step00_Download_Resources.py); nothing is downloaded.
    catalogue = load_eqtl_catalogue(root, False)
    contexts = select_catalogue_contexts(catalogue, args, allowed_types) if not catalogue.empty else catalogue
    contexts = resolve_catalogue_local_sources(root, contexts)
    if not contexts.empty:
        check_local_tabix_connectivity(tabix, contexts)

    custom_registry = ensure_custom_registry_files(
        root, load_custom_registry(root, args.registry), False
    )
    gwas_meta = load_gwas_metadata(paths["S06"], accession)
    gwas = load_gwas_locus(study_dir, locus_id)

    banner(f"STEP09 LOCUS TASK — {accession} | {locus_id}")
    print(f"Phenotype      : {phenotype}")
    print(f"Ancestry       : {ancestry_label} ({ancestry_code})")
    print(f"GWAS           : {gwas_meta['PHENOTYPE_TYPE']}, N {gwas_meta['COLOC_GWAS_N']} [{gwas_meta['GWAS_N_SOURCE']}]")
    print(f"QTL datasets   : {len(contexts)} eQTL Catalogue + "
          f"{int(custom_registry['ENABLED_BOOL'].sum()) if not custom_registry.empty else 0} custom")

    resource_status = resource_status_table(allowed_types, contexts, custom_registry)
    resource_status.to_csv(locus_out / "resource_status.tsv", sep="\t", index=False)

    candidates_frames = []
    audits = []
    if not contexts.empty:
        cand, audit = discover_catalogue_candidates_for_locus(
            tabix, contexts, locus_id, gwas, args
        )
        if not cand.empty:
            candidates_frames.append(cand)
        audits.extend(audit)
    if not custom_registry.empty:
        cand, audit = discover_custom_candidates_for_locus(
            tabix, custom_registry, locus_id, gwas, args, allowed_types
        )
        if not cand.empty:
            candidates_frames.append(cand)
        audits.extend(audit)

    candidates = (
        pd.concat(candidates_frames, ignore_index=True, sort=False)
        if candidates_frames else pd.DataFrame(columns=[
            "SOURCE","LOCUS_ID","DATASET","DATASET_ID","QTL_TYPE","TISSUE",
            "CONDITION","MOLECULAR_TRAIT_ID","GENE_ID","DISCOVERY_MIN_P",
            "N_DISCOVERY_ROWS","QTL_N","SOURCE_URL","ACCESS_PATH","ACCESS_MODE",
        ])
    )
    if not candidates.empty:
        candidates = (
            candidates.drop_duplicates(
                ["LOCUS_ID","QTL_TYPE","DATASET_ID","TISSUE","CONDITION","MOLECULAR_TRAIT_ID"],
                keep="first",
            )
            .sort_values(["DATASET_ID","QTL_TYPE","DISCOVERY_MIN_P"], kind="stable")
            .reset_index(drop=True)
        )

    candidates.to_csv(locus_out / "locus_candidates.tsv", sep="\t", index=False)
    audit_df = pd.DataFrame(audits)
    audit_df.to_csv(locus_out / "locus_resource_audit.tsv", sep="\t", index=False)
    print_candidates(candidates)

    r_worker = write_r_worker(locus_out)
    locus_resources = {"gwas": gwas}
    results = []

    for number, (_, candidate) in enumerate(candidates.iterrows(), start=1):
        banner(
            f"FORMAL COLOC {number}/{len(candidates)} — {locus_id} | "
            f"{candidate['QTL_TYPE']} | {candidate['TISSUE']} | {candidate['MOLECULAR_TRAIT_ID']}"
        )
        result = process_candidate(
            root, accession, study_dir, ancestry_code, gwas_meta, candidate,
            contexts, custom_registry, plink2, tabix, rscript,
            r_worker, locus_out, locus_resources, args,
        )
        results.append(result)
        print(
            f"RESULT tested={result['COLOC_TESTED']} status={result['STATUS']} "
            f"method={result.get('COLOC_METHOD','')} PP.H4={result['MAX_PP_H4']} "
            f"pair={result.get('BEST_SIGNAL_PAIR','')}"
        )

    summary = pd.DataFrame(results)
    summary.to_csv(locus_out / "locus_formal_coloc.tsv", sep="\t", index=False)
    print_formal_summary(summary)

    write_json(
        locus_meta,
        {
            "STEP09_VERSION": STEP09_VERSION,
            "STUDY_ACCESSION": accession,
            "LOCUS_ID": locus_id,
            "N_CANDIDATES": len(candidates),
            "N_TESTED": int((summary["COLOC_TESTED"] == "YES").sum()) if not summary.empty else 0,
            "N_STRONG": int((summary["COLOCALIZED_STRONG"] == "YES").sum()) if not summary.empty else 0,
            "N_RESOURCE_ERRORS": int((audit_df["STATUS"] == "ERROR").sum()) if not audit_df.empty else 0,
            "STATUS": (
                "COMPLETE"
                if not summary.empty
                and int((summary["COLOC_TESTED"] == "YES").sum()) > 0
                and int((summary["STATUS"] == "FAILED").sum()) == 0
                else "PARTIAL_FORMAL_ERROR"
                if not summary.empty
                and int((summary["COLOC_TESTED"] == "YES").sum()) > 0
                and int((summary["STATUS"] == "FAILED").sum()) > 0
                else "NOT_TESTED_FORMAL_ERROR"
                if not summary.empty
                and int((summary["STATUS"] == "FAILED").sum()) > 0
                else "NOT_TESTED_NO_FORMAL_RESULT"
            ),
            "COMPLETED_UTC": utcnow(),
        },
    )

    aggregate_study_outputs(root, phenotype, ancestry_label, accession, args)


def study_for_index(root: Path, phenotype: str, ancestry_label: str, index: int) -> str:
    """Accession of study N in the Step02 download manifest (the ARRAY_ID the
    Step02 planner prints) - the same N in Steps 02, 05, 06, 07, 08 and 09."""
    path = project_paths(root, phenotype, ancestry_label)["S02"] / "GWAS_download_manifest.tsv"
    if not path.is_file():
        raise SystemExit(f"No Step02 download manifest: {path}\nRun the Step02 planner first.")
    manifest = pd.read_csv(path, sep="\t", dtype=str)
    ids = pd.to_numeric(manifest["ARRAY_ID"], errors="coerce")
    selected = manifest[ids == index]
    if selected.empty:
        listing = "\n".join(f"  {int(i):>3}  {a}" for i, a in zip(ids, manifest["STUDY_ACCESSION"]) if pd.notna(i))
        raise SystemExit(f"--index {index} is not in {path}\nAvailable:\n{listing}")
    return str(selected.iloc[0]["STUDY_ACCESSION"]).strip()


def available_qtl_resources(root: Path, args) -> tuple[pd.DataFrame, pd.DataFrame]:
    allowed = {v.strip().lower() for v in args.qtl_types.split(",") if v.strip()}
    catalogue = load_eqtl_catalogue(root, False)
    contexts = select_catalogue_contexts(catalogue, args, allowed) if not catalogue.empty else catalogue
    contexts = resolve_catalogue_local_sources(root, contexts)
    registry = ensure_custom_registry_files(root, load_custom_registry(root, args.registry), False)
    return contexts, registry


def passthrough_options(args) -> list[str]:
    options = ["--registry", args.registry, "--qtl-types", args.qtl_types,
               "--candidate-p", str(args.candidate_p), "--min-common-snps", str(args.min_common_snps),
               "--max-common-snps", str(args.max_common_snps), "--strong-pp4", str(args.strong_pp4),
               "--suggestive-pp4", str(args.suggestive_pp4),
               "--max-candidates-per-context-locus", str(args.max_candidates_per_context_locus),
               "--ld-maf", str(args.ld_maf), "--ld-mismatch-s", str(args.ld_mismatch_s)]
    if args.catalog_all_studies:
        options.append("--catalog-all-studies")
    if args.catalog_study_regex:
        options += ["--catalog-study-regex", args.catalog_study_regex]
    if args.force:
        options.append("--force")
    return options


def planner_mode(root: Path, args) -> None:
    phenotype = args.phenotype.strip()
    ancestry_code, ancestry_label = canonical_ancestry(args.ancestry)
    paths = project_paths(root, phenotype, ancestry_label)
    studies, excluded = discover_step06_studies(paths["S06"])
    if args.study:
        studies = [x for x in studies if x == args.study]
    if not studies:
        raise SystemExit(f"No fine-mapped studies under {paths['S06']} - run Step06 first.")

    contexts, registry = available_qtl_resources(root, args)
    n_custom = int(registry["ENABLED_BOOL"].sum()) if not registry.empty else 0
    if contexts.empty and n_custom == 0:
        raise SystemExit(
            "No molecular-QTL data available locally. Download the eQTL Catalogue datasets with\n"
            "  envs/pipeline/bin/python Step00_Download_Resources.py --groups eqtl_catalogue qtl_catalogue "
            "--output Step00_Download_QTL.sh\n  sbatch Step00_Download_QTL.sh\n"
            f"or enable datasets in {args.registry}"
        )

    rows = []
    for study_index, accession in enumerate(studies, start=1):
        loci = list_loci(paths["S06"] / accession)
        if not loci:
            excluded.append({"STUDY_ACCESSION": accession, "STEP06_STATUS": "", "REASON": "no locus with a SuSiE fit"})
        for locus_id in loci:
            rows.append({"STUDY_INDEX": study_index, "STUDY_ACCESSION": accession, "LOCUS_ID": locus_id,
                         "PHENOTYPE": phenotype, "ANCESTRY_CODE": ancestry_code, "ANCESTRY_LABEL": ancestry_label})
    if not rows:
        raise SystemExit("No fine-mapped loci with a SuSiE fit (see Step06 --status).")

    n_tasks = min(len(rows), MAX_ARRAY_TASKS)
    manifest = pd.DataFrame(rows)
    manifest.insert(0, "COLOC_TASK_ID", [i % n_tasks + 1 for i in range(len(rows))])
    paths["OUT"].mkdir(parents=True, exist_ok=True)
    paths["LOG"].mkdir(parents=True, exist_ok=True)
    manifest_file = paths["OUT"] / "coloc_manifest.tsv"
    manifest.to_csv(manifest_file, sep="\t", index=False)
    pd.DataFrame(excluded, columns=["STUDY_ACCESSION", "STEP06_STATUS", "REASON"]).to_csv(
        paths["OUT"] / "coloc_excluded_studies.tsv", sep="\t", index=False)

    array = (f"1-{n_tasks}" if n_tasks > 1 else "1") + (f"%{args.max_parallel}" if args.max_parallel else "")
    directives = {**args.slurm_values, "job-name": f"COLOC_{slugify(phenotype)}_{ancestry_code.lower()}"[:100],
                  "array": array, "output": f"{paths['LOG']}/coloc.%A_%a.out",
                  "error": f"{paths['LOG']}/coloc.%A_%a.err"}
    order = ["job-name", "account", "partition", "qos", "nodes", "ntasks", "cpus-per-task", "mem", "time",
             "array", "output", "error"]
    header = "\n".join(f"#SBATCH --{k}={directives[k]}" for k in order if directives.get(k))
    header += "".join(f"\n#SBATCH --{k}={v}" for k, v in directives.items() if k not in order and v)
    python = root / "envs" / "pipeline" / "bin" / "python"
    if not python.exists():
        python = Path(sys.executable).resolve()
    command = [str(python), str(Path(__file__).resolve()), "--phenotype", phenotype, "--ancestry", ancestry_code,
               "--worker-manifest", str(manifest_file)] + passthrough_options(args)
    bash_file = root / f"Step09_Coloc_{slugify(phenotype)}_{slugify(ancestry_label)}.sh"
    bash_file.write_text(f"""#!/bin/bash --login
{header}

set -euo pipefail
cd {shlex.quote(str(root))}
export PATH={shlex.quote(str(root / 'envs' / 'pipeline' / 'bin'))}:"$PATH"
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-{args.cpus}}}"

{' '.join(shlex.quote(x) for x in command)}
""", encoding="utf-8", newline="\n")
    bash_file.chmod(0o755)

    banner("STEP 09 COLOCALISATION PLANNER - ONE ARRAY TASK PER LOCUS")
    print(f"Phenotype       : {phenotype}")
    print(f"Ancestry        : {ancestry_label} ({ancestry_code})")
    print(f"Studies         : {len(studies)}  ({len(excluded)} excluded, see coloc_excluded_studies.tsv)")
    print(f"Locus tasks     : {len(manifest)}  -> array {array}"
          + (f"  ({len(manifest) / n_tasks:.1f} loci per task)" if len(manifest) > n_tasks else ""))
    print(f"QTL datasets    : {len(contexts)} eQTL Catalogue (downloaded) + {n_custom} custom")
    print()
    per_study = manifest.groupby("STUDY_ACCESSION")["COLOC_TASK_ID"].agg(["min", "max", "size"])
    for accession, r in per_study.iterrows():
        print(f"  {accession:<15} {int(r['size']):>4} loci   tasks {int(r['min'])}-{int(r['max'])}")
    print(f"\nManifest: {manifest_file}\n\nRUN ALL LOCI:\n  sbatch {bash_file.name}")
    if len(manifest) <= n_tasks and len(per_study) > 1:
        accession, r = next(iter(per_study.iterrows()))
        print(f"ONLY {accession}:\n  sbatch --array={int(r['min'])}-{int(r['max'])} {bash_file.name}")
    base = f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {shlex.quote(phenotype)} --ancestry {ancestry_code}"
    print(f"THEN:\n  {base} --status\n  {base} --stats")


def run_locus(root: Path, row: pd.Series, args) -> bool:
    try:
        process_locus_task(root, row, args)
        return True
    except Exception as exc:
        paths = project_paths(root, str(row["PHENOTYPE"]), str(row["ANCESTRY_LABEL"]))
        out = paths["OUT"] / str(row["STUDY_ACCESSION"]) / "loci" / str(row["LOCUS_ID"])
        out.mkdir(parents=True, exist_ok=True)
        (out / "LOCUS_FAILED.txt").write_text(f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}",
                                              encoding="utf-8")
        print(f"[FAILED] {row['STUDY_ACCESSION']} {row['LOCUS_ID']}: {exc}", file=sys.stderr)
        return False


def worker_mode(root: Path, args) -> None:
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    args.cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", args.cpus))
    manifest = pd.read_csv(args.worker_manifest, sep="\t", dtype=str)
    selected = manifest[pd.to_numeric(manifest["COLOC_TASK_ID"], errors="coerce") == task]
    if selected.empty:
        raise SystemExit(f"No manifest rows for array task {task}")
    failed = [f"{r['STUDY_ACCESSION']} {r['LOCUS_ID']}" for _, r in selected.iterrows() if not run_locus(root, r, args)]
    if failed:
        raise SystemExit("Colocalisation failed for: " + ", ".join(failed))


def index_mode(root: Path, args) -> None:
    phenotype = args.phenotype.strip()
    ancestry_code, ancestry_label = canonical_ancestry(args.ancestry)
    paths = project_paths(root, phenotype, ancestry_label)
    accession = study_for_index(root, phenotype, ancestry_label, args.index)
    loci = list_loci(paths["S06"] / accession)
    if not loci:
        raise SystemExit(
            f"Study {args.index} ({accession}) has no fine-mapped locus with a SuSiE fit. Run:\n"
            f"  envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype {shlex.quote(phenotype)} "
            f"--ancestry {ancestry_code} --index {args.index}")
    if args.locus:
        if args.locus not in loci:
            raise SystemExit(f"{args.locus} is not a fine-mapped locus of {accession}: {', '.join(loci)}")
        loci = [args.locus]
    args.cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", args.cpus))
    failed = []
    for locus_id in loci:
        row = pd.Series({"STUDY_ACCESSION": accession, "LOCUS_ID": locus_id, "PHENOTYPE": phenotype,
                         "ANCESTRY_CODE": ancestry_code, "ANCESTRY_LABEL": ancestry_label})
        if not run_locus(root, row, args):
            failed.append(locus_id)
    if failed:
        raise SystemExit(f"{accession}: colocalisation failed for {', '.join(failed)}")


def check_resources_mode(root: Path, args) -> None:
    banner("STEP 09 - MOLECULAR-QTL RESOURCES (local)")
    contexts, registry = available_qtl_resources(root, args)
    if not contexts.empty:
        counts = contexts.groupby(["study_label", "QTL_TYPE"]).size().reset_index(name="DATASETS")
        print(counts.to_string(index=False))
    n_custom = int(registry["ENABLED_BOOL"].sum()) if not registry.empty else 0
    print(f"\nCustom registry ({args.registry}): {n_custom} enabled datasets")
    _, ancestry_label = canonical_ancestry(args.ancestry)
    code = canonical_ancestry(args.ancestry)[0]
    panels = all(Path(str(reference_prefix(root, code, c)) + ".pgen").exists() for c in range(1, 23))
    print(f"1000G {code} LD panels (Step04): {'OK' if panels else 'MISSING - run Step04'}")
    try:
        rscript = find_executable(root, args.rscript, "RSCRIPT", "Rscript")
        print("R coloc / susieR: " + ", ".join(f"{k} {v}" for k, v in verify_r_packages(rscript).items()))
    except Exception as exc:
        print(f"R coloc / susieR: MISSING ({exc})")


def status_mode(root: Path, args) -> None:
    phenotype = args.phenotype.strip()
    _, ancestry_label = canonical_ancestry(args.ancestry)
    paths = project_paths(root, phenotype, ancestry_label)
    studies, _ = discover_step06_studies(paths["S06"])
    banner(f"STEP 09 STATUS - {phenotype} / {ancestry_label}")
    print(f"{'STUDY':<15}{'LOCI':>6}{'DONE':>6}{'FAILED':>8}{'PENDING':>9}{'TESTS':>7}{'STRONG':>8}")
    totals = {"DONE": 0, "FAILED": 0, "PENDING": 0}
    for accession in studies:
        loci = list_loci(paths["S06"] / accession)
        done = failed = 0
        for locus in loci:
            d = paths["OUT"] / accession / "loci" / locus
            if (d / "locus_summary.json").exists():
                done += 1
            elif (d / "LOCUS_FAILED.txt").exists():
                failed += 1
        summary = safe_json(paths["OUT"] / accession / "coloc_summary.json")
        pending = len(loci) - done - failed
        totals["DONE"] += done
        totals["FAILED"] += failed
        totals["PENDING"] += pending
        print(f"{accession:<15}{len(loci):>6}{done:>6}{failed:>8}{pending:>9}"
              f"{summary.get('N_TESTS', 0):>7}{summary.get('N_STRONG', 0):>8}")
    print("\nLoci    : " + "  ".join(f"{k}={v}" for k, v in totals.items()))


def stats_mode(root: Path, args) -> None:
    phenotype = args.phenotype.strip()
    _, ancestry_label = canonical_ancestry(args.ancestry)
    paths = project_paths(root, phenotype, ancestry_label)
    banner(f"STEP 09 COLOCALISATION RESULTS - {phenotype} / {ancestry_label}")
    rows = []
    for d in sorted(p for p in paths["OUT"].iterdir() if p.is_dir()) if paths["OUT"].exists() else []:
        summary = read_table(d / f"{d.name}_formal_coloc_summary.tsv")
        if summary.empty:
            continue
        tested = summary[summary["COLOC_TESTED"].astype(str) == "YES"]
        strong = tested[tested["COLOCALIZED_STRONG"].astype(str) == "YES"]
        robust = tested[tested.get("PRIOR_ROBUST_STRONG", pd.Series(dtype=str)).astype(str) == "YES"]
        rows.append({"STUDY_ACCESSION": d.name, "N_LOCI": summary["LOCUS_ID"].nunique(),
                     "N_CANDIDATE_TRAITS": len(summary), "N_TESTED": len(tested), "N_STRONG": len(strong),
                     "N_PRIOR_ROBUST_STRONG": len(robust),
                     "N_LOCI_WITH_STRONG": strong["LOCUS_ID"].nunique() if len(strong) else 0,
                     "N_LD_MISMATCH_WARNINGS": int((tested.get("LD_MISMATCH_WARNING", pd.Series(dtype=str))
                                                    .astype(str) == "YES").sum())})
        print(f"\n{d.name}: {len(tested)} tests in {summary['LOCUS_ID'].nunique()} loci; "
              f"{len(strong)} strong (PP.H4>={args.strong_pp4}), {len(robust)} robust to the low p12 prior")
        top = tested.assign(_H4=numeric(tested["MAX_PP_H4"])).sort_values("_H4", ascending=False).head(10)
        for _, r in top.iterrows():
            print(f"   {r['LOCUS_ID']} {r['QTL_TYPE']:<7} {str(r['TISSUE'])[:30]:<30} {r['GENE_ID'] or r['MOLECULAR_TRAIT_ID']}"
                  f"  PP.H4 {r['_H4']:.3f}  {r['COLOC_CLASS']}"
                  + ("  [LD mismatch]" if str(r.get("LD_MISMATCH_WARNING", "")) == "YES" else ""))
    if not rows:
        print("No colocalisation results yet (see --status).")
        return
    output = paths["OUT"] / "coloc_stats.tsv"
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
    print(f"\nTable: {output}")


def main() -> None:
    args = arguments()
    root = PROJECT_ROOT
    write_custom_registry_template(root)
    if args.check_resources:
        check_resources_mode(root, args)
    elif args.status:
        status_mode(root, args)
    elif args.stats:
        stats_mode(root, args)
    elif args.worker_manifest:
        worker_mode(root, args)
    elif args.index is not None:
        index_mode(root, args)
    else:
        planner_mode(root, args)


if __name__ == "__main__":
    main()
