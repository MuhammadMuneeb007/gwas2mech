#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 06 - Regional extraction + ancestry-matched LD + SuSiE-RSS
===============================================================================

    # plan + SLURM array (one task per study; >1000 studies share tasks)
    envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR
    sbatch Step06_FineMap_GWAS_migraine_european.sh

    # progress (read-only), per-study statistics, plots
    envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR --status
    envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR --stats
    envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR --plots

    # study 1 (the number the Step02 planner prints) on this machine; only its
    # Step05 clumping has to be done
    envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype migraine --ancestry EUR --index 1

Inputs
    Step02  02_summary_stats/<p>/<a>/qc/<GCST>/<GCST>_GRCh38_QC.tsv.gz
            (CHR POS SNP EA NEA BETA SE P EAF MAF N INFO Z)
            01_gwas_catalog/<p>/<a>/GWAS_study_metadata.tsv
            (PHENOTYPE_TYPE N_CASES N_CONTROLS N_TOTAL N_EFF NEEDS_REVIEW ...)
    Step04  resources/1000G/<ANC>/chr<N>_<ANC>_GRCh38.pgen/.pvar/.psam
    Step05  05_ld_clumping/<p>/<a>/<GCST>/lead_variants.tsv

Sample size (N) given to SuSiE-RSS - from Step02, never re-guessed here:
    Binary      N_EFF = 4*cases*controls/(cases+controls)   (default; --susie-n total = N_TOTAL)
    Continuous  N_TOTAL
    Unknown or NEEDS_REVIEW studies are NOT fine-mapped (fill VERIFIED_PHENOTYPE_TYPE in
    Step02's phenotype_type_review.tsv, or --include-unreviewed). The N used and why
    is recorded per study (GWAS_N, GWAS_N_SOURCE).

Per study / locus
  1. lead variants +- --locus-window-kb (1000 kb) merged into loci,
  2. all GWAS variants of each locus extracted (not only significant ones),
  3. matched to the 1000G panel by rsID or position AND both alleles; Z is
     oriented to the 1000G REF allele (EA=ALT -> -Z),
  4. signed LD with plink2 --r-unphased ref-based on exactly those variants,
     in the same order (checked in R),
  5. susieR::susie_rss(z, R, n), estimate_residual_variance = FALSE (reference LD),
  6. diagnostics: estimate_s_rss (LD / z consistency), kriging_rss (possible
     allele flips), convergence, credible-set purity, and an LD-independent
     L = 1 fit (top variant agreement). Each locus gets LOCUS_QC = PASS / WARN
     with the reasons.

Outputs
    06_finemapping/<p>/<a>/
        finemapping_manifest.tsv, finemapping_excluded_studies.tsv, finemapping_stats.tsv
        plots/<GCST>/<L>_regional.png, <L>_ld.png, <GCST>_pip_genome.png   (--plots)
        <GCST>/loci_definition.tsv, regional_gwas/<L>.tsv.gz,
               loci/<L>/{matched_variants.tsv, <ANC>_LD.unphased.vcor1.bin(.vars),
                         susie/{susie_variants.tsv, credible_sets.tsv, susie_summary.tsv,
                                susie_diagnostics.tsv, kriging.tsv, susie_fit.rds}}
               <GCST>_finemapped_variants.tsv.gz, <GCST>_95pct_credible_sets.tsv,
               <GCST>_finemapping_summary.tsv, <GCST>_failed_loci.tsv

SLURM settings come from the #SBATCH header of HPC_CPU.sh; --partition --qos
--account --cpus --mem --time --max-parallel override them.
===============================================================================
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

PROJECT_ROOT = Path.cwd().resolve()
MAX_ARRAY_TASKS = 1000
DEFAULT_CPUS = 8
DEFAULT_LOCUS_WINDOW_KB = 1000
DEFAULT_REFERENCE_MAF = 0.01
DEFAULT_REFERENCE_GENO = 0.05
DEFAULT_MIN_VARIANTS = 10
DEFAULT_MAX_LD_VARIANTS = 15000
DEFAULT_COVERAGE = 0.95
DEFAULT_SUSIE_L = 10
DEFAULT_MAX_SUSIE_L = 20
DEFAULT_LD_RIDGE = 1e-4
DEFAULT_CHUNK_SIZE = 500_000
DEFAULT_ESTIMATE_S_MAX = DEFAULT_MAX_LD_VARIANTS
DEFAULT_MAX_ITER = 500

# Bumped whenever the R helper changes, so finished loci are re-run.
R_HELPER_VERSION = "2"

# Heuristic LD-mismatch flags on susieR::estimate_s_rss (0 = z and LD fully
# consistent). Not formal tests: they mark loci for inspection.
S_RSS_WARN = 0.1
S_RSS_HIGH = 0.3

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


# =============================================================================
# GENERIC HELPERS
# =============================================================================


def banner(text: str) -> None:
    print()
    print("=" * 88)
    print(text)
    print("=" * 88)


def normalize_text(value: Any) -> str:
    value = "" if value is None else str(value)
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def slugify(value: Any) -> str:
    return normalize_text(value).replace(" ", "_")


def resolve_ancestry(value: str) -> tuple[str, str]:
    key = normalize_text(value)
    if key not in ANCESTRY_ALIASES:
        allowed = ", ".join(sorted({x[0] for x in ANCESTRY_ALIASES.values()}))
        raise ValueError(f"Unsupported ancestry '{value}'. Supported codes: {allowed}")
    return ANCESTRY_ALIASES[key]


def clean_chr(value: Any) -> int | None:
    text = str(value).strip().lower().replace("chr", "")
    try:
        chrom = int(float(text))
    except Exception:
        return None
    return chrom if 1 <= chrom <= 22 else None


def clean_allele(value: Any) -> str | None:
    if pd.isna(value):
        return None
    value = str(value).strip().upper()
    if not value or value in {"NAN", "NA", ".", "<NA>"}:
        return None
    return value


def clean_rsid(value: Any) -> str | None:
    if pd.isna(value):
        return None
    value = str(value).strip()
    return value if re.fullmatch(r"rs\d+", value, flags=re.I) else None


def allele_key(a: Any, b: Any) -> tuple[str, str] | None:
    a = clean_allele(a)
    b = clean_allele(b)
    if a is None or b is None:
        return None
    return tuple(sorted((a, b)))


def find_column(columns, candidates):
    lookup = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def safe_int(value: Any) -> int | None:
    if value is None or pd.isna(value):
        return None
    if isinstance(value, (int, np.integer)):
        return int(value) if int(value) > 0 else None
    if isinstance(value, float):
        return int(round(value)) if np.isfinite(value) and value > 0 else None
    text = str(value).strip().replace(",", "")
    match = re.search(r"\d+", text)
    if not match:
        return None
    out = int(match.group())
    return out if out > 0 else None


def parse_numbers(text: Any) -> list[int]:
    if text is None or pd.isna(text):
        return []
    values = []
    for token in re.findall(r"\d[\d,]*", str(text)):
        try:
            number = int(token.replace(",", ""))
            if number >= 100:
                values.append(number)
        except ValueError:
            pass
    return values


def effective_n(n_cases: int | None, n_controls: int | None) -> int | None:
    if not n_cases or not n_controls or n_cases <= 0 or n_controls <= 0:
        return None
    return int(round(4.0 / ((1.0 / n_cases) + (1.0 / n_controls))))


def run_command(command, log_file: Path | None = None) -> str:
    print()
    print(" ".join(shlex.quote(str(x)) for x in command), flush=True)
    process = subprocess.Popen(
        [str(x) for x in command],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    captured = []
    handle = None
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handle = log_file.open("w", encoding="utf-8")
    try:
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            captured.append(line)
            if handle:
                handle.write(line)
                handle.flush()
    finally:
        if handle:
            handle.close()
    code = process.wait()
    output = "".join(captured)
    if code != 0:
        raise RuntimeError(
            f"Command failed with exit code {code}:\n"
            + " ".join(str(x) for x in command)
            + "\n\n"
            + output[-5000:]
        )
    return output


# =============================================================================
# SLURM SETTINGS (#SBATCH header of HPC_CPU.sh; command line overrides)
# =============================================================================

PER_JOB_KEYS = {"job-name", "output", "error", "array"}


def read_slurm_template(path):
    values = {}
    if not Path(path).is_file():
        return values
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        match = re.match(r"^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
        if match and (match.group(2) or "").strip():
            values[match.group(1)] = match.group(2).strip()
    return values


def validate_slurm_time(value):
    text = str(value).strip()
    if not re.fullmatch(r"(\d+-)?\d{1,3}(:\d{2}){0,2}", text):
        raise argparse.ArgumentTypeError(f"walltime {value!r} must look like 12:00:00 or 2-00:00:00")
    return text


def resolve_slurm(args):
    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": str(DEFAULT_CPUS),
              "mem": "32G", "time": "24:00:00"}
    values.update({k: v for k, v in read_slurm_template(PROJECT_ROOT / args.template).items()
                   if k not in PER_JOB_KEYS})
    for key, value in (("partition", args.partition), ("qos", args.qos), ("account", args.account),
                       ("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time)):
        if value not in (None, ""):
            values[key] = str(value)
    args.slurm_values = values
    args.cpus = int(values["cpus-per-task"])
    args.partition, args.time, args.mem = values.get("partition"), values.get("time"), values.get("mem")
    return args


def assign_task_ids(n_items):
    """At most MAX_ARRAY_TASKS array tasks; task k runs items k, k+T, k+2T, ..."""
    n_tasks = min(n_items, MAX_ARRAY_TASKS)
    return [i % n_tasks + 1 for i in range(n_items)], n_tasks


def sbatch_header_block(args, job_name, output, error, n_tasks):
    array = f"1-{n_tasks}" if n_tasks > 1 else "1"
    if args.max_parallel and int(args.max_parallel) > 0:
        array += f"%{int(args.max_parallel)}"
    directives = {**args.slurm_values, "job-name": str(job_name)[:100], "array": array,
                  "output": output, "error": error}
    order = ["job-name", "account", "partition", "qos", "nodes", "ntasks",
             "cpus-per-task", "mem", "time", "array", "output", "error"]
    lines = [f"#SBATCH --{k}={directives[k]}" for k in order if directives.get(k)]
    lines += [f"#SBATCH --{k}={v}" for k, v in directives.items() if k not in order and v]
    return "\n".join(lines)


# =============================================================================
# ARGUMENTS
# =============================================================================


def arguments(argv=None):
    parser = argparse.ArgumentParser(description="SuSiE-RSS fine-mapping of Step05 loci.")
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)

    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--index", type=int, default=None,
                       help="Study N of the Step02 download manifest (the number Step02 prints): "
                            "fine-map it on this machine; other studies need not be ready.")
    modes.add_argument("--status", action="store_true",
                       help="READ-ONLY: fine-mapping done / partial / failed / pending / excluded per study.")
    modes.add_argument("--stats", action="store_true",
                       help="Per-study loci, credible sets, PIPs, N used and locus QC (finemapping_stats.tsv).")
    modes.add_argument("--plots", action="store_true",
                       help="Regional association + PIP + genes, LD heatmaps, genome-wide PIP (plots/).")

    parser.add_argument("--susie-n", choices=["effective", "total"], default="effective",
                        help="N for BINARY studies: N_EFF (default) or N_TOTAL. Continuous studies use N_TOTAL.")
    parser.add_argument("--include-unreviewed", action="store_true",
                        help="Also fine-map studies with PHENOTYPE_TYPE=Unknown / NEEDS_REVIEW (not recommended).")
    parser.add_argument("--sample-size", type=int, default=None,
                        help="Last-resort N when Step02 has no sample size for a study.")
    parser.add_argument("--locus-window-kb", type=int, default=DEFAULT_LOCUS_WINDOW_KB)
    parser.add_argument("--reference-maf", type=float, default=DEFAULT_REFERENCE_MAF)
    parser.add_argument("--reference-geno", type=float, default=DEFAULT_REFERENCE_GENO)
    parser.add_argument("--min-variants", type=int, default=DEFAULT_MIN_VARIANTS)
    parser.add_argument("--max-ld-variants", type=int, default=DEFAULT_MAX_LD_VARIANTS)
    parser.add_argument("--coverage", type=float, default=DEFAULT_COVERAGE)
    parser.add_argument("--susie-L", type=int, default=DEFAULT_SUSIE_L, dest="susie_L")
    parser.add_argument("--max-susie-L", type=int, default=DEFAULT_MAX_SUSIE_L, dest="max_susie_L")
    parser.add_argument("--ld-ridge", type=float, default=DEFAULT_LD_RIDGE)
    parser.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER)
    parser.add_argument("--refine", action="store_true",
                        help="susie refine=TRUE (escapes local optima; slower).")
    parser.add_argument("--estimate-s-max", type=int, default=DEFAULT_ESTIMATE_S_MAX,
                        help="Largest locus (variants) for estimate_s_rss / kriging_rss diagnostics.")
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="Plan even when Step05 is unfinished for some studies.")

    slurm = parser.add_argument_group("SLURM (default: the #SBATCH header of HPC_CPU.sh)")
    slurm.add_argument("--template", default="HPC_CPU.sh")
    slurm.add_argument("--partition", default=None)
    slurm.add_argument("--qos", default=None)
    slurm.add_argument("--account", default=None)
    slurm.add_argument("--cpus", type=int, default=None)
    slurm.add_argument("--mem", "--memory", dest="mem", default=None)
    slurm.add_argument("--time", type=validate_slurm_time, default=None)
    slurm.add_argument("--max-parallel", type=int, default=0,
                       help="Maximum simultaneously running array tasks (0 = no throttle).")

    args = parser.parse_args(argv)
    resolve_slurm(args)

    if args.cpus < 1 or args.max_parallel < 0:
        parser.error("--cpus must be >= 1 and --max-parallel >= 0 (0 = no throttle)")
    if args.locus_window_kb < 1:
        parser.error("--locus-window-kb must be >= 1")
    if not (0 < args.reference_maf < 0.5):
        parser.error("--reference-maf must be between 0 and 0.5")
    if not (0 <= args.reference_geno < 1):
        parser.error("--reference-geno must be in [0,1)")
    if not (0 < args.coverage < 1):
        parser.error("--coverage must be between 0 and 1")
    if args.min_variants < 2 or args.max_ld_variants < args.min_variants:
        parser.error("Invalid fine-mapping variant limits")
    if args.susie_L < 1 or args.max_susie_L < args.susie_L:
        parser.error("Invalid SuSiE L settings")
    if args.index is not None and args.index < 1:
        parser.error("--index must be >= 1")
    return args


# =============================================================================
# STEP02 STUDY METADATA -> N FOR SUSIE-RSS
# =============================================================================


def study_metadata_file(root: Path, phenotype_slug: str, ancestry_slug: str) -> Path:
    return root / "01_gwas_catalog" / phenotype_slug / ancestry_slug / "GWAS_study_metadata.tsv"


def read_study_metadata(path: Path) -> dict[str, dict]:
    if not path.is_file():
        raise FileNotFoundError(
            f"Step02 study metadata not found:\n  {path}\n"
            "It is written by every Step02 run; refresh it with:\n"
            "  envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype ... --ancestry ... --status"
        )
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return {str(r["GCST_ID"]).strip(): r.to_dict() for _, r in table.iterrows()}


def choose_sample_size(meta: dict | None, args) -> tuple[int | None, str, str]:
    """(N, source, problem). problem non-empty -> the study is not fine-mapped."""
    if not meta:
        return None, "", "not in Step02 GWAS_study_metadata.tsv"

    phenotype_type = str(meta.get("PHENOTYPE_TYPE", "")).strip()
    review = str(meta.get("NEEDS_REVIEW", "")).strip().lower() == "true"

    if (phenotype_type not in ("Binary", "Continuous") or review) and not args.include_unreviewed:
        return None, "", (
            f"PHENOTYPE_TYPE={phenotype_type or 'Unknown'}"
            + (" NEEDS_REVIEW" if review else "")
            + " - verify it in Step02 phenotype_type_review.tsv (VERIFIED_PHENOTYPE_TYPE)"
        )

    n_total, n_eff = safe_int(meta.get("N_TOTAL")), safe_int(meta.get("N_EFF"))

    if phenotype_type == "Binary":
        if args.susie_n == "effective" and n_eff:
            return n_eff, "N_EFF = 4*cases*controls/(cases+controls) (Step02)", ""
        if n_total:
            why = "--susie-n total" if args.susie_n == "total" else "no case/control counts for N_EFF"
            return n_total, f"N_TOTAL (binary; {why})", ""
    elif n_total:
        return n_total, f"N_TOTAL ({phenotype_type.lower() or 'unknown type'})", ""

    if args.sample_size:
        return int(args.sample_size), "--sample-size", ""

    return None, "", "no sample size in Step02 metadata (use --sample-size)"


# =============================================================================
# REFERENCE HELPERS
# =============================================================================


def reference_prefix(root: Path, ancestry_code: str, chromosome: int) -> Path:
    return root / "resources" / "1000G" / ancestry_code / f"chr{chromosome}_{ancestry_code}_GRCh38"


def verify_reference(root: Path, ancestry_code: str) -> None:
    missing = []
    for chrom in range(1, 23):
        prefix = reference_prefix(root, ancestry_code, chrom)
        for ext in (".pgen", ".pvar", ".psam"):
            path = Path(str(prefix) + ext)
            if not path.exists() or path.stat().st_size == 0:
                missing.append(str(path))
    if missing:
        raise FileNotFoundError(
            f"Ancestry reference {ancestry_code} is incomplete. Missing {len(missing)} files.\n"
            + "\n".join(missing[:20])
            + f"\n\nBuild it with:\n  envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry {ancestry_code}"
        )


def count_reference_samples(root: Path, ancestry_code: str) -> int:
    psam = Path(str(reference_prefix(root, ancestry_code, 1)) + ".psam")
    count = 0
    with psam.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            count += 1
    if count <= 0:
        raise RuntimeError(f"Could not count samples in {psam}")
    return count


# =============================================================================
# TOOLS + ONE FINE-MAPPING MANIFEST ROW
# =============================================================================


def resolve_finemap_tools(root: Path) -> tuple[Path, Path, Path]:
    """Return (pipeline_python, plink2, rscript) after verifying susieR."""
    pipeline_python = root / "envs" / "pipeline" / "bin" / "python"
    plink2 = root / "envs" / "pipeline" / "bin" / "plink2"
    rscript = root / "envs" / "pipeline" / "bin" / "Rscript"
    if not rscript.exists():
        located = shutil.which("Rscript")
        rscript = Path(located).resolve() if located else rscript

    for label, path in [
        ("pipeline Python", pipeline_python),
        ("PLINK2", plink2),
        ("Rscript", rscript),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Required {label} not found: {path}\nRun: python Step00_Installation.py")

    susie_check = subprocess.run(
        [str(rscript), "-e", 'quit(status=ifelse(requireNamespace("susieR", quietly=TRUE),0,1))'],
        capture_output=True,
        text=True,
    )
    if susie_check.returncode != 0:
        raise RuntimeError("R package susieR is not available in the configured R environment")

    return pipeline_python, plink2, rscript


def build_finemap_row(
    args,
    accession: str,
    qc_file: Path,
    clump_root: Path,
    output_root: Path,
    meta: dict | None,
    reference_n: int,
    rscript: Path,
    plink2: Path,
    phenotype: str,
    ancestry_code: str,
    ancestry_label: str,
) -> tuple[dict | None, list[str], str]:
    """Return (manifest_row, missing_inputs, exclusion_reason)."""
    study_clump_dir = clump_root / accession
    lead_file = study_clump_dir / "lead_variants.tsv"
    clump_summary = study_clump_dir / "clumping_summary.tsv"

    missing = [
        str(p) for p in (lead_file, clump_summary, qc_file)
        if not p.exists() or p.stat().st_size == 0
    ]
    if missing:
        return None, missing, ""

    gwas_n, n_source, problem = choose_sample_size(meta, args)
    if problem:
        return None, [], problem

    meta = meta or {}
    return (
        {
            "FINEMAP_TASK_ID": 0,
            "STUDY_ACCESSION": accession,
            "PHENOTYPE": phenotype,
            "ANCESTRY": ancestry_label,
            "ANCESTRY_CODE": ancestry_code,
            "ANCESTRY_LABEL": ancestry_label,
            "QC_FILE": str(qc_file),
            "LEAD_VARIANTS_FILE": str(lead_file.resolve()),
            "CLUMPING_SUMMARY_FILE": str(clump_summary.resolve()),
            "OUTPUT_DIR": str((output_root / accession).resolve()),
            "PHENOTYPE_TYPE": meta.get("PHENOTYPE_TYPE", ""),
            "EFFECT_TYPE": meta.get("EFFECT_TYPE", ""),
            "N_CASES": meta.get("N_CASES", ""),
            "N_CONTROLS": meta.get("N_CONTROLS", ""),
            "N_TOTAL": meta.get("N_TOTAL", ""),
            "N_EFF": meta.get("N_EFF", ""),
            "GWAS_N": int(gwas_n),
            "GWAS_N_SOURCE": n_source,
            "REFERENCE_N": reference_n,
            "LOCUS_WINDOW_KB": args.locus_window_kb,
            "REFERENCE_MAF": args.reference_maf,
            "REFERENCE_GENO": args.reference_geno,
            "MIN_VARIANTS": args.min_variants,
            "MAX_LD_VARIANTS": args.max_ld_variants,
            "COVERAGE": args.coverage,
            "SUSIE_L": args.susie_L,
            "MAX_SUSIE_L": args.max_susie_L,
            "LD_RIDGE": args.ld_ridge,
            "MAX_ITER": args.max_iter,
            "REFINE": bool(args.refine),
            "ESTIMATE_S_MAX": args.estimate_s_max,
            "CHUNK_SIZE": args.chunk_size,
            "R_SCRIPT": str(rscript),
            "PLINK2": str(plink2),
        },
        [],
        "",
    )


# =============================================================================
# PLANNER
# =============================================================================


def planner_mode(args) -> None:
    root = PROJECT_ROOT
    phenotype = args.phenotype.strip()
    phenotype_slug = slugify(phenotype)
    ancestry_code, ancestry_label = resolve_ancestry(args.ancestry)
    ancestry_slug = slugify(ancestry_label)

    pipeline_python, plink2, rscript = resolve_finemap_tools(root)

    verify_reference(root, ancestry_code)
    reference_n = count_reference_samples(root, ancestry_code)

    clump_root = root / "05_ld_clumping" / phenotype_slug / ancestry_slug
    clump_manifest_file = clump_root / "clumping_manifest.tsv"
    if not clump_manifest_file.exists():
        raise FileNotFoundError(
            f"Step05 manifest not found:\n{clump_manifest_file}\nRun Step05 first."
        )
    clump_manifest = pd.read_csv(clump_manifest_file, sep="\t", dtype=str)
    study_meta = read_study_metadata(study_metadata_file(root, phenotype_slug, ancestry_slug))

    output_root = root / "06_finemapping" / phenotype_slug / ancestry_slug
    log_root = root / "logs" / "step06_finemapping" / phenotype_slug / ancestry_slug
    output_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    rows = []
    incomplete = []
    excluded = []

    banner("STEP 06 - FINEMAPPING PLANNER")
    print(f"Phenotype    : {phenotype}")
    print(f"Ancestry     : {ancestry_label} ({ancestry_code})")
    print(f"Reference N  : {reference_n}  (1000G {ancestry_code} LD panel)")
    print(f"Binary N     : {'N_EFF' if args.susie_n == 'effective' else 'N_TOTAL'} (--susie-n {args.susie_n})")
    print(f"Step05 rows  : {len(clump_manifest)}")
    print()

    for _, crow in clump_manifest.iterrows():
        accession = str(crow["STUDY_ACCESSION"]).strip()
        qc_file = Path(str(crow["QC_FILE"])).resolve()

        row, missing, problem = build_finemap_row(
            args,
            accession=accession,
            qc_file=qc_file,
            clump_root=clump_root,
            output_root=output_root,
            meta=study_meta.get(accession),
            reference_n=reference_n,
            rscript=rscript,
            plink2=plink2,
            phenotype=phenotype,
            ancestry_code=ancestry_code,
            ancestry_label=ancestry_label,
        )

        if missing:
            incomplete.append((accession, missing))
            print(f"INCOMPLETE {accession}: Step05 outputs missing")
            continue
        if problem:
            excluded.append({"STUDY_ACCESSION": accession, "REASON": problem})
            print(f"EXCLUDED   {accession}: {problem}")
            continue

        rows.append(row)
        print(f"READY      {accession}: {row['PHENOTYPE_TYPE']}  N={row['GWAS_N']:,}  [{row['GWAS_N_SOURCE']}]")

    pd.DataFrame(excluded, columns=["STUDY_ACCESSION", "REASON"]).to_csv(
        output_root / "finemapping_excluded_studies.tsv", sep="\t", index=False
    )

    if incomplete:
        report = output_root / "finemapping_incomplete_inputs.tsv"
        pd.DataFrame(
            [{"STUDY_ACCESSION": acc, "PROBLEM": " | ".join(problems)} for acc, problems in incomplete]
        ).to_csv(report, sep="\t", index=False)
        print()
        print(f"Incomplete studies: {len(incomplete)}  (report: {report})")
        if not args.allow_incomplete:
            raise RuntimeError(
                "Step05 is not complete for all expected studies. "
                "Finish Step05 and rerun. Use --allow-incomplete only for testing."
            )

    if not rows:
        raise RuntimeError("No studies are ready for fine-mapping (see finemapping_excluded_studies.tsv)")

    manifest = pd.DataFrame(rows)
    task_ids, n_tasks = assign_task_ids(len(manifest))
    manifest["FINEMAP_TASK_ID"] = task_ids
    manifest_file = (output_root / "finemapping_manifest.tsv").resolve()
    manifest.to_csv(manifest_file, sep="\t", index=False)

    bash_file = root / f"Step06_FineMap_GWAS_{phenotype_slug}_{ancestry_slug}.sh"
    job_name = f"FM_{phenotype_slug}_{ancestry_code.lower()}"[:100]
    script_path = Path(__file__).resolve()

    sbatch_header = sbatch_header_block(
        args,
        job_name=job_name,
        output=f"{log_root}/finemap.%A_%a.out",
        error=f"{log_root}/finemap.%A_%a.err",
        n_tasks=n_tasks,
    )

    bash_text = f'''#!/bin/bash --login
{sbatch_header}

set -euo pipefail
cd {shlex.quote(str(root))}
export PATH={shlex.quote(str(root / "envs" / "pipeline" / "bin"))}:"$PATH"
export GWAS_FINEMAP_MANIFEST={shlex.quote(str(manifest_file))}
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-{args.cpus}}}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"

printf '%s\\n' "STEP 06 - SUSIE-RSS FINEMAPPING | job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID on $(hostname)"

{shlex.quote(str(pipeline_python))} {shlex.quote(str(script_path))} --phenotype {shlex.quote(phenotype)} --ancestry {shlex.quote(ancestry_code)}
'''
    bash_file.write_text(bash_text, encoding="utf-8", newline="\n")
    bash_file.chmod(0o755)

    base = f"envs/pipeline/bin/python {script_path.name} --phenotype {shlex.quote(phenotype)} --ancestry {ancestry_code}"
    banner("PLANNING COMPLETE")
    print(f"Studies ready       : {len(manifest)}")
    print(f"Studies excluded    : {len(excluded)}")
    print(f"Studies incomplete  : {len(incomplete)}")
    print(f"Array               : 1-{n_tasks}")
    print(f"Partition / time    : {args.partition or 'scheduler default'} / {args.time}")
    print(f"Memory / CPUs       : {args.mem} / {args.cpus}")
    print(f"Fine-map manifest   : {manifest_file}")
    print(f"Generated batch     : {bash_file}")
    print()
    print("NEXT COMMAND:")
    print(f"  sbatch {bash_file.name}")
    print("THEN:")
    print(f"  {base} --status")
    print(f"  {base} --stats")
    print(f"  {base} --plots")


# =============================================================================
# LOCUS DEFINITION
# =============================================================================


def load_leads_and_build_loci(lead_file: Path, window_kb: int) -> pd.DataFrame:
    leads = pd.read_csv(lead_file, sep="\t", low_memory=False)
    if leads.empty:
        return pd.DataFrame(
            columns=["LOCUS_ID", "CHR", "LOCUS_START", "LOCUS_END", "N_INDEPENDENT_SIGNALS", "LEAD_IDS", "LEAD_POSITIONS", "MIN_LEAD_P"]
        )

    chr_col = find_column(leads.columns, ["CHROM", "#CHROM", "CHR", "SOURCE_CHROMOSOME"])
    pos_col = find_column(leads.columns, ["POS", "BP"])
    p_col = find_column(leads.columns, ["P", "P_VALUE"])
    id_col = find_column(leads.columns, ["ID", "SNP", "SNPID", "rsID", "RSID"])
    if chr_col is None or pos_col is None:
        raise RuntimeError(f"Lead file lacks chromosome/position columns: {lead_file}")

    work = pd.DataFrame()
    work["CHR"] = leads[chr_col].map(clean_chr)
    work["POS"] = pd.to_numeric(leads[pos_col], errors="coerce")
    work["P"] = pd.to_numeric(leads[p_col], errors="coerce") if p_col else np.nan
    work["ID"] = leads[id_col].astype(str) if id_col else ""
    work = work.dropna(subset=["CHR", "POS"]).copy()
    work["CHR"] = work["CHR"].astype(int)
    work["POS"] = work["POS"].astype(int)
    work = work.sort_values(["CHR", "POS", "P"], kind="stable")

    window = window_kb * 1000
    merged = []
    locus_counter = 0
    for chrom, group in work.groupby("CHR", sort=True):
        current = None
        for _, row in group.iterrows():
            start = max(1, int(row["POS"]) - window)
            end = int(row["POS"]) + window
            signal = {
                "id": str(row["ID"]),
                "pos": int(row["POS"]),
                "p": float(row["P"]) if pd.notna(row["P"]) else np.nan,
            }
            if current is None or start > current["end"]:
                if current is not None:
                    merged.append(current)
                current = {"chr": int(chrom), "start": start, "end": end, "signals": [signal]}
            else:
                current["end"] = max(current["end"], end)
                current["signals"].append(signal)
        if current is not None:
            merged.append(current)

    records = []
    for item in merged:
        locus_counter += 1
        pvals = [x["p"] for x in item["signals"] if np.isfinite(x["p"])]
        records.append(
            {
                "LOCUS_ID": f"L{locus_counter:03d}",
                "CHR": item["chr"],
                "LOCUS_START": item["start"],
                "LOCUS_END": item["end"],
                "N_INDEPENDENT_SIGNALS": len(item["signals"]),
                "LEAD_IDS": ";".join(x["id"] for x in item["signals"]),
                "LEAD_POSITIONS": ";".join(str(x["pos"]) for x in item["signals"]),
                "MIN_LEAD_P": min(pvals) if pvals else np.nan,
            }
        )
    return pd.DataFrame(records)


# =============================================================================
# REGIONAL GWAS EXTRACTION
# =============================================================================


def inspect_qc_columns(path: Path) -> dict[str, str | None]:
    header = pd.read_csv(path, sep="\t", compression="infer", nrows=0)
    cols = list(header.columns)
    mapping = {
        "CHR": find_column(cols, ["CHR"]),
        "POS": find_column(cols, ["POS"]),
        "P": find_column(cols, ["P"]),
        "EA": find_column(cols, ["EA"]),
        "NEA": find_column(cols, ["NEA"]),
        "BETA": find_column(cols, ["BETA"]),
        "SE": find_column(cols, ["SE"]),
        "Z": find_column(cols, ["Z"]),
        "OR": find_column(cols, ["OR"]),
        "OR_95U": find_column(cols, ["OR_95U"]),
        "OR_95L": find_column(cols, ["OR_95L"]),
        "SNPID": find_column(cols, ["SNP", "SNPID", "rsID", "RSID"]),
        "N": find_column(cols, ["N"]),
    }
    required = ["CHR", "POS", "P", "EA", "NEA"]
    missing = [x for x in required if mapping[x] is None]
    if missing:
        raise RuntimeError(f"QC GWAS missing required columns {missing}: {path}")

    effect_ready = bool(
        mapping["Z"]
        or (mapping["BETA"] and mapping["SE"])
        or (
            mapping["OR"]
            and (
                mapping["SE"]
                or (mapping["OR_95U"] and mapping["OR_95L"])
            )
        )
    )
    if not effect_ready:
        raise RuntimeError(
            "QC GWAS needs Z, BETA+SE, OR+SE, or OR+95% CI "
            f"for SuSiE-RSS: {path}"
        )
    return mapping


def extract_regional_gwas_once(
    qc_file: Path,
    loci: pd.DataFrame,
    regional_dir: Path,
    chunk_size: int,
) -> dict[str, Path]:
    mapping = inspect_qc_columns(qc_file)
    usecols = [x for x in mapping.values() if x is not None]
    usecols = list(dict.fromkeys(usecols))
    by_chr = defaultdict(list)
    for _, row in loci.iterrows():
        by_chr[int(row["CHR"])].append(row.to_dict())

    buffers: dict[str, list[pd.DataFrame]] = {str(x): [] for x in loci["LOCUS_ID"]}
    reader = pd.read_csv(
        qc_file,
        sep="\t",
        compression="infer",
        usecols=usecols,
        chunksize=chunk_size,
        low_memory=False,
    )

    total = 0
    extracted = 0
    for chunk in reader:
        total += len(chunk)
        canon = pd.DataFrame(index=chunk.index)
        canon["CHR"] = chunk[mapping["CHR"]].map(clean_chr)
        canon["POS"] = pd.to_numeric(chunk[mapping["POS"]], errors="coerce")
        canon["P"] = pd.to_numeric(chunk[mapping["P"]], errors="coerce")
        canon["EA"] = chunk[mapping["EA"]].map(clean_allele)
        canon["NEA"] = chunk[mapping["NEA"]].map(clean_allele)
        canon["SNPID"] = chunk[mapping["SNPID"]].astype(str) if mapping["SNPID"] else ""
        canon["N"] = pd.to_numeric(chunk[mapping["N"]], errors="coerce") if mapping["N"] else np.nan
        # -------------------------------------------------------------
        # Effect statistic harmonisation for SuSiE-RSS.
        #
        # Preferred inputs:
        #   1. Z directly
        #   2. BETA + SE
        #   3. OR + SE       -> BETA = log(OR), Z = BETA / SE
        #   4. OR + 95% CI   -> reconstruct SE on the log-OR scale
        # -------------------------------------------------------------
        if mapping["BETA"]:
            canon["BETA"] = pd.to_numeric(
                chunk[mapping["BETA"]], errors="coerce"
            )
        elif mapping["OR"]:
            odds_ratio = pd.to_numeric(
                chunk[mapping["OR"]], errors="coerce"
            )
            odds_ratio = odds_ratio.where(odds_ratio > 0)
            canon["BETA"] = np.log(odds_ratio)
        else:
            canon["BETA"] = np.nan

        if mapping["SE"]:
            canon["SE"] = pd.to_numeric(
                chunk[mapping["SE"]], errors="coerce"
            )
        elif mapping["OR_95U"] and mapping["OR_95L"]:
            upper = pd.to_numeric(
                chunk[mapping["OR_95U"]], errors="coerce"
            )
            lower = pd.to_numeric(
                chunk[mapping["OR_95L"]], errors="coerce"
            )
            upper = upper.where(upper > 0)
            lower = lower.where(lower > 0)
            canon["SE"] = (
                np.log(upper) - np.log(lower)
            ) / (2.0 * 1.959963984540054)
        else:
            canon["SE"] = np.nan

        canon.loc[
            (~np.isfinite(canon["SE"])) | (canon["SE"] <= 0),
            "SE",
        ] = np.nan

        if mapping["Z"]:
            canon["Z"] = pd.to_numeric(
                chunk[mapping["Z"]], errors="coerce"
            )
        else:
            canon["Z"] = canon["BETA"] / canon["SE"]

        canon = canon.dropna(subset=["CHR", "POS", "P", "EA", "NEA", "Z"])
        if canon.empty:
            continue
        canon["CHR"] = canon["CHR"].astype(int)
        canon["POS"] = canon["POS"].astype(int)
        canon = canon[np.isfinite(canon["Z"])].copy()

        for chrom, locus_list in by_chr.items():
            chrom_data = canon[canon["CHR"] == chrom]
            if chrom_data.empty:
                continue
            for locus in locus_list:
                mask = chrom_data["POS"].between(
                    int(locus["LOCUS_START"]), int(locus["LOCUS_END"]), inclusive="both"
                )
                selected = chrom_data.loc[mask].copy()
                if not selected.empty:
                    buffers[str(locus["LOCUS_ID"])].append(selected)
                    extracted += len(selected)

    regional_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for _, locus in loci.iterrows():
        locus_id = str(locus["LOCUS_ID"])
        path = regional_dir / f"{locus_id}.tsv.gz"
        if buffers[locus_id]:
            frame = pd.concat(buffers[locus_id], ignore_index=True)
            frame = frame.sort_values(["CHR", "POS", "P"], kind="stable")
            frame = frame.drop_duplicates(subset=["CHR", "POS", "EA", "NEA"], keep="first")
        else:
            frame = pd.DataFrame(columns=["CHR", "POS", "P", "EA", "NEA", "SNPID", "N", "BETA", "SE", "Z"])
        frame.to_csv(path, sep="\t", index=False, compression="gzip")
        paths[locus_id] = path

    print(f"GWAS rows scanned : {total:,}")
    print(f"Regional rows hit : {extracted:,}")
    return paths


# =============================================================================
# REFERENCE REGION + MATCHING
# =============================================================================


def create_regional_pvar(
    plink2: Path,
    root: Path,
    ancestry_code: str,
    chromosome: int,
    start: int,
    end: int,
    locus_dir: Path,
) -> Path:
    prefix = locus_dir / "reference_region"
    pvar = Path(str(prefix) + ".pvar")
    if pvar.exists() and pvar.stat().st_size > 0:
        return pvar
    run_command(
        [
            plink2,
            "--pfile", reference_prefix(root, ancestry_code, chromosome),
            "--chr", chromosome,
            "--from-bp", start,
            "--to-bp", end,
            "--make-just-pvar",
            "--out", prefix,
        ],
        locus_dir / "reference_region.command.log",
    )
    if not pvar.exists() or pvar.stat().st_size == 0:
        raise RuntimeError(f"Regional reference PVAR not generated for chr{chromosome}:{start}-{end}")
    return pvar


def read_pvar(path: Path) -> pd.DataFrame:
    records = []
    header = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            parts = line.rstrip("\n").split("\t")
            if header and len(parts) >= len(header):
                row = dict(zip(header, parts))
                chrom = clean_chr(row.get("CHROM", row.get("#CHROM", "")))
                pos = safe_int(row.get("POS"))
                if chrom and pos:
                    records.append(
                        {
                            "REFERENCE_CHR": chrom,
                            "REFERENCE_POS": pos,
                            "REFERENCE_ID": row.get("ID", ""),
                            "REFERENCE_REF": clean_allele(row.get("REF")),
                            "REFERENCE_ALT": clean_allele(row.get("ALT")),
                        }
                    )
    return pd.DataFrame(records)


def match_gwas_to_reference(regional: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    by_id = defaultdict(list)
    by_position = defaultdict(list)
    for _, row in reference.iterrows():
        record = row.to_dict()
        by_id[str(row["REFERENCE_ID"])].append(record)
        by_position[int(row["REFERENCE_POS"])].append(record)

    results = []
    for _, row in regional.iterrows():
        ea = clean_allele(row.get("EA"))
        nea = clean_allele(row.get("NEA"))
        if ea is None or nea is None:
            continue
        key = allele_key(ea, nea)
        rsid = clean_rsid(row.get("SNPID"))
        pos = int(row["POS"])
        selected = None
        method = None

        if rsid and rsid in by_id:
            candidates = [c for c in by_id[rsid] if allele_key(c["REFERENCE_REF"], c["REFERENCE_ALT"]) == key]
            if len(candidates) == 1:
                selected = candidates[0]
                method = "rsID+alleles"
        if selected is None:
            candidates = [c for c in by_position.get(pos, []) if allele_key(c["REFERENCE_REF"], c["REFERENCE_ALT"]) == key]
            if len(candidates) == 1:
                selected = candidates[0]
                method = "POS+alleles"
        if selected is None:
            continue

        ref = clean_allele(selected["REFERENCE_REF"])
        alt = clean_allele(selected["REFERENCE_ALT"])
        z = float(row["Z"])
        if not np.isfinite(z):
            continue

        # PLINK --r-unphased ... ref-based uses the reference allele orientation.
        if ea == ref and nea == alt:
            z_ref = z
            orientation = "EA=REF"
        elif ea == alt and nea == ref:
            z_ref = -z
            orientation = "EA=ALT_FLIPPED"
        else:
            continue

        record = row.to_dict()
        record.update(
            {
                "REFERENCE_ID": selected["REFERENCE_ID"],
                "REFERENCE_POS": selected["REFERENCE_POS"],
                "REFERENCE_REF": ref,
                "REFERENCE_ALT": alt,
                "MATCH_METHOD": method,
                "ORIENTATION": orientation,
                "Z_REF": z_ref,
            }
        )
        if pd.notna(row.get("BETA", np.nan)):
            beta = float(row["BETA"])
            record["BETA_REF"] = beta if orientation == "EA=REF" else -beta
        results.append(record)

    if not results:
        return pd.DataFrame()
    matched = pd.DataFrame(results)
    matched = matched.sort_values("P", kind="stable").drop_duplicates("REFERENCE_ID", keep="first")
    return matched.reset_index(drop=True)


# =============================================================================
# LD MATRIX
# =============================================================================


def _parse_lead_positions(value: Any) -> list[int]:
    positions = []
    for token in str(value or "").split(";"):
        token = token.strip()
        if not token:
            continue
        try:
            pos = int(float(token))
        except Exception:
            continue
        if pos > 0:
            positions.append(pos)
    return sorted(set(positions))


def _balanced_ld_subset(
    matched: pd.DataFrame,
    max_variants: int,
    lead_positions: list[int],
) -> pd.DataFrame:
    """Deterministically cap an oversized locus while preserving lead neighborhoods.

    Variants are first assigned to their nearest independent lead signal.  An
    approximately equal quota is selected around each lead by physical distance,
    then any remaining slots are filled globally by distance and association P.
    This is only invoked AFTER the 1000G MAF/GENO filters have been applied.
    """
    if len(matched) <= max_variants:
        return matched.copy()

    work = matched.copy().reset_index(drop=True)
    work["P"] = pd.to_numeric(work["P"], errors="coerce")
    work["POS"] = pd.to_numeric(work["POS"], errors="coerce")

    leads = [int(x) for x in lead_positions if int(x) > 0]
    leads = sorted(set(leads))

    if not leads:
        # Fallback for an unexpected missing LEAD_POSITIONS field: retain the
        # strongest associations deterministically rather than failing outright.
        return (
            work.sort_values(["P", "POS", "REFERENCE_ID"], kind="stable")
            .head(max_variants)
            .reset_index(drop=True)
        )

    pos = work["POS"].to_numpy(dtype=float)
    lead_array = np.asarray(leads, dtype=float)
    distances = np.abs(pos[:, None] - lead_array[None, :])
    nearest_index = np.argmin(distances, axis=1)
    nearest_distance = distances[np.arange(len(work)), nearest_index]

    work["_NEAREST_LEAD_INDEX"] = nearest_index
    work["_NEAREST_LEAD_POS"] = [leads[i] for i in nearest_index]
    work["_LEAD_DISTANCE"] = nearest_distance

    quota = max(1, max_variants // len(leads))
    selected_indices: list[int] = []

    for lead_index in range(len(leads)):
        group = work[work["_NEAREST_LEAD_INDEX"] == lead_index]
        if group.empty:
            continue
        group = group.sort_values(
            ["_LEAD_DISTANCE", "P", "POS", "REFERENCE_ID"],
            kind="stable",
        )
        selected_indices.extend(group.head(quota).index.tolist())

    # De-duplicate and respect the exact cap.
    selected_indices = list(dict.fromkeys(selected_indices))[:max_variants]

    remaining_slots = max_variants - len(selected_indices)
    if remaining_slots > 0:
        remainder = work.drop(index=selected_indices, errors="ignore")
        remainder = remainder.sort_values(
            ["_LEAD_DISTANCE", "P", "POS", "REFERENCE_ID"],
            kind="stable",
        )
        selected_indices.extend(remainder.head(remaining_slots).index.tolist())

    selected = work.loc[selected_indices].copy()
    selected = selected.sort_values(["P", "POS", "REFERENCE_ID"], kind="stable")
    return selected.reset_index(drop=True)


def calculate_ld(
    plink2: Path,
    root: Path,
    ancestry_code: str,
    chromosome: int,
    matched: pd.DataFrame,
    locus_dir: Path,
    maf: float,
    geno: float,
    threads: int,
    max_variants: int,
    lead_positions: list[int] | None = None,
) -> tuple[Path, Path, int]:
    """Build a validated LD matrix with post-reference-QC size protection."""
    lead_positions = lead_positions or []

    # ------------------------------------------------------------------
    # 1. Start from all GWAS/reference matched variants.
    # ------------------------------------------------------------------
    all_ids_file = locus_dir / "matched_reference_ids_all.txt"
    matched["REFERENCE_ID"].to_csv(all_ids_file, index=False, header=False)

    # ------------------------------------------------------------------
    # 2. Apply 1000G reference QC BEFORE enforcing the LD size cap.
    #    --write-snplist gives the exact variants surviving MAF/GENO and
    #    biallelic filters without constructing the dense LD matrix yet.
    # ------------------------------------------------------------------
    qc_prefix = locus_dir / f"{ancestry_code}_reference_qc"
    run_command(
        [
            plink2,
            "--pfile", reference_prefix(root, ancestry_code, chromosome),
            "--extract", all_ids_file,
            "--maf", maf,
            "--geno", geno,
            "--min-alleles", 2,
            "--max-alleles", 2,
            "--write-snplist",
            "--threads", threads,
            "--out", qc_prefix,
        ],
        locus_dir / f"{ancestry_code}_reference_qc.command.log",
    )

    snplist = Path(str(qc_prefix) + ".snplist")
    if not snplist.exists():
        raise RuntimeError("PLINK2 reference-QC SNP list is missing")

    surviving_ids = [
        x.strip() for x in snplist.read_text().splitlines() if x.strip()
    ]
    surviving_set = set(surviving_ids)

    qc_matched = matched[
        matched["REFERENCE_ID"].astype(str).isin(surviving_set)
    ].copy()

    # Preserve PLINK/reference order where possible for deterministic auditing.
    order = {variant_id: i for i, variant_id in enumerate(surviving_ids)}
    qc_matched["_REFERENCE_QC_ORDER"] = (
        qc_matched["REFERENCE_ID"].astype(str).map(order)
    )
    qc_matched = qc_matched.sort_values(
        ["_REFERENCE_QC_ORDER", "P"], kind="stable"
    ).drop(columns=["_REFERENCE_QC_ORDER"])
    qc_matched = qc_matched.reset_index(drop=True)

    qc_matched.to_csv(
        locus_dir / "reference_qc_matched_variants.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    if len(qc_matched) < 2:
        raise RuntimeError(
            f"Only {len(qc_matched)} matched variants survive reference MAF/GENO QC"
        )

    # ------------------------------------------------------------------
    # 3. If still oversized, deterministically retain balanced local
    #    neighborhoods around the independent lead signals.
    # ------------------------------------------------------------------
    selected = qc_matched
    selection_applied = len(qc_matched) > max_variants

    if selection_applied:
        print(
            f"[WARNING] {len(qc_matched):,} variants survive reference QC; "
            f"LD safety cap is {max_variants:,}."
        )
        print(
            "[WARNING] Applying deterministic balanced selection around "
            "independent lead positions."
        )
        selected = _balanced_ld_subset(
            qc_matched,
            max_variants=max_variants,
            lead_positions=lead_positions,
        )

    selected.to_csv(
        locus_dir / "ld_selected_variants.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    selection_summary = {
        "N_MATCHED_BEFORE_REFERENCE_QC": int(len(matched)),
        "N_AFTER_REFERENCE_QC": int(len(qc_matched)),
        "MAX_LD_VARIANTS": int(max_variants),
        "SELECTION_APPLIED": bool(selection_applied),
        "N_SELECTED_FOR_LD": int(len(selected)),
        "LEAD_POSITIONS": [int(x) for x in lead_positions],
        "SELECTION_METHOD": (
            "balanced_nearest_lead"
            if selection_applied and lead_positions
            else "association_rank_fallback"
            if selection_applied
            else "all_reference_qc_variants"
        ),
    }
    (locus_dir / "ld_variant_selection_summary.json").write_text(
        json.dumps(selection_summary, indent=2), encoding="utf-8"
    )

    ids_file = locus_dir / "matched_reference_ids.txt"
    selected["REFERENCE_ID"].to_csv(ids_file, index=False, header=False)

    # ------------------------------------------------------------------
    # 4. Build the dense signed LD matrix from the final selected IDs.
    # ------------------------------------------------------------------
    ld_prefix = locus_dir / f"{ancestry_code}_LD"
    run_command(
        [
            plink2,
            "--pfile", reference_prefix(root, ancestry_code, chromosome),
            "--extract", ids_file,
            "--maf", maf,
            "--geno", geno,
            "--min-alleles", 2,
            "--max-alleles", 2,
            "--r-unphased", "square", "bin4", "ref-based", "yes-really",
            "--threads", threads,
            "--out", ld_prefix,
        ],
        locus_dir / f"{ancestry_code}_LD.command.log",
    )

    ld_file = locus_dir / f"{ancestry_code}_LD.unphased.vcor1.bin"
    vars_file = locus_dir / f"{ancestry_code}_LD.unphased.vcor1.bin.vars"
    if not ld_file.exists() or not vars_file.exists():
        raise RuntimeError("PLINK2 LD output is missing")

    variants = [x.strip() for x in vars_file.read_text().splitlines() if x.strip()]
    p = len(variants)
    expected_bytes = p * p * 4
    actual_bytes = ld_file.stat().st_size
    if p < 2 or actual_bytes != expected_bytes:
        raise RuntimeError(
            f"LD matrix validation failed: p={p}, "
            f"expected_bytes={expected_bytes}, actual_bytes={actual_bytes}"
        )

    if p > max_variants:
        raise RuntimeError(
            f"Internal error: final LD matrix has {p:,} variants, "
            f"above safety cap {max_variants:,}"
        )

    return ld_file, vars_file, p


# =============================================================================
# SUSIE R HELPER
#
# susie_rss(z, R, n) with estimate_residual_variance = FALSE (reference-panel LD;
# Zou et al. 2022, PLoS Genet), plus diagnostics:
#   estimate_s_rss  - consistency of z with R (0 = consistent)
#   kriging_rss     - variants whose z disagrees with LD neighbours
#                     (logLR > 2 and |z| > 2: possible allele flips / mismatches)
#   convergence     - fit$converged within max_iter
#   purity          - min / mean / median |r| inside each credible set
#   L = 1 fit       - single-effect posterior does not use the LD matrix; its top
#                     variant is an LD-independent cross-check of the main fit
# =============================================================================

SUSIE_R_CODE = r'''options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 15) stop("Expected 15 arguments")
stats_file <- args[1]; ld_file <- args[2]; vars_file <- args[3]; outdir <- args[4]
n_gwas <- as.numeric(args[5]); n_ref <- as.numeric(args[6]); L <- as.integer(args[7])
coverage <- as.numeric(args[8]); ld_ridge <- as.numeric(args[9])
estimate_s_max_p <- as.integer(args[10]); locus_id <- args[11]
max_iter <- as.integer(args[12]); refine <- toupper(args[13]) %in% c("TRUE", "1")
s_warn <- as.numeric(args[14]); s_high <- as.numeric(args[15])

suppressPackageStartupMessages(library(susieR))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
logcon <- file(file.path(outdir, "susie_rss.log"), open="wt")
sink(logcon, type="output"); sink(logcon, type="message")
on.exit({sink(type="message"); sink(type="output"); close(logcon)}, add=TRUE)

cat("SuSiE-RSS locus:", locus_id, " susieR", as.character(packageVersion("susieR")), "\n")
cat("GWAS N:", n_gwas, " Reference N:", n_ref, " L:", L, " coverage:", coverage, "\n")

vars <- scan(vars_file, what="character", quiet=TRUE)
p <- length(vars)
if (p < 2) stop("Fewer than two LD variants")
con <- file(ld_file, "rb")
values <- readBin(con, what="numeric", n=p*p, size=4, endian="little")
close(con)
if (length(values) != p*p) stop("LD matrix length mismatch")
R <- matrix(values, nrow=p, ncol=p, byrow=TRUE)
rm(values); invisible(gc())
if (any(!is.finite(R))) stop("LD matrix contains non-finite values")
R <- (R + t(R))/2
diag(R) <- 1
if (is.finite(ld_ridge) && ld_ridge > 0) { R <- R * (1-ld_ridge); diag(R) <- 1 }

stats <- read.delim(stats_file, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c("REFERENCE_ID","Z_REF") %in% names(stats))) stop("Matched stats missing REFERENCE_ID/Z_REF")
idx <- match(vars, stats$REFERENCE_ID)
if (any(is.na(idx))) stop("LD variants missing from matched GWAS")
stats <- stats[idx,,drop=FALSE]
if (!all(stats$REFERENCE_ID == vars)) stop("GWAS and LD order mismatch")
z <- as.numeric(stats$Z_REF)
if (any(!is.finite(z))) stop("Non-finite Z scores")

# ---- LD / z consistency -------------------------------------------------------
s_est <- NA_real_; s_status <- "NOT_RUN"; s_error <- NA_character_; n_flips <- NA_integer_
if (p <= estimate_s_max_p) {
  r <- tryCatch(list(ok=TRUE, v=susieR::estimate_s_rss(z=z, R=R, n=n_gwas)),
                error=function(e) list(ok=FALSE, v=conditionMessage(e)))
  if (r$ok) { s_est <- r$v; s_status <- "SUCCESS" } else { s_status <- "FAILED"; s_error <- r$v }
  if (r$ok) {
    k <- tryCatch(susieR::kriging_rss(z=z, R=R, n=n_gwas, s=s_est)$conditional_dist,
                  error=function(e) NULL)
    if (!is.null(k)) {
      k <- cbind(REFERENCE_ID=stats$REFERENCE_ID, POS=stats$POS, k)
      flagged <- k$logLR > 2 & abs(k$z) > 2
      n_flips <- sum(flagged, na.rm=TRUE)
      k$POSSIBLE_ALLELE_FLIP <- flagged
      write.table(k, file.path(outdir, "kriging.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    }
  }
} else {
  s_status <- paste0("SKIPPED_P_GT_", estimate_s_max_p)
}

# ---- main fit -------------------------------------------------------------------
fit_args <- list(z=z, R=R, n=n_gwas, L=min(L,p), coverage=coverage, min_abs_corr=0.5,
                 estimate_residual_variance=FALSE, estimate_prior_variance=TRUE,
                 refine=refine, max_iter=max_iter, verbose=FALSE)
fit <- do.call(susieR::susie_rss, fit_args)
saveRDS(fit, file.path(outdir, "susie_fit.rds"))
converged <- isTRUE(fit$converged)

# ---- LD-independent single-effect fit (L = 1) -------------------------------------
fit1 <- tryCatch(do.call(susieR::susie_rss, modifyList(fit_args, list(L=1, refine=FALSE))),
                 error=function(e) NULL)

pip <- fit$pip
if (length(pip) != nrow(stats)) stop("PIP length mismatch")
variant <- stats
variant$PIP <- as.numeric(pip)
variant$PIP_L1 <- if (is.null(fit1)) NA_real_ else as.numeric(fit1$pip)
variant$IN_95_CREDIBLE_SET <- FALSE
variant$CREDIBLE_SET <- ""

cs_rows <- list(); purity_min <- c()
if (!is.null(fit$sets) && !is.null(fit$sets$cs)) {
  cs <- fit$sets$cs
  purity <- fit$sets$purity
  for (i in seq_along(cs)) {
    members <- as.integer(cs[[i]])
    members <- members[members >= 1 & members <= nrow(variant)]
    if (length(members) == 0) next
    cs_name <- names(cs)[i]
    if (is.null(cs_name) || is.na(cs_name) || cs_name == "") cs_name <- paste0("L", i)
    cs_name <- paste0(locus_id, "_", cs_name)
    pmin <- if (!is.null(purity)) purity[i, 1] else NA_real_
    pmean <- if (!is.null(purity)) purity[i, 2] else NA_real_
    pmed <- if (!is.null(purity)) purity[i, 3] else NA_real_
    purity_min <- c(purity_min, pmin)
    cov_i <- if (!is.null(fit$sets$coverage)) fit$sets$coverage[i] else NA_real_
    l_index <- if (!is.null(fit$sets$cs_index)) fit$sets$cs_index[i] else NA_integer_
    log10bf <- if (!is.na(l_index) && !is.null(fit$lbf)) fit$lbf[l_index] / log(10) else NA_real_
    variant$IN_95_CREDIBLE_SET[members] <- TRUE
    for (j in members) {
      existing <- variant$CREDIBLE_SET[j]
      variant$CREDIBLE_SET[j] <- ifelse(existing == "", cs_name, paste(existing, cs_name, sep=";"))
      cs_rows[[length(cs_rows)+1]] <- data.frame(
        LOCUS_ID=locus_id, CREDIBLE_SET=cs_name, REFERENCE_ID=variant$REFERENCE_ID[j],
        CHR=variant$CHR[j], POS=variant$POS[j], PIP=variant$PIP[j],
        CS_SIZE=length(members), CS_COVERAGE=cov_i, CS_LOG10BF=log10bf,
        PURITY_MIN_ABS_R=pmin, PURITY_MEAN_ABS_R=pmean, PURITY_MEDIAN_ABS_R=pmed,
        stringsAsFactors=FALSE)
    }
  }
}
variant <- variant[order(-variant$PIP),,drop=FALSE]
write.table(variant, file.path(outdir,"susie_variants.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (length(cs_rows) > 0) cs_df <- do.call(rbind, cs_rows) else cs_df <- data.frame(
  LOCUS_ID=character(), CREDIBLE_SET=character(), REFERENCE_ID=character(), CHR=integer(),
  POS=integer(), PIP=numeric(), CS_SIZE=integer(), CS_COVERAGE=numeric(), CS_LOG10BF=numeric(),
  PURITY_MIN_ABS_R=numeric(), PURITY_MEAN_ABS_R=numeric(), PURITY_MEDIAN_ABS_R=numeric())
write.table(cs_df, file.path(outdir,"credible_sets.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ---- locus QC -----------------------------------------------------------------------
top <- variant$REFERENCE_ID[1]
top_l1 <- if (all(is.na(variant$PIP_L1))) NA_character_ else variant$REFERENCE_ID[which.max(variant$PIP_L1)]
notes <- c()
if (!converged) notes <- c(notes, paste0("not converged in ", max_iter, " iterations"))
if (s_status != "SUCCESS") notes <- c(notes, paste0("LD consistency not checked (", s_status, ")"))
if (!is.na(s_est) && s_est > s_high) notes <- c(notes, sprintf("strong LD/z mismatch (s=%.2f)", s_est))
if (!is.na(s_est) && s_est > s_warn && s_est <= s_high) notes <- c(notes, sprintf("LD/z mismatch (s=%.2f)", s_est))
if (!is.na(n_flips) && n_flips > 0) notes <- c(notes, paste0(n_flips, " possible allele flips (kriging_rss)"))
if (length(purity_min) > 0 && any(purity_min < 0.5, na.rm=TRUE)) notes <- c(notes, "credible set purity < 0.5")

summary <- data.frame(
  LOCUS_ID=locus_id, N_VARIANTS=p,
  N_CREDIBLE_SETS=ifelse(is.null(fit$sets$cs), 0, length(fit$sets$cs)),
  MAX_PIP=max(variant$PIP, na.rm=TRUE), TOP_VARIANT=top, TOP_VARIANT_L1=top_l1,
  TOP_AGREES_L1=ifelse(is.na(top_l1), NA, top == top_l1),
  CONVERGED=converged, N_ITER=ifelse(is.null(fit$niter), NA, fit$niter),
  MIN_CS_PURITY=ifelse(length(purity_min) > 0, min(purity_min, na.rm=TRUE), NA),
  ESTIMATE_S_RSS=s_est, ESTIMATE_S_STATUS=s_status, N_POSSIBLE_ALLELE_FLIPS=n_flips,
  MAX_ABS_Z=max(abs(z)), SUSIE_L=min(L,p), GWAS_N=n_gwas, REFERENCE_N=n_ref,
  COVERAGE=coverage, LD_RIDGE=ld_ridge, SUSIE_VERSION=as.character(packageVersion("susieR")),
  LOCUS_QC=ifelse(length(notes) == 0, "PASS", "WARN"),
  LOCUS_QC_NOTES=paste(notes, collapse="; "),
  stringsAsFactors=FALSE)
write.table(summary, file.path(outdir,"susie_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

diagnostics <- data.frame(LOCUS_ID=locus_id, ESTIMATE_S_RSS=s_est, ESTIMATE_S_STATUS=s_status,
                          ESTIMATE_S_ERROR=s_error, N_POSSIBLE_ALLELE_FLIPS=n_flips,
                          CONVERGED=converged, MAX_ABS_Z=max(abs(z)), stringsAsFactors=FALSE)
write.table(diagnostics, file.path(outdir,"susie_diagnostics.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
'''


def write_susie_r_script(path: Path) -> Path:
    """Always (re)written so an older helper is never reused."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text(encoding="utf-8") != SUSIE_R_CODE:
        path.write_text(SUSIE_R_CODE, encoding="utf-8", newline="\n")
    return path


def augment_locus_summary(locus_dir: Path, regional: pd.DataFrame, lead_positions: list[int]) -> None:
    """Add Python-side locus checks to susie/susie_summary.tsv (LD cap, N spread,
    leads absent from the LD panel) and update LOCUS_QC."""
    summary_file = locus_dir / "susie" / "susie_summary.tsv"
    summary = pd.read_csv(summary_file, sep="\t")
    notes = [x for x in str(summary.at[0, "LOCUS_QC_NOTES"]).split("; ") if x and x != "nan"]

    selection_file = locus_dir / "ld_variant_selection_summary.json"
    if selection_file.exists():
        selection = json.loads(selection_file.read_text(encoding="utf-8"))
        summary["N_MATCHED_BEFORE_REFERENCE_QC"] = selection.get("N_MATCHED_BEFORE_REFERENCE_QC")
        summary["N_AFTER_REFERENCE_QC"] = selection.get("N_AFTER_REFERENCE_QC")
        summary["LD_CAP_APPLIED"] = bool(selection.get("SELECTION_APPLIED"))
        if selection.get("SELECTION_APPLIED"):
            notes.append(f"LD cap: {selection.get('N_SELECTED_FOR_LD')} of "
                         f"{selection.get('N_AFTER_REFERENCE_QC')} variants used")

    n = pd.to_numeric(regional.get("N", pd.Series(dtype=float)), errors="coerce").dropna()
    if len(n):
        ratio = float(n.min() / n.max()) if n.max() > 0 else np.nan
        summary["N_VARIANT_MIN"], summary["N_VARIANT_MAX"] = int(n.min()), int(n.max())
        summary["N_VARIANT_MIN_MAX_RATIO"] = round(ratio, 4)
        if np.isfinite(ratio) and ratio < 0.5:
            notes.append(f"per-variant N varies (min/max {ratio:.2f}); SuSiE-RSS assumes one N")

    variants_file = locus_dir / "susie" / "susie_variants.tsv"
    if variants_file.exists():
        fm = pd.read_csv(variants_file, sep="\t", usecols=["POS", "PIP", "IN_95_CREDIBLE_SET"])
        # A high PIP outside every credible set is typical of LD-reference noise
        # (weak z, strong "explained" signal): flag, do not interpret.
        outside = int(((fm["PIP"] > 0.5) & (fm["IN_95_CREDIBLE_SET"].astype(str) != "True")).sum())
        summary["N_HIGH_PIP_OUTSIDE_CS"] = outside
        if outside:
            notes.append(f"{outside} variant(s) with PIP>0.5 outside credible sets (likely LD artefact)")
        if lead_positions:
            positions = set(fm["POS"].astype(int))
            absent = [p for p in lead_positions if p not in positions]
            summary["N_LEADS_NOT_IN_LD"] = len(absent)
            if absent:
                notes.append(f"{len(absent)} Step05 lead position(s) not in the LD variant set")

    summary["LOCUS_QC"] = "PASS" if not notes else "WARN"
    summary["LOCUS_QC_NOTES"] = "; ".join(notes)
    summary.to_csv(summary_file, sep="\t", index=False)


# =============================================================================
# PER-LOCUS RESUME
#
# A locus is skipped on rerun only if loci/<L>/LOCUS_COMPLETE.json exists,
# its coordinates, input-file signatures and fine-mapping parameters all match
# the current run, and the SuSiE outputs are present. Otherwise it is re-run.
# =============================================================================

LOCUS_COMPLETE_MARKER = "LOCUS_COMPLETE.json"
LOCUS_REQUIRED_OUTPUTS = (
    "susie_variants.tsv",
    "susie_summary.tsv",
    "susie_fit.rds",
)


def file_signature(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {"PATH": str(path), "SIZE": int(stat.st_size), "MTIME_NS": int(stat.st_mtime_ns)}


def locus_identity(locus) -> dict[str, Any]:
    return {
        "LOCUS_ID": str(locus["LOCUS_ID"]),
        "CHR": int(locus["CHR"]),
        "LOCUS_START": int(locus["LOCUS_START"]),
        "LOCUS_END": int(locus["LOCUS_END"]),
        "N_INDEPENDENT_SIGNALS": int(locus["N_INDEPENDENT_SIGNALS"]),
        "LEAD_POSITIONS": str(locus.get("LEAD_POSITIONS", "")),
    }


def read_valid_locus_marker(locus_dir: Path, locus, run_signature: dict) -> dict | None:
    marker = locus_dir / LOCUS_COMPLETE_MARKER
    if not marker.exists() or (locus_dir / "LOCUS_FAILED.txt").exists():
        return None
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
    except Exception:
        return None
    if payload.get("LOCUS") != locus_identity(locus):
        return None
    if payload.get("RUN_SIGNATURE") != run_signature:
        return None
    susie_dir = locus_dir / "susie"
    for name in LOCUS_REQUIRED_OUTPUTS:
        path = susie_dir / name
        if not path.exists() or path.stat().st_size == 0:
            return None
    return payload


def write_locus_marker(
    locus_dir: Path,
    locus,
    run_signature: dict,
    n_regional: int,
    n_matched: int,
) -> None:
    payload = {
        "LOCUS": locus_identity(locus),
        "RUN_SIGNATURE": run_signature,
        "N_REGIONAL_GWAS_VARIANTS": int(n_regional),
        "N_MATCHED_VARIANTS": int(n_matched),
        "COMPLETED_UTC": datetime.now(timezone.utc).isoformat(),
    }
    temporary = locus_dir / (LOCUS_COMPLETE_MARKER + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(temporary, locus_dir / LOCUS_COMPLETE_MARKER)


def read_locus_outputs(
    susie_dir: Path,
    accession: str,
    phenotype: str,
    ancestry_code: str,
    locus_id: str,
    chromosome: int,
    start: int,
    end: int,
    n_signals: int,
    n_regional: int,
    n_matched: int,
) -> tuple[pd.DataFrame, pd.DataFrame | None, pd.DataFrame]:
    """Load one locus' SuSiE outputs exactly as the worker always has."""
    variants = pd.read_csv(susie_dir / "susie_variants.tsv", sep="\t", low_memory=False)
    variants.insert(0, "STUDY_ACCESSION", accession)
    variants.insert(1, "PHENOTYPE", phenotype)
    variants.insert(2, "ANCESTRY_CODE", ancestry_code)
    variants.insert(3, "LOCUS_ID", locus_id)

    cs = None
    cs_path = susie_dir / "credible_sets.tsv"
    if cs_path.exists() and cs_path.stat().st_size > 0:
        cs = pd.read_csv(cs_path, sep="\t", low_memory=False)
        if not cs.empty:
            cs.insert(0, "STUDY_ACCESSION", accession)
            cs.insert(1, "PHENOTYPE", phenotype)
            cs.insert(2, "ANCESTRY_CODE", ancestry_code)
        else:
            cs = None

    s = pd.read_csv(susie_dir / "susie_summary.tsv", sep="\t", low_memory=False)
    s.insert(0, "STUDY_ACCESSION", accession)
    s["CHR"] = chromosome
    s["LOCUS_START"] = start
    s["LOCUS_END"] = end
    s["N_INDEPENDENT_SIGNALS"] = n_signals
    s["N_REGIONAL_GWAS_VARIANTS"] = n_regional
    s["N_MATCHED_VARIANTS"] = n_matched
    return variants, cs, s


# =============================================================================
# WORKER
# =============================================================================


def worker_mode() -> None:
    manifest_value = os.environ.get("GWAS_FINEMAP_MANIFEST")
    task_value = os.environ.get("SLURM_ARRAY_TASK_ID")
    if not manifest_value or not task_value:
        raise RuntimeError("GWAS_FINEMAP_MANIFEST and SLURM_ARRAY_TASK_ID are required in worker mode")

    manifest_file = Path(manifest_value).resolve()
    manifest = pd.read_csv(manifest_file, sep="\t", dtype=str)
    task_id = int(task_value)
    task_numbers = pd.to_numeric(manifest["FINEMAP_TASK_ID"], errors="coerce")
    selected = manifest[task_numbers == task_id]
    if selected.empty:
        raise RuntimeError(f"No manifest rows for array task {task_id}")

    failed = []
    for _, row in selected.iterrows():
        try:
            finemap_one_study(row, task_id)
        except Exception as error:
            failed.append(f"{row['STUDY_ACCESSION']}: {error}")

    if failed:
        raise SystemExit("Fine-mapping incomplete for:\n  " + "\n  ".join(failed))


def finemap_one_study(row, task_id: int) -> None:
    root = PROJECT_ROOT
    accession = str(row["STUDY_ACCESSION"])
    phenotype = str(row["PHENOTYPE"])
    ancestry_code = str(row["ANCESTRY_CODE"]).upper()
    ancestry_label = str(row["ANCESTRY_LABEL"])
    qc_file = Path(row["QC_FILE"]).resolve()
    lead_file = Path(row["LEAD_VARIANTS_FILE"]).resolve()
    output_dir = Path(row["OUTPUT_DIR"]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    gwas_n = int(float(row["GWAS_N"]))
    gwas_n_source = str(row["GWAS_N_SOURCE"])
    reference_n = int(float(row["REFERENCE_N"]))
    locus_window_kb = int(float(row["LOCUS_WINDOW_KB"]))
    reference_maf = float(row["REFERENCE_MAF"])
    reference_geno = float(row["REFERENCE_GENO"])
    min_variants = int(float(row["MIN_VARIANTS"]))
    max_ld_variants = int(float(row["MAX_LD_VARIANTS"]))
    coverage = float(row["COVERAGE"])
    susie_l = int(float(row["SUSIE_L"]))
    max_susie_l = int(float(row["MAX_SUSIE_L"]))
    ld_ridge = float(row["LD_RIDGE"])
    chunk_size = int(float(row["CHUNK_SIZE"]))
    max_iter = int(float(row.get("MAX_ITER", DEFAULT_MAX_ITER) or DEFAULT_MAX_ITER))
    refine = str(row.get("REFINE", "False")).strip().lower() in {"true", "1", "yes"}
    estimate_s_max = int(float(row.get("ESTIMATE_S_MAX", DEFAULT_ESTIMATE_S_MAX) or DEFAULT_ESTIMATE_S_MAX))
    phenotype_type = str(row.get("PHENOTYPE_TYPE", ""))
    plink2 = Path(row["PLINK2"]).resolve()
    rscript = Path(row["R_SCRIPT"]).resolve()
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", DEFAULT_CPUS))
    # Optional FORCE column: re-run even when outputs exist.
    force = str(row.get("FORCE", "False")).strip().lower() in {"true", "1", "yes", "y"}

    final_variants = output_dir / f"{accession}_finemapped_variants.tsv.gz"
    final_cs = output_dir / f"{accession}_95pct_credible_sets.tsv"
    final_summary = output_dir / f"{accession}_finemapping_summary.tsv"
    failed_file = output_dir / f"{accession}_failed_loci.tsv"
    summary_json = output_dir / "finemapping_summary.json"
    hard_failure = output_dir / "FINEMAPPING_FAILED.txt"

    if not force and final_variants.exists() and final_summary.exists() and summary_json.exists() and not hard_failure.exists():
        banner("FINEMAPPING ALREADY COMPLETE")
        print(accession)
        return

    banner("STEP 06 - SUSIE-RSS FINEMAPPING WORKER")
    print(f"Study          : {accession}")
    print(f"Phenotype      : {phenotype}")
    print(f"Ancestry       : {ancestry_label} ({ancestry_code})")
    print(f"Phenotype type : {phenotype_type}")
    print(f"GWAS N         : {gwas_n} [{gwas_n_source}]")
    print(f"Reference N    : {reference_n}")
    print(f"QC GWAS        : {qc_file}")
    print(f"Lead variants  : {lead_file}")
    print(f"Output         : {output_dir}")

    try:
        loci = load_leads_and_build_loci(lead_file, locus_window_kb)
        loci_file = output_dir / "loci_definition.tsv"
        loci.to_csv(loci_file, sep="\t", index=False)
        if loci.empty:
            # A study can complete Step 05 successfully yet have zero independent
            # genome-wide-significant lead variants.  That is a valid scientific
            # outcome, not a pipeline failure.  Write explicit empty deliverables
            # and a machine-readable completion status so downstream stages can
            # skip this study cleanly.
            pd.DataFrame(
                columns=[
                    "STUDY_ACCESSION", "PHENOTYPE", "ANCESTRY_CODE",
                    "LOCUS_ID", "REFERENCE_ID", "CHR", "POS", "PIP",
                ]
            ).to_csv(final_variants, sep="\t", index=False, compression="gzip")

            pd.DataFrame(
                columns=[
                    "STUDY_ACCESSION", "LOCUS_ID", "CREDIBLE_SET",
                    "REFERENCE_ID", "CHR", "POS", "PIP",
                ]
            ).to_csv(final_cs, sep="\t", index=False)

            no_loci_row = {
                "STUDY_ACCESSION": accession,
                "PHENOTYPE": phenotype,
                "ANCESTRY_CODE": ancestry_code,
                "ANCESTRY_LABEL": ancestry_label,
                "GWAS_N": gwas_n,
                "GWAS_N_SOURCE": gwas_n_source,
                "REFERENCE_N": reference_n,
                "STATUS": "NO_SIGNIFICANT_LOCI",
                "N_LOCI": 0,
                "N_LOCI_SUCCESS": 0,
                "N_LOCI_FAILED": 0,
                "MESSAGE": "Step 05 produced zero independent lead variants; no fine-mapping was required.",
            }
            pd.DataFrame([no_loci_row]).to_csv(
                final_summary, sep="\t", index=False
            )

            pd.DataFrame(
                columns=[
                    "STUDY_ACCESSION", "LOCUS_ID", "CHR",
                    "LOCUS_START", "LOCUS_END", "ERROR",
                ]
            ).to_csv(failed_file, sep="\t", index=False)

            summary = {
                **no_loci_row,
                "FINEMAPPED_VARIANTS_FILE": str(final_variants),
                "CREDIBLE_SETS_FILE": str(final_cs),
                "SUMMARY_FILE": str(final_summary),
                "FAILED_LOCI_FILE": str(failed_file),
                "COMPLETED_UTC": datetime.now(timezone.utc).isoformat(),
            }
            summary_json.write_text(
                json.dumps(summary, indent=2), encoding="utf-8"
            )

            if hard_failure.exists():
                hard_failure.unlink()

            banner("FINEMAPPING SKIPPED - NO SIGNIFICANT LOCI")
            print(f"Study  : {accession}")
            print("Status : NO_SIGNIFICANT_LOCI")
            print("Step 05 produced zero independent lead variants.")
            print("This is a valid completed outcome; there is nothing to fine-map.")
            return

        regional_dir = output_dir / "regional_gwas"
        locus_root = output_dir / "loci"
        locus_root.mkdir(parents=True, exist_ok=True)

        # Per-locus resume: everything that can change a locus result.
        run_signature = {
            "QC_FILE": file_signature(qc_file),
            "LEAD_VARIANTS_FILE": file_signature(lead_file),
            "GWAS_N": gwas_n,
            "REFERENCE_N": reference_n,
            "LOCUS_WINDOW_KB": locus_window_kb,
            "REFERENCE_MAF": reference_maf,
            "REFERENCE_GENO": reference_geno,
            "MIN_VARIANTS": min_variants,
            "MAX_LD_VARIANTS": max_ld_variants,
            "COVERAGE": coverage,
            "SUSIE_L": susie_l,
            "MAX_SUSIE_L": max_susie_l,
            "LD_RIDGE": ld_ridge,
            "ESTIMATE_S_MAX": estimate_s_max,
            "MAX_ITER": max_iter,
            "REFINE": refine,
            "R_HELPER_VERSION": R_HELPER_VERSION,
        }
        completed_loci = {}
        for _, locus in loci.iterrows():
            marker = read_valid_locus_marker(
                locus_root / str(locus["LOCUS_ID"]), locus, run_signature
            )
            if marker is not None:
                completed_loci[str(locus["LOCUS_ID"])] = marker
        pending_loci = loci[~loci["LOCUS_ID"].astype(str).isin(completed_loci)]
        print(f"Loci validated  : {len(completed_loci)} (skipped)")
        print(f"Loci to run     : {len(pending_loci)}")

        # Regions are extracted independently per locus, so extracting only
        # the pending loci gives identical regional files for those loci.
        regional_paths = (
            extract_regional_gwas_once(qc_file, pending_loci, regional_dir, chunk_size)
            if not pending_loci.empty
            else {}
        )
        r_helper = write_susie_r_script(output_dir / "run_susie_rss.R")

        fine_frames = []
        cs_frames = []
        summary_frames = []
        failures = []

        for _, locus in loci.iterrows():
            locus_id = str(locus["LOCUS_ID"])
            chromosome = int(locus["CHR"])
            start = int(locus["LOCUS_START"])
            end = int(locus["LOCUS_END"])
            n_signals = int(locus["N_INDEPENDENT_SIGNALS"])
            locus_dir = locus_root / locus_id
            locus_dir.mkdir(parents=True, exist_ok=True)

            if locus_id in completed_loci:
                marker = completed_loci[locus_id]
                print(f"[SKIP] {locus_id}: output already validated")
                variants, cs, s = read_locus_outputs(
                    locus_dir / "susie", accession, phenotype, ancestry_code,
                    locus_id, chromosome, start, end, n_signals,
                    marker["N_REGIONAL_GWAS_VARIANTS"], marker["N_MATCHED_VARIANTS"],
                )
                fine_frames.append(variants)
                if cs is not None:
                    cs_frames.append(cs)
                summary_frames.append(s)
                continue

            stale_marker = locus_dir / LOCUS_COMPLETE_MARKER
            if stale_marker.exists():
                stale_marker.unlink()

            try:
                regional = pd.read_csv(regional_paths[locus_id], sep="\t", compression="infer", low_memory=False)
                if len(regional) < min_variants:
                    raise RuntimeError(f"Only {len(regional)} regional GWAS variants; minimum is {min_variants}")

                pvar = create_regional_pvar(
                    plink2, root, ancestry_code, chromosome, start, end, locus_dir
                )
                reference = read_pvar(pvar)
                matched = match_gwas_to_reference(regional, reference)
                if len(matched) < min_variants:
                    raise RuntimeError(
                        f"Only {len(matched)} GWAS/reference matched variants; minimum is {min_variants}"
                    )
                matched_file = locus_dir / "matched_variants.tsv"
                matched.to_csv(matched_file, sep="\t", index=False)

                lead_positions = _parse_lead_positions(
                    locus.get("LEAD_POSITIONS", "")
                )
                ld_file, vars_file, p = calculate_ld(
                    plink2=plink2,
                    root=root,
                    ancestry_code=ancestry_code,
                    chromosome=chromosome,
                    matched=matched,
                    locus_dir=locus_dir,
                    maf=reference_maf,
                    geno=reference_geno,
                    threads=threads,
                    max_variants=max_ld_variants,
                    lead_positions=lead_positions,
                )
                if p < min_variants:
                    raise RuntimeError(f"Only {p} variants remain in LD matrix; minimum is {min_variants}")

                susie_dir = locus_dir / "susie"
                L = min(max(susie_l, n_signals), max_susie_l, p)
                run_command(
                    [
                        rscript,
                        r_helper,
                        matched_file,
                        ld_file,
                        vars_file,
                        susie_dir,
                        gwas_n,
                        reference_n,
                        L,
                        coverage,
                        ld_ridge,
                        estimate_s_max,
                        locus_id,
                        max_iter,
                        refine,
                        S_RSS_WARN,
                        S_RSS_HIGH,
                    ],
                    susie_dir / "r_command.log",
                )

                augment_locus_summary(locus_dir, regional, lead_positions)

                variants, cs, s = read_locus_outputs(
                    susie_dir, accession, phenotype, ancestry_code,
                    locus_id, chromosome, start, end, n_signals,
                    len(regional), len(matched),
                )
                fine_frames.append(variants)
                if cs is not None:
                    cs_frames.append(cs)
                summary_frames.append(s)

                locus_failure_marker = locus_dir / "LOCUS_FAILED.txt"
                if locus_failure_marker.exists():
                    locus_failure_marker.unlink()

                write_locus_marker(
                    locus_dir, locus, run_signature, len(regional), len(matched)
                )

            except Exception as exc:
                failure_text = f"{type(exc).__name__}: {exc}"
                failures.append(
                    {
                        "STUDY_ACCESSION": accession,
                        "LOCUS_ID": locus_id,
                        "CHR": chromosome,
                        "LOCUS_START": start,
                        "LOCUS_END": end,
                        "ERROR": failure_text,
                    }
                )
                (locus_dir / "LOCUS_FAILED.txt").write_text(
                    failure_text + "\n\n" + traceback.format_exc(), encoding="utf-8"
                )
                print(f"[FAILED] {locus_id}: {failure_text}", file=sys.stderr)

        if fine_frames:
            all_fine = pd.concat(fine_frames, ignore_index=True)
            all_fine.to_csv(final_variants, sep="\t", index=False, compression="gzip")
        else:
            pd.DataFrame().to_csv(final_variants, sep="\t", index=False, compression="gzip")

        if cs_frames:
            pd.concat(cs_frames, ignore_index=True).to_csv(final_cs, sep="\t", index=False)
        else:
            pd.DataFrame(columns=["STUDY_ACCESSION", "LOCUS_ID", "CREDIBLE_SET", "REFERENCE_ID", "CHR", "POS", "PIP"]).to_csv(
                final_cs, sep="\t", index=False
            )

        if summary_frames:
            pd.concat(summary_frames, ignore_index=True).to_csv(final_summary, sep="\t", index=False)
        else:
            pd.DataFrame().to_csv(final_summary, sep="\t", index=False)

        pd.DataFrame(failures).to_csv(failed_file, sep="\t", index=False)
        all_summary = pd.concat(summary_frames, ignore_index=True) if summary_frames else pd.DataFrame()
        summary = {
            "STUDY_ACCESSION": accession,
            "PHENOTYPE": phenotype,
            "ANCESTRY_CODE": ancestry_code,
            "PHENOTYPE_TYPE": phenotype_type,
            "N_CASES": row.get("N_CASES", ""),
            "N_CONTROLS": row.get("N_CONTROLS", ""),
            "N_TOTAL": row.get("N_TOTAL", ""),
            "N_EFF": row.get("N_EFF", ""),
            "N_LOCI_QC_PASS": int((all_summary.get("LOCUS_QC", pd.Series(dtype=str)) == "PASS").sum()),
            "N_LOCI_QC_WARN": int((all_summary.get("LOCUS_QC", pd.Series(dtype=str)) == "WARN").sum()),
            "GWAS_N": gwas_n,
            "GWAS_N_SOURCE": gwas_n_source,
            "REFERENCE_N": reference_n,
            "N_LOCI": int(len(loci)),
            "N_LOCI_SUCCESS": int(len(summary_frames)),
            "N_LOCI_FAILED": int(len(failures)),
            "FINEMAPPED_VARIANTS_FILE": str(final_variants),
            "CREDIBLE_SETS_FILE": str(final_cs),
            "SUMMARY_FILE": str(final_summary),
            "FAILED_LOCI_FILE": str(failed_file),
            "COMPLETED_UTC": datetime.now(timezone.utc).isoformat(),
        }
        summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        if hard_failure.exists():
            hard_failure.unlink()

        banner("FINEMAPPING FINISHED")
        print(f"Loci total   : {len(loci)}")
        print(f"Loci success : {len(summary_frames)}")
        print(f"Loci failed  : {len(failures)}")
        print(f"Variants     : {final_variants}")
        print(f"Credible sets: {final_cs}")

        if failures:
            raise RuntimeError(
                f"Fine-mapping completed partially: {len(failures)} of {len(loci)} loci failed. "
                f"See {failed_file}"
            )

    except Exception as exc:
        hard_failure.write_text(
            f"STEP 06 FINEMAPPING FAILED\n\n{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}",
            encoding="utf-8",
        )
        raise


# =============================================================================
# DIRECT SINGLE-INDEX MODE
# =============================================================================


def project_manifest_path(args) -> Path:
    ancestry_code, ancestry_label = resolve_ancestry(args.ancestry)
    phenotype_slug = slugify(args.phenotype)
    ancestry_slug = slugify(ancestry_label)
    return (
        Path.cwd().resolve()
        / "06_finemapping"
        / phenotype_slug
        / ancestry_slug
        / "finemapping_manifest.tsv"
    )


def study_for_index(phenotype_slug: str, ancestry_slug: str, index: int) -> str:
    """Accession of study N in the Step02 download manifest (the ARRAY_ID the
    Step02 planner prints) - the same N in Step02, Step05, Step06 and Step07."""
    path = PROJECT_ROOT / "01_gwas_catalog" / phenotype_slug / ancestry_slug / "GWAS_download_manifest.tsv"
    if not path.is_file():
        raise SystemExit(f"No Step02 download manifest: {path}\nRun the Step02 planner first.")
    manifest = pd.read_csv(path, sep="\t", dtype=str)
    ids = pd.to_numeric(manifest["ARRAY_ID"], errors="coerce")
    selected = manifest[ids == index]
    if selected.empty:
        listing = "\n".join(f"  {int(i):>3}  {a}" for i, a in zip(ids, manifest["STUDY_ACCESSION"]) if pd.notna(i))
        raise SystemExit(f"--index {index} is not in {path}\nAvailable:\n{listing}")
    return str(selected.iloc[0]["STUDY_ACCESSION"]).strip()


def direct_index_mode(args) -> None:
    """--index N: study N of the Step02 download manifest, fine-mapped on this
    machine with the normal worker. Only this study has to be through Step05."""
    root = PROJECT_ROOT
    phenotype = args.phenotype.strip()
    ancestry_code, ancestry_label = resolve_ancestry(args.ancestry)
    p_slug, a_slug = slugify(phenotype), slugify(ancestry_label)
    accession = study_for_index(p_slug, a_slug, args.index)
    quoted = shlex.quote(phenotype)

    clump_root = root / "05_ld_clumping" / p_slug / a_slug
    if not (clump_root / accession / "lead_variants.tsv").is_file():
        raise SystemExit(
            f"Step05 clumping of study {args.index} ({accession}) is not done yet. Run:\n"
            f"  envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype {quoted} "
            f"--ancestry {ancestry_code} --index {args.index}")

    pipeline_python, plink2, rscript = resolve_finemap_tools(root)
    verify_reference(root, ancestry_code)
    reference_n = count_reference_samples(root, ancestry_code)
    meta = read_study_metadata(study_metadata_file(root, p_slug, a_slug)).get(accession)
    output_root = root / "06_finemapping" / p_slug / a_slug
    output_root.mkdir(parents=True, exist_ok=True)
    qc_file = root / "02_summary_stats" / p_slug / a_slug / "qc" / accession / f"{accession}_GRCh38_QC.tsv.gz"

    row, missing, problem = build_finemap_row(
        args, accession=accession, qc_file=qc_file, clump_root=clump_root, output_root=output_root,
        meta=meta, reference_n=reference_n, rscript=rscript, plink2=plink2, phenotype=phenotype,
        ancestry_code=ancestry_code, ancestry_label=ancestry_label)
    if missing:
        raise SystemExit(f"Missing inputs for {accession}:\n  " + "\n  ".join(missing))
    if problem:
        raise SystemExit(f"{accession} is not fine-mapped: {problem}\n"
                         "(--include-unreviewed fine-maps it anyway; not recommended)")

    banner(f"STEP 06 - INDEX {args.index}: {accession} (SuSiE-RSS on this machine)")
    print(f"{row['PHENOTYPE_TYPE']}  N={row['GWAS_N']:,}  [{row['GWAS_N_SOURCE']}]")
    row["FINEMAP_TASK_ID"] = 1
    single = output_root / f".index_{args.index}_finemap_manifest.tsv"
    pd.DataFrame([row]).to_csv(single, sep="\t", index=False)

    saved = {k: os.environ.get(k) for k in ("SLURM_ARRAY_TASK_ID", "GWAS_FINEMAP_MANIFEST", "SLURM_CPUS_PER_TASK")}
    os.environ["SLURM_ARRAY_TASK_ID"] = "1"
    os.environ["GWAS_FINEMAP_MANIFEST"] = str(single)
    os.environ.setdefault("SLURM_CPUS_PER_TASK", str(args.cpus))
    try:
        worker_mode()
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        single.unlink(missing_ok=True)

    base = f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {quoted} --ancestry {ancestry_code}"
    print(f"\nResults: {base} --stats\nFigures: envs/pipeline/bin/python Step07_Plots.py --phenotype {quoted} "
          f"--ancestry {ancestry_code} --index {args.index}")


# =============================================================================
# --status / --stats / --plots
# =============================================================================

GRCH38_LENGTHS = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}


def project_paths(args) -> dict:
    ancestry_code, ancestry_label = resolve_ancestry(args.ancestry)
    phenotype_slug = slugify(args.phenotype.strip())
    ancestry_slug = slugify(ancestry_label)
    output_root = PROJECT_ROOT / "06_finemapping" / phenotype_slug / ancestry_slug
    return {
        "phenotype": args.phenotype.strip(),
        "code": ancestry_code,
        "label": ancestry_label,
        "p_slug": phenotype_slug,
        "a_slug": ancestry_slug,
        "output_root": output_root,
        "manifest": output_root / "finemapping_manifest.tsv",
        "excluded": output_root / "finemapping_excluded_studies.tsv",
        "clump_manifest": PROJECT_ROOT / "05_ld_clumping" / phenotype_slug / ancestry_slug / "clumping_manifest.tsv",
        "plots": output_root / "plots",
        "base": f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {shlex.quote(args.phenotype.strip())} "
                f"--ancestry {ancestry_code}",
    }


def finished_studies(P) -> list[str]:
    """Studies with Step06 output, including ones run only with --index."""
    if not P["output_root"].exists():
        return []
    return sorted(d.name for d in P["output_root"].iterdir() if (d / "finemapping_summary.json").exists())


def study_records(P) -> dict:
    """accession -> manifest row (or the study's finemapping_summary.json)."""
    records = {}
    if P["manifest"].is_file():
        for _, m in pd.read_csv(P["manifest"], sep="\t", dtype=str).iterrows():
            records[m["STUDY_ACCESSION"]] = m.to_dict()
    for accession in finished_studies(P):
        if accession not in records:
            records[accession] = json.loads(
                (P["output_root"] / accession / "finemapping_summary.json").read_text(encoding="utf-8"))
    return records


def finemap_state(output_dir: Path, accession: str) -> tuple[str, dict]:
    summary_json = output_dir / "finemapping_summary.json"
    if (output_dir / "FINEMAPPING_FAILED.txt").exists() and not summary_json.exists():
        return "FAILED", {}
    if summary_json.exists():
        info = json.loads(summary_json.read_text(encoding="utf-8"))
        if info.get("STATUS") == "NO_SIGNIFICANT_LOCI":
            return "NO_LOCI", info
        if int(info.get("N_LOCI_FAILED", 0)) > 0:
            return "PARTIAL", info
        return "DONE", info
    return "PENDING", {}


def status_mode(args) -> None:
    P = project_paths(args)
    banner(f"STEP 06 STATUS - {P['phenotype']} / {P['label']} ({P['code']})")

    planned = pd.read_csv(P["manifest"], sep="\t", dtype=str) if P["manifest"].is_file() else pd.DataFrame()
    excluded = pd.read_csv(P["excluded"], sep="\t", dtype=str) if P["excluded"].is_file() else pd.DataFrame()
    clumped = pd.read_csv(P["clump_manifest"], sep="\t", dtype=str) if P["clump_manifest"].is_file() else pd.DataFrame()

    accessions = list(dict.fromkeys(
        list(clumped.get("STUDY_ACCESSION", [])) + list(planned.get("STUDY_ACCESSION", []))
        + finished_studies(P)))
    if not accessions:
        print("No Step05 / Step06 manifest yet. Run Step05, then the Step06 planner.")
        return

    reasons = dict(zip(excluded.get("STUDY_ACCESSION", []), excluded.get("REASON", [])))
    in_plan = set(planned.get("STUDY_ACCESSION", []))
    counts = defaultdict(int)

    print(f"{'STUDY':<15}{'STATE':<11}{'LOCI':>6}{'OK':>5}{'FAIL':>6}{'WARN':>6}  DETAIL")
    for accession in accessions:
        out = P["output_root"] / accession
        state, info = finemap_state(out, accession)
        if state == "PENDING" and accession in reasons:
            state = "EXCLUDED"
        elif state == "PENDING" and accession not in in_plan and not (P["output_root"] / accession).exists():
            state = "NOT_PLANNED"
        counts[state] += 1

        loci = ok = fail = warn = ""
        detail = reasons.get(accession, "")[:70] if state == "EXCLUDED" else ""
        if info:
            loci, ok, fail = info.get("N_LOCI", ""), info.get("N_LOCI_SUCCESS", ""), info.get("N_LOCI_FAILED", "")
            summary = out / f"{accession}_finemapping_summary.tsv"
            if summary.exists() and summary.stat().st_size > 2:
                s = pd.read_csv(summary, sep="\t")
                if "LOCUS_QC" in s.columns:
                    warn = int((s["LOCUS_QC"] == "WARN").sum())
            detail = f"N={info.get('GWAS_N', '')} [{info.get('GWAS_N_SOURCE', '')}]"[:70]
        if state == "FAILED":
            detail = (out / "FINEMAPPING_FAILED.txt").read_text(encoding="utf-8", errors="replace").splitlines()[2][:70]
        print(f"{accession:<15}{state:<11}{loci!s:>6}{ok!s:>5}{fail!s:>6}{warn!s:>6}  {detail}")

    print()
    print("Studies : " + "  ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print()
    if counts["PENDING"] or counts["FAILED"] or counts["PARTIAL"] or counts["NOT_PLANNED"]:
        print("Plan / re-run (finished loci are reused):")
        print(f"  {P['base']}")
        print(f"  sbatch Step06_FineMap_GWAS_{P['p_slug']}_{P['a_slug']}.sh")
    if counts["EXCLUDED"]:
        print("Excluded studies need a verified phenotype type (Step02 phenotype_type_review.tsv),")
        print("then re-run Step02 --qc and this planner.")
    if counts["DONE"] or counts["PARTIAL"]:
        print(f"Results: {P['base']} --stats   |   plots: {P['base']} --plots")


def stats_mode(args) -> None:
    P = project_paths(args)
    banner(f"STEP 06 STATISTICS - {P['phenotype']} / {P['label']} ({P['code']})")
    rows = []
    for accession, m in study_records(P).items():
        out = P["output_root"] / accession
        state, info = finemap_state(out, accession)
        if state in ("PENDING", "FAILED"):
            continue

        summary_file = out / f"{accession}_finemapping_summary.tsv"
        s = pd.read_csv(summary_file, sep="\t") if summary_file.stat().st_size > 2 else pd.DataFrame()
        cs_file = out / f"{accession}_95pct_credible_sets.tsv"
        cs = pd.read_csv(cs_file, sep="\t") if cs_file.exists() and cs_file.stat().st_size > 2 else pd.DataFrame()
        fm_file = out / f"{accession}_finemapped_variants.tsv.gz"
        fm = pd.read_csv(fm_file, sep="\t", low_memory=False) if fm_file.exists() else pd.DataFrame()

        cs_sizes = cs.groupby("CREDIBLE_SET").size() if not cs.empty else pd.Series(dtype=int)
        row = {
            "STUDY_ACCESSION": accession,
            "PHENOTYPE_TYPE": m.get("PHENOTYPE_TYPE", ""),
            "N_CASES": m.get("N_CASES", ""), "N_CONTROLS": m.get("N_CONTROLS", ""),
            "N_TOTAL": m.get("N_TOTAL", ""), "N_EFF": m.get("N_EFF", ""),
            "GWAS_N_USED": m.get("GWAS_N", ""), "GWAS_N_SOURCE": m.get("GWAS_N_SOURCE", ""),
            "LD_REFERENCE": f"1000G {P['code']} (n={m.get('REFERENCE_N', '')})",
            "N_LOCI": info.get("N_LOCI", 0),
            "N_LOCI_SUCCESS": info.get("N_LOCI_SUCCESS", 0),
            "N_LOCI_FAILED": info.get("N_LOCI_FAILED", 0),
            "N_LOCI_QC_PASS": int((s.get("LOCUS_QC", pd.Series(dtype=str)) == "PASS").sum()),
            "N_LOCI_QC_WARN": int((s.get("LOCUS_QC", pd.Series(dtype=str)) == "WARN").sum()),
            "N_LOCI_NOT_CONVERGED": int((s.get("CONVERGED", pd.Series(dtype=bool)).astype(str) == "False").sum()),
            "MEDIAN_ESTIMATE_S_RSS": round(float(pd.to_numeric(s.get("ESTIMATE_S_RSS"), errors="coerce").median()), 4)
            if "ESTIMATE_S_RSS" in s else "",
            "N_POSSIBLE_ALLELE_FLIPS": int(pd.to_numeric(s.get("N_POSSIBLE_ALLELE_FLIPS"), errors="coerce").fillna(0).sum())
            if "N_POSSIBLE_ALLELE_FLIPS" in s else "",
            "N_TOP_AGREES_L1": int((s.get("TOP_AGREES_L1", pd.Series(dtype=str)).astype(str) == "True").sum()),
            "N_CREDIBLE_SETS": int(cs_sizes.size),
            "MEDIAN_CS_SIZE": float(cs_sizes.median()) if cs_sizes.size else "",
            "N_CS_SINGLE_VARIANT": int((cs_sizes == 1).sum()),
            "N_CS_SIZE_LE_5": int((cs_sizes <= 5).sum()),
            "N_VARIANTS_PIP_GT_0.5": int((fm.get("PIP", pd.Series(dtype=float)) > 0.5).sum()),
            "N_VARIANTS_PIP_GT_0.9": int((fm.get("PIP", pd.Series(dtype=float)) > 0.9).sum()),
            "N_HIGH_PIP_OUTSIDE_CS": int(pd.to_numeric(s.get("N_HIGH_PIP_OUTSIDE_CS"), errors="coerce").fillna(0).sum())
            if "N_HIGH_PIP_OUTSIDE_CS" in s else "",
        }
        rows.append(row)

        print()
        print(f"{accession}: {row['PHENOTYPE_TYPE']} | N used {row['GWAS_N_USED']} [{row['GWAS_N_SOURCE']}]")
        print(f"   loci {row['N_LOCI']} (ok {row['N_LOCI_SUCCESS']}, failed {row['N_LOCI_FAILED']}) | "
              f"QC PASS {row['N_LOCI_QC_PASS']} / WARN {row['N_LOCI_QC_WARN']} | not converged "
              f"{row['N_LOCI_NOT_CONVERGED']} | median s {row['MEDIAN_ESTIMATE_S_RSS']} | allele-flip flags "
              f"{row['N_POSSIBLE_ALLELE_FLIPS']}")
        print(f"   credible sets {row['N_CREDIBLE_SETS']} (median size {row['MEDIAN_CS_SIZE']}, single-variant "
              f"{row['N_CS_SINGLE_VARIANT']}) | PIP>0.5 {row['N_VARIANTS_PIP_GT_0.5']} | PIP>0.9 "
              f"{row['N_VARIANTS_PIP_GT_0.9']} | PIP>0.5 outside CS (artefact flag) {row['N_HIGH_PIP_OUTSIDE_CS']} | "
              f"top variant = L1 top in {row['N_TOP_AGREES_L1']} loci")
        if not fm.empty and "PIP" in fm:
            id_col = "SNPID" if "SNPID" in fm.columns else "REFERENCE_ID"
            top = fm.sort_values("PIP", ascending=False).head(5)
            print("   top PIP: " + "; ".join(
                f"{r[id_col]} ({r['LOCUS_ID']}) {r['PIP']:.2f}" for _, r in top.iterrows()))
        if not s.empty and "LOCUS_QC_NOTES" in s:
            for _, r in s[s["LOCUS_QC"] == "WARN"].head(10).iterrows():
                print(f"   WARN {r['LOCUS_ID']} chr{r['CHR']}:{r['LOCUS_START']}-{r['LOCUS_END']}: {r['LOCUS_QC_NOTES']}")

    if not rows:
        print("No finished studies yet (see --status).")
        return

    output = P["output_root"] / "finemapping_stats.tsv"
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
    print(f"\nTable: {output}")


# ---- plots ----------------------------------------------------------------------

def gene_table() -> pd.DataFrame:
    """Protein-coding + lncRNA genes from the Step00 GENCODE GTF (cached as TSV)."""
    gtfs = sorted((PROJECT_ROOT / "resources" / "gencode").glob("*/gencode.v*.annotation.gtf.gz"))
    if not gtfs:
        return pd.DataFrame(columns=["CHR", "START", "END", "STRAND", "GENE", "TYPE"])
    gtf = gtfs[-1]
    cache = gtf.with_name(gtf.name.replace(".gtf.gz", ".genes.tsv"))
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    rows = []
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = clean_chr(f[0])
            gtype = re.search(r'gene_type "([^"]+)"', f[8])
            name = re.search(r'gene_name "([^"]+)"', f[8])
            if chrom and gtype and name and gtype.group(1) in ("protein_coding", "lncRNA"):
                rows.append((chrom, int(f[3]), int(f[4]), f[6], name.group(1), gtype.group(1)))
    table = pd.DataFrame(rows, columns=["CHR", "START", "END", "STRAND", "GENE", "TYPE"])
    table.to_csv(cache, sep="\t", index=False)
    return table


def lead_r(locus_dir: Path, ancestry_code: str, lead_id: str):
    """Signed r of every LD variant with lead_id (memory-mapped row of the PLINK matrix)."""
    ld_file = locus_dir / f"{ancestry_code}_LD.unphased.vcor1.bin"
    vars_file = locus_dir / f"{ancestry_code}_LD.unphased.vcor1.bin.vars"
    if not ld_file.exists() or not vars_file.exists():
        return None, None
    ids = [x.strip() for x in vars_file.read_text().splitlines() if x.strip()]
    if lead_id not in ids:
        return ids, None
    matrix = np.memmap(ld_file, dtype="<f4", mode="r", shape=(len(ids), len(ids)))
    return ids, np.asarray(matrix[ids.index(lead_id)], dtype=float)


R2_BINS = [(0.8, "#d7191c", "0.8-1"), (0.6, "#fdae61", "0.6-0.8"), (0.4, "#4dac26", "0.4-0.6"),
           (0.2, "#8fbfe0", "0.2-0.4"), (0.0, "#2c3e91", "0-0.2")]


def plot_locus(plt, P, accession: str, locus: pd.Series, genes: pd.DataFrame, out_dir: Path) -> list[Path]:
    study_dir = P["output_root"] / accession
    locus_id = str(locus["LOCUS_ID"])
    locus_dir = study_dir / "loci" / locus_id
    variants_file = locus_dir / "susie" / "susie_variants.tsv"
    if not variants_file.exists():
        return []

    fm = pd.read_csv(variants_file, sep="\t", low_memory=False)
    regional_file = study_dir / "regional_gwas" / f"{locus_id}.tsv.gz"
    regional = pd.read_csv(regional_file, sep="\t") if regional_file.exists() else pd.DataFrame()
    summary_file = locus_dir / "susie" / "susie_summary.tsv"
    summary = pd.read_csv(summary_file, sep="\t").iloc[0] if summary_file.exists() else pd.Series(dtype=object)

    chrom, start, end = int(locus["CHR"]), int(locus["LOCUS_START"]), int(locus["LOCUS_END"])
    fm["LOGP"] = -np.log10(pd.to_numeric(fm["P"], errors="coerce").clip(lower=1e-300))
    index_variant = fm.loc[fm["LOGP"].idxmax()]
    ids, r = lead_r(locus_dir, P["code"], str(index_variant["REFERENCE_ID"]))
    fm["R2"] = np.nan
    if r is not None:
        fm["R2"] = fm["REFERENCE_ID"].map(dict(zip(ids, r ** 2)))

    label_col = "SNPID" if "SNPID" in fm.columns else "REFERENCE_ID"
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True,
                             gridspec_kw={"height_ratios": [3, 2, 1.3]})

    # (a) association, coloured by r2 with the index variant
    ax = axes[0]
    if not regional.empty:
        in_ld = set(fm["POS"].astype(int))
        other = regional[~regional["POS"].astype(int).isin(in_ld)]
        ax.scatter(other["POS"], -np.log10(pd.to_numeric(other["P"], errors="coerce").clip(lower=1e-300)),
                   s=6, color="#c8c8c8", label="not in LD panel / filtered", zorder=1)
    for low, color, label in R2_BINS:
        upper = 1.01 if low == 0.8 else low + 0.2
        sel = fm[(fm["R2"] >= low) & (fm["R2"] < upper)]
        ax.scatter(sel["POS"], sel["LOGP"], s=12, color=color, label=f"r2 {label}", zorder=2)
    na = fm[fm["R2"].isna()]
    ax.scatter(na["POS"], na["LOGP"], s=10, color="#7f7f7f", zorder=2)
    ax.scatter([index_variant["POS"]], [index_variant["LOGP"]], s=90, marker="D", color="#7b1fa2",
               edgecolors="black", zorder=3, label=f"index {index_variant[label_col]}")
    ax.axhline(-np.log10(5e-8), color="#c0392b", lw=0.8, ls="--")
    ax.set_ylabel(r"$-\log_{10}(P)$")
    top = max(float(fm["LOGP"].max()), -np.log10(5e-8))
    ax.set_ylim(0, top * 1.35)  # headroom: the legend never hides a variant
    ax.legend(fontsize=7, ncol=7, loc="upper left", frameon=False, handletextpad=0.2, columnspacing=0.8)
    qc = f"{summary.get('LOCUS_QC', '')} {summary.get('LOCUS_QC_NOTES', '')}".strip()
    s_value = pd.to_numeric(pd.Series([summary.get("ESTIMATE_S_RSS")]), errors="coerce").iloc[0]
    s_text = f"{s_value:.2f}" if pd.notna(s_value) else "not computed"
    ax.set_title(f"{accession} {locus_id} chr{chrom}:{start:,}-{end:,}   N={summary.get('GWAS_N', '')}   "
                 f"LD/z mismatch s={s_text}   index r2 vs {index_variant[label_col]}\n{qc[:140]}", fontsize=9)

    # (b) PIP, credible sets coloured
    ax = axes[1]
    ax.scatter(fm["POS"], fm["PIP"], s=10, color="#9e9e9e", zorder=1)
    palette = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#17becf"]
    sets = [x for x in dict.fromkeys(";".join(fm["CREDIBLE_SET"].fillna("").astype(str)).split(";")) if x]
    for i, cs in enumerate(sets):
        sel = fm[fm["CREDIBLE_SET"].fillna("").astype(str).str.split(";").apply(lambda v: cs in v)]
        ax.scatter(sel["POS"], sel["PIP"], s=26, color=palette[i % len(palette)], zorder=2,
                   label=f"{cs.split('_')[-1]} ({len(sel)} var)")
        best = sel.loc[sel["PIP"].idxmax()]
        ax.annotate(str(best[label_col]), (best["POS"], best["PIP"]), fontsize=7, xytext=(3, 3),
                    textcoords="offset points")
    if "PIP_L1" in fm:
        ax.scatter(fm["POS"], fm["PIP_L1"], s=8, marker="x", color="#555555", alpha=0.5, zorder=1,
                   label="PIP L=1 (LD-free)")
    ax.set_ylim(-0.03, 1.3)
    ax.set_ylabel("PIP")
    ax.legend(fontsize=7, ncol=6, loc="upper left", frameon=False, handletextpad=0.2, columnspacing=0.8)

    # (c) genes
    ax = axes[2]
    window = genes[(genes["CHR"] == chrom) & (genes["END"] >= start) & (genes["START"] <= end)] \
        if not genes.empty else genes
    window = window.sort_values("START")
    rows_end = []
    for _, g in window.iterrows():
        row = next((i for i, e in enumerate(rows_end) if g["START"] > e + (end - start) * 0.03), None)
        if row is None:
            rows_end.append(g["END"])
            row = len(rows_end) - 1
        else:
            rows_end[row] = g["END"]
        y = -row
        color = "#1f4e79" if g["TYPE"] == "protein_coding" else "#9e9e9e"
        ax.plot([max(g["START"], start), min(g["END"], end)], [y, y], color=color, lw=3,
                solid_capstyle="butt")
        if g["TYPE"] == "protein_coding":
            ax.text(max(g["START"], start), y + 0.25, g["GENE"] + ("" if g["STRAND"] == "+" else " <"),
                    fontsize=6, style="italic")
    ax.set_ylim(-max(len(rows_end), 1) - 0.2, 0.9)
    ax.set_yticks([])
    ax.set_ylabel("genes" if not window.empty else "genes\n(no GTF)")
    ax.set_xlim(start, end)
    ax.set_xlabel(f"chr{chrom} position (Mb, GRCh38)")
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1e6:.2f}"))
    fig.tight_layout()
    files = [out_dir / f"{locus_id}_regional.png"]
    fig.savefig(files[0], dpi=150)
    plt.close(fig)

    # LD heatmap of up to 250 variants (credible sets + strongest associations), by position
    if ids is not None:
        chosen = pd.concat([fm[fm["IN_95_CREDIBLE_SET"].astype(str) == "True"],
                            fm.sort_values("LOGP", ascending=False).head(250)]).drop_duplicates("REFERENCE_ID")
        chosen = chosen.head(250).sort_values("POS")
        index = [ids.index(x) for x in chosen["REFERENCE_ID"] if x in ids]
        if len(index) >= 2:
            matrix = np.memmap(locus_dir / f"{P['code']}_LD.unphased.vcor1.bin", dtype="<f4", mode="r",
                               shape=(len(ids), len(ids)))
            sub = np.asarray(matrix[np.ix_(index, index)], dtype=float)
            fig, ax = plt.subplots(figsize=(7, 6))
            image = ax.imshow(sub, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
            fig.colorbar(image, ax=ax, label="signed r (REF alleles, 1000G)")
            ax.set_title(f"{accession} {locus_id}: LD of {len(index)} variants (credible sets + top P)",
                         fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            fig.tight_layout()
            files.append(out_dir / f"{locus_id}_ld.png")
            fig.savefig(files[-1], dpi=150)
            plt.close(fig)
    return files


def plot_genome_pip(plt, P, accession: str, out_dir: Path) -> Path | None:
    fm_file = P["output_root"] / accession / f"{accession}_finemapped_variants.tsv.gz"
    if not fm_file.exists():
        return None
    fm = pd.read_csv(fm_file, sep="\t", low_memory=False)
    if fm.empty or "PIP" not in fm:
        return None
    offsets, start = {}, 0
    for chrom in range(1, 23):
        offsets[chrom] = start
        start += GRCH38_LENGTHS[chrom]
    fm = fm.dropna(subset=["CHR", "POS"])
    x = fm["POS"].astype(int) + fm["CHR"].astype(int).map(offsets)
    in_cs = fm["IN_95_CREDIBLE_SET"].astype(str) == "True"
    fig, ax = plt.subplots(figsize=(16, 4.5))
    ax.scatter(x[~in_cs], fm.loc[~in_cs, "PIP"], s=4, color="#b0bec5")
    ax.scatter(x[in_cs], fm.loc[in_cs, "PIP"], s=12, color="#c0392b", label="in 95% credible set")
    label_col = "SNPID" if "SNPID" in fm.columns else "REFERENCE_ID"
    for _, r in fm[fm["PIP"] > 0.5].iterrows():
        ax.annotate(str(r[label_col]), (int(r["POS"]) + offsets[int(r["CHR"])], r["PIP"]), fontsize=6,
                    rotation=30, xytext=(2, 2), textcoords="offset points")
    ax.set_xticks([offsets[c] + GRCH38_LENGTHS[c] / 2 for c in range(1, 23)])
    ax.set_xticklabels([str(c) for c in range(1, 23)], fontsize=8)
    ax.set_xlim(0, start)
    ax.set_ylim(-0.02, 1.05)
    ax.set_ylabel("PIP (SuSiE-RSS)")
    ax.set_xlabel("Chromosome (GRCh38)")
    ax.set_title(f"{accession} - fine-mapped variants across the genome")
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    path = out_dir / f"{accession}_pip_genome.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plots_mode(args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    P = project_paths(args)
    banner(f"STEP 06 PLOTS - {P['phenotype']} / {P['label']} ({P['code']})")
    genes = gene_table()
    if genes.empty:
        print("No GENCODE GTF under resources/gencode/ - regional plots have no gene track "
              "(Step00_Download_Resources.py --groups genome).")

    for accession in study_records(P):
        state, _ = finemap_state(P["output_root"] / accession, accession)
        if state not in ("DONE", "PARTIAL"):
            print(f"{accession}: {state} - skipped")
            continue
        out_dir = P["plots"] / accession
        out_dir.mkdir(parents=True, exist_ok=True)
        loci = pd.read_csv(P["output_root"] / accession / "loci_definition.tsv", sep="\t")
        made = []
        for _, locus in loci.iterrows():
            made += plot_locus(plt, P, accession, locus, genes, out_dir)
        genome = plot_genome_pip(plt, P, accession, out_dir)
        if genome:
            made.append(genome)
        print(f"{accession}: {len(made)} plots in {out_dir}")

    print("\nGenome-wide Manhattan / QQ / locus maps: Step05_Clump_GWAS_Loci.py ... --plots")


# =============================================================================
# MAIN
# =============================================================================


def main():
    # The generated batch repeats --phenotype/--ancestry for provenance; the
    # worker itself uses the parameters stored in the manifest.
    args = arguments()

    if args.status:
        status_mode(args)
    elif args.stats:
        stats_mode(args)
    elif args.plots:
        plots_mode(args)
    elif args.index is not None:
        direct_index_mode(args)
    elif os.environ.get("SLURM_ARRAY_TASK_ID") and os.environ.get("GWAS_FINEMAP_MANIFEST"):
        worker_mode()
    else:
        planner_mode(args)


if __name__ == "__main__":
    main()
