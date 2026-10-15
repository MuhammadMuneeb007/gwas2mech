#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 05 - Ancestry-matched LD clumping / independent lead variants
===============================================================================

    # plan + SLURM array (one task per QC'd GWAS; >1000 studies share tasks)
    envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype "Parkinson's disease" --ancestry EUR
    sbatch Step05_Clump_GWAS_parkinson_s_disease_european.sh

    # progress (read-only), per-study statistics, plots
    envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype "Parkinson's disease" --ancestry EUR --stats
    envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype "Parkinson's disease" --ancestry EUR --plots

    # study 2 (the number the Step02 planner prints) on this machine; only its
    # Step02 QC has to be done
    envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype "Parkinson's disease" --ancestry EUR --index 2

Inputs
    Step02 QC files   02_summary_stats/<phenotype>/<ancestry>/qc/<GCST>/<GCST>_GRCh38_QC.tsv.gz
                      (+ <GCST>_QC_summary.tsv: CLUMPING_READY)
    Step04 panels     resources/1000G/<ANC>/chr<N>_<ANC>_GRCh38.pgen/.pvar/.psam

For each GWAS the worker
  1. streams the QC file once and keeps variants with P <= p2,
  2. maps them to the 1000G panel by chromosome + position + BOTH alleles
     (either orientation). A variant whose alleles do not match the panel
     variant at that position is NOT used (never position-only),
  3. runs PLINK2 --clump per chromosome (p1 5e-8, p2 1e-2, r2 0.1, 1000 kb),
  4. writes the independent lead variants used by Step06.

Outputs
    05_ld_clumping/<phenotype>/<ancestry>/
        clumping_manifest.tsv, clumping_excluded_studies.tsv, clumping_stats.tsv
        plots/<GCST>_manhattan.png, <GCST>_qq.png, <GCST>_loci_map.png   (--plots)
        <GCST>/lead_variants.tsv, mapped_candidates.tsv.gz,
               chromosome_clumping_summary.tsv, clumping_summary.tsv/.json,
               chr01..chr22/ (clump_input.tsv, plink.clumps, plink.log)

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
import subprocess
import sys
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULTS
# =============================================================================

PROJECT_ROOT = Path.cwd().resolve()
MAX_ARRAY_TASKS = 1000
DEFAULT_CPUS = 4
DEFAULT_P1 = 5e-8
DEFAULT_P2 = 1e-2
DEFAULT_R2 = 0.1
DEFAULT_KB = 1000
DEFAULT_CHUNK_SIZE = 500_000

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
# DISPLAY / TEXT HELPERS
# =============================================================================


def banner(text: str) -> None:
    print()
    print("=" * 88)
    print(text)
    print("=" * 88)


def normalize_text(value) -> str:
    value = str(value).lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    value = re.sub(r"\s+", " ", value)
    return value.strip()


def slugify(value) -> str:
    return normalize_text(value).replace(" ", "_")


def canonical_ancestry(value: str) -> tuple[str, str]:
    key = normalize_text(value)
    if key not in ANCESTRY_ALIASES:
        valid = "EUR, AFR, EAS, SAS, AMR"
        raise ValueError(
            f"Unsupported ancestry: {value!r}. "
            f"Use one of {valid} (or the corresponding full label)."
        )
    return ANCESTRY_ALIASES[key]


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


def add_slurm_arguments(parser):
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


def resolve_slurm(args):
    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": str(DEFAULT_CPUS),
              "mem": "16G", "time": "12:00:00"}
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


def find_column(columns, candidates):
    lookup = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def clean_chr_series(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip().str.replace(
        r"^chr", "", regex=True, case=False
    )
    numeric = pd.to_numeric(text, errors="coerce")
    return numeric.astype("Int64")


def clean_allele(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    x = str(value).strip().upper()
    if not x or x in {"NA", "NAN", "<NA>", "NONE"}:
        return None
    return x


# =============================================================================
# CLI
# =============================================================================


def arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Ancestry-matched LD clumping of Step02 QC'd GWAS (PLINK2)."
    )

    parser.add_argument("--phenotype", required=True)
    parser.add_argument(
        "--ancestry",
        required=True,
        help="EUR, AFR, EAS, SAS, AMR or corresponding full ancestry label.",
    )

    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--index", type=int, default=None,
                       help="Study N of the Step02 download manifest (the number Step02 prints): "
                            "clump it on this machine; other studies need not be ready.")
    modes.add_argument("--status", action="store_true",
                       help="READ-ONLY: clumping done / failed / pending per study.")
    modes.add_argument("--stats", action="store_true",
                       help="Per-study lead variants, mapping rate, unmatched variants (clumping_stats.tsv).")
    modes.add_argument("--plots", action="store_true",
                       help="Manhattan, QQ and chromosome locus map per study (plots/).")

    parser.add_argument("--p1", type=float, default=DEFAULT_P1)
    parser.add_argument("--p2", type=float, default=DEFAULT_P2)
    parser.add_argument("--r2", type=float, default=DEFAULT_R2)
    parser.add_argument("--kb", type=int, default=DEFAULT_KB)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Allow clumping to proceed when some expected Step02 QC files are still missing.",
    )

    add_slurm_arguments(parser)

    args = parser.parse_args(argv)
    return resolve_slurm(args)


def validate_arguments(args) -> None:
    if not args.phenotype.strip():
        raise ValueError("--phenotype cannot be empty.")

    if args.cpus < 1:
        raise ValueError("--cpus must be >= 1")
    if args.max_parallel < 0:
        raise ValueError("--max-parallel must be >= 0 (0 = no throttle)")
    if args.index is not None and args.index < 1:
        raise ValueError("--index must be >= 1")
    if args.chunk_size < 10_000:
        raise ValueError("--chunk-size must be >= 10000")

    if not (0 < args.p1 <= 1):
        raise ValueError("--p1 must be in (0, 1]")
    if not (0 < args.p2 <= 1):
        raise ValueError("--p2 must be in (0, 1]")
    if args.p1 > args.p2:
        raise ValueError("--p1 must be <= --p2")
    if not (0 < args.r2 <= 1):
        raise ValueError("--r2 must be in (0, 1]")
    if args.kb < 1:
        raise ValueError("--kb must be >= 1")


# =============================================================================
# REFERENCE VALIDATION
# =============================================================================


def reference_prefix(root: Path, ancestry_code: str, chrom: int) -> Path:
    return (
        root
        / "resources"
        / "1000G"
        / ancestry_code
        / f"chr{chrom}_{ancestry_code}_GRCh38"
    )


def validate_reference(root: Path, ancestry_code: str) -> None:
    missing = []
    for chrom in range(1, 23):
        prefix = reference_prefix(root, ancestry_code, chrom)
        for suffix in (".pgen", ".pvar", ".psam"):
            path = Path(str(prefix) + suffix)
            if not path.exists() or path.stat().st_size <= 0:
                missing.append(str(path))

    if missing:
        preview = "\n".join(f"  {x}" for x in missing[:20])
        if len(missing) > 20:
            preview += f"\n  ... and {len(missing) - 20} more"
        raise RuntimeError(
            f"1000G {ancestry_code} reference is incomplete. Missing:\n{preview}\n\n"
            f"Build it with:\n  envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry {ancestry_code}"
        )


# =============================================================================
# ONE CLUMPING MANIFEST ROW
# =============================================================================


def build_clump_row(
    args,
    task_id: int,
    accession: str,
    qc_file: Path,
    output_root: Path,
    phenotype: str,
    ancestry_code: str,
    ancestry_label: str,
) -> tuple[dict | None, dict | None]:
    """Return (manifest_row, None) or (None, exclusion_record)."""
    summary_file = qc_file.parent / f"{accession}_QC_summary.tsv"

    clumping_ready = True
    qc_status = "UNKNOWN"
    reason = ""

    if summary_file.exists() and summary_file.stat().st_size > 0:
        try:
            summary_df = pd.read_csv(summary_file, sep="\t", dtype=str)
            if not summary_df.empty:
                summary_row = summary_df.iloc[0]
                qc_status = str(summary_row.get("QC_STATUS", "UNKNOWN"))
                ready_text = str(summary_row.get("CLUMPING_READY", "True")).strip().lower()
                clumping_ready = ready_text in {"true", "1", "yes", "y"}
        except Exception as exc:
            reason = f"Could not read QC summary: {exc}"

    if not clumping_ready:
        return None, {
            "STUDY_ACCESSION": accession,
            "QC_FILE": str(qc_file.resolve()),
            "QC_STATUS": qc_status,
            "REASON": reason or "Step02 marked CLUMPING_READY=False",
        }

    study_output = output_root / accession
    return {
        "CLUMP_TASK_ID": task_id,
        "STUDY_ACCESSION": accession,
        "QC_FILE": str(qc_file.resolve()),
        "QC_SUMMARY_FILE": str(summary_file.resolve()) if summary_file.exists() else "",
        "QC_STATUS": qc_status,
        "OUTPUT_DIR": str(study_output.resolve()),
        "PHENOTYPE": phenotype,
        "ANCESTRY": ancestry_label,
        "ANCESTRY_CODE": ancestry_code,
        "ANCESTRY_LABEL": ancestry_label,
        "P1": args.p1,
        "P2": args.p2,
        "R2": args.r2,
        "KB": args.kb,
        "CHUNK_SIZE": args.chunk_size,
    }, None


# =============================================================================
# PLANNER
# =============================================================================


def planner_mode(args) -> None:
    validate_arguments(args)

    root = Path.cwd().resolve()
    phenotype = args.phenotype.strip()
    phenotype_slug = slugify(phenotype)
    ancestry_code, ancestry_label = canonical_ancestry(args.ancestry)
    ancestry_slug = slugify(ancestry_label)

    pipeline_python = root / "envs" / "pipeline" / "bin" / "python"
    plink2 = root / "envs" / "pipeline" / "bin" / "plink2"
    script_path = Path(__file__).resolve()

    if not pipeline_python.exists():
        raise FileNotFoundError(f"Pipeline Python not found: {pipeline_python}")
    if not plink2.exists():
        raise FileNotFoundError(f"PLINK2 not found: {plink2}")

    validate_reference(root, ancestry_code)

    qc_root = (
        root
        / "02_summary_stats"
        / phenotype_slug
        / ancestry_slug
        / "qc"
    )

    if not qc_root.exists():
        raise FileNotFoundError(
            "Step02 QC directory not found:\n"
            f"  {qc_root}\n\n"
            "Run Step02_Plan_QC_GWAS.py --qc for this phenotype/ancestry first."
        )

    # Resolve expected studies from the Step02 QC manifest when available.
    qc_manifest_file = (
        root
        / "01_gwas_catalog"
        / phenotype_slug
        / ancestry_slug
        / "GWAS_QC_manifest.tsv"
    )

    expected_accessions = []
    if qc_manifest_file.exists() and qc_manifest_file.stat().st_size > 0:
        qc_manifest = pd.read_csv(
            qc_manifest_file, sep="\t", dtype=str, low_memory=False
        )
        if "STUDY_ACCESSION" in qc_manifest.columns:
            expected_accessions = [
                str(x).strip()
                for x in qc_manifest["STUDY_ACCESSION"].dropna().tolist()
                if str(x).strip()
            ]

    qc_files = sorted(qc_root.glob("*/*_GRCh38_QC.tsv.gz"))
    by_accession = {path.parent.name: path for path in qc_files}

    if not qc_files:
        raise RuntimeError(
            "No completed Step02 QC files were found under:\n"
            f"  {qc_root}"
        )

    missing_qc_accessions = [
        accession for accession in expected_accessions if accession not in by_accession
    ]
    if missing_qc_accessions and not args.allow_partial:
        preview = ", ".join(missing_qc_accessions[:20])
        if len(missing_qc_accessions) > 20:
            preview += f", ... (+{len(missing_qc_accessions)-20} more)"
        raise RuntimeError(
            "Step02 QC is not complete for all expected studies.\n"
            f"Missing QC outputs for {len(missing_qc_accessions)} accession(s): {preview}\n"
            "Wait for Step02 to finish, or use --allow-partial intentionally."
        )

    output_root = (
        root
        / "05_ld_clumping"
        / phenotype_slug
        / ancestry_slug
    )
    log_root = (
        root
        / "logs"
        / "step05_clumping"
        / phenotype_slug
        / ancestry_slug
    )

    output_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    rows = []
    exclusions = []

    for qc_file in qc_files:
        accession = qc_file.parent.name

        row, exclusion = build_clump_row(
            args,
            task_id=len(rows) + 1,
            accession=accession,
            qc_file=qc_file,
            output_root=output_root,
            phenotype=phenotype,
            ancestry_code=ancestry_code,
            ancestry_label=ancestry_label,
        )

        if exclusion is not None:
            exclusions.append(exclusion)
            continue

        rows.append(row)

    if exclusions:
        pd.DataFrame(exclusions).to_csv(
            output_root / "clumping_excluded_studies.tsv", sep="\t", index=False
        )

    if not rows:
        raise RuntimeError(
            "No Step02 studies are marked ready for clumping. "
            "Inspect the QC summaries and clumping_excluded_studies.tsv."
        )

    manifest = pd.DataFrame(rows)
    task_ids, n_tasks = assign_task_ids(len(manifest))
    manifest["CLUMP_TASK_ID"] = task_ids
    manifest_file = (output_root / "clumping_manifest.tsv").resolve()
    manifest.to_csv(manifest_file, sep="\t", index=False)

    bash_file = (
        root
        / f"Step05_Clump_GWAS_{phenotype_slug}_{ancestry_slug}.sh"
    )

    job_name = f"CLUMP_{phenotype_slug}_{ancestry_code.lower()}"[:100]
    array_spec = f"1-{n_tasks}" + (f"%{args.max_parallel}" if args.max_parallel else "")

    sbatch_header = sbatch_header_block(
        args,
        job_name=job_name,
        output=f"{log_root}/clump.%A_%a.out",
        error=f"{log_root}/clump.%A_%a.err",
        n_tasks=n_tasks,
    )

    bash_text = f'''#!/bin/bash --login
{sbatch_header}

set -euo pipefail

cd {shlex.quote(str(root))}
export PATH={shlex.quote(str(root / "envs" / "pipeline" / "bin"))}:"$PATH"

export GWAS_CLUMP_MANIFEST={shlex.quote(str(manifest_file))}

printf '%s\n' "=============================================================="
printf '%s\n' "STEP 05 - LD CLUMPING"
printf '%s\n' "=============================================================="
printf 'Job ID        : %s\n' "$SLURM_JOB_ID"
printf 'Array task ID : %s\n' "$SLURM_ARRAY_TASK_ID"
printf 'Hostname      : %s\n' "$(hostname)"
printf 'Python        : %s\n' {shlex.quote(str(pipeline_python))}
printf 'Phenotype     : %s\n' {shlex.quote(phenotype)}
printf 'Ancestry      : %s (%s)\n' {shlex.quote(ancestry_label)} {shlex.quote(ancestry_code)}
printf '%s\n' "=============================================================="

{shlex.quote(str(pipeline_python))} {shlex.quote(str(script_path))}

printf '%s\n' "=============================================================="
printf 'STEP 05 TASK COMPLETE: %s\n' "$SLURM_ARRAY_TASK_ID"
printf '%s\n' "=============================================================="
'''

    bash_file.write_text(bash_text, encoding="utf-8")
    bash_file.chmod(0o755)

    banner("STEP 05 - LD CLUMPING PLANNER")
    print(f"Phenotype          : {phenotype}")
    print(f"Ancestry           : {ancestry_label}")
    print(f"Ancestry code      : {ancestry_code}")
    print(f"QC GWAS files      : {len(qc_files)}")
    print(f"Clumping-ready     : {len(manifest)}")
    print(f"Excluded           : {len(exclusions)}")
    print(f"Missing QC outputs : {len(missing_qc_accessions)}")
    print(f"Reference          : resources/1000G/{ancestry_code}/")
    print(f"Partition          : {args.partition or 'scheduler default'}")
    print(f"Walltime           : {args.time}")
    print(f"Memory/task        : {args.mem}")
    print(f"CPUs/task          : {args.cpus}")
    print(f"Array              : {array_spec}")
    print()
    print("Clumping parameters:")
    print(f"  p1               : {args.p1}")
    print(f"  p2               : {args.p2}")
    print(f"  r2               : {args.r2}")
    print(f"  kb               : {args.kb}")
    print()
    print(f"Manifest:\n  {manifest_file}")
    print(f"Output root:\n  {output_root}")
    print(f"Generated SLURM script:\n  {bash_file}")
    print()
    print("NEXT COMMAND:")
    print()
    print(f"  sbatch {bash_file.name}")
    print()
    print("THEN:")
    base = f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {shlex.quote(phenotype)} --ancestry {ancestry_code}"
    print(f"  {base} --status")
    print(f"  {base} --stats")
    print(f"  {base} --plots")
    print()


# =============================================================================
# GWAS STREAMING
# =============================================================================


def inspect_gwas_columns(path: Path):
    header = pd.read_csv(path, sep="\t", compression="infer", nrows=0)
    columns = list(header.columns)

    mapping = {
        "CHR": find_column(columns, ["CHR"]),
        "POS": find_column(columns, ["POS"]),
        "P": find_column(columns, ["P"]),
        "EA": find_column(columns, ["EA"]),
        "NEA": find_column(columns, ["NEA"]),
        "SNPID": find_column(columns, ["SNP", "rsID", "RSID", "SNPID"]),
    }

    required = ["CHR", "POS", "P"]
    missing = [x for x in required if mapping[x] is None]
    if missing:
        raise RuntimeError(
            f"GWAS file is missing required columns {missing}: {path}"
        )

    return mapping


def stream_candidate_variants(
    qc_file: Path,
    p2: float,
    chunk_size: int,
):
    mapping = inspect_gwas_columns(qc_file)

    usecols = [mapping["CHR"], mapping["POS"], mapping["P"]]
    for optional in ("EA", "NEA", "SNPID"):
        if mapping[optional] is not None:
            usecols.append(mapping[optional])
    usecols = list(dict.fromkeys(usecols))

    per_chrom = defaultdict(list)
    n_total = 0
    n_valid_p = 0
    n_p2 = 0

    reader = pd.read_csv(
        qc_file,
        sep="\t",
        compression="infer",
        usecols=usecols,
        chunksize=chunk_size,
        low_memory=False,
    )

    for chunk in reader:
        n_total += len(chunk)

        renamed = pd.DataFrame(index=chunk.index)
        renamed["CHR"] = clean_chr_series(chunk[mapping["CHR"]])
        renamed["POS"] = pd.to_numeric(chunk[mapping["POS"]], errors="coerce").astype("Int64")
        renamed["P"] = pd.to_numeric(chunk[mapping["P"]], errors="coerce")

        if mapping["EA"] is not None:
            renamed["EA"] = chunk[mapping["EA"]].map(clean_allele)
        else:
            renamed["EA"] = None

        if mapping["NEA"] is not None:
            renamed["NEA"] = chunk[mapping["NEA"]].map(clean_allele)
        else:
            renamed["NEA"] = None

        if mapping["SNPID"] is not None:
            renamed["SOURCE_ID"] = chunk[mapping["SNPID"]].astype(str)
        else:
            renamed["SOURCE_ID"] = ""

        valid = (
            renamed["CHR"].between(1, 22, inclusive="both")
            & renamed["POS"].notna()
            & renamed["P"].notna()
            & renamed["P"].between(0, 1, inclusive="both")
        )
        n_valid_p += int(valid.sum())

        renamed = renamed.loc[valid]
        selected = renamed.loc[renamed["P"] <= p2].copy()
        n_p2 += len(selected)

        if selected.empty:
            continue

        for chrom, frame in selected.groupby("CHR", sort=False):
            per_chrom[int(chrom)].append(frame)

    combined = {}
    for chrom, frames in per_chrom.items():
        df = pd.concat(frames, ignore_index=True)
        df = df.sort_values("P", kind="stable")
        # Exact duplicate GWAS records are unnecessary for clumping.
        subset = ["POS", "P", "EA", "NEA"]
        df = df.drop_duplicates(subset=subset, keep="first")
        combined[chrom] = df

    stats = {
        "N_GWAS_ROWS": n_total,
        "N_VALID_ROWS": n_valid_p,
        "N_P_LE_P2": n_p2,
    }
    return combined, stats


# =============================================================================
# PVAR MAPPING
# =============================================================================


def load_reference_matches(pvar_path: Path, target_positions: set[int]):
    """
    Stream a PLINK2 .pvar and retain only target positions.

    Returns mapping:
      position -> list of dicts with ID, REF, ALT
    """
    result = defaultdict(list)

    with open(pvar_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.rstrip("\n").lstrip("#").split("\t")
                break

        if header is None:
            raise RuntimeError(f"Could not find PVAR header in {pvar_path}")

        index = {name: i for i, name in enumerate(header)}
        required = ["POS", "ID", "REF", "ALT"]
        missing = [x for x in required if x not in index]
        if missing:
            raise RuntimeError(f"PVAR missing columns {missing}: {pvar_path}")

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            try:
                pos = int(fields[index["POS"]])
            except (ValueError, IndexError):
                continue

            if pos not in target_positions:
                continue

            try:
                variant_id = fields[index["ID"]]
                ref = fields[index["REF"]].upper()
                alt = fields[index["ALT"]].upper()
            except IndexError:
                continue

            result[pos].append(
                {
                    "ID": variant_id,
                    "REF": ref,
                    "ALT": alt,
                }
            )

    return result


def map_gwas_to_reference(gwas: pd.DataFrame, pvar_path: Path) -> tuple[pd.DataFrame, dict]:
    if gwas.empty:
        return pd.DataFrame(columns=["ID", "P", "CHR", "POS", "EA", "NEA", "SOURCE_ID"]), {
            "N_INPUT": 0,
            "N_MAPPED": 0,
            "N_UNMAPPED": 0,
            "N_AMBIGUOUS": 0,
            "N_ALLELE_MISMATCH": 0,
            "N_NO_ALLELES": 0,
        }

    positions = set(int(x) for x in gwas["POS"].dropna().astype(int).unique())
    reference = load_reference_matches(pvar_path, positions)

    mapped_rows = []
    n_unmapped = 0
    n_ambiguous = 0
    n_allele_mismatch = 0
    n_no_alleles = 0

    for row in gwas.itertuples(index=False):
        pos = int(row.POS)
        options = reference.get(pos, [])

        if not options:
            n_unmapped += 1
            continue

        ea = clean_allele(row.EA)
        nea = clean_allele(row.NEA)

        # Same chromosome (per-chromosome PVAR) + same GRCh38 position + the
        # same two alleles, in either orientation. Both files are GRCh38 and
        # forward strand (GWAS Catalog harmonised + 1000G 2019), so no strand
        # flipping. A different variant at the same position is NOT a match.
        if not (ea and nea):
            n_no_alleles += 1
            continue

        allele_matches = [
            opt
            for opt in options
            if {opt["REF"], opt["ALT"]} == {ea, nea}
        ]

        if len(allele_matches) == 1:
            chosen = allele_matches[0]
        elif len(allele_matches) > 1:
            n_ambiguous += 1
            continue
        else:
            n_allele_mismatch += 1
            continue

        mapped_rows.append(
            {
                "ID": chosen["ID"],
                "P": float(row.P),
                "CHR": int(row.CHR),
                "POS": pos,
                "REF": chosen["REF"],
                "ALT": chosen["ALT"],
                "EA": ea,
                "NEA": nea,
                "SOURCE_ID": getattr(row, "SOURCE_ID", ""),
            }
        )

    mapped = pd.DataFrame(mapped_rows)
    if not mapped.empty:
        mapped = (
            mapped.sort_values("P", kind="stable")
            .drop_duplicates(subset=["ID"], keep="first")
            .reset_index(drop=True)
        )

    stats = {
        "N_INPUT": int(len(gwas)),
        "N_MAPPED": int(len(mapped)),
        "N_UNMAPPED": int(n_unmapped),
        "N_AMBIGUOUS": int(n_ambiguous),
        "N_ALLELE_MISMATCH": int(n_allele_mismatch),
        "N_NO_ALLELES": int(n_no_alleles),
    }
    return mapped, stats


# =============================================================================
# PLINK CLUMPING
# =============================================================================


def run_plink_clump(
    plink2: Path,
    pfile_prefix: Path,
    clump_input: Path,
    out_prefix: Path,
    p1: float,
    p2: float,
    r2: float,
    kb: int,
    threads: int,
):
    cmd = [
        str(plink2),
        "--pfile", str(pfile_prefix),
        "--clump", str(clump_input),
        "--clump-id-field", "ID",
        "--clump-p-field", "P",
        "--clump-p1", str(p1),
        "--clump-p2", str(p2),
        "--clump-r2", str(r2),
        "--clump-kb", str(kb),
        "--threads", str(threads),
        "--out", str(out_prefix),
    ]

    print("$ " + " ".join(shlex.quote(x) for x in cmd))
    subprocess.run(cmd, check=True)


# =============================================================================
# WORKER
# =============================================================================


def worker_mode() -> None:
    task_value = os.environ.get("SLURM_ARRAY_TASK_ID")
    manifest_value = os.environ.get("GWAS_CLUMP_MANIFEST")

    if not task_value:
        raise RuntimeError("SLURM_ARRAY_TASK_ID is not defined.")
    if not manifest_value:
        raise RuntimeError("GWAS_CLUMP_MANIFEST is not defined.")

    task_id = int(task_value)
    manifest_file = Path(manifest_value).resolve()
    root = Path.cwd().resolve()
    plink2 = root / "envs" / "pipeline" / "bin" / "plink2"

    manifest = pd.read_csv(manifest_file, sep="\t", dtype=str, low_memory=False)
    required = [
        "CLUMP_TASK_ID",
        "STUDY_ACCESSION",
        "QC_FILE",
        "OUTPUT_DIR",
        "PHENOTYPE",
        "ANCESTRY_CODE",
        "ANCESTRY_LABEL",
        "P1",
        "P2",
        "R2",
        "KB",
        "CHUNK_SIZE",
    ]
    missing = [x for x in required if x not in manifest.columns]
    if missing:
        raise RuntimeError(f"Clumping manifest missing columns: {missing}")

    task_numbers = pd.to_numeric(manifest["CLUMP_TASK_ID"], errors="coerce")
    selected = manifest.loc[task_numbers == task_id]
    if selected.empty:
        raise RuntimeError(f"No manifest rows for array task {task_id}")

    failed = []
    for _, row in selected.iterrows():
        try:
            clump_one_study(row, task_id, root, plink2)
        except Exception as error:
            failed.append(f"{row['STUDY_ACCESSION']}: {error}")

    if failed:
        raise SystemExit("Clumping failed for:\n  " + "\n  ".join(failed))


def clump_one_study(row, task_id, root, plink2) -> None:
    accession = str(row["STUDY_ACCESSION"]).strip()
    qc_file = Path(row["QC_FILE"]).resolve()
    output_dir = Path(row["OUTPUT_DIR"]).resolve()
    phenotype = str(row["PHENOTYPE"]).strip()
    ancestry_code = str(row["ANCESTRY_CODE"]).strip().upper()
    ancestry_label = str(row["ANCESTRY_LABEL"]).strip()
    p1 = float(row["P1"])
    p2 = float(row["P2"])
    r2 = float(row["R2"])
    kb = int(float(row["KB"]))
    chunk_size = int(float(row["CHUNK_SIZE"]))
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", DEFAULT_CPUS))
    # Optional FORCE column: re-run even when outputs exist.
    force = str(row.get("FORCE", "False")).strip().lower() in {"true", "1", "yes", "y"}

    output_dir.mkdir(parents=True, exist_ok=True)

    lead_output = output_dir / "lead_variants.tsv"
    mapped_output = output_dir / "mapped_candidates.tsv.gz"
    summary_tsv = output_dir / "clumping_summary.tsv"
    summary_json = output_dir / "clumping_summary.json"
    failed_output = output_dir / "CLUMPING_FAILED.txt"

    # Resume only when all key deliverables are present.
    if (
        not force
        and lead_output.exists()
        and summary_tsv.exists()
        and summary_json.exists()
        and mapped_output.exists()
    ):
        banner("CLUMPING ALREADY COMPLETE")
        print(f"Study: {accession}")
        print(f"Output: {output_dir}")
        return

    if not qc_file.exists() or qc_file.stat().st_size <= 0:
        raise FileNotFoundError(f"QC file not found or empty: {qc_file}")

    banner("STEP 05 - LD CLUMPING WORKER")
    print(f"Task          : {task_id}")
    print(f"Study         : {accession}")
    print(f"Phenotype     : {phenotype}")
    print(f"Ancestry      : {ancestry_label} ({ancestry_code})")
    print(f"QC file       : {qc_file}")
    print(f"Output        : {output_dir}")
    print(f"p1 / p2       : {p1} / {p2}")
    print(f"r2 / kb       : {r2} / {kb}")
    print(f"Threads       : {threads}")

    try:
        banner("STREAMING QC GWAS")
        per_chrom, stream_stats = stream_candidate_variants(
            qc_file=qc_file,
            p2=p2,
            chunk_size=chunk_size,
        )

        print(f"GWAS rows      : {stream_stats['N_GWAS_ROWS']:,}")
        print(f"Valid rows     : {stream_stats['N_VALID_ROWS']:,}")
        print(f"P <= p2        : {stream_stats['N_P_LE_P2']:,}")

        all_mapped = []
        all_leads = []
        chromosome_summaries = []

        for chrom in range(1, 23):
            chrom_dir = output_dir / f"chr{chrom:02d}"
            chrom_dir.mkdir(parents=True, exist_ok=True)

            gwas_chr = per_chrom.get(
                chrom,
                pd.DataFrame(columns=["CHR", "POS", "P", "EA", "NEA", "SOURCE_ID"]),
            )

            prefix = reference_prefix(root, ancestry_code, chrom)
            pvar = Path(str(prefix) + ".pvar")

            mapped, map_stats = map_gwas_to_reference(gwas_chr, pvar)

            if not mapped.empty:
                all_mapped.append(mapped)

            clump_input = chrom_dir / "clump_input.tsv"
            clump_report_prefix = chrom_dir / "plink"

            # Remove stale chromosome-level reports from an interrupted prior run.
            for suffix in (
                ".clumps",
                ".clumps.zst",
                ".clumps.missing_id",
                ".clumps.missing_id.zst",
                ".clumps.missing_allele",
                ".clumps.missing_allele.zst",
            ):
                stale = Path(str(clump_report_prefix) + suffix)
                if stale.exists():
                    stale.unlink()

            # PLINK only needs ID and P for this biallelic reference.
            mapped[["ID", "P"]].to_csv(
                clump_input,
                sep="\t",
                index=False,
            )

            n_p1 = int((mapped["P"] <= p1).sum()) if not mapped.empty else 0
            clumps_file = Path(str(clump_report_prefix) + ".clumps")

            if not mapped.empty and n_p1 > 0:
                banner(f"CHR {chrom} CLUMPING")
                run_plink_clump(
                    plink2=plink2,
                    pfile_prefix=prefix,
                    clump_input=clump_input,
                    out_prefix=clump_report_prefix,
                    p1=p1,
                    p2=p2,
                    r2=r2,
                    kb=kb,
                    threads=threads,
                )

            n_leads = 0
            if clumps_file.exists() and clumps_file.stat().st_size > 0:
                try:
                    clumps = pd.read_csv(
                        clumps_file,
                        sep=r"\s+",
                        engine="python",
                    )
                except pd.errors.EmptyDataError:
                    clumps = pd.DataFrame()

                if not clumps.empty:
                    clumps.insert(0, "STUDY_ACCESSION", accession)
                    clumps.insert(1, "PHENOTYPE", phenotype)
                    clumps.insert(2, "ANCESTRY_CODE", ancestry_code)
                    clumps.insert(3, "ANCESTRY_LABEL", ancestry_label)
                    clumps.insert(4, "SOURCE_CHROMOSOME", chrom)
                    all_leads.append(clumps)
                    n_leads = len(clumps)

            chromosome_summaries.append(
                {
                    "CHR": chrom,
                    "N_GWAS_P_LE_P2": map_stats["N_INPUT"],
                    "N_MAPPED": map_stats["N_MAPPED"],
                    "N_UNMAPPED": map_stats["N_UNMAPPED"],
                    "N_AMBIGUOUS": map_stats["N_AMBIGUOUS"],
                    "N_ALLELE_MISMATCH": map_stats["N_ALLELE_MISMATCH"],
                    "N_NO_ALLELES": map_stats["N_NO_ALLELES"],
                    "N_P_LE_P1_MAPPED": n_p1,
                    "N_LEAD_VARIANTS": n_leads,
                }
            )

        if all_mapped:
            mapped_all = pd.concat(all_mapped, ignore_index=True)
            mapped_all = mapped_all.sort_values(["CHR", "POS", "P"], kind="stable")
        else:
            mapped_all = pd.DataFrame(
                columns=["ID", "P", "CHR", "POS", "REF", "ALT", "EA", "NEA", "SOURCE_ID"]
            )

        mapped_all.to_csv(
            mapped_output,
            sep="\t",
            index=False,
            compression="gzip",
        )

        if all_leads:
            lead_all = pd.concat(all_leads, ignore_index=True)
            p_column = "P" if "P" in lead_all.columns else None
            if p_column:
                lead_all = lead_all.sort_values(p_column, kind="stable")
        else:
            lead_all = pd.DataFrame(
                columns=[
                    "STUDY_ACCESSION",
                    "PHENOTYPE",
                    "ANCESTRY_CODE",
                    "ANCESTRY_LABEL",
                    "SOURCE_CHROMOSOME",
                    "CHROM",
                    "POS",
                    "ID",
                    "P",
                ]
            )

        lead_all.to_csv(lead_output, sep="\t", index=False)

        chromosome_summary_file = output_dir / "chromosome_clumping_summary.tsv"
        chromosome_df = pd.DataFrame(chromosome_summaries)
        chromosome_df.to_csv(chromosome_summary_file, sep="\t", index=False)

        total_mapped = int(chromosome_df["N_MAPPED"].sum())
        total_input_for_mapping = int(chromosome_df["N_GWAS_P_LE_P2"].sum())
        mapping_rate = (
            total_mapped / total_input_for_mapping
            if total_input_for_mapping > 0
            else None
        )

        summary = {
            "CLUMP_TASK_ID": task_id,
            "STUDY_ACCESSION": accession,
            "PHENOTYPE": phenotype,
            "ANCESTRY": ancestry_label,
            "ANCESTRY_CODE": ancestry_code,
            "ANCESTRY_LABEL": ancestry_label,
            "BUILD": "GRCh38",
            "QC_FILE": str(qc_file),
            "OUTPUT_DIR": str(output_dir),
            "REFERENCE_ROOT": str(root / "resources" / "1000G" / ancestry_code),
            "P1": p1,
            "P2": p2,
            "R2": r2,
            "KB": kb,
            "N_GWAS_ROWS": stream_stats["N_GWAS_ROWS"],
            "N_VALID_ROWS": stream_stats["N_VALID_ROWS"],
            "N_P_LE_P2": stream_stats["N_P_LE_P2"],
            "N_REFERENCE_MAPPED": total_mapped,
            "N_NOT_IN_REFERENCE": int(chromosome_df["N_UNMAPPED"].sum()),
            "N_ALLELE_MISMATCH": int(chromosome_df["N_ALLELE_MISMATCH"].sum()),
            "N_AMBIGUOUS": int(chromosome_df["N_AMBIGUOUS"].sum()),
            "N_NO_ALLELES": int(chromosome_df["N_NO_ALLELES"].sum()),
            "MAPPING_RATE": mapping_rate,
            "N_INDEX_CANDIDATES_P_LE_P1": int(chromosome_df["N_P_LE_P1_MAPPED"].sum()),
            "N_LEAD_VARIANTS": int(len(lead_all)),
            "LEAD_VARIANTS_FILE": str(lead_output),
            "MAPPED_CANDIDATES_FILE": str(mapped_output),
            "CHROMOSOME_SUMMARY_FILE": str(chromosome_summary_file),
            "THREADS": threads,
            "COMPLETED_UTC": datetime.now(timezone.utc).isoformat(),
        }

        pd.DataFrame([summary]).to_csv(summary_tsv, sep="\t", index=False)
        summary_json.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

        if failed_output.exists():
            failed_output.unlink()

        banner("LD CLUMPING COMPLETE")
        print(f"Study                  : {accession}")
        print(f"P <= p2 candidates     : {stream_stats['N_P_LE_P2']:,}")
        print(f"Reference mapped       : {total_mapped:,}")
        print(f"Mapping rate           : {mapping_rate}")
        print(f"Index candidates p<=p1 : {summary['N_INDEX_CANDIDATES_P_LE_P1']:,}")
        print(f"Lead variants          : {len(lead_all):,}")
        print()
        print(f"Lead variants:\n  {lead_output}")
        print(f"Mapped candidates:\n  {mapped_output}")
        print(f"Summary:\n  {summary_tsv}")

    except Exception as error:
        failure_text = (
            "STEP 05 LD CLUMPING FAILED\n\n"
            f"Task ID: {task_id}\n"
            f"Accession: {accession}\n"
            f"QC file: {qc_file}\n\n"
            f"Error:\n{type(error).__name__}: {error}\n\n"
            f"{traceback.format_exc()}"
        )
        failed_output.write_text(failure_text, encoding="utf-8")
        print(failure_text, file=sys.stderr)
        raise


# =============================================================================
# DIRECT SINGLE-INDEX MODE
# =============================================================================


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
    """--index N: study N of the Step02 download manifest, clumped on this
    machine with the normal worker. Only this study has to be QC'd."""
    validate_arguments(args)

    root = PROJECT_ROOT
    phenotype = args.phenotype.strip()
    ancestry_code, ancestry_label = canonical_ancestry(args.ancestry)
    p_slug, a_slug = slugify(phenotype), slugify(ancestry_label)
    accession = study_for_index(p_slug, a_slug, args.index)
    quoted = shlex.quote(phenotype)

    qc_file = root / "02_summary_stats" / p_slug / a_slug / "qc" / accession / f"{accession}_GRCh38_QC.tsv.gz"
    if not qc_file.is_file():
        raise SystemExit(
            f"Step02 QC of study {args.index} ({accession}) is not done yet. Run:\n"
            f"  envs/pipeline/bin/python Step02_Plan_QC_GWAS.py --phenotype {quoted} "
            f"--ancestry {ancestry_code} --index {args.index}")

    if not (root / "envs" / "pipeline" / "bin" / "plink2").exists():
        raise SystemExit("plink2 missing in envs/pipeline. Run: python Step00_Installation.py")
    validate_reference(root, ancestry_code)

    output_root = root / "05_ld_clumping" / p_slug / a_slug
    output_root.mkdir(parents=True, exist_ok=True)
    row, exclusion = build_clump_row(args, task_id=1, accession=accession, qc_file=qc_file,
                                     output_root=output_root, phenotype=phenotype,
                                     ancestry_code=ancestry_code, ancestry_label=ancestry_label)
    if exclusion is not None:
        raise SystemExit(f"{accession} is not clumped: {exclusion['REASON']}")

    banner(f"STEP 05 - INDEX {args.index}: {accession} (clumping on this machine)")
    single = output_root / f".index_{args.index}_clump_manifest.tsv"
    pd.DataFrame([row]).to_csv(single, sep="\t", index=False)

    saved = {k: os.environ.get(k) for k in ("SLURM_ARRAY_TASK_ID", "GWAS_CLUMP_MANIFEST", "SLURM_CPUS_PER_TASK")}
    os.environ["SLURM_ARRAY_TASK_ID"] = "1"
    os.environ["GWAS_CLUMP_MANIFEST"] = str(single)
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

    print(f"\nNext: envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype {quoted} "
          f"--ancestry {ancestry_code} --index {args.index}")


# =============================================================================
# --status / --stats / --plots
# =============================================================================

# GRCh38 autosome lengths (bp), for Manhattan x-axis and chromosome maps.
GRCH38_LENGTHS = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}


def project_paths(args) -> dict:
    root = PROJECT_ROOT
    ancestry_code, ancestry_label = canonical_ancestry(args.ancestry)
    phenotype_slug = slugify(args.phenotype.strip())
    ancestry_slug = slugify(ancestry_label)
    output_root = root / "05_ld_clumping" / phenotype_slug / ancestry_slug
    return {
        "phenotype": args.phenotype.strip(),
        "code": ancestry_code,
        "label": ancestry_label,
        "qc_root": root / "02_summary_stats" / phenotype_slug / ancestry_slug / "qc",
        "output_root": output_root,
        "manifest": output_root / "clumping_manifest.tsv",
        "plots": output_root / "plots",
        "base": f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {shlex.quote(args.phenotype.strip())} "
                f"--ancestry {ancestry_code}",
    }


def study_list(P) -> list[tuple[str, Path]]:
    """(accession, QC file) for every QC'd study of this phenotype/ancestry."""
    return [(f.parent.name, f) for f in sorted(P["qc_root"].glob("*/*_GRCh38_QC.tsv.gz"))]


def clump_state(P, accession: str) -> str:
    out = P["output_root"] / accession
    if (out / "lead_variants.tsv").is_file() and (out / "clumping_summary.tsv").is_file():
        return "DONE"
    if (out / "CLUMPING_FAILED.txt").is_file():
        return "FAILED"
    return "PENDING"


def status_mode(args) -> None:
    P = project_paths(args)
    planned = set()
    if P["manifest"].is_file():
        planned = set(pd.read_csv(P["manifest"], sep="\t", dtype=str)["STUDY_ACCESSION"].astype(str))

    banner(f"STEP 05 STATUS - {P['phenotype']} / {P['label']} ({P['code']})")
    studies = study_list(P)
    if not studies:
        print(f"No Step02 QC files under {P['qc_root']}")
        return

    counts = defaultdict(int)
    print(f"{'STUDY':<15}{'CLUMPING':<10}{'LEADS':>7}  DETAIL")
    for accession, _ in studies:
        state = clump_state(P, accession)
        if state == "PENDING" and accession not in planned:
            state = "NOT_PLANNED"
        counts[state] += 1
        leads, detail = "", ""
        out = P["output_root"] / accession
        if state == "DONE":
            summary = pd.read_csv(out / "clumping_summary.tsv", sep="\t").iloc[0]
            leads = str(int(summary.get("N_LEAD_VARIANTS", 0)))
        elif state == "FAILED":
            lines = (out / "CLUMPING_FAILED.txt").read_text(encoding="utf-8", errors="replace").splitlines()
            detail = next((x for x in lines if x.startswith(("RuntimeError", "FileNotFound", "ValueError",
                                                              "CalledProcess", "KeyError"))), "")[:90]
        print(f"{accession:<15}{state:<10}{leads:>7}  {detail}")

    print()
    print("Studies : " + "  ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print()
    if counts["NOT_PLANNED"] or counts["FAILED"] or counts["PENDING"]:
        print("Plan / re-run (finished studies are skipped):")
        print(f"  {P['base']}")
        print(f"  sbatch Step05_Clump_GWAS_{slugify(P['phenotype'])}_{slugify(P['label'])}.sh")
    else:
        print("Every QC'd study is clumped. Next:")
        print(f"  {P['base']} --stats")
        print(f"  envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype {shlex.quote(P['phenotype'])} "
              f"--ancestry {P['code']}")


def count_loci(leads: pd.DataFrame, window_kb: int) -> int:
    """Leads within +-window_kb merged into one locus (as Step06 does)."""
    if leads.empty:
        return 0
    chr_col = find_column(leads.columns, ["#CHROM", "CHROM", "CHR"])
    pos_col = find_column(leads.columns, ["POS", "BP"])
    work = pd.DataFrame({"CHR": pd.to_numeric(leads[chr_col], errors="coerce"),
                         "POS": pd.to_numeric(leads[pos_col], errors="coerce")}).dropna()
    n = 0
    for _, group in work.sort_values(["CHR", "POS"]).groupby("CHR"):
        end = -1
        for pos in group["POS"]:
            if pos - window_kb * 1000 > end:
                n += 1
            end = max(end, pos + window_kb * 1000)
    return n


def stats_mode(args) -> None:
    P = project_paths(args)
    banner(f"STEP 05 STATISTICS - {P['phenotype']} / {P['label']} ({P['code']})")
    rows = []

    for accession, qc_file in study_list(P):
        if clump_state(P, accession) != "DONE":
            continue
        out = P["output_root"] / accession
        s = pd.read_csv(out / "clumping_summary.tsv", sep="\t").iloc[0].to_dict()
        leads = pd.read_csv(out / "lead_variants.tsv", sep="\t", low_memory=False)

        qc_summary = qc_file.parent / f"{accession}_QC_summary.tsv"
        gws = None
        if qc_summary.is_file():
            q = pd.read_csv(qc_summary, sep="\t").iloc[0]
            gws = q.get("N_GENOME_WIDE_SIGNIFICANT")

        mapped_gws = int(s.get("N_INDEX_CANDIDATES_P_LE_P1", 0))
        row = {
            "STUDY_ACCESSION": accession,
            "N_GWAS_VARIANTS": s.get("N_GWAS_ROWS"),
            "N_P_LE_P2": s.get("N_P_LE_P2"),
            "N_MAPPED_TO_1000G": s.get("N_REFERENCE_MAPPED"),
            "MAPPING_RATE": round(float(s["MAPPING_RATE"]), 4) if pd.notna(s.get("MAPPING_RATE")) else "",
            "N_NOT_IN_1000G": s.get("N_NOT_IN_REFERENCE", ""),
            "N_ALLELE_MISMATCH": s.get("N_ALLELE_MISMATCH", ""),
            "N_GENOME_WIDE_SIGNIFICANT": gws if gws is not None else "",
            "N_GWS_MAPPED": mapped_gws,
            "N_GWS_NOT_IN_PANEL": int(gws) - mapped_gws if gws is not None and pd.notna(gws) else "",
            "N_LEAD_VARIANTS": s.get("N_LEAD_VARIANTS"),
            "N_LOCI_1MB": count_loci(leads, 1000),
            "P1": s.get("P1"), "P2": s.get("P2"), "R2": s.get("R2"), "KB": s.get("KB"),
        }
        rows.append(row)

        print()
        print(f"{accession}: {int(row['N_LEAD_VARIANTS']):,} independent leads in {row['N_LOCI_1MB']} loci "
              f"(+-1 Mb merged)")
        print(f"   P<=p2 variants {int(row['N_P_LE_P2']):,} | mapped to 1000G {int(row['N_MAPPED_TO_1000G']):,} "
              f"({row['MAPPING_RATE']}) | not in panel {row['N_NOT_IN_1000G']} | allele mismatch "
              f"{row['N_ALLELE_MISMATCH']}")
        if row["N_GWS_NOT_IN_PANEL"] not in ("", 0):
            print(f"   ! {row['N_GWS_NOT_IN_PANEL']} genome-wide significant variants are not in the "
                  f"1000G {P['code']} panel (rare / absent / allele mismatch) and cannot lead a clump")
        if not leads.empty:
            id_col = find_column(leads.columns, ["ID"])
            chr_col = find_column(leads.columns, ["#CHROM", "CHROM", "CHR"])
            pos_col = find_column(leads.columns, ["POS"])
            top = leads.sort_values("P").head(5)
            print("   top leads: " + "; ".join(
                f"{r[id_col]} chr{r[chr_col]}:{r[pos_col]} P={r['P']:.2g}" for _, r in top.iterrows()))

    if not rows:
        print("No finished studies yet (see --status).")
        return

    output = P["output_root"] / "clumping_stats.tsv"
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
    print(f"\nTable: {output}")


def read_association(qc_file: Path) -> pd.DataFrame:
    frames = []
    usecols = lambda c: c in ("CHR", "POS", "P")  # noqa: E731
    for chunk in pd.read_csv(qc_file, sep="\t", usecols=usecols, chunksize=2_000_000, low_memory=False):
        chunk = chunk.copy()
        chunk["CHR"] = clean_chr_series(chunk["CHR"])
        chunk["POS"] = pd.to_numeric(chunk["POS"], errors="coerce")
        chunk["P"] = pd.to_numeric(chunk["P"], errors="coerce")
        chunk = chunk.dropna()
        chunk = chunk[chunk["CHR"].between(1, 22) & chunk["P"].between(0, 1)]
        frames.append(chunk.astype({"CHR": "int8", "POS": "int64", "P": "float64"}))
    data = pd.concat(frames, ignore_index=True)
    data["LOGP"] = -np.log10(data["P"].clip(lower=1e-300))
    return data


def plot_study(P, accession: str, qc_file: Path, plt) -> list[Path]:
    out = P["output_root"] / accession
    data = read_association(qc_file)
    leads = pd.read_csv(out / "lead_variants.tsv", sep="\t", low_memory=False) \
        if (out / "lead_variants.tsv").is_file() else pd.DataFrame()
    P["plots"].mkdir(parents=True, exist_ok=True)
    files = []

    offsets, start = {}, 0
    for chrom in range(1, 23):
        offsets[chrom] = start
        start += GRCH38_LENGTHS[chrom]

    chi2 = (np.asarray(_norm_isf(data["P"].to_numpy() / 2)) ** 2)
    lambda_gc = float(np.median(chi2) / 0.4549364)

    # ---- Manhattan (all P < 1e-3 plus a fixed 10% sample of the rest) ----
    rng = np.random.default_rng(1)
    keep = (data["P"] < 1e-3) | (rng.random(len(data)) < 0.10)
    shown = data[keep]
    x = shown["POS"].to_numpy() + shown["CHR"].map(offsets).to_numpy()

    fig, ax = plt.subplots(figsize=(16, 5.5))
    colors = np.where(shown["CHR"].to_numpy() % 2 == 0, "#2b6c9e", "#8fb8d9")
    ax.scatter(x, shown["LOGP"], c=colors, s=3, linewidths=0, rasterized=True)
    ax.axhline(-np.log10(5e-8), color="#c0392b", lw=1, ls="--", label="P = 5e-8")
    ax.axhline(-np.log10(1e-5), color="#7f8c8d", lw=0.8, ls=":", label="P = 1e-5")

    if not leads.empty:
        chr_col = find_column(leads.columns, ["#CHROM", "CHROM", "CHR"])
        pos_col = find_column(leads.columns, ["POS"])
        id_col = find_column(leads.columns, ["ID"])
        lead_chr = pd.to_numeric(leads[chr_col], errors="coerce").astype(int)
        lead_x = leads[pos_col].astype(int) + lead_chr.map(offsets)
        lead_y = -np.log10(leads["P"].clip(lower=1e-300))
        ax.scatter(lead_x, lead_y, s=26, facecolors="none", edgecolors="#c0392b", linewidths=1.1,
                   label=f"independent leads ({len(leads)})")
        for _, r in leads.sort_values("P").head(25).iterrows():
            ax.annotate(str(r[id_col]), (int(r[pos_col]) + offsets[int(r[chr_col])],
                                         -np.log10(max(r["P"], 1e-300))),
                        fontsize=6, xytext=(2, 3), textcoords="offset points", rotation=30)

    ax.set_xticks([offsets[c] + GRCH38_LENGTHS[c] / 2 for c in range(1, 23)])
    ax.set_xticklabels([str(c) for c in range(1, 23)], fontsize=8)
    ax.set_xlim(0, start)
    ax.set_ylim(0, max(10.0, float(data["LOGP"].max()) * 1.08))
    ax.set_xlabel("Chromosome (GRCh38)")
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_title(f"{accession} - {P['phenotype']} ({P['code']})   variants {len(data):,}   "
                 f"lambda GC {lambda_gc:.3f}")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    files.append(P["plots"] / f"{accession}_manhattan.png")
    fig.savefig(files[-1], dpi=150)
    plt.close(fig)

    # ---- QQ ----
    observed = np.sort(data["LOGP"].to_numpy())[::-1]
    n = len(observed)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    idx = np.unique(np.concatenate([np.arange(min(n, 20000)),
                                    np.linspace(0, n - 1, min(n, 20000)).astype(int)]))
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(expected[idx], observed[idx], s=4, color="#2b6c9e", linewidths=0, rasterized=True)
    lim = max(expected.max(), 1.0)
    ax.plot([0, lim], [0, lim], color="#c0392b", lw=1)
    ax.set_xlabel(r"Expected $-\log_{10}(P)$")
    ax.set_ylabel(r"Observed $-\log_{10}(P)$")
    ax.set_title(f"{accession} QQ   lambda GC = {lambda_gc:.3f}")
    fig.tight_layout()
    files.append(P["plots"] / f"{accession}_qq.png")
    fig.savefig(files[-1], dpi=150)
    plt.close(fig)

    # ---- chromosome map of lead variants / +-1 Mb loci ----
    fig, ax = plt.subplots(figsize=(12, 7))
    for chrom in range(1, 23):
        y = 23 - chrom
        ax.add_patch(plt.Rectangle((0, y - 0.3), GRCH38_LENGTHS[chrom], 0.6, color="#e5e8ec", lw=0))
    if not leads.empty:
        for _, r in leads.iterrows():
            chrom, pos = int(r[chr_col]), int(r[pos_col])
            y = 23 - chrom
            ax.add_patch(plt.Rectangle((max(0, pos - 1_000_000), y - 0.3), 2_000_000, 0.6,
                                       color="#f5b7b1", lw=0))
            ax.plot([pos, pos], [y - 0.36, y + 0.36], color="#c0392b", lw=1.2)
        for _, r in leads.sort_values("P").head(30).iterrows():
            ax.annotate(str(r[id_col]), (int(r[pos_col]), 23 - int(r[chr_col]) + 0.38),
                        fontsize=6, ha="center")
    ax.set_yticks([23 - c for c in range(1, 23)])
    ax.set_yticklabels([f"chr{c}" for c in range(1, 23)], fontsize=8)
    ax.set_xlim(0, max(GRCH38_LENGTHS.values()) * 1.02)
    ax.set_ylim(0.3, 22.8)
    ax.set_xlabel("Position (Mb, GRCh38)")
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1e6:.0f}"))
    ax.set_title(f"{accession} - independent leads (red) and +-1 Mb fine-mapping windows")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    files.append(P["plots"] / f"{accession}_loci_map.png")
    fig.savefig(files[-1], dpi=150)
    plt.close(fig)

    return files


def _norm_isf(values):
    from scipy.stats import norm
    return norm.isf(np.clip(values, 1e-300, 1.0))


def plots_mode(args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    P = project_paths(args)
    banner(f"STEP 05 PLOTS - {P['phenotype']} / {P['label']} ({P['code']})")

    made = 0
    for accession, qc_file in study_list(P):
        if clump_state(P, accession) != "DONE":
            print(f"{accession}: clumping not finished - skipped")
            continue
        for path in plot_study(P, accession, qc_file, plt):
            print(f"  {path}")
        made += 1

    print(f"\n{made} studies plotted in {P['plots']}")


# =============================================================================
# MAIN
# =============================================================================


def main() -> None:
    # Generated SLURM-array worker mode takes priority and requires no CLI args.
    if os.environ.get("SLURM_ARRAY_TASK_ID") and os.environ.get("GWAS_CLUMP_MANIFEST"):
        worker_mode()
        return

    args = arguments()

    if args.status:
        status_mode(args)
        return
    if args.stats:
        stats_mode(args)
        return
    if args.plots:
        plots_mode(args)
        return

    # Direct single-study execution using the same manifest/worker implementation.
    if args.index is not None:
        direct_index_mode(args)
        return

    # Standard planner mode: create/refresh manifest + SLURM array script.
    planner_mode(args)


if __name__ == "__main__":
    main()
