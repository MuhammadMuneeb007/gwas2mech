#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 04 - Ancestry-specific 1000 Genomes GRCh38 LD reference panels
===============================================================================

Builds the PLINK2 panels that Step05 (clumping) and Step06 (fine-mapping LD)
read. Panels are shared by every phenotype, so they are built once per
ancestry, inside the project:

    resources/1000G/<ANC>/chr<N>_<ANC>_GRCh38.pgen / .pvar / .psam
    resources/1000G/<ANC>/samples.txt            (1000G sample IDs used)
    resources/1000G/<ANC>/<ANC>_reference_manifest.tsv

Inputs come from Step00_Download_Resources.py --groups 1000g (never downloaded
here):

    resources/1000G/raw/ALL.chr<N>.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz (+ .tbi)
    resources/1000G/metadata/integrated_call_samples_v3.20130502.ALL.panel

Usage:

    # plan + SLURM array (one task per ancestry x chromosome)
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR
    sbatch Step04_Prepare_1000G_EUR_GRCh38.sh

    # several ancestries / all five 1000G super-populations
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR EAS
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry all

    # progress (read-only) and panel statistics
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --status
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --stats

    # chromosome 22 directly on this machine (no planner needed)
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --index 22

Each chromosome: bcftools subsets the ancestry samples, plink2 converts to
PGEN, keeps rsIDs, names variants without an ID CHR:POS:REF:ALT and removes
duplicate IDs. REF/ALT are the 1000G GRCh38 alleles; Step05/06 match GWAS
alleles to them explicitly and orient Z to REF.

Note: the 1000G EUR panel has ~503 people (CEU, FIN, GBR, IBS, TSI). That is a
standard, public LD reference but a small one; for very large GWAS a larger
ancestry-matched panel (or in-sample LD) gives more reliable fine-mapping.

SLURM settings come from the #SBATCH header of HPC_CPU.sh; --partition --qos
--account --cpus --mem --time --max-parallel override them.
===============================================================================
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path.cwd().resolve()
CPUS_PER_TASK = 4
MAX_ARRAY_TASKS = 1000

REFERENCE_ROOT = PROJECT_ROOT / "resources" / "1000G"
RAW_DIR = REFERENCE_ROOT / "raw"
METADATA_DIR = REFERENCE_ROOT / "metadata"
PANEL_FILE = METADATA_DIR / "integrated_call_samples_v3.20130502.ALL.panel"

BASE_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
)

# 1000 Genomes super-population -> sub-populations.
POPULATIONS = {
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "AFR": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "EAS": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "SAS": ["BEB", "GIH", "ITU", "PJL", "STU"],
    "AMR": ["CLM", "MXL", "PEL", "PUR"],
}

ANCESTRY_ALIASES = {
    "eur": "EUR", "european": "EUR",
    "afr": "AFR", "african": "AFR",
    "eas": "EAS", "east asian": "EAS",
    "sas": "SAS", "south asian": "SAS",
    "amr": "AMR", "american": "AMR", "hispanic": "AMR", "latino": "AMR",
    "hispanic or latin american": "AMR",
}


# =============================================================================
# HELPERS
# =============================================================================

def banner(text):

    print()
    print("=" * 80)
    print(text)
    print("=" * 80)


def normalize_text(value):

    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", str(value).lower())).strip()


def ancestry_code(value):

    code = ANCESTRY_ALIASES.get(normalize_text(value))

    if code is None:
        raise SystemExit(
            f"Ancestry {value!r} has no 1000 Genomes panel. "
            "Use EUR, AFR, EAS, SAS or AMR (or 'all')."
        )

    return code


def chromosome_filename(chromosome):

    return (
        f"ALL.chr{chromosome}."
        "shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    )


def panel_prefix(code, chromosome):

    return REFERENCE_ROOT / code / f"chr{chromosome}_{code}_GRCh38"


def panel_files(code, chromosome):

    prefix = panel_prefix(code, chromosome)

    return [Path(str(prefix) + ext) for ext in (".pgen", ".pvar", ".psam")]


def panel_complete(code, chromosome):

    return all(p.is_file() and p.stat().st_size > 0 for p in panel_files(code, chromosome))


def failure_file(code, chromosome):

    return REFERENCE_ROOT / code / f"chr{chromosome}_{code}_FAILED.txt"


def require_program(program):

    local = PROJECT_ROOT / "envs" / "pipeline" / "bin" / program
    path = str(local) if local.exists() else shutil.which(program)

    if path is None:
        raise RuntimeError(
            f"Required program is not available: {program}\n"
            "Run: python Step00_Installation.py"
        )

    return path


def require_input(path):

    path = Path(path)

    if path.is_file() and path.stat().st_size > 0:
        return path

    raise FileNotFoundError(
        f"Missing 1000 Genomes input: {path}\n"
        "Step04 does not download shared data. Run:\n"
        "  envs/pipeline/bin/python Step00_Download_Resources.py --groups 1000g "
        "--output Step00_Download_1000G.sh\n"
        "  sbatch Step00_Download_1000G.sh"
    )


def run_command(command):

    print()
    print(" ".join(str(x) for x in command), flush=True)
    subprocess.run([str(x) for x in command], check=True)


def count_nonheader_lines(path):

    count = 0

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                count += 1

    return count


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

    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": str(CPUS_PER_TASK),
              "mem": "16G", "time": "12:00:00"}
    values.update({k: v for k, v in read_slurm_template(PROJECT_ROOT / args.template).items()
                   if k not in PER_JOB_KEYS})

    for key, value in (("partition", args.partition), ("qos", args.qos), ("account", args.account),
                       ("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time)):
        if value not in (None, ""):
            values[key] = str(value)

    args.slurm_values = values
    args.cpus = int(values["cpus-per-task"])

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
# SAMPLE LIST
# =============================================================================

def create_sample_list(code):

    panel = pd.read_csv(require_input(PANEL_FILE), sep="\t", dtype=str, low_memory=False)

    for column in ("sample", "pop"):
        if column not in panel.columns:
            raise RuntimeError(f"1000 Genomes panel missing column: {column}")

    selected = panel[panel["pop"].isin(POPULATIONS[code])].copy()

    if selected.empty:
        raise RuntimeError(f"No samples for populations {', '.join(POPULATIONS[code])}")

    output = REFERENCE_ROOT / code / "samples.txt"
    output.parent.mkdir(parents=True, exist_ok=True)
    selected[["sample"]].to_csv(output, index=False, header=False)

    return output, selected


def collect_reference_summaries(code):

    files = sorted((REFERENCE_ROOT / code).glob(f"chr*_{code}_summary.tsv"))
    frames = []

    for file in files:
        try:
            frames.append(pd.read_csv(file, sep="\t"))
        except Exception:
            continue

    if not frames:
        return None

    combined = pd.concat(frames, ignore_index=True)
    combined["CHR"] = pd.to_numeric(combined["CHR"], errors="coerce")
    combined = combined.sort_values("CHR")

    output = REFERENCE_ROOT / code / f"{code}_reference_manifest.tsv"
    combined.to_csv(output, sep="\t", index=False)

    return combined


# =============================================================================
# PLANNER
# =============================================================================

def panel_label(codes):

    return "_".join(codes) if len(codes) < len(POPULATIONS) else "ALL"


def manifest_file(codes):
    """One manifest per ancestry set, so a queued job is never overwritten."""

    return REFERENCE_ROOT / f"step04_manifest_{panel_label(codes)}.tsv"


def planner_mode(args):

    banner("STEP 04 - LD REFERENCE PLANNER")

    script_path = Path(__file__).resolve()
    python_path = PROJECT_ROOT / "envs" / "pipeline" / "bin" / "python"

    if not python_path.exists():
        python_path = Path(sys.executable).resolve()

    print(f"bcftools : {require_program('bcftools')}")
    print(f"plink2   : {require_program('plink2')}")

    rows = []

    for code in args.codes:

        sample_file, selected = create_sample_list(code)

        banner(f"{code} SAMPLE PANEL: {len(selected):,} people")
        print(selected["pop"].value_counts().sort_index().to_string())
        print(f"Sample list: {sample_file}")

        for chromosome in range(1, 23):

            if panel_complete(code, chromosome) and not args.force:
                continue

            filename = chromosome_filename(chromosome)
            rows.append({
                "CHR": chromosome,
                "ANCESTRY_CODE": code,
                "POPULATIONS": ",".join(POPULATIONS[code]),
                "SAMPLE_FILE": str(sample_file),
                "VCF_URL": f"{BASE_URL}/{filename}",
                "RAW_VCF": str(RAW_DIR / filename),
                "RAW_TBI": str(RAW_DIR / (filename + ".tbi")),
                "OUTPUT_PREFIX": str(panel_prefix(code, chromosome)),
                "TMP_DIR": str(REFERENCE_ROOT / "tmp" / code),
            })

        collect_reference_summaries(code)

    if not rows:
        banner("ALL REQUESTED PANELS ARE COMPLETE")
        print("Nothing to run. Panel statistics: add --stats")
        return

    # Inputs must exist before submitting (Step04 never downloads).
    for row in rows:
        require_input(row["RAW_VCF"])
        require_input(row["RAW_TBI"])

    task_ids, n_tasks = assign_task_ids(len(rows))
    manifest = pd.DataFrame(rows)
    manifest.insert(0, "TASK_ID", task_ids)
    manifest_file(args.codes).parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_file(args.codes), sep="\t", index=False)

    log_dir = PROJECT_ROOT / "logs" / "step04_ld_reference"
    log_dir.mkdir(parents=True, exist_ok=True)

    label = panel_label(args.codes)
    bash_file = PROJECT_ROOT / f"Step04_Prepare_1000G_{label}_GRCh38.sh"
    header = sbatch_header_block(args, f"1000G_{label}", f"{log_dir}/ld.%A_%a.out",
                                 f"{log_dir}/ld.%A_%a.err", n_tasks)

    bash_file.write_text(f"""#!/bin/bash --login
{header}

set -euo pipefail

ROOT="{PROJECT_ROOT}"
cd "$ROOT"
export PATH="$ROOT/envs/pipeline/bin:$PATH"
export LD_REFERENCE_MANIFEST="{manifest_file(args.codes)}"
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-{args.cpus}}}"

echo "STEP 04 - 1000G LD reference | job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID on $(hostname)"

"{python_path}" "{script_path}"
""", encoding="utf-8", newline="\n")
    bash_file.chmod(0o755)

    banner("STEP 04 BASH GENERATED")
    print(f"Ancestries          : {', '.join(args.codes)}")
    print(f"Chromosome tasks    : {len(rows)} (array 1-{n_tasks})")
    print(f"Already complete    : {22 * len(args.codes) - len(rows)}")
    print(f"Panels              : {REFERENCE_ROOT}/<ANC>/chr<N>_<ANC>_GRCh38.pgen/.pvar/.psam")
    print(f"Manifest            : {manifest_file(args.codes)}")
    print()
    print("RUN:")
    print(f"  sbatch {bash_file.name}")
    print("THEN:")
    print(f"  envs/pipeline/bin/python {script_path.name} --ancestry {' '.join(args.codes)} --status")


# =============================================================================
# WORKER (one ancestry x chromosome per manifest row)
# =============================================================================

def process_chromosome(row, threads):

    chromosome = int(row["CHR"])
    code = str(row["ANCESTRY_CODE"])
    sample_file = Path(row["SAMPLE_FILE"])
    raw_vcf = Path(row["RAW_VCF"])
    raw_tbi = Path(row["RAW_TBI"])
    output_prefix = Path(row["OUTPUT_PREFIX"])
    temporary_directory = Path(row["TMP_DIR"])

    temporary_directory.mkdir(parents=True, exist_ok=True)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    pgen, pvar, psam = panel_files(code, chromosome)
    summary_file = output_prefix.parent / f"chr{chromosome}_{code}_summary.tsv"

    def write_summary():
        pd.DataFrame([{
            "CHR": chromosome,
            "ANCESTRY": code,
            "SAMPLES": count_nonheader_lines(psam),
            "VARIANTS": count_nonheader_lines(pvar),
            "RAW_VCF": str(raw_vcf),
            "PGEN": str(pgen),
            "PVAR": str(pvar),
            "PSAM": str(psam),
            "COMPLETED_UTC": datetime.now(timezone.utc).isoformat(),
        }]).to_csv(summary_file, sep="\t", index=False)

    if panel_complete(code, chromosome):
        banner(f"REFERENCE ALREADY COMPLETE: chr{chromosome} {code}")
        if not summary_file.exists():
            write_summary()
        return

    bcftools = require_program("bcftools")
    plink2 = require_program("plink2")
    require_input(raw_vcf)
    require_input(raw_tbi)

    banner(f"STEP 04 WORKER - chr{chromosome} {code} ({threads} threads)")

    # 1. ancestry subset
    ancestry_vcf = temporary_directory / f"chr{chromosome}.{code}.GRCh38.vcf.gz"

    if not (ancestry_vcf.exists() and ancestry_vcf.stat().st_size > 0):
        temporary_vcf = Path(str(ancestry_vcf) + ".part")
        temporary_vcf.unlink(missing_ok=True)
        run_command([bcftools, "view", "--threads", threads, "-S", sample_file,
                     "--force-samples", "-Oz", "-o", temporary_vcf, raw_vcf])
        if not temporary_vcf.exists() or temporary_vcf.stat().st_size == 0:
            raise RuntimeError(f"Failed to create ancestry VCF for chr{chromosome}.")
        temporary_vcf.replace(ancestry_vcf)

    # 2. PLINK2 (incomplete previous output removed first)
    for path in (pgen, pvar, psam):
        path.unlink(missing_ok=True)

    run_command([plink2, "--vcf", ancestry_vcf,
                 "--set-missing-var-ids", "@:#:$r:$a",
                 "--new-id-max-allele-len", "1000",
                 "--rm-dup", "force-first",
                 "--make-pgen", "--threads", threads, "--out", output_prefix])

    for path in (pgen, pvar, psam):
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"PLINK2 output missing: {path}")

    write_summary()

    ancestry_vcf.unlink(missing_ok=True)
    Path(str(ancestry_vcf) + ".tbi").unlink(missing_ok=True)
    failure_file(code, chromosome).unlink(missing_ok=True)

    banner("CHROMOSOME COMPLETE")
    print(f"chr{chromosome} {code}: {count_nonheader_lines(psam):,} samples, "
          f"{count_nonheader_lines(pvar):,} variants")


def worker_mode():

    task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    manifest = pd.read_csv(os.environ["LD_REFERENCE_MANIFEST"], sep="\t", dtype=str)
    selected = manifest[pd.to_numeric(manifest["TASK_ID"], errors="coerce") == task_id]

    if selected.empty:
        raise RuntimeError(f"No manifest rows for array task {task_id}")

    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", CPUS_PER_TASK))
    failed = []

    for _, row in selected.iterrows():
        try:
            process_chromosome(row, threads)
        except Exception as error:
            code, chromosome = row["ANCESTRY_CODE"], int(row["CHR"])
            failure_file(code, chromosome).write_text(
                f"{type(error).__name__}: {error}\n\n{traceback.format_exc()}", encoding="utf-8")
            print(f"[FAILED] chr{chromosome} {code}: {error}", file=sys.stderr)
            failed.append(f"chr{chromosome} {code}")

    for code in selected["ANCESTRY_CODE"].unique():
        collect_reference_summaries(code)

    if failed:
        raise SystemExit(f"Failed: {', '.join(failed)}")


def direct_index_mode(args):
    """--index N: chromosome N of every requested ancestry, built on this machine
    (no planner needed; a finished panel is skipped)."""

    if not 1 <= args.index <= 22:
        raise SystemExit("For Step04, --index is a chromosome number (1-22).")

    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", args.cpus))

    for code in args.codes:
        sample_file, _ = create_sample_list(code)
        filename = chromosome_filename(args.index)
        process_chromosome({
            "CHR": args.index,
            "ANCESTRY_CODE": code,
            "SAMPLE_FILE": str(sample_file),
            "RAW_VCF": str(RAW_DIR / filename),
            "RAW_TBI": str(RAW_DIR / (filename + ".tbi")),
            "OUTPUT_PREFIX": str(panel_prefix(code, args.index)),
            "TMP_DIR": str(REFERENCE_ROOT / "tmp" / code),
        }, threads)
        collect_reference_summaries(code)

    print(f"\nchr{args.index} done for {', '.join(args.codes)}. Progress: "
          f"envs/pipeline/bin/python {Path(__file__).name} --ancestry {' '.join(args.codes)} --status")


# =============================================================================
# --status / --stats
# =============================================================================

def status_mode(args):

    banner("STEP 04 STATUS - 1000G LD reference panels")
    print(f"{'ANC':<5}{'DONE':>6}{'FAILED':>8}{'PENDING':>9}  DETAIL")

    for code in args.codes:
        done = [c for c in range(1, 23) if panel_complete(code, c)]
        failed = [c for c in range(1, 23) if c not in done and failure_file(code, c).is_file()]
        pending = [c for c in range(1, 23) if c not in done and c not in failed]
        detail = ""
        if failed:
            detail = "failed: " + ",".join(f"chr{c}" for c in failed)
        elif pending:
            detail = "pending: " + ",".join(f"chr{c}" for c in pending[:8]) + (" ..." if len(pending) > 8 else "")
        print(f"{code:<5}{len(done):>6}{len(failed):>8}{len(pending):>9}  {detail}")

    print()
    print(f"Panels: {REFERENCE_ROOT}/<ANC>/   (Step05 / Step06 need all 22 chromosomes)")


def stats_mode(args):

    banner("STEP 04 STATISTICS - 1000G LD reference panels")

    rows = []

    for code in args.codes:

        table = collect_reference_summaries(code)
        sample_file = REFERENCE_ROOT / code / "samples.txt"

        if PANEL_FILE.is_file():
            panel = pd.read_csv(PANEL_FILE, sep="\t", dtype=str)
            pops = panel[panel["pop"].isin(POPULATIONS[code])]["pop"].value_counts().sort_index()
            print(f"\n{code}: {int(pops.sum()):,} people  ("
                  + ", ".join(f"{p} {n}" for p, n in pops.items()) + ")")

        if table is None:
            print("  no completed chromosomes yet")
            continue

        print(f"  chromosomes complete : {len(table)}/22")
        print(f"  samples in panel     : {int(table['SAMPLES'].max()):,}")
        print(f"  variants (total)     : {int(table['VARIANTS'].sum()):,}")
        print("  variants per chromosome:")
        print("   " + "  ".join(f"chr{int(r.CHR)} {int(r.VARIANTS):,}" for r in table.itertuples()))

        rows.append({"ANCESTRY": code, "CHROMOSOMES": len(table),
                     "SAMPLES": int(table["SAMPLES"].max()),
                     "VARIANTS": int(table["VARIANTS"].sum()),
                     "SAMPLE_FILE": str(sample_file)})

    if rows:
        output = REFERENCE_ROOT / "step04_reference_stats.tsv"
        pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
        print(f"\nTable: {output}")


# =============================================================================
# MAIN
# =============================================================================

def arguments(argv=None):

    parser = argparse.ArgumentParser(description="Build ancestry-specific 1000G GRCh38 PLINK2 LD panels.")
    parser.add_argument("--ancestry", nargs="+", required=True,
                        help="EUR AFR EAS SAS AMR (or full names), or 'all'.")
    parser.add_argument("--phenotype", default=None,
                        help="Accepted for symmetry with Step05/06; panels do not depend on it.")
    parser.add_argument("--force", action="store_true", help="Rebuild panels that already exist.")

    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--status", action="store_true", help="READ-ONLY: built / failed / pending panels.")
    modes.add_argument("--stats", action="store_true", help="Samples per population and variants per chromosome.")
    modes.add_argument("--index", type=int, default=None,
                       help="Build chromosome N (1-22) on this machine; no planner needed.")

    add_slurm_arguments(parser)

    args = parser.parse_args(argv)

    if any(x.lower() == "all" for x in args.ancestry):
        args.codes = list(POPULATIONS)
    else:
        args.codes = list(dict.fromkeys(ancestry_code(x) for x in args.ancestry))

    return resolve_slurm(args)


def main():

    if os.environ.get("SLURM_ARRAY_TASK_ID") and os.environ.get("LD_REFERENCE_MANIFEST"):
        worker_mode()
        return

    args = arguments()

    if args.status:
        status_mode(args)
    elif args.stats:
        stats_mode(args)
    elif args.index is not None:
        direct_index_mode(args)
    else:
        planner_mode(args)


if __name__ == "__main__":

    main()
