#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 04 - Ancestry-specific GRCh38 LD reference panels
===============================================================================

Builds the PLINK2 panels that Step05 (clumping) and Step06 (fine-mapping LD)
read. Panels are shared by every phenotype, so they are built once per
reference panel x ancestry, inside the project. Choose the panel with
--reference-panel (the same name is then given to Step05 / Step06):

    1000G      (default) 1000 Genomes phase 3 on GRCh38, biallelic SNV/INDEL
               (2019), 2504 unrelated samples
                 inputs : Step00_Download_Resources.py --groups 1000g
                 panels : resources/1000G/<ANC>/chr<N>_<ANC>_GRCh38.pgen/.pvar/.psam
    1000G30x   1000 Genomes 30x on GRCh38 (NYGC; Byrska-Bishop et al. 2022),
               filtered phased SNV/INDEL/SV panel, the 2504 unrelated of the
               3202 samples (related samples would distort LD)
                 inputs : Step00_Download_Resources.py --groups 1000g30x
                 panels : resources/LD_panels/1000G30x/<ANC>/chr<N>_<ANC>_GRCh38.*
    <NAME>     your own panel, resources/LD_panels/<NAME>/<ANC>/, either
                 - built here from per-chromosome VCFs (bgzipped + .tbi/.csi):
                   --panel-vcf "/path/chr{chr}.vcf.gz" --panel-samples samples.tsv
                   (samples.tsv: sample ID, then ancestry code EUR/AFR/EAS/SAS/AMR or a
                   1000G population code; header optional). Both are saved in
                   resources/LD_panels/<NAME>/panel_config.json, so later runs only
                   need --reference-panel <NAME>.
                 - or ready-made: put chr<N>_<ANC>_GRCh38.pgen/.pvar/.psam there.

Per panel and ancestry:
    <panel>/<ANC>/samples.txt                  (sample IDs used)
    <panel>/<ANC>/<ANC>_reference_manifest.tsv (samples + variants per chromosome)

Usage:

    # plan + SLURM array (one task per ancestry x chromosome)
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR
    sbatch Step04_Prepare_1000G_EUR_GRCh38.sh

    # ---- 1000 Genomes 30x panel (inputs: Step00 --groups 1000g30x) ----
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --reference-panel 1000G30x
    sbatch Step04_Prepare_1000G30x_EUR_GRCh38.sh
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --reference-panel 1000G30x --status
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --reference-panel 1000G30x --stats
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --reference-panel 1000G30x --verify
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --reference-panel 1000G30x --index 22

    # ---- a custom panel ----
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --reference-panel MyPanel \\
        --panel-vcf "/data/mypanel/chr{chr}.vcf.gz" --panel-samples /data/mypanel/samples.tsv

    # several ancestries / all five super-populations
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR EAS
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry all

    # progress (read-only), panel statistics, full PLINK2 validation, panel list
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --status
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --stats
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --verify
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --list-panels

    # chromosome 22 directly on this machine (no planner needed)
    envs/pipeline/bin/python Step04_Prepare_LD_Reference.py --ancestry EUR --index 22

Each chromosome: bcftools subsets the ancestry samples (a missing sample ID is
an error), plink2 converts to PGEN, keeps the VCF IDs, names variants without an
ID CHR:POS:REF:ALT and removes duplicate IDs. REF/ALT are the panel's GRCh38
alleles; Step05/06 match GWAS alleles to them explicitly and orient Z to REF.

Note: the 1000G EUR panel has ~503 people (CEU, FIN, GBR, IBS, TSI) in both
1000G and 1000G30x. That is a standard, public LD reference but a small one; for
very large GWAS a larger ancestry-matched panel (or in-sample LD) gives more
reliable fine-mapping - pass it as a custom panel.

SLURM settings come from the #SBATCH header of HPC_CPU.sh; --partition --qos
--account --cpus --mem --time --max-parallel override them.
===============================================================================
"""

from __future__ import annotations

import argparse
import json
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

LD_PANELS_ROOT = PROJECT_ROOT / "resources" / "LD_panels"
THIRTY_X_ROOT = LD_PANELS_ROOT / "1000G30x"
THIRTY_X_PED = THIRTY_X_ROOT / "metadata" / "20130606_g1k_3202_samples_ped_population.txt"

DEFAULT_PANEL = "1000G"

# Built-in reference panels. vcf: per-chromosome input, {chr} = 1..22.
BUILTIN_PANELS = {
    "1000G": {
        "description": "1000 Genomes phase 3, GRCh38 biallelic SNV/INDEL (2019), 2504 unrelated",
        "root": REFERENCE_ROOT,
        "vcf": str(RAW_DIR / "ALL.chr{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"),
        "samples": "phase3",
        "download_group": "1000g",
        # Biallelic release: a duplicate ID with different genotypes is an error.
        "rm_dup": "error",
    },
    "1000G30x": {
        "description": "1000 Genomes 30x (NYGC 2022), phased SNV/INDEL/SV, 2504 unrelated of 3202",
        "root": THIRTY_X_ROOT,
        "vcf": str(THIRTY_X_ROOT / "raw"
                   / "1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        "samples": "30x",
        "download_group": "1000g30x",
        # SNV/INDEL and SV calls can share an ID: conflicting copies are dropped
        # (listed in .rmdup.mismatch) instead of failing the chromosome.
        "rm_dup": "exclude-mismatch",
    },
}

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
            f"Ancestry {value!r} has no LD panel. "
            "Use EUR, AFR, EAS, SAS or AMR (or 'all')."
        )

    return code


# =============================================================================
# REFERENCE PANELS
# =============================================================================

def panel_config_file(name):

    return LD_PANELS_ROOT / name / "panel_config.json"


def resolve_panel(name, vcf=None, samples=None):
    """Definition of a reference panel: built-in, or custom under
    resources/LD_panels/<name>/ (build inputs from --panel-vcf / --panel-samples,
    saved in panel_config.json; without them, ready-made PLINK2 files)."""

    if name in BUILTIN_PANELS:
        if vcf or samples:
            raise SystemExit(f"--panel-vcf / --panel-samples are for custom panels, not {name}.")
        return {"name": name, **BUILTIN_PANELS[name]}

    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", str(name)):
        raise SystemExit(f"Panel name {name!r}: use letters, digits, '_', '-' or '.'")

    config_file = panel_config_file(name)
    config = json.loads(config_file.read_text(encoding="utf-8")) if config_file.is_file() else {}

    if vcf or samples:
        if not (vcf and samples):
            raise SystemExit("A custom panel built here needs both --panel-vcf and --panel-samples.")
        if "{chr}" not in vcf:
            raise SystemExit('--panel-vcf must contain {chr}, e.g. "/data/mypanel/chr{chr}.vcf.gz"')
        config = {
            "vcf": str(Path(vcf).expanduser().resolve()),
            "samples": str(Path(samples).expanduser().resolve()),
            "created_utc": datetime.now(timezone.utc).isoformat(),
        }
        config_file.parent.mkdir(parents=True, exist_ok=True)
        config_file.write_text(json.dumps(config, indent=2), encoding="utf-8")

    return {
        "name": name,
        "description": f"custom panel ({config.get('vcf') or 'ready-made PLINK2 files'})",
        "root": LD_PANELS_ROOT / name,
        "vcf": config.get("vcf"),
        "samples": config.get("samples"),
        "download_group": None,
        "rm_dup": "exclude-mismatch",
    }


def panel_vcf(panel, chromosome):

    return Path(str(panel["vcf"]).replace("{chr}", str(chromosome)))


def panel_index(vcf):
    """The .tbi (or .csi) next to a VCF."""

    for suffix in (".tbi", ".csi"):
        if Path(str(vcf) + suffix).is_file():
            return Path(str(vcf) + suffix)

    return Path(str(vcf) + ".tbi")


def panel_prefix(panel, code, chromosome):

    return Path(panel["root"]) / code / f"chr{chromosome}_{code}_GRCh38"


def panel_files(panel, code, chromosome):

    prefix = panel_prefix(panel, code, chromosome)

    return [Path(str(prefix) + ext) for ext in (".pgen", ".pvar", ".psam")]


def panel_complete(panel, code, chromosome):
    """Lightweight structural check; publication release also requires --verify."""
    files = panel_files(panel, code, chromosome)
    if not all(p.is_file() and p.stat().st_size > 0 for p in files):
        return False
    sample_file = Path(panel["root"]) / code / "samples.txt"
    if sample_file.is_file():
        try:
            with sample_file.open(encoding="utf-8") as fh:
                expected = sum(bool(line.strip()) for line in fh)
            if expected == 0 or count_nonheader_lines(files[2]) != expected:
                return False
        except OSError:
            return False
    return True


def failure_file(panel, code, chromosome):

    return Path(panel["root"]) / code / f"chr{chromosome}_{code}_FAILED.txt"


def require_program(program):

    local = PROJECT_ROOT / "envs" / "pipeline" / "bin" / program
    path = str(local) if local.exists() else shutil.which(program)

    if path is None:
        raise RuntimeError(
            f"Required program is not available: {program}\n"
            "Run: python Step00_Installation.py"
        )

    return path


def require_input(path, panel=None):

    path = Path(path)

    if path.is_file() and path.stat().st_size > 0:
        return path

    group = (panel or {}).get("download_group", "1000g") if panel is not None else "1000g"
    script = f"Step00_Download_{(panel or {}).get('name', DEFAULT_PANEL)}.sh"

    if group:
        raise FileNotFoundError(
            f"Missing reference-panel input: {path}\n"
            "Step04 does not download shared data. Run:\n"
            f"  envs/pipeline/bin/python Step00_Download_Resources.py --groups {group} "
            f"--output {script}\n"
            f"  sbatch {script}"
        )

    raise FileNotFoundError(
        f"Missing input of custom panel {panel['name']}: {path}\n"
        f"Check --panel-vcf / --panel-samples (saved in {panel_config_file(panel['name'])})."
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

def ancestry_directory(panel, code):
    """<panel root>/<code>/ as a usable directory (a symlink to a directory is
    fine). A dangling symlink, e.g. left by an older setup, holds no data and is
    replaced; anything else that is not a directory stops with instructions."""

    path = Path(panel["root"]) / code

    if path.is_symlink() and not path.exists():
        print(f"Replacing broken symlink {path} -> {os.readlink(path)} with a directory")
        path.unlink()
    elif path.exists() and not path.is_dir():
        raise SystemExit(
            f"{path} exists but is not a directory, so the {code} panels cannot be written there.\n"
            f"Look at it with:  ls -la {path}\n"
            f"then move it away, e.g.:  mv {path} {path}.old"
        )

    path.mkdir(parents=True, exist_ok=True)

    return path


def phase3_samples():
    """The 2504 unrelated 1000 Genomes phase 3 samples (sample, pop)."""

    panel = pd.read_csv(require_input(PANEL_FILE), sep="\t", dtype=str, low_memory=False)

    for column in ("sample", "pop"):
        if column not in panel.columns:
            raise RuntimeError(f"1000 Genomes panel missing column: {column}")

    return panel


def custom_samples(path):
    """Custom sample table: sample ID, then ancestry code or population (header optional)."""

    table = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, engine="python", comment="#")

    if table.shape[1] < 2:
        raise SystemExit(f"{path}: need two columns (sample ID, ancestry code or population)")

    first = [str(x).strip().upper() for x in table.iloc[0, :2]]
    if first[0] in ("SAMPLE", "IID", "ID", "SAMPLE_ID", "SAMPLEID") or \
            first[1] in ("ANCESTRY", "SUPERPOPULATION", "SUPER_POP", "POP", "POPULATION", "GROUP"):
        table = table.iloc[1:]

    return pd.DataFrame({"sample": table.iloc[:, 0].str.strip(), "pop": table.iloc[:, 1].str.strip()})


def panel_samples(panel, code):

    if panel["samples"] == "phase3":
        samples = phase3_samples()
        return samples[samples["pop"].isin(POPULATIONS[code])].copy()

    if panel["samples"] == "30x":
        ped = pd.read_csv(require_input(THIRTY_X_PED, panel), sep=r"\s+", dtype=str, engine="python")
        unrelated = set(phase3_samples()["sample"])
        selected = ped[ped["SampleID"].isin(unrelated) & ped["Population"].isin(POPULATIONS[code])]
        return pd.DataFrame({"sample": selected["SampleID"], "pop": selected["Population"]})

    if not panel["samples"]:
        raise SystemExit(
            f"Custom panel {panel['name']} has no build inputs. Pass --panel-vcf and --panel-samples, "
            f"or put chr<N>_{code}_GRCh38.pgen/.pvar/.psam in {Path(panel['root']) / code}/"
        )

    samples = custom_samples(require_input(panel["samples"], panel))
    label = samples["pop"].str.upper()

    return samples[(label == code) | samples["pop"].isin(POPULATIONS[code])].copy()


def create_sample_list(panel, code):

    selected = panel_samples(panel, code)

    if selected.empty:
        raise RuntimeError(f"No {code} samples in panel {panel['name']}")

    output = ancestry_directory(panel, code) / "samples.txt"
    selected[["sample"]].to_csv(output, index=False, header=False)

    return output, selected


def collect_reference_summaries(panel, code):

    directory = Path(panel["root"]) / code
    files = sorted(directory.glob(f"chr*_{code}_summary.tsv"))
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

    output = directory / f"{code}_reference_manifest.tsv"
    combined.to_csv(output, sep="\t", index=False)

    return combined


# =============================================================================
# PLANNER
# =============================================================================

def panel_label(codes):

    return "_".join(codes) if len(codes) < len(POPULATIONS) else "ALL"


def manifest_file(panel, codes):
    """Unique snapshot per planning run: running jobs never read a rewritten manifest."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return Path(panel["root"]) / f"step04_manifest_{panel_label(codes)}_{stamp}_{os.getpid()}.tsv"


def chromosome_row(panel, code, chromosome, sample_file):

    vcf = panel_vcf(panel, chromosome)

    return {
        "CHR": chromosome,
        "ANCESTRY_CODE": code,
        "REFERENCE_PANEL": panel["name"],
        "POPULATIONS": ",".join(POPULATIONS[code]),
        "SAMPLE_FILE": str(sample_file),
        "RAW_VCF": str(vcf),
        "RAW_TBI": str(panel_index(vcf)),
        "OUTPUT_PREFIX": str(panel_prefix(panel, code, chromosome)),
        "TMP_DIR": str(Path(panel["root"]) / "tmp" / code),
    }


def planner_mode(args):

    panel = args.panel

    banner(f"STEP 04 - LD REFERENCE PLANNER ({panel['name']})")

    script_path = Path(__file__).resolve()
    python_path = PROJECT_ROOT / "envs" / "pipeline" / "bin" / "python"

    if not python_path.exists():
        python_path = Path(sys.executable).resolve()

    print(f"Panel    : {panel['name']} - {panel['description']}")
    print(f"Location : {panel['root']}/<ANC>/")
    print(f"bcftools : {require_program('bcftools')}")
    print(f"plink2   : {require_program('plink2')}")

    rows = []

    for code in args.codes:

        if not panel["vcf"]:
            # Ready-made custom panel: nothing to build.
            missing = [c for c in range(1, 23) if not panel_complete(panel, code, c)]
            if missing:
                panel_samples(panel, code)   # stops with instructions
            print(f"\n{code}: ready-made panel, all 22 chromosomes present")
            continue

        sample_file, selected = create_sample_list(panel, code)

        banner(f"{code} SAMPLE PANEL ({panel['name']}): {len(selected):,} people")
        print(selected["pop"].value_counts().sort_index().to_string())
        print(f"Sample list: {sample_file}")

        for chromosome in range(1, 23):

            if panel_complete(panel, code, chromosome) and not args.force:
                continue

            rows.append(chromosome_row(panel, code, chromosome, sample_file))

        collect_reference_summaries(panel, code)

    if not rows:
        banner("ALL REQUESTED PANELS ARE COMPLETE")
        print("Nothing to run. Panel statistics: add --stats")
        return

    # Inputs must exist before submitting (Step04 never downloads).
    for row in rows:
        require_input(row["RAW_VCF"], panel)
        require_input(row["RAW_TBI"], panel)

    task_ids, n_tasks = assign_task_ids(len(rows))
    manifest = pd.DataFrame(rows)
    manifest.insert(0, "TASK_ID", task_ids)
    manifest_path = manifest_file(panel, args.codes)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, sep="\t", index=False)

    log_dir = PROJECT_ROOT / "logs" / "step04_ld_reference" / panel["name"]
    log_dir.mkdir(parents=True, exist_ok=True)

    label = panel_label(args.codes)
    bash_file = PROJECT_ROOT / f"Step04_Prepare_{panel['name']}_{label}_GRCh38.sh"
    header = sbatch_header_block(args, f"{panel['name']}_{label}", f"{log_dir}/ld.%A_%a.out",
                                 f"{log_dir}/ld.%A_%a.err", n_tasks)

    bash_file.write_text(f"""#!/bin/bash --login
{header}

set -euo pipefail

ROOT="{PROJECT_ROOT}"
cd "$ROOT"
export PATH="$ROOT/envs/pipeline/bin:$PATH"
export LD_REFERENCE_MANIFEST="{manifest_path}"
export LD_FORCE_REBUILD="{int(args.force)}"
export OMP_NUM_THREADS="${{SLURM_CPUS_PER_TASK:-{args.cpus}}}"

echo "STEP 04 - {panel['name']} LD reference | job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID on $(hostname)"

"{python_path}" "{script_path}"
""", encoding="utf-8", newline="\n")
    bash_file.chmod(0o755)

    option = "" if panel["name"] == DEFAULT_PANEL else f" --reference-panel {panel['name']}"
    banner("STEP 04 BASH GENERATED")
    print(f"Reference panel     : {panel['name']}")
    print(f"Ancestries          : {', '.join(args.codes)}")
    print(f"Chromosome tasks    : {len(rows)} (array 1-{n_tasks})")
    print(f"Already complete    : {22 * len(args.codes) - len(rows)}")
    print(f"Panels              : {panel['root']}/<ANC>/chr<N>_<ANC>_GRCh38.pgen/.pvar/.psam")
    print(f"Manifest            : {manifest_path}")
    print()
    print("RUN:")
    print(f"  sbatch {bash_file.name}")
    print("THEN:")
    print(f"  envs/pipeline/bin/python {script_path.name} --ancestry {' '.join(args.codes)}{option} --status")
    if option:
        print(f"Step05 / Step06 use this panel with{option}")


# =============================================================================
# WORKER (one ancestry x chromosome per manifest row)
# =============================================================================

def process_chromosome(row, threads, force=False):

    panel = resolve_panel(str(row.get("REFERENCE_PANEL") or DEFAULT_PANEL))
    chromosome = int(row["CHR"])
    code = str(row["ANCESTRY_CODE"])
    sample_file = Path(row["SAMPLE_FILE"])
    raw_vcf = Path(row["RAW_VCF"])
    raw_tbi = Path(row["RAW_TBI"])
    output_prefix = Path(row["OUTPUT_PREFIX"])
    temporary_directory = Path(row["TMP_DIR"])

    temporary_directory.mkdir(parents=True, exist_ok=True)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    pgen, pvar, psam = panel_files(panel, code, chromosome)
    summary_file = output_prefix.parent / f"chr{chromosome}_{code}_summary.tsv"

    def write_summary():
        pd.DataFrame([{
            "CHR": chromosome,
            "ANCESTRY": code,
            "REFERENCE_PANEL": panel["name"],
            "SAMPLES": count_nonheader_lines(psam),
            "VARIANTS": count_nonheader_lines(pvar),
            "RAW_VCF": str(raw_vcf),
            "PGEN": str(pgen),
            "PVAR": str(pvar),
            "PSAM": str(psam),
            "COMPLETED_UTC": datetime.now(timezone.utc).isoformat(),
        }]).to_csv(summary_file, sep="\t", index=False)

    if panel_complete(panel, code, chromosome) and not force:
        banner(f"REFERENCE ALREADY COMPLETE: chr{chromosome} {code} ({panel['name']})")
        if not summary_file.exists():
            write_summary()
        return

    bcftools = require_program("bcftools")
    plink2 = require_program("plink2")
    require_input(raw_vcf, panel)
    require_input(raw_tbi, panel)

    banner(f"STEP 04 WORKER - chr{chromosome} {code} {panel['name']} ({threads} threads)")

    # 1. ancestry subset. Do not silently ignore missing sample IDs.
    ancestry_vcf = temporary_directory / f"chr{chromosome}.{code}.GRCh38.vcf.gz"

    if force:
        ancestry_vcf.unlink(missing_ok=True)
    if not (ancestry_vcf.exists() and ancestry_vcf.stat().st_size > 0):
        temporary_vcf = Path(str(ancestry_vcf) + ".part")
        temporary_vcf.unlink(missing_ok=True)
        run_command([bcftools, "view", "--threads", threads, "-S", sample_file,
                     "-Oz", "-o", temporary_vcf, raw_vcf])
        if not temporary_vcf.exists() or temporary_vcf.stat().st_size == 0:
            raise RuntimeError(f"Failed to create ancestry VCF for chr{chromosome}.")
        temporary_vcf.replace(ancestry_vcf)

    # 2. PLINK2 (incomplete previous output removed first). Long (SV) alleles are
    #    truncated in generated IDs only, never in the alleles.
    for path in (pgen, pvar, psam):
        path.unlink(missing_ok=True)

    run_command([plink2, "--vcf", ancestry_vcf,
                 "--set-missing-var-ids", "@:#:$r:$a",
                 "--new-id-max-allele-len", "1000", "truncate",
                 "--rm-dup", panel["rm_dup"], "list",
                 "--make-pgen", "--threads", threads, "--out", output_prefix])

    for path in (pgen, pvar, psam):
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"PLINK2 output missing: {path}")

    expected_samples = sum(bool(line.strip()) for line in sample_file.open(encoding="utf-8"))
    observed_samples = count_nonheader_lines(psam)
    if expected_samples != observed_samples:
        raise RuntimeError(f"Sample mismatch: expected {expected_samples}, observed {observed_samples}")
    if count_nonheader_lines(pvar) == 0:
        raise RuntimeError("PLINK2 output has zero variants")
    # Complete genotype file validation, not merely a nonempty-file check.
    run_command([plink2, "--pfile", output_prefix, "--validate", "--out", str(output_prefix) + "_validation"])
    write_summary()

    ancestry_vcf.unlink(missing_ok=True)
    Path(str(ancestry_vcf) + ".tbi").unlink(missing_ok=True)
    failure_file(panel, code, chromosome).unlink(missing_ok=True)

    banner("CHROMOSOME COMPLETE")
    print(f"chr{chromosome} {code} ({panel['name']}): {count_nonheader_lines(psam):,} samples, "
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
        panel = resolve_panel(str(row.get("REFERENCE_PANEL") or DEFAULT_PANEL))
        try:
            process_chromosome(row, threads, force=os.environ.get("LD_FORCE_REBUILD") == "1")
        except Exception as error:
            code, chromosome = row["ANCESTRY_CODE"], int(row["CHR"])
            failure_file(panel, code, chromosome).write_text(
                f"{type(error).__name__}: {error}\n\n{traceback.format_exc()}", encoding="utf-8")
            print(f"[FAILED] chr{chromosome} {code} ({panel['name']}): {error}", file=sys.stderr)
            failed.append(f"chr{chromosome} {code}")

    for code in selected["ANCESTRY_CODE"].unique():
        panel_name = selected.loc[selected["ANCESTRY_CODE"] == code, "REFERENCE_PANEL"].iloc[0] \
            if "REFERENCE_PANEL" in selected.columns else DEFAULT_PANEL
        collect_reference_summaries(resolve_panel(str(panel_name or DEFAULT_PANEL)), code)

    if failed:
        raise SystemExit(f"Failed: {', '.join(failed)}")


def direct_index_mode(args):
    """--index N: chromosome N of every requested ancestry, built on this machine
    (no planner needed; a finished panel is skipped)."""

    if not 1 <= args.index <= 22:
        raise SystemExit("For Step04, --index is a chromosome number (1-22).")

    panel = args.panel
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", args.cpus))

    for code in args.codes:
        sample_file, _ = create_sample_list(panel, code)
        process_chromosome(chromosome_row(panel, code, args.index, sample_file), threads, force=args.force)
        collect_reference_summaries(panel, code)

    option = "" if panel["name"] == DEFAULT_PANEL else f" --reference-panel {panel['name']}"
    print(f"\nchr{args.index} done for {', '.join(args.codes)} ({panel['name']}). Progress: "
          f"envs/pipeline/bin/python {Path(__file__).name} --ancestry {' '.join(args.codes)}{option} --status")


# =============================================================================
# --status / --stats / --verify / --list-panels
# =============================================================================

def status_mode(args):

    panel = args.panel
    banner(f"STEP 04 STATUS - {panel['name']} LD reference panels")
    print(f"{panel['description']}")
    print(f"{'ANC':<5}{'DONE':>6}{'FAILED':>8}{'PENDING':>9}  DETAIL")

    for code in args.codes:
        done = [c for c in range(1, 23) if panel_complete(panel, code, c)]
        failed = [c for c in range(1, 23) if c not in done and failure_file(panel, code, c).is_file()]
        pending = [c for c in range(1, 23) if c not in done and c not in failed]
        detail = ""
        if failed:
            detail = "failed: " + ",".join(f"chr{c}" for c in failed)
        elif pending:
            detail = "pending: " + ",".join(f"chr{c}" for c in pending[:8]) + (" ..." if len(pending) > 8 else "")
        print(f"{code:<5}{len(done):>6}{len(failed):>8}{len(pending):>9}  {detail}")

    print()
    print(f"Panels: {panel['root']}/<ANC>/   (Step05 / Step06 need all 22 chromosomes)")


def verify_mode(args):
    """Full read-only PLINK2 validation of all 22 filesets per ancestry."""
    panel = args.panel
    plink2 = require_program("plink2")
    for code in args.codes:
        sample_list = Path(panel["root"]) / code / "samples.txt"
        expected = sum(bool(line.strip()) for line in sample_list.open(encoding="utf-8")) \
            if sample_list.exists() else None
        for chromosome in range(1, 23):
            if not panel_complete(panel, code, chromosome):
                raise SystemExit(f"Missing or inconsistent panel chr{chromosome} {code} ({panel['name']})")
            prefix = panel_prefix(panel, code, chromosome)
            observed = count_nonheader_lines(Path(str(prefix) + ".psam"))
            if expected is not None and observed != expected:
                raise SystemExit(f"Sample mismatch chr{chromosome}: {observed} vs {expected}")
            run_command([plink2, "--pfile", prefix, "--validate", "--out", str(prefix) + "_validation"])
            print(f"[VERIFIED] {panel['name']} {code} chromosome {chromosome}: {observed} samples", flush=True)
    print(f"All requested {panel['name']} ancestry panels passed structural checks and PLINK2 validation.")


def stats_mode(args):

    panel = args.panel
    banner(f"STEP 04 STATISTICS - {panel['name']} LD reference panels")
    print(panel["description"])

    rows = []

    for code in args.codes:

        table = collect_reference_summaries(panel, code)
        sample_file = Path(panel["root"]) / code / "samples.txt"

        try:
            pops = panel_samples(panel, code)["pop"].value_counts().sort_index()
            print(f"\n{code}: {int(pops.sum()):,} people  ("
                  + ", ".join(f"{p} {n}" for p, n in pops.items()) + ")")
        except (FileNotFoundError, SystemExit, RuntimeError):
            print(f"\n{code}:")

        if table is None:
            print("  no completed chromosomes yet")
            continue

        print(f"  chromosomes complete : {len(table)}/22")
        print(f"  samples in panel     : {int(table['SAMPLES'].max()):,}")
        print(f"  variants (total)     : {int(table['VARIANTS'].sum()):,}")
        print("  variants per chromosome:")
        print("   " + "  ".join(f"chr{int(r.CHR)} {int(r.VARIANTS):,}" for r in table.itertuples()))

        rows.append({"REFERENCE_PANEL": panel["name"], "ANCESTRY": code, "CHROMOSOMES": len(table),
                     "SAMPLES": int(table["SAMPLES"].max()),
                     "VARIANTS": int(table["VARIANTS"].sum()),
                     "SAMPLE_FILE": str(sample_file)})

    if rows:
        output = Path(panel["root"]) / "step04_reference_stats.tsv"
        pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
        print(f"\nTable: {output}")


def list_panels_mode(args):

    banner("STEP 04 - REFERENCE PANELS")
    names = list(BUILTIN_PANELS)
    if LD_PANELS_ROOT.is_dir():
        names += sorted(d.name for d in LD_PANELS_ROOT.iterdir() if d.is_dir() and d.name not in BUILTIN_PANELS)

    print(f"{'PANEL':<14}" + "".join(f"{c:>6}" for c in POPULATIONS) + "  DESCRIPTION")
    for name in names:
        panel = resolve_panel(name)
        done = [sum(panel_complete(panel, c, x) for x in range(1, 23)) for c in POPULATIONS]
        print(f"{name:<14}" + "".join(f"{d:>4}/22"[-6:] for d in done) + f"  {panel['description']}")
    print("\nchromosomes built per ancestry. Use one with --reference-panel <PANEL> in Step04 / Step05 / Step06.")


# =============================================================================
# MAIN
# =============================================================================

def arguments(argv=None):

    parser = argparse.ArgumentParser(description="Build ancestry-specific GRCh38 PLINK2 LD reference panels.")
    parser.add_argument("--ancestry", nargs="+", required=True,
                        help="EUR AFR EAS SAS AMR (or full names), or 'all'.")
    parser.add_argument("--reference-panel", default=DEFAULT_PANEL,
                        help="1000G (default), 1000G30x, or the name of a custom panel under resources/LD_panels/.")
    parser.add_argument("--panel-vcf", default=None,
                        help='Custom panel: per-chromosome VCFs with {chr}, e.g. "/data/panel/chr{chr}.vcf.gz".')
    parser.add_argument("--panel-samples", default=None,
                        help="Custom panel: sample ID + ancestry code (EUR, ...) or 1000G population per line.")
    parser.add_argument("--phenotype", default=None,
                        help="Accepted for symmetry with Step05/06; panels do not depend on it.")
    parser.add_argument("--force", action="store_true", help="Rebuild panels that already exist.")

    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--status", action="store_true", help="READ-ONLY: built / failed / pending panels.")
    modes.add_argument("--stats", action="store_true", help="Samples per population and variants per chromosome.")
    modes.add_argument("--verify", action="store_true", help="Validate all 22 panels with PLINK2 (use a compute node).")
    modes.add_argument("--list-panels", action="store_true", help="Every known reference panel and what is built.")
    modes.add_argument("--index", type=int, default=None,
                       help="Build chromosome N (1-22) on this machine; no planner needed.")

    add_slurm_arguments(parser)

    args = parser.parse_args(argv)

    if any(x.lower() == "all" for x in args.ancestry):
        args.codes = list(POPULATIONS)
    else:
        args.codes = list(dict.fromkeys(ancestry_code(x) for x in args.ancestry))

    args.panel = resolve_panel(args.reference_panel, args.panel_vcf, args.panel_samples)

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
    elif args.verify:
        verify_mode(args)
    elif args.list_panels:
        list_panels_mode(args)
    elif args.index is not None:
        direct_index_mode(args)
    else:
        planner_mode(args)


if __name__ == "__main__":

    main()
