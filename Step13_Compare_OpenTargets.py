#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2m - STEP 13
COMPARE STEP12 MASTER EVIDENCE WITH OPEN TARGETS GENETICS / PLATFORM BULK DATA
===============================================================================

GWAS2Mech PIPELINE USAGE (same flags in Steps 12-17)
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR                       # run on this machine
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR --slurm               # or write a one-job SLURM script
    sbatch Step13_OpenTargets_parkinson_s_disease_european.sh
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR --stats
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR --index 2             # study 2 only -> <ancestry>/<GCST>/
    # 1000 Genomes 30x LD panel: every input / output folder is <ancestry>_1000G30x
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --slurm
    sbatch Step13_OpenTargets_parkinson_s_disease_european_1000G30x.sh
    envs/pipeline/bin/python Step13_Compare_OpenTargets.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

PURPOSE
-------
Step13 is the external benchmarking / novelty-screen layer after Step12.
It does NOT alter Step12 evidence and does NOT call absence from Open Targets
"novel" by itself.

Inputs from Step12 v2:
  12_master_table/<phenotype>/<ancestry>/
    Step12_Locus_Gene_Master.tsv          (required)
    Step12_Locus_Master.tsv               (required)
    Step12_Variant_Master.tsv             (recommended)
    Step12_Signal_Master.tsv              (recommended)
    Step12_Molecular_Evidence_Long.tsv    (optional audit layer)
    Step12_Strong_Shared_Signals.tsv      (optional audit layer)
    Step12_Pathway_Input_Genes.tsv        (optional downstream link)

It also reads Step06 fine-mapping files when available to reconstruct complete
GWAS locus intervals / variant sets.

Open Targets bulk data (local cache):
  resources/opentargets/26.09/
    credible_set/       (required)
    l2g_prediction/     (required)
    colocalisation/     (optional but strongly recommended)
    study/              (optional; resolves molQTL study metadata)

By default, corrupt/truncated Open Targets shards are automatically retried via rsync.

The script compares:
  1. GWAS2m loci vs Open Targets credible sets from the SAME GWAS study
  2. method-aware lead-variant / credible-set / PIP agreement
     (SuSiE-style and PICS-style loci are interpreted differently)
  3. Step12 effector genes vs Open Targets L2G predictions
  4. Step12 formal coloc genes vs Open Targets molecular-QTL colocalisation
  5. Step12 bulk-tissue vs cell-type mechanisms vs Open Targets evidence
  6. individual Step12 formal signal rows vs the external Open Targets benchmark
  7. Open Targets coverage gaps vs genuine scientific disagreement

Missing Open Targets evidence is never treated as a biological negative.

Main output:
  Step13_Gene_Comparison.tsv

Conservative interpretation labels:
  A_ESTABLISHED_OT_GENE_AND_MOLQTL
  B_OT_SUPPORTED_GENE_POTENTIAL_NEW_MECHANISM
  C_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE
  D_POTENTIAL_NEW_LOCUS_CANDIDATE_REQUIRES_EXTERNAL_VALIDATION
  COVERAGE_OR_MAPPING_UNRESOLVED

IMPORTANT
---------
"Open Targets did not contain it" is NOT proof of novelty. Tier D is a
candidate for GWAS Catalog / literature / independent-replication review.

INSTALL
-------
  pip install duckdb pandas numpy scipy

RUN
---
  python Step13_Compare_OpenTargets.py \
      --phenotype "parkinson's disease" \
      --ancestry EUR

Use an existing cache elsewhere:
  python Step13_Compare_OpenTargets.py \
      --phenotype "parkinson's disease" \
      --ancestry EUR \
      --ot-cache /path/to/opentargets/26.09

Single study:
  python Step13_Compare_OpenTargets.py \
      --phenotype "parkinson's disease" \
      --ancestry EUR \
      --study GCST90319903
===============================================================================
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

try:
    import duckdb
except ImportError:
    duckdb = None

try:
    from scipy.stats import pearsonr, spearmanr
except ImportError:
    pearsonr = None
    spearmanr = None



# =============================================================================
# GWAS2Mech PIPELINE FRAME (same block in Steps 12-17)
#   --reference-panel  LD panel of Steps 04-06: results in <ancestry>_<PANEL>/ folders
#   --index N          only study N of the Step02 list: results in <ancestry>[_<PANEL>]/<GCST>/
#   --status           which inputs are there, are the results current, what to run next
#   --stats            results summary from the step's own output (nothing is recomputed)
#   --check-resources  inputs, Python modules and external resources; exit 1 if not ready
#   --slurm            write a one-job SLURM script (HPC_CPU.sh header) instead of running here
# =============================================================================

import argparse as _fargparse
import json as _fjson
import re as _fre
import shlex as _fshlex
import sys as _fsys
from pathlib import Path as _FPath

STEP_FRAME = {
    "NAME": "Step13", "TITLE": "comparison with Open Targets", "FOLDER": "13_opentargets_comparison",
    "SCRIPT_PREFIX": "Step13_OpenTargets", "SUMMARY": "Step13_summary.json",
    "PREVIOUS": "Step12_Master_Table.py",
    "INPUTS": [
        ("Step12", "12_master_table", ["Step12_Locus_Gene_Master.tsv"], True, True),
        ("Step12", "12_master_table", ["Step12_Locus_Master.tsv"], True, True),
        ("Step06", "06_finemapping", ["GCST*/*_finemapped_variants.tsv.gz"], False, False),
    ],
    "MODULES": [("duckdb", True, "reads the Open Targets parquet files"), ("scipy", False, "correlations")],
    "STATS_TABLES": [
        ("Step13_Gene_Comparison.tsv", ["STUDY_ACCESSION", "LOCUS_ID", "EFFECTOR_GENE_SYMBOL", "EFFECTOR_GENE_ID",
                                        "BEST_H4", "OPEN_TARGETS_EVIDENCE_CLASS", "NOVELTY_TIER_STEP13"], "BEST_H4",
         "GWAS2Mech GENES vs OPEN TARGETS"),
        ("Step13_Novelty_Tier_Summary.tsv", ["NOVELTY_TIER_STEP13", "N_GENES", "N_ROWS", "COUNT", "N"], "",
         "TIER SUMMARY"),
    ],
}


def frame_extra_checks(ns, P) -> list:
    """Open Targets bulk data (Step00_Download_Resources.py --groups opentargets)."""
    release = DEFAULT_OT_RELEASE
    for i, arg in enumerate(_fsys.argv):
        if arg == "--ot-release" and i + 1 < len(_fsys.argv):
            release = _fsys.argv[i + 1]
    cache = P["root"] / "resources" / "opentargets" / release
    hint = f"Step00_Download_Resources.py --groups opentargets --opentargets-release {release}"
    out = []
    for logical in DATASET_ALIASES:
        path = resolve_ot_dataset(cache, logical) if cache.exists() else None
        required = logical in {"credible_set", "l2g_prediction"}
        out.append((f"Open Targets {logical}", path is not None,
                    str(path) if path else f"missing under {cache} ({hint})", required))
    return out


DEFAULT_PANEL = "1000G"
PIPELINE: dict[str, str] = {"panel": DEFAULT_PANEL, "study": "", "index": ""}
_FRAME_ANCESTRY = {
    "eur": ("EUR", "European"), "european": ("EUR", "European"), "afr": ("AFR", "African"),
    "african": ("AFR", "African"), "eas": ("EAS", "East Asian"), "east asian": ("EAS", "East Asian"),
    "sas": ("SAS", "South Asian"), "south asian": ("SAS", "South Asian"),
    "amr": ("AMR", "Hispanic or Latin American"), "hispanic": ("AMR", "Hispanic or Latin American"),
    "latino": ("AMR", "Hispanic or Latin American"),
    "hispanic or latin american": ("AMR", "Hispanic or Latin American"),
}


def _frame_slug(value) -> str:
    return _fre.sub(r"\s+", "_", _fre.sub(r"[^a-z0-9]+", " ", str(value).lower()).strip())


def panel_slug(ancestry_slug: str) -> str:
    """Folder of Steps 04-10 results for the LD panel: <ancestry> (1000G) or <ancestry>_<PANEL>."""
    return ancestry_slug if PIPELINE["panel"] == DEFAULT_PANEL else f"{ancestry_slug}_{PIPELINE['panel']}"


def results_slug(ancestry_slug: str) -> str:
    """Folder of Steps 12-17 results: <ancestry>[_<PANEL>], plus /<GCST> with --index."""
    base = panel_slug(ancestry_slug)
    return f"{base}/{PIPELINE['study']}" if PIPELINE["study"] else base


FRAME_HELP = f"""
GWAS2Mech pipeline options ({STEP_FRAME['NAME']}):
  --reference-panel PANEL   LD panel of Steps 04-06 (default 1000G; e.g. 1000G30x): every input and
                            output folder is <ancestry>_<PANEL>
  --index N                 only study N of the Step02 list (results in <ancestry>[_<PANEL>]/<GCST>/)
  --status                  inputs present? results current? what to run next
  --stats                   summary of this step's results (from its output files)
  --check-resources         inputs, Python modules and external resources (exit 1 if not ready)
  --slurm                   write a one-job SLURM script (header of HPC_CPU.sh) instead of running
  --template --partition --qos --account --cpus --mem --time   SLURM overrides for --slurm
"""


def frame_arguments():
    pre = _fargparse.ArgumentParser(add_help=False)
    pre.add_argument("--phenotype")
    pre.add_argument("--ancestry")
    pre.add_argument("--reference-panel", default=DEFAULT_PANEL)
    pre.add_argument("--index", type=int, default=None)
    for flag in ("--status", "--stats", "--check-resources", "--slurm"):
        pre.add_argument(flag, action="store_true")
    pre.add_argument("--template", default="HPC_CPU.sh")
    for flag in ("--partition", "--qos", "--account", "--cpus", "--mem", "--time"):
        pre.add_argument(flag, default=None)
    ns, rest = pre.parse_known_args(_fsys.argv[1:])
    if ns.phenotype is not None:
        rest = ["--phenotype", ns.phenotype] + rest
    if ns.ancestry is not None:
        rest = ["--ancestry", ns.ancestry] + rest
    return ns, rest


def frame_paths(ns) -> dict:
    root = _FPath.cwd().resolve()
    key = _frame_slug(ns.ancestry).replace("_", " ")
    if key not in _FRAME_ANCESTRY:
        raise SystemExit(f"Unsupported ancestry {ns.ancestry!r}; use EUR, AFR, EAS, SAS or AMR.")
    code, label = _FRAME_ANCESTRY[key]
    p, a = _frame_slug(ns.phenotype), _frame_slug(label)
    if ns.index is not None:
        manifest = root / "01_gwas_catalog" / p / a / "GWAS_download_manifest.tsv"
        if not manifest.is_file():
            raise SystemExit(f"--index needs the Step02 download manifest: {manifest}")
        import pandas as _fpd
        m = _fpd.read_csv(manifest, sep="\t", dtype=str)
        hit = m[_fpd.to_numeric(m["ARRAY_ID"], errors="coerce") == ns.index]
        if hit.empty:
            raise SystemExit(f"--index {ns.index} is not in {manifest}")
        PIPELINE["study"], PIPELINE["index"] = str(hit.iloc[0]["STUDY_ACCESSION"]).strip(), str(ns.index)
    inputs = []
    for step, folder, patterns, required, upstream in STEP_FRAME["INPUTS"]:
        base = root / folder / p / (results_slug(a) if upstream else panel_slug(a))
        inputs.append({"STEP": step, "DIR": base, "PATTERNS": patterns, "REQUIRED": required, "UPSTREAM": upstream})
    return {"root": root, "p": p, "a": a, "code": code, "label": label,
            "OUT": root / STEP_FRAME["FOLDER"] / p / results_slug(a), "INPUTS": inputs,
            "LOG": root / "logs" / STEP_FRAME["FOLDER"] / p / results_slug(a).replace("/", "_")}


def _frame_matches(entry: dict) -> list:
    files = []
    for pattern in entry["PATTERNS"]:
        if PIPELINE["study"] and not entry["UPSTREAM"]:
            pattern = pattern.replace("GCST*", PIPELINE["study"])
        files += [f for f in entry["DIR"].glob(pattern) if f.is_file() and f.stat().st_size > 0]
    return files


def frame_command(ns) -> str:
    extra = (f" --reference-panel {ns.reference_panel}" if ns.reference_panel != DEFAULT_PANEL else "") + \
            (f" --index {ns.index}" if ns.index is not None else "")
    return (f"envs/pipeline/bin/python {_FPath(_fsys.argv[0]).name} --phenotype {_fshlex.quote(ns.phenotype)} "
            f"--ancestry {ns.ancestry}{extra}")


def frame_inputs_report(P: dict) -> list[str]:
    """Print the input inventory; return the missing required inputs."""
    missing = []
    print(f"{'INPUT':<9}{'STATE':<10}{'FILES':>7}  FOLDER")
    for entry in P["INPUTS"]:
        files = _frame_matches(entry)
        state = "OK" if files else ("MISSING" if entry["REQUIRED"] else "absent")
        print(f"{entry['STEP']:<9}{state:<10}{len(files):>7}  {entry['DIR']}"
              + ("" if entry["REQUIRED"] else "   (optional)"))
        if entry["REQUIRED"] and not files:
            missing.append(f"{entry['STEP']}: {entry['DIR']} ({', '.join(entry['PATTERNS'])})")
    return missing


def frame_status(ns, P: dict) -> None:
    print("=" * 112)
    print(f"{STEP_FRAME['NAME']} STATUS - {STEP_FRAME['TITLE']} - {ns.phenotype} / {P['label']} - LD panel "
          f"{PIPELINE['panel']}" + (f" - study {ns.index} ({PIPELINE['study']})" if PIPELINE["study"] else ""))
    print("=" * 112)
    missing = frame_inputs_report(P)
    summary = P["OUT"] / STEP_FRAME["SUMMARY"]
    newest_input = max((f.stat().st_mtime for e in P["INPUTS"] for f in _frame_matches(e)), default=0)
    if summary.exists():
        when = summary.stat().st_mtime
        import datetime as _fdt
        stamp = _fdt.datetime.fromtimestamp(when).strftime("%Y-%m-%d %H:%M")
        state = "STALE (an input is newer - run again)" if newest_input > when else "DONE (current)"
        print(f"\nRESULTS  : {state}, written {stamp}\n           {P['OUT']}")
        outputs = sorted(f.name for f in P["OUT"].glob(f"{STEP_FRAME['NAME']}_*") if f.is_file())
        print(f"           {len(outputs)} files: {', '.join(outputs[:6])}{' ...' if len(outputs) > 6 else ''}")
    else:
        print(f"\nRESULTS  : NOT RUN YET ({P['OUT']})")
    if missing:
        print("\nMissing required input:\n  " + "\n  ".join(missing))
        print(f"Run the previous step first: {STEP_FRAME['PREVIOUS']}")
    else:
        print(f"\nRUN      : {frame_command(ns)}\n           {frame_command(ns)} --slurm   (one SLURM job)")
    print(f"THEN     : {frame_command(ns)} --stats")


def frame_stats(ns, P: dict) -> None:
    import pandas as _fpd
    print("=" * 112)
    print(f"{STEP_FRAME['NAME']} RESULTS - {STEP_FRAME['TITLE']} - {ns.phenotype} / {P['label']} - LD panel "
          f"{PIPELINE['panel']}" + (f" - study {ns.index} ({PIPELINE['study']})" if PIPELINE["study"] else ""))
    print("=" * 112)
    summary = P["OUT"] / STEP_FRAME["SUMMARY"]
    if not summary.exists():
        raise SystemExit(f"No results yet ({summary}). Run: {frame_command(ns)}")
    info = _fjson.loads(summary.read_text(encoding="utf-8"))
    for key, value in info.items():
        if isinstance(value, (int, float, str)) and not str(key).upper().endswith(("FILE", "DIR", "PATH")):
            if isinstance(value, str) and ("/" in value or "\\" in value):
                continue
            print(f"  {key:<44} {value:,}" if isinstance(value, int) else f"  {key:<44} {value}")
        elif isinstance(value, dict) and value and all(isinstance(v, (int, float)) for v in value.values()):
            print(f"  {key:<44} " + ", ".join(f"{k}={v}" for k, v in list(value.items())[:12]))
    for table, columns, sort_col, title in STEP_FRAME["STATS_TABLES"]:
        path = P["OUT"] / table
        if not path.exists() or path.stat().st_size == 0:
            continue
        x = _fpd.read_csv(path, sep="\t", low_memory=False)
        cols = [c for c in columns if c in x.columns]
        if not cols or x.empty:
            continue
        if sort_col in x.columns:
            x = x.assign(_S=_fpd.to_numeric(x[sort_col], errors="coerce")).sort_values("_S", ascending=False)
        print(f"\n{title} ({len(x):,} rows; top 25)  [{path.name}]")
        with _fpd.option_context("display.width", 250, "display.max_columns", None, "display.max_colwidth", 32):
            print(x[cols].head(25).to_string(index=False))
    print(f"\nAll tables: {P['OUT']}")


def frame_check(ns, P: dict) -> bool:
    import importlib.util as _fimp
    print("=" * 112)
    print(f"{STEP_FRAME['NAME']} RESOURCE CHECK - {STEP_FRAME['TITLE']} - LD panel {PIPELINE['panel']}")
    print("=" * 112)
    problems = [f"input {m}" for m in frame_inputs_report(P)]
    print()
    for module, required, use in STEP_FRAME["MODULES"]:
        ok = _fimp.find_spec(module) is not None
        print(f"  {'OK     ' if ok else ('MISSING' if required else 'absent ')}  python {module:<22} {use}")
        if required and not ok:
            problems.append(f"python module {module} (Step00_Installation.py or pip install {module})")
    for label, ok, detail, required in frame_extra_checks(ns, P):
        print(f"  {'OK     ' if ok else ('MISSING' if required else 'absent ')}  {label:<29} {detail}")
        if required and not ok:
            problems.append(f"{label}: {detail}")
    print()
    if problems:
        print(f"NOT READY - {len(problems)} problem(s):")
        for problem in problems:
            print(f"  - {problem}")
        return False
    print(f"READY.\nRUN: {frame_command(ns)}\n     {frame_command(ns)} --slurm")
    return True


def frame_slurm(ns, P: dict, rest: list[str]) -> None:
    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": "4", "mem": "32G", "time": "12:00:00"}
    template = P["root"] / ns.template
    if template.is_file():
        for line in template.read_text(encoding="utf-8").splitlines():
            match = _fre.match(r"^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
            if match and (match.group(2) or "").strip() and match.group(1) not in {"job-name", "output", "error", "array"}:
                values[match.group(1)] = match.group(2).strip()
    for key, value in (("partition", ns.partition), ("qos", ns.qos), ("account", ns.account),
                       ("cpus-per-task", ns.cpus), ("mem", ns.mem), ("time", ns.time)):
        if value not in (None, ""):
            values[key] = str(value)
    P["LOG"].mkdir(parents=True, exist_ok=True)
    tag = results_slug(P["a"]).replace("/", "_")
    name = f"{STEP_FRAME['SCRIPT_PREFIX']}_{P['p']}_{tag}"
    values.update({"job-name": name[:100], "output": f"{P['LOG']}/{name}.%j.out", "error": f"{P['LOG']}/{name}.%j.err"})
    order = ["job-name", "account", "partition", "qos", "nodes", "ntasks", "cpus-per-task", "mem", "time",
             "output", "error"]
    header = "\n".join(f"#SBATCH --{k}={values[k]}" for k in order if values.get(k))
    header += "".join(f"\n#SBATCH --{k}={v}" for k, v in values.items() if k not in order and v)
    python = P["root"] / "envs" / "pipeline" / "bin" / "python"
    args = rest + (["--reference-panel", ns.reference_panel] if ns.reference_panel != DEFAULT_PANEL else []) + \
        (["--index", str(ns.index)] if ns.index is not None else [])
    command = [str(python), str(_FPath(_fsys.argv[0]).resolve())] + args
    script = P["root"] / f"{name}.sh"
    script.write_text(f"""#!/bin/bash --login
{header}

set -euo pipefail
cd {_fshlex.quote(str(P['root']))}
export PATH={_fshlex.quote(str(P['root'] / 'envs' / 'pipeline' / 'bin'))}:"$PATH"

{' '.join(_fshlex.quote(x) for x in command)}
""", encoding="utf-8", newline="\n")
    script.chmod(0o755)
    print(f"Job script: {script.name}  (one job; {values.get('cpus-per-task')} CPUs, {values.get('mem')}, "
          f"{values.get('time')})\nRUN:\n  sbatch {script.name}\nTHEN:\n  {frame_command(ns)} --status\n"
          f"  {frame_command(ns)} --stats")


def pipeline_main(run) -> None:
    ns, rest = frame_arguments()
    if any(h in _fsys.argv[1:] for h in ("-h", "--help")):
        print(FRAME_HELP)
        _fsys.argv = [_fsys.argv[0]] + rest
        run()
        return
    if not ns.phenotype or not ns.ancestry:
        _fsys.argv = [_fsys.argv[0]] + rest
        run()          # the step's own parser reports the missing --phenotype / --ancestry
        return
    if not _fre.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", ns.reference_panel):
        raise SystemExit(f"--reference-panel {ns.reference_panel!r}: use letters, digits, '_', '-' or '.'")
    PIPELINE["panel"] = ns.reference_panel
    P = frame_paths(ns)
    if ns.status:
        frame_status(ns, P)
    elif ns.stats:
        frame_stats(ns, P)
    elif ns.check_resources:
        raise SystemExit(0 if frame_check(ns, P) else 1)
    elif ns.slurm:
        frame_slurm(ns, P, rest)
    else:
        _fsys.argv = [_fsys.argv[0]] + rest
        run()

# ============================ end of pipeline frame ===========================


VERSION = "2.0.0"
DEFAULT_OT_RELEASE = "26.09"
DEFAULT_STRONG_H4 = 0.80
DEFAULT_SUGGESTIVE_H4 = 0.50
DEFAULT_L2G_TOP_K = 5
DEFAULT_LOCUS_PADDING = 250_000
DEFAULT_REPAIR_RETRIES = 2

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

DATASET_ALIASES = {
    "credible_set": [
        "credible_set", "credible_sets", "credibleSet", "credibleSets",
    ],
    "l2g_prediction": [
        "l2g_prediction", "l2g_predictions", "l2gPrediction", "l2gPredictions",
    ],
    "colocalisation": [
        "colocalisation", "colocalization", "colocalisation_coloc",
        "colocalisationColoc", "coloc",
    ],
    "study": ["study", "studies"],
}


# =============================================================================
# GENERIC HELPERS
# =============================================================================

def banner(text: str) -> None:
    print()
    print("=" * 130)
    print(text)
    print("=" * 130)


def norm_text(value: Any) -> str:
    text = re.sub(r"[^a-z0-9]+", " ", str(value).lower())
    return re.sub(r"\s+", " ", text).strip()


def slugify(value: Any) -> str:
    return norm_text(value).replace(" ", "_")


def canonical_ancestry(value: str) -> tuple[str, str]:
    key = norm_text(value)
    if key not in ANCESTRY_ALIASES:
        raise SystemExit(
            f"Unsupported ancestry {value!r}. Use EUR, AFR, EAS, SAS, or AMR."
        )
    return ANCESTRY_ALIASES[key]


def clean_str(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    x = str(value).strip()
    return "" if x.lower() in {"", "nan", "none", "null", "na", "."} else x


def canonical_chr(value: Any) -> str:
    x = clean_str(value)
    x = re.sub(r"^chr", "", x, flags=re.I)
    if x.endswith(".0"):
        x = x[:-2]
    return x.upper()


def canonical_variant(value: Any) -> str:
    x = clean_str(value)
    if not x:
        return ""
    x = re.sub(r"^chr", "", x, flags=re.I)
    parts = re.split(r"[:_]", x)
    if len(parts) < 4:
        return ""
    chrom = canonical_chr(parts[0])
    try:
        pos = str(int(float(parts[1])))
    except Exception:
        return ""
    ref = clean_str(parts[2]).upper()
    alt = clean_str(parts[3]).upper()
    if not chrom or not ref or not alt:
        return ""
    return f"{chrom}:{pos}:{ref}:{alt}"


def variant_from_parts(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    return canonical_variant(f"{chrom}:{pos}:{ref}:{alt}")


def variant_pos(value: Any) -> float:
    x = canonical_variant(value)
    if not x:
        return np.nan
    try:
        return float(x.split(":")[1])
    except Exception:
        return np.nan


def gene_base(value: Any) -> str:
    x = clean_str(value)
    if not x:
        return ""
    m = re.search(r"(ENSG\d+)", x, flags=re.I)
    if m:
        return m.group(1).upper()
    return x.split(".")[0]


def boolish(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return clean_str(value).lower() in {"1", "true", "t", "yes", "y"}


def first_existing(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    cols = list(columns)
    lower = {str(c).lower(): c for c in cols}
    for candidate in candidates:
        if candidate in cols:
            return candidate
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def dictionary_get(obj: Any, candidates: Iterable[str], default: Any = None) -> Any:
    if not isinstance(obj, dict):
        return default
    lower = {str(k).lower(): k for k in obj}
    for c in candidates:
        if c in obj:
            return obj[c]
        k = lower.get(c.lower())
        if k is not None:
            return obj[k]
    return default


def as_list(value: Any) -> list:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return []


def safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def finite_float(value: Any) -> float:
    try:
        x = float(value)
    except Exception:
        return np.nan
    return x if math.isfinite(x) else np.nan


def semicolon_unique(values: Iterable[Any], limit: int | None = None) -> str:
    out: list[str] = []
    for value in values:
        x = clean_str(value)
        if not x:
            continue
        for token in re.split(r"[;,]", x):
            token = token.strip()
            if token and token not in out:
                out.append(token)
                if limit is not None and len(out) >= limit:
                    return ";".join(out)
    return ";".join(out)


def set_from_semicolon(value: Any) -> set[str]:
    x = clean_str(value)
    if not x:
        return set()
    return {z.strip() for z in re.split(r"[;,]", x) if z.strip()}


def jaccard(a: Iterable[Any], b: Iterable[Any]) -> float:
    aa = {clean_str(x) for x in a if clean_str(x)}
    bb = {clean_str(x) for x in b if clean_str(x)}
    if not aa and not bb:
        return np.nan
    u = aa | bb
    return len(aa & bb) / len(u) if u else np.nan


def interval_overlap(start1: Any, end1: Any, start2: Any, end2: Any) -> int:
    vals = [finite_float(x) for x in [start1, end1, start2, end2]]
    if not all(math.isfinite(x) for x in vals):
        return 0
    a, b, c, d = [int(x) for x in vals]
    return max(0, min(b, d) - max(a, c) + 1)


def interval_jaccard(start1: Any, end1: Any, start2: Any, end2: Any) -> float:
    ov = interval_overlap(start1, end1, start2, end2)
    vals = [finite_float(x) for x in [start1, end1, start2, end2]]
    if not all(math.isfinite(x) for x in vals):
        return np.nan
    a, b, c, d = [int(x) for x in vals]
    union = max(b, d) - min(a, c) + 1
    return ov / union if union > 0 else np.nan


def interval_gap(start1: Any, end1: Any, start2: Any, end2: Any) -> float:
    vals = [finite_float(x) for x in [start1, end1, start2, end2]]
    if not all(math.isfinite(x) for x in vals):
        return np.nan
    a, b, c, d = [int(x) for x in vals]
    if max(a, c) <= min(b, d):
        return 0.0
    if b < c:
        return float(c - b)
    return float(a - d)


def safe_corr(x: pd.Series, y: pd.Series, method: str) -> float:
    a = pd.to_numeric(x, errors="coerce")
    b = pd.to_numeric(y, errors="coerce")
    good = a.notna() & b.notna()
    a = a[good].to_numpy()
    b = b[good].to_numpy()
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    try:
        if method == "spearman" and spearmanr is not None:
            return float(spearmanr(a, b).statistic)
        if method == "pearson" and pearsonr is not None:
            return float(pearsonr(a, b).statistic)
    except Exception:
        pass
    return np.nan


def json_safe(obj: Any) -> Any:
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return None if np.isnan(obj) else float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    raise TypeError(f"Not JSON serializable: {type(obj).__name__}")


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", low_memory=False)


# =============================================================================
# PATHS / INPUTS
# =============================================================================

def pipeline_paths(root: Path, phenotype: str, ancestry_label: str, ot_cache: str | None, release: str) -> dict[str, Path]:
    p = slugify(phenotype)
    a = slugify(ancestry_label)
    cache = Path(ot_cache).resolve() if ot_cache else root / "resources" / "opentargets" / release
    return {
        "S06": root / "06_finemapping" / p / panel_slug(a),
        "S12": root / "12_master_table" / p / results_slug(a),
        "OUT": root / "13_opentargets_comparison" / p / results_slug(a),
        "OT_CACHE": cache,
    }


def arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="GWAS2m Step13: compare Step12 master evidence with Open Targets")
    p.add_argument("--phenotype", required=True)
    p.add_argument("--ancestry", required=True)
    p.add_argument("--root", default=".")
    p.add_argument("--study", default="", help="Optional GCST study filter")
    p.add_argument("--ot-cache", default="", help="Local Open Targets bulk-data cache")
    p.add_argument("--ot-release", default=DEFAULT_OT_RELEASE)
    p.add_argument("--strong-h4", type=float, default=DEFAULT_STRONG_H4)
    p.add_argument("--suggestive-h4", type=float, default=DEFAULT_SUGGESTIVE_H4)
    p.add_argument("--ot-coloc-h4", type=float, default=0.80)
    p.add_argument("--l2g-top-k", type=int, default=DEFAULT_L2G_TOP_K)
    p.add_argument("--locus-padding", type=int, default=DEFAULT_LOCUS_PADDING)
    p.add_argument("--allow-no-coloc", action="store_true", help="Continue when OT colocalisation dataset is absent")
    p.add_argument(
        "--download-ot",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Download missing Open Targets datasets with rsync during the run (default: off - the data "
            "come from Step00_Download_Resources.py --groups opentargets; needs internet)."
        ),
    )
    p.add_argument(
        "--compact-screen",
        action="store_true",
        help="Print summary counts only instead of large detailed terminal tables.",
    )
    p.add_argument(
        "--repair-ot",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Re-download unreadable/truncated Open Targets parquet shards with rsync (default: off; "
            "unreadable shards are reported and skipped; needs internet)."
        ),
    )
    p.add_argument(
        "--repair-retries",
        type=int,
        default=DEFAULT_REPAIR_RETRIES,
        help="Number of rsync retries per corrupt Open Targets shard (default: 2).",
    )
    p.add_argument(
        "--rsync-base",
        default="",
        help=(
            "Optional Open Targets rsync output base. Default is derived from "
            "--ot-release, e.g. rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.09/output"
        ),
    )
    return p.parse_args()


def load_step12(
    paths: dict[str, Path],
    study: str = "",
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    """Load the integrated Step12 v2 evidence model.

    The first two tables are mandatory.  The remaining v2 tables are optional
    so Step13 can still audit older/reduced Step12 runs, but mechanism-level
    outputs are richest when all files are present.
    """
    base = paths["S12"]
    files = {
        "genes": base / "Step12_Locus_Gene_Master.tsv",
        "loci": base / "Step12_Locus_Master.tsv",
        "variants": base / "Step12_Variant_Master.tsv",
        "signals": base / "Step12_Signal_Master.tsv",
        "molecular": base / "Step12_Molecular_Evidence_Long.tsv",
        "strong": base / "Step12_Strong_Shared_Signals.tsv",
        "pathway": base / "Step12_Pathway_Input_Genes.tsv",
    }

    missing = [str(files[k]) for k in ["genes", "loci"] if not files[k].exists()]
    if missing:
        raise SystemExit(
            "Missing required Step12 v2 output(s):\n  " + "\n  ".join(missing) +
            "\nRun Step12_Master_Table.py first."
        )

    data = {k: read_tsv(v) for k, v in files.items()}

    if study:
        for key in ["genes", "loci", "variants", "signals", "molecular", "strong"]:
            df = data[key]
            if not df.empty and "STUDY_ACCESSION" in df.columns:
                data[key] = df[df["STUDY_ACCESSION"].astype(str).eq(study)].copy()

    if data["genes"].empty:
        raise SystemExit("Step12 locus-gene master is empty after filtering.")

    banner("STEP13 - STEP12 v2 INPUT INVENTORY")
    for key in ["genes", "loci", "variants", "signals", "molecular", "strong", "pathway"]:
        path = files[key]
        print(f"{key:12s}: {len(data[key]):>10,} rows  |  {path if path.exists() else 'NOT FOUND'}")

    return (
        data["genes"], data["loci"], data["variants"], data["signals"],
        data["molecular"], data["strong"], data["pathway"],
    )


# =============================================================================
# GWAS2m LOCUS / VARIANT REPRESENTATION
# =============================================================================

def identify_variant_columns(df: pd.DataFrame) -> dict[str, str | None]:
    return {
        "chr": first_existing(df.columns, ["CHR", "REFERENCE_CHR", "CHROM", "chromosome"]),
        "pos": first_existing(df.columns, ["REFERENCE_POS", "POS", "position"]),
        "ref": first_existing(df.columns, ["REFERENCE_REF", "REF", "referenceAllele"]),
        "alt": first_existing(df.columns, ["REFERENCE_ALT", "ALT", "alternateAllele"]),
        "locus": first_existing(df.columns, ["LOCUS_ID"]),
        "pip": first_existing(df.columns, ["PIP", "GWAS_PIP"]),
        "cs95": first_existing(df.columns, [
            "IN_95_CREDIBLE_SET", "SELECTED_CREDIBLE_SET", "CREDIBLE_SET", "IS_95_CREDIBLE_SET",
        ]),
    }


def normalize_gwas_variants(df: pd.DataFrame, study: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    cols = identify_variant_columns(df)
    if not all(cols[k] for k in ["chr", "pos", "ref", "alt", "locus"]):
        return pd.DataFrame()

    x = df.copy()
    x["STUDY_ACCESSION"] = study
    x["GWAS2M_LOCUS_ID"] = x[cols["locus"]].astype(str)
    x["VARIANT_KEY"] = [
        variant_from_parts(c, p, r, a)
        for c, p, r, a in zip(x[cols["chr"]], x[cols["pos"]], x[cols["ref"]], x[cols["alt"]])
    ]
    x["GWAS2M_PIP"] = safe_numeric(x[cols["pip"]]) if cols["pip"] else np.nan
    if cols["cs95"]:
        raw = x[cols["cs95"]]
        # Numeric credible-set IDs should count as membership when non-empty.
        x["GWAS2M_IS_95_CS"] = raw.map(lambda v: boolish(v) or clean_str(v) not in {"", "0", "False", "false", "NO", "no"})
    else:
        x["GWAS2M_IS_95_CS"] = False
    x["GWAS2M_CHR"] = x[cols["chr"]].map(canonical_chr)
    x["GWAS2M_POS"] = safe_numeric(x[cols["pos"]])
    return x[[
        "STUDY_ACCESSION", "GWAS2M_LOCUS_ID", "GWAS2M_CHR", "GWAS2M_POS",
        "VARIANT_KEY", "GWAS2M_PIP", "GWAS2M_IS_95_CS",
    ]].drop_duplicates(["STUDY_ACCESSION", "GWAS2M_LOCUS_ID", "VARIANT_KEY"])


def load_gwas_variants(paths: dict[str, Path], step12_variants: pd.DataFrame, studies: Iterable[str]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for study in sorted({clean_str(x) for x in studies if clean_str(x)}):
        fine = paths["S06"] / study / f"{study}_finemapped_variants.tsv.gz"
        if fine.exists() and fine.stat().st_size > 0:
            try:
                z = pd.read_csv(fine, sep="\t", compression="infer", low_memory=False)
                n = normalize_gwas_variants(z, study)
                if not n.empty:
                    frames.append(n)
                    continue
            except Exception as exc:
                print(f"[WARN] Could not read Step06 variants for {study}: {exc}")

        # Fallback to Step12 variant master.
        if not step12_variants.empty:
            z = step12_variants.copy()
            if "STUDY_ACCESSION" in z.columns:
                z = z[z["STUDY_ACCESSION"].astype(str).eq(study)].copy()
            n = normalize_gwas_variants(z, study)
            if not n.empty:
                frames.append(n)

    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def build_gwas_loci(locus_master: pd.DataFrame, gwas_variants: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, base in locus_master.iterrows():
        study = clean_str(base.get("STUDY_ACCESSION"))
        locus = clean_str(base.get("LOCUS_ID"))
        g = (
            gwas_variants[
                gwas_variants["STUDY_ACCESSION"].astype(str).eq(study)
                & gwas_variants["GWAS2M_LOCUS_ID"].astype(str).eq(locus)
            ].copy()
            if not gwas_variants.empty else pd.DataFrame()
        )

        row = base.to_dict()
        row["STUDY_ID"] = study
        row["GWAS2M_LOCUS_ID"] = locus
        if not g.empty:
            pip = safe_numeric(g["GWAS2M_PIP"])
            lead = g.loc[pip.idxmax()] if pip.notna().any() else g.iloc[0]
            row["GWAS2M_CHR"] = clean_str(lead.get("GWAS2M_CHR"))
            row["GWAS2M_START"] = int(g["GWAS2M_POS"].min()) if g["GWAS2M_POS"].notna().any() else np.nan
            row["GWAS2M_END"] = int(g["GWAS2M_POS"].max()) if g["GWAS2M_POS"].notna().any() else np.nan
            row["GWAS2M_LEAD_VARIANT"] = clean_str(lead.get("VARIANT_KEY"))
            row["GWAS2M_LEAD_PIP"] = finite_float(lead.get("GWAS2M_PIP"))
            row["GWAS2M_N_VARIANTS"] = len(g)
            row["GWAS2M_N_CS95"] = int(g["GWAS2M_IS_95_CS"].map(boolish).sum())
        else:
            top = canonical_variant(base.get("TOP_PIP_VARIANT", ""))
            row["GWAS2M_CHR"] = top.split(":")[0] if top else ""
            pos = variant_pos(top)
            row["GWAS2M_START"] = pos
            row["GWAS2M_END"] = pos
            row["GWAS2M_LEAD_VARIANT"] = top
            row["GWAS2M_LEAD_PIP"] = finite_float(base.get("ANNOT_MAX_PIP"))
            row["GWAS2M_N_VARIANTS"] = 0
            row["GWAS2M_N_CS95"] = 0
        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# OPEN TARGETS PARQUET ACCESS
# =============================================================================

def require_duckdb() -> None:
    if duckdb is None:
        raise RuntimeError(
            "duckdb is required. Install in the active environment with:\n"
            "  pip install duckdb"
        )


def dataset_has_parquet(path: Path) -> bool:
    return path.exists() and any(path.rglob("*.parquet"))


# Cache Parquet health so each shard is validated only once per Step13 run.
_PARQUET_HEALTH_CACHE: dict[str, tuple[list[str], list[dict[str, str]]]] = {}
_PARQUET_HEALTH_ROWS: list[dict[str, str]] = []
_PARQUET_EXPECTED_COUNTS: dict[str, int] = {}
_OT_REPAIR_ROWS: list[dict[str, Any]] = []


def healthy_parquet_files(path: Path, *, allow_all_bad: bool = False) -> list[str]:
    """Return readable Parquet shards and audit corrupt/truncated files.

    Open Targets bulk downloads are sharded. A single interrupted shard should
    not destroy locus/L2G comparisons from otherwise healthy datasets. We test
    each shard independently with DuckDB, exclude unreadable shards, and record
    them for an explicit QC report. Downstream code MUST treat an incomplete
    colocalisation dataset conservatively (not as evidence of absence).
    """
    require_duckdb()
    key = str(path.resolve())
    cached = _PARQUET_HEALTH_CACHE.get(key)
    if cached is not None:
        return cached[0]

    raw = sorted(path.rglob("*.parquet")) if path.exists() else []
    # Preserve the number of shards seen before any repair/quarantine operation.
    # This prevents a failed repair from being mistaken for a complete dataset
    # merely because the corrupt shard was moved out of the way.
    _PARQUET_EXPECTED_COUNTS.setdefault(key, len(raw))
    good: list[str] = []
    bad: list[dict[str, str]] = []

    con = duckdb.connect(database=":memory:")
    try:
        for fp in raw:
            try:
                q = f"SELECT * FROM read_parquet({sql_quote(str(fp))}) LIMIT 0"
                con.execute(q)
                good.append(str(fp))
                _PARQUET_HEALTH_ROWS.append({
                    "DATASET": path.name,
                    "FILE": str(fp),
                    "STATUS": "OK",
                    "ERROR": "",
                })
            except Exception as exc:
                row = {
                    "DATASET": path.name,
                    "FILE": str(fp),
                    "STATUS": "CORRUPT_OR_TRUNCATED",
                    "ERROR": f"{type(exc).__name__}: {exc}",
                }
                bad.append(row)
                _PARQUET_HEALTH_ROWS.append(row)
    finally:
        con.close()

    _PARQUET_HEALTH_CACHE[key] = (good, bad)

    if bad:
        print()
        print(f"[WARN] Open Targets dataset {path.name!r} contains {len(bad)} unreadable Parquet shard(s).")
        print(f"       Healthy shards retained: {len(good)}/{len(raw)}")
        for row in bad[:20]:
            print(f"       BAD: {row['FILE']}")
        if len(bad) > 20:
            print(f"       ... plus {len(bad)-20} additional bad shard(s)")

    if raw and not good and not allow_all_bad:
        raise RuntimeError(
            f"All Parquet shards are unreadable for Open Targets dataset: {path}"
        )

    return good


def _clear_parquet_health(path: Path) -> None:
    key = str(path.resolve())
    _PARQUET_HEALTH_CACHE.pop(key, None)
    # Keep only the latest health state for this dataset in the final QC table.
    _PARQUET_HEALTH_ROWS[:] = [
        row for row in _PARQUET_HEALTH_ROWS
        if clean_str(row.get("DATASET")) != path.name
    ]


def _validate_single_parquet(path: Path) -> tuple[bool, str]:
    require_duckdb()
    if not path.exists() or path.stat().st_size == 0:
        return False, "FILE_MISSING_OR_EMPTY"
    con = duckdb.connect(database=":memory:")
    try:
        con.execute(
            f"SELECT * FROM read_parquet({sql_quote(str(path))}) LIMIT 0"
        )
        return True, ""
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"
    finally:
        con.close()


def official_ot_rsync_base(release: str, override: str = "") -> str:
    if clean_str(override):
        return clean_str(override).rstrip("/")
    return (
        "rsync.ebi.ac.uk::pub/databases/opentargets/"
        f"platform/{release}/output"
    )


def repair_bad_parquet_shards(
    dataset_path: Path,
    release: str,
    retries: int,
    rsync_base: str = "",
) -> dict[str, Any]:
    """Retry only unreadable Open Targets shards from the official release.

    A corrupt local shard is first renamed to *.corrupt.<timestamp>, then the
    exact same shard is requested from the official rsync release. If exact-file
    transfer fails, an incremental dataset-directory rsync is attempted as a
    fallback. Every attempt is audited.
    """
    if shutil.which("rsync") is None:
        return {
            "DATASET": dataset_path.name,
            "ATTEMPTED": False,
            "SUCCESS": False,
            "MESSAGE": "rsync executable not found",
        }

    # Populate initial health/cache and preserve expected shard count.
    healthy_parquet_files(dataset_path, allow_all_bad=True)
    key = str(dataset_path.resolve())
    good, bad = _PARQUET_HEALTH_CACHE.get(key, ([], []))
    expected = _PARQUET_EXPECTED_COUNTS.get(key, len(good) + len(bad))
    if not bad:
        return {
            "DATASET": dataset_path.name,
            "ATTEMPTED": False,
            "SUCCESS": True,
            "MESSAGE": "No corrupt shards detected",
        }

    base = official_ot_rsync_base(release, rsync_base)
    dataset_remote_name = dataset_path.name
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    print()
    print(f"[REPAIR] {dataset_path.name}: retrying {len(bad)} corrupt shard(s) from Open Targets {release}")

    for badrow in bad:
        fp = Path(badrow["FILE"])
        try:
            rel = fp.relative_to(dataset_path)
        except Exception:
            rel = Path(fp.name)
        backup = fp.with_name(fp.name + f".corrupt.{stamp}")
        if fp.exists():
            fp.rename(backup)

        repaired = False
        last_error = ""
        attempts = max(1, int(retries))
        for attempt in range(1, attempts + 1):
            remote = f"{base}/{dataset_remote_name}/{rel.as_posix()}"
            fp.parent.mkdir(parents=True, exist_ok=True)
            command = [
                "rsync", "-av", "--partial", "--info=progress2",
                remote, str(fp.parent) + "/",
            ]
            print(f"[REPAIR] attempt {attempt}/{attempts}: {' '.join(command)}")
            proc = subprocess.run(command, text=True, capture_output=True, check=False)
            ok, err = _validate_single_parquet(fp)
            _OT_REPAIR_ROWS.append({
                "DATASET": dataset_path.name,
                "FILE": str(fp),
                "BACKUP": str(backup),
                "ATTEMPT": attempt,
                "RSYNC_RETURN_CODE": proc.returncode,
                "VALID_AFTER_DOWNLOAD": "YES" if ok else "NO",
                "ERROR": err or proc.stderr.strip(),
            })
            if proc.returncode == 0 and ok:
                repaired = True
                print(f"[REPAIR] OK: {fp.name}")
                break
            last_error = err or proc.stderr.strip() or "unknown rsync/Parquet validation error"
            time.sleep(min(5, attempt))

        if not repaired:
            # Fallback: ask rsync to reconcile the whole dataset directory.
            remote_dir = f"{base}/{dataset_remote_name}/"
            command = [
                "rsync", "-av", "--partial", "--info=progress2",
                remote_dir, str(dataset_path) + "/",
            ]
            print(f"[REPAIR] exact shard retry failed; dataset-level incremental retry: {' '.join(command)}")
            proc = subprocess.run(command, text=True, capture_output=True, check=False)
            ok, err = _validate_single_parquet(fp)
            _OT_REPAIR_ROWS.append({
                "DATASET": dataset_path.name,
                "FILE": str(fp),
                "BACKUP": str(backup),
                "ATTEMPT": "DATASET_FALLBACK",
                "RSYNC_RETURN_CODE": proc.returncode,
                "VALID_AFTER_DOWNLOAD": "YES" if ok else "NO",
                "ERROR": err or proc.stderr.strip() or last_error,
            })

    _clear_parquet_health(dataset_path)
    good2 = healthy_parquet_files(dataset_path, allow_all_bad=True)
    key = str(dataset_path.resolve())
    _, bad2 = _PARQUET_HEALTH_CACHE.get(key, ([], []))
    success = len(bad2) == 0 and len(good2) >= expected
    return {
        "DATASET": dataset_path.name,
        "ATTEMPTED": True,
        "SUCCESS": success,
        "EXPECTED_SHARDS": expected,
        "HEALTHY_SHARDS_AFTER": len(good2),
        "BAD_SHARDS_AFTER": len(bad2),
        "MESSAGE": "Repair complete" if success else "Repair incomplete; inspect repair/QC audit",
    }


def parquet_dataset_complete(path: Path | None) -> bool:
    if path is None:
        return False
    healthy_parquet_files(path, allow_all_bad=True)
    key = str(path.resolve())
    good, bad = _PARQUET_HEALTH_CACHE.get(key, ([], []))
    expected = _PARQUET_EXPECTED_COUNTS.get(key, len(good) + len(bad))
    return bool(good) and len(bad) == 0 and len(good) >= expected


def parquet_dataset_has_healthy_shards(path: Path | None) -> bool:
    if path is None:
        return False
    good = healthy_parquet_files(path, allow_all_bad=True)
    return bool(good)



def remote_ot_dataset_names(release: str, rsync_base: str = "") -> set[str]:
    """List dataset folders in the official Open Targets release."""
    if shutil.which("rsync") is None:
        raise RuntimeError(
            "rsync is required to discover/download missing Open Targets datasets. "
            "Load/install rsync or rerun with --no-download-ot and provide --ot-cache."
        )
    base = official_ot_rsync_base(release, rsync_base)
    result = subprocess.run(
        ["rsync", "--list-only", base.rstrip("/") + "/"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Could not list the official Open Targets rsync release.\n"
            f"Base: {base}\n"
            f"stderr:\n{result.stderr}"
        )
    names: set[str] = set()
    for line in result.stdout.splitlines():
        parts = line.split()
        if not parts:
            continue
        name = parts[-1].rstrip("/")
        if name not in {"", ".", ".."}:
            names.add(name)
    return names


def sync_ot_dataset(remote_name: str, cache: Path, release: str, rsync_base: str = "") -> Path:
    if shutil.which("rsync") is None:
        raise RuntimeError("rsync is not available; cannot download Open Targets datasets.")
    base = official_ot_rsync_base(release, rsync_base)
    destination = cache / remote_name
    destination.mkdir(parents=True, exist_ok=True)
    remote = f"{base.rstrip('/')}/{remote_name}/"
    command = [
        "rsync", "-rltz", "--partial", "--info=progress2",
        remote, str(destination) + "/",
    ]
    print("\n[Open Targets download]")
    print(" ".join(command))
    subprocess.run(command, check=True)
    return destination


def ensure_ot_datasets(
    cache: Path,
    release: str,
    download: bool,
    rsync_base: str = "",
) -> dict[str, Path | None]:
    """Resolve cached OT datasets and optionally download missing ones.

    credible_set and l2g_prediction are required. colocalisation and study are
    strongly recommended because absence-based mechanism interpretation is
    suppressed when those datasets are incomplete/unavailable.
    """
    cache.mkdir(parents=True, exist_ok=True)
    resolved = {k: resolve_ot_dataset(cache, k) for k in DATASET_ALIASES}
    missing = [k for k, v in resolved.items() if v is None]
    if not missing or not download:
        return resolved

    banner("OPEN TARGETS - DISCOVERING / DOWNLOADING MISSING DATASETS")
    print("Missing:", ", ".join(missing))
    available = remote_ot_dataset_names(release, rsync_base)

    for logical in missing:
        aliases = DATASET_ALIASES[logical]
        remote_name = next((a for a in aliases if a in available), None)
        if remote_name is None:
            if logical in {"credible_set", "l2g_prediction"}:
                raise RuntimeError(
                    f"Required Open Targets dataset {logical!r} was not found in release {release}. "
                    f"Tried aliases: {aliases}"
                )
            print(f"[WARN] Optional/recommended OT dataset {logical!r} was not found in release {release}.")
            continue
        sync_ot_dataset(remote_name, cache, release, rsync_base)
        resolved[logical] = resolve_ot_dataset(cache, logical)

    return resolved


def resolve_ot_dataset(cache: Path, logical_name: str) -> Path | None:
    aliases = DATASET_ALIASES[logical_name]
    for alias in aliases:
        p = cache / alias
        if dataset_has_parquet(p):
            return p

    # More permissive discovery for already-downloaded folders with changed names.
    if cache.exists():
        wanted = {norm_text(a).replace(" ", "") for a in aliases}
        for p in cache.iterdir():
            if not p.is_dir():
                continue
            n = norm_text(p.name).replace(" ", "")
            if n in wanted and dataset_has_parquet(p):
                return p
    return None


def parquet_files(path: Path) -> list[str]:
    # Only return shards that DuckDB can actually open.
    return healthy_parquet_files(path)


def sql_quote(value: Any) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def duckdb_columns(path: Path) -> list[str]:
    require_duckdb()
    files = parquet_files(path)
    if not files:
        return []
    con = duckdb.connect(database=":memory:")
    try:
        flist = "[" + ",".join(sql_quote(x) for x in files) + "]"
        q = f"DESCRIBE SELECT * FROM read_parquet({flist}, union_by_name=true)"
        return con.execute(q).fetchdf()["column_name"].astype(str).tolist()
    finally:
        con.close()


def subset_parquet(path: Path, column_candidates: list[str], values: Iterable[Any]) -> tuple[pd.DataFrame, str | None]:
    require_duckdb()
    files = parquet_files(path)
    if not files:
        return pd.DataFrame(), None
    columns = duckdb_columns(path)
    column = first_existing(columns, column_candidates)
    if column is None:
        return pd.DataFrame(), None

    wanted = sorted({clean_str(x) for x in values if clean_str(x)})
    if not wanted:
        return pd.DataFrame(), column

    con = duckdb.connect(database=":memory:")
    try:
        flist = "[" + ",".join(sql_quote(x) for x in files) + "]"
        vals = "(" + ",".join(sql_quote(x) for x in wanted) + ")"
        q = f'''SELECT * FROM read_parquet({flist}, union_by_name=true)
                WHERE CAST("{column}" AS VARCHAR) IN {vals}'''
        return con.execute(q).fetchdf(), column
    finally:
        con.close()


# =============================================================================
# OPEN TARGETS CREDIBLE SETS
# =============================================================================

def extract_variant_from_item(item: Any) -> str:
    if not isinstance(item, dict):
        return ""
    direct = dictionary_get(item, ["variantId", "variant_id", "variant"])
    if isinstance(direct, dict):
        direct = dictionary_get(direct, ["id", "variantId", "variant_id"])
    if direct:
        z = canonical_variant(direct)
        if z:
            return z
    chrom = dictionary_get(item, ["chromosome", "chr"])
    pos = dictionary_get(item, ["position", "pos"])
    ref = dictionary_get(item, ["referenceAllele", "reference_allele", "ref"])
    alt = dictionary_get(item, ["alternateAllele", "alternate_allele", "alt"])
    if all(v is not None for v in [chrom, pos, ref, alt]):
        return variant_from_parts(chrom, pos, ref, alt)
    return ""


def flatten_ot_credible_sets(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()

    cols = df.columns
    sid_col = first_existing(cols, ["studyLocusId", "study_locus_id"])
    study_col = first_existing(cols, ["studyId", "study_id"])
    method_col = first_existing(cols, ["finemappingMethod", "finemapping_method"])
    chr_col = first_existing(cols, ["chromosome", "chr"])
    pos_col = first_existing(cols, ["position", "pos"])
    qtl_gene_col = first_existing(cols, ["qtlGeneId", "qtl_gene_id", "geneId", "gene_id"])
    lead_col = first_existing(cols, ["variantId", "variant_id", "leadVariantId", "lead_variant_id"])
    locus_col = first_existing(cols, ["locus", "variants", "credibleSetVariants", "credible_set_variants"])

    if sid_col is None:
        raise RuntimeError("Open Targets credible_set dataset has no studyLocusId column")

    meta_rows: list[dict[str, Any]] = []
    var_rows: list[dict[str, Any]] = []

    for _, rec in df.iterrows():
        otid = clean_str(rec.get(sid_col))
        study = clean_str(rec.get(study_col)) if study_col else ""
        method = clean_str(rec.get(method_col)) if method_col else ""
        qtl_gene = gene_base(rec.get(qtl_gene_col)) if qtl_gene_col else ""
        members = as_list(rec.get(locus_col)) if locus_col else []

        local: list[dict[str, Any]] = []
        for item in members:
            if not isinstance(item, dict):
                continue
            variant = extract_variant_from_item(item)
            if not variant:
                continue
            pip = dictionary_get(item, ["posteriorProbability", "posterior_probability", "pip", "PIP"])
            is95 = dictionary_get(item, ["is95CredibleSet", "is95_credible_set", "in95CredibleSet"], default=True)
            local.append({
                "OT_STUDY_LOCUS_ID": otid,
                "OT_STUDY_ID": study,
                "VARIANT_KEY": variant,
                "OT_PIP": finite_float(pip),
                "OT_IS_95_CS": boolish(is95),
            })

        var_rows.extend(local)
        lead = canonical_variant(rec.get(lead_col)) if lead_col else ""
        if not lead and local:
            locdf = pd.DataFrame(local)
            pp = safe_numeric(locdf["OT_PIP"])
            lead = clean_str(locdf.loc[pp.idxmax(), "VARIANT_KEY"]) if pp.notna().any() else clean_str(locdf.iloc[0]["VARIANT_KEY"])

        chrom = canonical_chr(rec.get(chr_col)) if chr_col else ""
        pos = finite_float(rec.get(pos_col)) if pos_col else np.nan
        if lead:
            chrom = lead.split(":")[0]
            if not math.isfinite(pos):
                pos = variant_pos(lead)

        positions = [variant_pos(x["VARIANT_KEY"]) for x in local]
        positions = [x for x in positions if math.isfinite(x)]
        start = min(positions) if positions else pos
        end = max(positions) if positions else pos

        meta_rows.append({
            "OT_STUDY_LOCUS_ID": otid,
            "OT_STUDY_ID": study,
            "OT_METHOD": method,
            "OT_QTL_GENE_ID": qtl_gene,
            "OT_CHR": chrom,
            "OT_START": start,
            "OT_END": end,
            "OT_LEAD_VARIANT": lead,
            "OT_LEAD_POSITION": variant_pos(lead) if lead else pos,
            "OT_N_VARIANTS": len(local),
            "OT_N_CS95": sum(boolish(x["OT_IS_95_CS"]) for x in local),
        })

    return pd.DataFrame(meta_rows).drop_duplicates("OT_STUDY_LOCUS_ID"), pd.DataFrame(var_rows)


# =============================================================================
# LOCUS MATCHING
# =============================================================================

def candidate_locus_pairs(gwas_loci: pd.DataFrame, ot_loci: pd.DataFrame, gwas_variants: pd.DataFrame, ot_variants: pd.DataFrame, padding: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    gvar_map: dict[tuple[str, str], set[str]] = {}
    if not gwas_variants.empty:
        for (study, locus), g in gwas_variants.groupby(["STUDY_ACCESSION", "GWAS2M_LOCUS_ID"], sort=False):
            gvar_map[(str(study), str(locus))] = set(g["VARIANT_KEY"].dropna().astype(str))

    ovar_map: dict[str, set[str]] = {}
    if not ot_variants.empty:
        for otid, g in ot_variants.groupby("OT_STUDY_LOCUS_ID", sort=False):
            ovar_map[str(otid)] = set(g["VARIANT_KEY"].dropna().astype(str))

    for _, g in gwas_loci.iterrows():
        study = clean_str(g.get("STUDY_ID"))
        locus = clean_str(g.get("GWAS2M_LOCUS_ID"))
        chrom = canonical_chr(g.get("GWAS2M_CHR"))
        candidates = ot_loci[ot_loci["OT_STUDY_ID"].astype(str).eq(study)].copy()
        if chrom:
            candidates = candidates[candidates["OT_CHR"].astype(str).map(canonical_chr).eq(chrom)].copy()

        gv = gvar_map.get((study, locus), set())
        for _, o in candidates.iterrows():
            otid = clean_str(o.get("OT_STUDY_LOCUS_ID"))
            ov = ovar_map.get(otid, set())
            shared = gv & ov
            lead_exact = bool(clean_str(g.get("GWAS2M_LEAD_VARIANT")) and clean_str(g.get("GWAS2M_LEAD_VARIANT")) == clean_str(o.get("OT_LEAD_VARIANT")))
            ov_bp = interval_overlap(g.get("GWAS2M_START"), g.get("GWAS2M_END"), o.get("OT_START"), o.get("OT_END"))
            ij = interval_jaccard(g.get("GWAS2M_START"), g.get("GWAS2M_END"), o.get("OT_START"), o.get("OT_END"))
            gap = interval_gap(g.get("GWAS2M_START"), g.get("GWAS2M_END"), o.get("OT_START"), o.get("OT_END"))
            lead_pos = finite_float(o.get("OT_LEAD_POSITION"))
            gs = finite_float(g.get("GWAS2M_START")); ge = finite_float(g.get("GWAS2M_END"))
            lead_inside = math.isfinite(lead_pos) and math.isfinite(gs) and math.isfinite(ge) and gs <= lead_pos <= ge

            method = clean_str(o.get("OT_METHOD"))
            method_lower = method.lower()
            is_pics = "pics" in method_lower
            is_susie = "susie" in method_lower

            # Method-aware matching. PICS is usually represented primarily by a
            # lead and variant probabilities, so a one-base PICS locus must not be
            # interpreted as a biological interval comparable to a SuSiE CS.
            if is_pics:
                eligible = (
                    lead_exact or bool(shared) or lead_inside
                    or (math.isfinite(gap) and gap <= padding)
                )
            else:
                eligible = (
                    lead_exact or bool(shared) or ov_bp > 0 or lead_inside
                    or (math.isfinite(gap) and gap <= padding)
                )
            if not eligible:
                continue

            if lead_exact:
                basis = "EXACT_LEAD"
            elif shared:
                basis = "SHARED_VARIANT"
            elif is_pics and lead_inside:
                basis = "PICS_LEAD_INSIDE_GWAS2M_LOCUS"
            elif is_pics:
                basis = "PICS_LEAD_WITHIN_PADDING"
            elif is_susie and ov_bp > 0:
                basis = "SUSIE_INTERVAL_OVERLAP"
            elif ov_bp > 0:
                basis = "INTERVAL_OVERLAP"
            elif lead_inside:
                basis = "OT_LEAD_INSIDE_GWAS2M_LOCUS"
            else:
                basis = "WITHIN_PADDING"

            # Deterministic score for one-to-one primary matching. Interval
            # overlap contributes strongly for non-PICS methods but not for PICS.
            score = (
                (1_000_000 if lead_exact else 0)
                + len(shared) * 10_000
                + (0 if is_pics else int(ov_bp > 0) * 5_000)
                + int(lead_inside) * 2_000
                + (0 if not math.isfinite(gap) else max(0, padding - gap) / max(1, padding))
                + (0 if is_pics or not math.isfinite(ij) else ij)
            )

            rows.append({
                "STUDY_ID": study,
                "GWAS2M_LOCUS_ID": locus,
                "OT_STUDY_LOCUS_ID": otid,
                "GWAS2M_CHR": chrom,
                "GWAS2M_START": g.get("GWAS2M_START"),
                "GWAS2M_END": g.get("GWAS2M_END"),
                "GWAS2M_LEAD_VARIANT": clean_str(g.get("GWAS2M_LEAD_VARIANT")),
                "OT_METHOD": method,
                "OT_START": o.get("OT_START"),
                "OT_END": o.get("OT_END"),
                "OT_LEAD_VARIANT": clean_str(o.get("OT_LEAD_VARIANT")),
                "EXACT_LEAD": lead_exact,
                "N_SHARED_VARIANTS": len(shared),
                "VARIANT_JACCARD": jaccard(gv, ov),
                "INTERVAL_OVERLAP_BP": ov_bp,
                "INTERVAL_JACCARD": ij,
                "METHOD_APPROPRIATE_INTERVAL_JACCARD": np.nan if is_pics else ij,
                "LOCUS_GAP_BP": gap,
                "OT_LEAD_INSIDE_GWAS2M_LOCUS": lead_inside,
                "MATCH_BASIS": basis,
                "MATCH_SCORE": score,
            })

    return pd.DataFrame(rows)


def select_primary_matches(pair_candidates: pd.DataFrame) -> pd.DataFrame:
    if pair_candidates.empty:
        return pd.DataFrame()
    rows: list[pd.Series] = []
    for study, x in pair_candidates.groupby("STUDY_ID", sort=True):
        x = x.sort_values(
            ["MATCH_SCORE", "N_SHARED_VARIANTS", "INTERVAL_OVERLAP_BP"],
            ascending=[False, False, False],
            kind="stable",
        )
        used_g: set[str] = set()
        used_o: set[str] = set()
        for _, r in x.iterrows():
            g = clean_str(r["GWAS2M_LOCUS_ID"])
            o = clean_str(r["OT_STUDY_LOCUS_ID"])
            if g in used_g or o in used_o:
                continue
            used_g.add(g); used_o.add(o)
            rows.append(r)
    return pd.DataFrame(rows).reset_index(drop=True) if rows else pd.DataFrame()


# =============================================================================
# VARIANT / PIP COMPARISON
# =============================================================================

def compare_variants(matches: pd.DataFrame, gwas_variants: pd.DataFrame, ot_variants: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, m in matches.iterrows():
        study = clean_str(m.get("STUDY_ID")); gl = clean_str(m.get("GWAS2M_LOCUS_ID")); ol = clean_str(m.get("OT_STUDY_LOCUS_ID"))
        g = gwas_variants[
            gwas_variants["STUDY_ACCESSION"].astype(str).eq(study)
            & gwas_variants["GWAS2M_LOCUS_ID"].astype(str).eq(gl)
        ].copy() if not gwas_variants.empty else pd.DataFrame()
        o = ot_variants[ot_variants["OT_STUDY_LOCUS_ID"].astype(str).eq(ol)].copy() if not ot_variants.empty else pd.DataFrame()

        base = {k: m.get(k) for k in ["STUDY_ID", "GWAS2M_LOCUS_ID", "OT_STUDY_LOCUS_ID", "OT_METHOD", "MATCH_BASIS"]}
        if g.empty or o.empty:
            rows.append({**base, "DETAIL_AVAILABLE": False, "N_GWAS2M_VARIANTS": len(g), "N_OT_VARIANTS": len(o)})
            continue

        gp = safe_numeric(g["GWAS2M_PIP"]); op = safe_numeric(o["OT_PIP"])
        gtop = clean_str(g.loc[gp.idxmax(), "VARIANT_KEY"]) if gp.notna().any() else clean_str(g.iloc[0]["VARIANT_KEY"])
        otop = clean_str(o.loc[op.idxmax(), "VARIANT_KEY"]) if op.notna().any() else clean_str(o.iloc[0]["VARIANT_KEY"])
        gall = set(g["VARIANT_KEY"].dropna().astype(str)); oall = set(o["VARIANT_KEY"].dropna().astype(str))
        g95 = set(g.loc[g["GWAS2M_IS_95_CS"].map(boolish), "VARIANT_KEY"].dropna().astype(str))
        o95 = set(o.loc[o["OT_IS_95_CS"].map(boolish), "VARIANT_KEY"].dropna().astype(str))

        merged = g[["VARIANT_KEY", "GWAS2M_PIP"]].merge(o[["VARIANT_KEY", "OT_PIP"]], on="VARIANT_KEY", how="inner")
        mae = float((safe_numeric(merged["GWAS2M_PIP"]) - safe_numeric(merged["OT_PIP"])).abs().mean()) if not merged.empty else np.nan
        rows.append({
            **base,
            "DETAIL_AVAILABLE": True,
            "N_GWAS2M_VARIANTS": len(gall),
            "N_OT_VARIANTS": len(oall),
            "N_SHARED_VARIANTS": len(gall & oall),
            "VARIANT_JACCARD": jaccard(gall, oall),
            "GWAS2M_TOP_VARIANT": gtop,
            "OT_TOP_VARIANT": otop,
            "TOP_VARIANT_AGREEMENT": gtop == otop and bool(gtop),
            "N_GWAS2M_CS95": len(g95),
            "N_OT_CS95": len(o95),
            "N_SHARED_CS95": len(g95 & o95),
            "CS95_JACCARD": jaccard(g95, o95),
            "N_SHARED_FOR_PIP": len(merged),
            "PIP_SPEARMAN": safe_corr(merged["GWAS2M_PIP"], merged["OT_PIP"], "spearman") if not merged.empty else np.nan,
            "PIP_PEARSON": safe_corr(merged["GWAS2M_PIP"], merged["OT_PIP"], "pearson") if not merged.empty else np.nan,
            "PIP_MAE": mae,
        })
    return pd.DataFrame(rows)


# =============================================================================
# VARIANT-LEVEL DETAIL + PERFORMANCE / MISS DIAGNOSTICS
# =============================================================================

def build_variant_level_detail(
    matches: pd.DataFrame,
    gwas_variants: pd.DataFrame,
    ot_variants: pd.DataFrame,
) -> pd.DataFrame:
    """One row per variant in the union of each matched locus pair."""
    rows: list[dict[str, Any]] = []
    for _, m in matches.iterrows():
        study = clean_str(m.get("STUDY_ID"))
        gl = clean_str(m.get("GWAS2M_LOCUS_ID"))
        ol = clean_str(m.get("OT_STUDY_LOCUS_ID"))
        g = (
            gwas_variants[
                gwas_variants["STUDY_ACCESSION"].astype(str).eq(study)
                & gwas_variants["GWAS2M_LOCUS_ID"].astype(str).eq(gl)
            ].copy()
            if not gwas_variants.empty else pd.DataFrame()
        )
        o = (
            ot_variants[ot_variants["OT_STUDY_LOCUS_ID"].astype(str).eq(ol)].copy()
            if not ot_variants.empty else pd.DataFrame()
        )
        if g.empty and o.empty:
            continue
        gg = g[["VARIANT_KEY", "GWAS2M_PIP", "GWAS2M_IS_95_CS"]].drop_duplicates("VARIANT_KEY") if not g.empty else pd.DataFrame(columns=["VARIANT_KEY", "GWAS2M_PIP", "GWAS2M_IS_95_CS"])
        oo = o[["VARIANT_KEY", "OT_PIP", "OT_IS_95_CS"]].drop_duplicates("VARIANT_KEY") if not o.empty else pd.DataFrame(columns=["VARIANT_KEY", "OT_PIP", "OT_IS_95_CS"])
        z = gg.merge(oo, on="VARIANT_KEY", how="outer")
        for _, r in z.iterrows():
            gp = finite_float(r.get("GWAS2M_PIP"))
            op = finite_float(r.get("OT_PIP"))
            g_present = math.isfinite(gp) or clean_str(r.get("GWAS2M_IS_95_CS")) != ""
            o_present = math.isfinite(op) or clean_str(r.get("OT_IS_95_CS")) != ""
            if g_present and o_present:
                status = "BOTH"
            elif g_present:
                status = "GWAS2M_ONLY_VARIANT"
            else:
                status = "OPEN_TARGETS_ONLY_VARIANT"
            rows.append({
                "STUDY_ID": study,
                "GWAS2M_LOCUS_ID": gl,
                "OT_STUDY_LOCUS_ID": ol,
                "MATCH_BASIS": clean_str(m.get("MATCH_BASIS")),
                "VARIANT_KEY": clean_str(r.get("VARIANT_KEY")),
                "VARIANT_STATUS": status,
                "GWAS2M_PIP": gp,
                "OT_PIP": op,
                "ABS_PIP_DIFFERENCE": abs(gp - op) if math.isfinite(gp) and math.isfinite(op) else np.nan,
                "GWAS2M_IN_95_CS": boolish(r.get("GWAS2M_IS_95_CS")),
                "OT_IN_95_CS": boolish(r.get("OT_IS_95_CS")),
                "BOTH_IN_95_CS": boolish(r.get("GWAS2M_IS_95_CS")) and boolish(r.get("OT_IS_95_CS")),
            })
    return pd.DataFrame(rows)


def build_missed_loci_diagnostic(
    gcomp: pd.DataFrame,
    ocomp: pd.DataFrame,
    pair_candidates: pd.DataFrame,
    gwas_loci: pd.DataFrame,
    ot_loci: pd.DataFrame,
    padding: int,
) -> pd.DataFrame:
    """Explain algorithmically why each unmatched locus was not a primary match.

    These are matching/coverage explanations, not claims about biological truth.
    """
    rows: list[dict[str, Any]] = []

    def candidate_rows_for_g(study: str, locus: str) -> pd.DataFrame:
        if pair_candidates.empty:
            return pd.DataFrame()
        return pair_candidates[
            pair_candidates["STUDY_ID"].astype(str).eq(study)
            & pair_candidates["GWAS2M_LOCUS_ID"].astype(str).eq(locus)
        ].copy()

    def candidate_rows_for_o(study: str, otid: str) -> pd.DataFrame:
        if pair_candidates.empty:
            return pd.DataFrame()
        return pair_candidates[
            pair_candidates["STUDY_ID"].astype(str).eq(study)
            & pair_candidates["OT_STUDY_LOCUS_ID"].astype(str).eq(otid)
        ].copy()

    for _, r in gcomp.iterrows():
        status = clean_str(r.get("COMPARISON_STATUS"))
        if status == "MATCHED":
            continue
        study = clean_str(r.get("STUDY_ID")); locus = clean_str(r.get("GWAS2M_LOCUS_ID"))
        c = candidate_rows_for_g(study, locus)
        nearest_gap = np.nan; nearest_ot = ""
        same = ot_loci[ot_loci["OT_STUDY_ID"].astype(str).eq(study)].copy() if not ot_loci.empty else pd.DataFrame()
        chrom = canonical_chr(r.get("GWAS2M_CHR"))
        if chrom and not same.empty:
            same = same[same["OT_CHR"].map(canonical_chr).eq(chrom)].copy()
        if not same.empty:
            gaps = same.apply(lambda x: interval_gap(r.get("GWAS2M_START"), r.get("GWAS2M_END"), x.get("OT_START"), x.get("OT_END")), axis=1)
            if gaps.notna().any():
                idx = gaps.astype(float).idxmin(); nearest_gap = finite_float(gaps.loc[idx]); nearest_ot = clean_str(same.loc[idx, "OT_STUDY_LOCUS_ID"])
        if status == "GWAS2M_ONLY_NO_OT_STUDY_COVERAGE":
            why = "NO_OPEN_TARGETS_CREDIBLE_SET_COVERAGE_FOR_STUDY"
        elif not c.empty:
            why = "CANDIDATE_OVERLAP_EXISTED_BUT_LOST_ONE_TO_ONE_PRIMARY_MATCH"
        elif same.empty:
            why = "NO_OPEN_TARGETS_LOCUS_ON_SAME_CHROMOSOME_FOR_STUDY"
        elif math.isfinite(nearest_gap) and nearest_gap > padding:
            why = "NEAREST_OPEN_TARGETS_LOCUS_OUTSIDE_MATCH_WINDOW"
        else:
            why = "NO_ELIGIBLE_VARIANT_OR_INTERVAL_MATCH_UNDER_CURRENT_RULES"
        rows.append({
            "SOURCE": "GWAS2M_ONLY",
            "COMPARISON_STATUS": status,
            "STUDY_ID": study,
            "GWAS2M_LOCUS_ID": locus,
            "OT_STUDY_LOCUS_ID": "",
            "CHR": canonical_chr(r.get("GWAS2M_CHR")),
            "LEAD_VARIANT": clean_str(r.get("GWAS2M_LEAD_VARIANT")),
            "HAS_NONPRIMARY_CANDIDATE_PAIR": "YES" if not c.empty else "NO",
            "N_NONPRIMARY_CANDIDATE_PAIRS": len(c),
            "NEAREST_OTHER_LOCUS_ID": nearest_ot,
            "NEAREST_GAP_BP": nearest_gap,
            "WHY_NOT_PRIMARY_MATCH": why,
        })

    for _, r in ocomp.iterrows():
        status = clean_str(r.get("COMPARISON_STATUS"))
        if status == "MATCHED":
            continue
        study = clean_str(r.get("OT_STUDY_ID")); otid = clean_str(r.get("OT_STUDY_LOCUS_ID"))
        c = candidate_rows_for_o(study, otid)
        nearest_gap = np.nan; nearest_g = ""
        same = gwas_loci[gwas_loci["STUDY_ID"].astype(str).eq(study)].copy() if not gwas_loci.empty else pd.DataFrame()
        chrom = canonical_chr(r.get("OT_CHR"))
        if chrom and not same.empty:
            same = same[same["GWAS2M_CHR"].map(canonical_chr).eq(chrom)].copy()
        if not same.empty:
            gaps = same.apply(lambda x: interval_gap(x.get("GWAS2M_START"), x.get("GWAS2M_END"), r.get("OT_START"), r.get("OT_END")), axis=1)
            if gaps.notna().any():
                idx = gaps.astype(float).idxmin(); nearest_gap = finite_float(gaps.loc[idx]); nearest_g = clean_str(same.loc[idx, "GWAS2M_LOCUS_ID"])
        if status == "OPENTARGETS_ONLY_NO_GWAS2M_STUDY":
            why = "NO_GWAS2M_LOCUS_COVERAGE_FOR_STUDY"
        elif not c.empty:
            why = "CANDIDATE_OVERLAP_EXISTED_BUT_LOST_ONE_TO_ONE_PRIMARY_MATCH"
        elif same.empty:
            why = "NO_GWAS2M_LOCUS_ON_SAME_CHROMOSOME_FOR_STUDY"
        elif math.isfinite(nearest_gap) and nearest_gap > padding:
            why = "NEAREST_GWAS2M_LOCUS_OUTSIDE_MATCH_WINDOW"
        else:
            why = "NO_ELIGIBLE_VARIANT_OR_INTERVAL_MATCH_UNDER_CURRENT_RULES"
        rows.append({
            "SOURCE": "OPEN_TARGETS_ONLY",
            "COMPARISON_STATUS": status,
            "STUDY_ID": study,
            "GWAS2M_LOCUS_ID": "",
            "OT_STUDY_LOCUS_ID": otid,
            "CHR": canonical_chr(r.get("OT_CHR")),
            "LEAD_VARIANT": clean_str(r.get("OT_LEAD_VARIANT")),
            "HAS_NONPRIMARY_CANDIDATE_PAIR": "YES" if not c.empty else "NO",
            "N_NONPRIMARY_CANDIDATE_PAIRS": len(c),
            "NEAREST_OTHER_LOCUS_ID": nearest_g,
            "NEAREST_GAP_BP": nearest_gap,
            "WHY_NOT_PRIMARY_MATCH": why,
        })
    return pd.DataFrame(rows)


def build_performance_metrics(
    gcomp: pd.DataFrame,
    ocomp: pd.DataFrame,
    matches: pd.DataFrame,
    variant_cmp: pd.DataFrame,
    gene_cmp: pd.DataFrame,
    strong_h4: float,
    ot_coloc_complete: bool,
) -> pd.DataFrame:
    """External benchmark metrics using OT as a reference, not a gold standard."""
    rows: list[dict[str, Any]] = []

    def add(level, metric, value, numerator=None, denominator=None, interpretation=""):
        rows.append({
            "LEVEL": level,
            "METRIC": metric,
            "VALUE": value,
            "NUMERATOR": numerator,
            "DENOMINATOR": denominator,
            "PERCENT": (100.0 * numerator / denominator) if numerator is not None and denominator not in (None, 0) else np.nan,
            "INTERPRETATION": interpretation,
        })

    tp = int((gcomp["COMPARISON_STATUS"].astype(str) == "MATCHED").sum()) if not gcomp.empty else 0
    fp = int((gcomp["COMPARISON_STATUS"].astype(str) == "GWAS2M_ONLY_COMPARABLE").sum()) if not gcomp.empty else 0
    fn = int((ocomp["COMPARISON_STATUS"].astype(str) == "OPENTARGETS_ONLY_COMPARABLE").sum()) if not ocomp.empty else 0
    coverage_gap = int((gcomp["COMPARISON_STATUS"].astype(str) == "GWAS2M_ONLY_NO_OT_STUDY_COVERAGE").sum()) if not gcomp.empty else 0
    precision = tp / (tp + fp) if tp + fp else np.nan
    recall = tp / (tp + fn) if tp + fn else np.nan
    f1 = 2 * precision * recall / (precision + recall) if math.isfinite(precision) and math.isfinite(recall) and precision + recall else np.nan
    add("LOCUS", "OT_AS_REFERENCE_TRUE_POSITIVE_MATCHES", tp, tp, tp + fp, "Matched GWAS2m loci among OT-comparable GWAS2m loci")
    add("LOCUS", "OT_AS_REFERENCE_FALSE_POSITIVE_GWAS2M_ONLY", fp, fp, tp + fp, "GWAS2m-only loci in studies where OT has credible-set coverage")
    add("LOCUS", "OT_AS_REFERENCE_FALSE_NEGATIVE_OT_ONLY", fn, fn, tp + fn, "OT credible sets in represented studies not recovered as primary GWAS2m matches")
    add("LOCUS", "OT_REFERENCE_PRECISION", precision, tp, tp + fp, "External concordance only; Open Targets is not biological ground truth")
    add("LOCUS", "OT_REFERENCE_RECALL", recall, tp, tp + fn, "External concordance only; Open Targets is not biological ground truth")
    add("LOCUS", "OT_REFERENCE_F1", f1, None, None, "Harmonic mean of OT-reference precision and recall")
    add("LOCUS", "GWAS2M_LOCI_EXCLUDED_NO_OT_STUDY_COVERAGE", coverage_gap, coverage_gap, len(gcomp), "Coverage gaps excluded from precision/recall")

    nmatch = len(matches)
    exact = int(matches.get("EXACT_LEAD", pd.Series(False, index=matches.index)).map(boolish).sum()) if nmatch else 0
    add("LOCUS", "EXACT_LEAD_AGREEMENT", exact / nmatch if nmatch else np.nan, exact, nmatch, "Same lead variant at matched loci")
    if not variant_cmp.empty:
        detail = variant_cmp[variant_cmp.get("DETAIL_AVAILABLE", pd.Series(False, index=variant_cmp.index)).map(boolish)].copy()
        nd = len(detail)
        top = int(detail.get("TOP_VARIANT_AGREEMENT", pd.Series(False, index=detail.index)).map(boolish).sum()) if nd else 0
        cs = int((safe_numeric(detail.get("N_SHARED_CS95", pd.Series(0, index=detail.index))).fillna(0) > 0).sum()) if nd else 0
        add("VARIANT", "TOP_PIP_VARIANT_AGREEMENT", top / nd if nd else np.nan, top, nd, "Same highest-PIP variant")
        add("VARIANT", "CS95_ANY_OVERLAP", cs / nd if nd else np.nan, cs, nd, "At least one shared 95% credible-set variant")
        for col, metric in [
            ("VARIANT_JACCARD", "MEDIAN_VARIANT_SET_JACCARD"),
            ("CS95_JACCARD", "MEDIAN_CS95_JACCARD"),
            ("PIP_SPEARMAN", "MEDIAN_PIP_SPEARMAN"),
            ("PIP_PEARSON", "MEDIAN_PIP_PEARSON"),
            ("PIP_MAE", "MEDIAN_PIP_MAE"),
        ]:
            vals = safe_numeric(detail.get(col, pd.Series(dtype=float))).dropna()
            add("VARIANT", metric, float(vals.median()) if len(vals) else np.nan, len(vals), nd, f"Median across {len(vals)} matched loci with finite {col}")

    if not gene_cmp.empty:
        h4 = safe_numeric(gene_cmp.get("BEST_H4", pd.Series(np.nan, index=gene_cmp.index)))
        strong = gene_cmp[h4 >= strong_h4].copy()
        strong_match = strong[strong.get("OT_LOCUS_MATCH", pd.Series("NO", index=strong.index)).astype(str).eq("YES")].copy()
        denom = len(strong_match)
        l2g_yes = strong_match.get("OT_L2G_GENE_PRESENT", pd.Series("NO", index=strong_match.index)).astype(str).eq("YES")
        rank = safe_numeric(strong_match.get("OT_L2G_RANK", pd.Series(np.nan, index=strong_match.index)))
        add("GENE", "STRONG_H4_ROWS", len(strong), len(strong), len(gene_cmp), f"BEST_H4 >= {strong_h4:g}")
        add("GENE", "STRONG_H4_AT_MATCHED_OT_LOCUS", denom, denom, len(strong), "Strong candidate rows at matched OT loci")
        add("GENE", "STRONG_GENE_FOUND_ANYWHERE_IN_OT_L2G", (int(l2g_yes.sum()) / denom if denom else np.nan), int(l2g_yes.sum()), denom, "Same strong GWAS2m effector appears anywhere in OT L2G")
        add("GENE", "STRONG_GENE_OT_L2G_TOP1", (int((rank == 1).sum()) / denom if denom else np.nan), int((rank == 1).sum()), denom, "Same strong effector is OT L2G rank 1")
        add("GENE", "STRONG_GENE_OT_L2G_TOP5", (int((rank <= 5).fillna(False).sum()) / denom if denom else np.nan), int((rank <= 5).fillna(False).sum()), denom, "Same strong effector is within OT L2G top 5")
        if ot_coloc_complete:
            coloc_yes = strong_match.get("OT_COLOC_GENE_PRESENT_H4_THRESHOLD", pd.Series("NO", index=strong_match.index)).astype(str).eq("YES")
            add("GENE", "STRONG_GENE_OT_MOLQTL_COLOC_SUPPORTED", (int(coloc_yes.sum()) / denom if denom else np.nan), int(coloc_yes.sum()), denom, "Same strong effector has OT molecular-QTL colocalisation above threshold")
        else:
            add("GENE", "STRONG_GENE_OT_MOLQTL_COLOC_SUPPORTED", np.nan, None, denom, "NOT INTERPRETABLE: OT colocalisation/study data incomplete")

    return pd.DataFrame(rows)


def save_tsv_and_csv(df: pd.DataFrame, outdir: Path, stem: str) -> None:
    df.to_csv(outdir / f"{stem}.tsv", sep="\t", index=False)
    df.to_csv(outdir / f"{stem}.csv", index=False)


# =============================================================================
# OPEN TARGETS L2G
# =============================================================================

def flatten_l2g(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    sid = first_existing(df.columns, ["studyLocusId", "study_locus_id"])
    gene = first_existing(df.columns, ["geneId", "gene_id"])
    score = first_existing(df.columns, ["score", "l2gScore", "l2g_score"])
    if not sid or not gene:
        return pd.DataFrame()
    x = pd.DataFrame({
        "OT_STUDY_LOCUS_ID": df[sid].astype(str),
        "OT_GENE_ID": df[gene].map(gene_base),
        "OT_L2G_SCORE": safe_numeric(df[score]) if score else np.nan,
    })
    x = x[x["OT_GENE_ID"].astype(str).str.len() > 0].copy()
    x = x.sort_values(["OT_STUDY_LOCUS_ID", "OT_L2G_SCORE"], ascending=[True, False], na_position="last")
    x["OT_L2G_RANK"] = x.groupby("OT_STUDY_LOCUS_ID").cumcount() + 1
    return x


# =============================================================================
# OPEN TARGETS COLOCALISATION / STUDY METADATA
# =============================================================================

def flatten_colocalisation(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    cols = df.columns
    sid = first_existing(cols, ["studyLocusId", "study_locus_id"])
    other = first_existing(cols, ["otherStudyLocusId", "other_study_locus_id"])
    nested = first_existing(cols, ["colocalisation", "colocalization", "rows"])
    rows: list[dict[str, Any]] = []

    def row_from_item(left: str, item: dict) -> dict[str, Any]:
        return {
            "OT_STUDY_LOCUS_ID": left,
            "OTHER_STUDY_LOCUS_ID": clean_str(dictionary_get(item, ["otherStudyLocusId", "other_study_locus_id"])),
            "OT_COLOC_METHOD": clean_str(dictionary_get(item, ["colocalisationMethod", "colocalizationMethod", "method"])),
            "OT_H3": finite_float(dictionary_get(item, ["h3", "H3"])),
            "OT_H4": finite_float(dictionary_get(item, ["h4", "H4"])),
            "OT_CLPP": finite_float(dictionary_get(item, ["clpp", "CLPP"])),
        }

    if sid and other:
        for _, rec in df.iterrows():
            rows.append(row_from_item(clean_str(rec.get(sid)), rec.to_dict()))
    elif sid and nested:
        for _, rec in df.iterrows():
            left = clean_str(rec.get(sid))
            for item in as_list(rec.get(nested)):
                if isinstance(item, dict):
                    rows.append(row_from_item(left, item))
    return pd.DataFrame(rows)


def study_metadata(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    sid = first_existing(df.columns, ["studyId", "study_id", "id"])
    if not sid:
        return pd.DataFrame()
    mapping = {
        "OT_OTHER_STUDY_ID": sid,
        "OT_STUDY_TYPE": first_existing(df.columns, ["studyType", "study_type"]),
        "OT_PROJECT_ID": first_existing(df.columns, ["projectId", "project_id"]),
        "OT_BIOSAMPLE_ID": first_existing(df.columns, ["biosampleFromSourceId", "biosample_from_source_id", "biosampleId", "biosample_id"]),
        "OT_TRAIT": first_existing(df.columns, ["traitFromSource", "trait_from_source", "trait"]),
        "OT_GENE_ID_FROM_STUDY": first_existing(df.columns, ["geneId", "gene_id"]),
    }
    out = pd.DataFrame()
    for outcol, source in mapping.items():
        out[outcol] = df[source] if source else ""
    out["OT_OTHER_STUDY_ID"] = out["OT_OTHER_STUDY_ID"].astype(str)
    out["OT_GENE_ID_FROM_STUDY"] = out["OT_GENE_ID_FROM_STUDY"].map(gene_base)
    return out.drop_duplicates("OT_OTHER_STUDY_ID")


def resolve_ot_coloc(coloc: pd.DataFrame, other_credible_meta: pd.DataFrame, study_meta: pd.DataFrame) -> pd.DataFrame:
    if coloc.empty:
        return pd.DataFrame()
    x = coloc.copy()
    if not other_credible_meta.empty:
        meta = other_credible_meta.rename(columns={
            "OT_STUDY_LOCUS_ID": "OTHER_STUDY_LOCUS_ID",
            "OT_STUDY_ID": "OT_OTHER_STUDY_ID",
            "OT_QTL_GENE_ID": "OT_COLOC_GENE_ID",
            "OT_METHOD": "OT_OTHER_FINEMAPPING_METHOD",
        })
        keep = [c for c in ["OTHER_STUDY_LOCUS_ID", "OT_OTHER_STUDY_ID", "OT_COLOC_GENE_ID", "OT_OTHER_FINEMAPPING_METHOD"] if c in meta.columns]
        x = x.merge(meta[keep].drop_duplicates("OTHER_STUDY_LOCUS_ID"), on="OTHER_STUDY_LOCUS_ID", how="left")
    if not study_meta.empty and "OT_OTHER_STUDY_ID" in x.columns:
        x = x.merge(study_meta, on="OT_OTHER_STUDY_ID", how="left")
    if "OT_COLOC_GENE_ID" not in x.columns:
        x["OT_COLOC_GENE_ID"] = ""
    if "OT_GENE_ID_FROM_STUDY" in x.columns:
        x["OT_COLOC_GENE_ID"] = [
            gene_base(a) or gene_base(b)
            for a, b in zip(x["OT_COLOC_GENE_ID"], x["OT_GENE_ID_FROM_STUDY"])
        ]
    return x


# =============================================================================
# MAIN GENE COMPARISON / INTERPRETATION
# =============================================================================

def classify_gene_row(
    locus_matched: bool,
    study_has_ot: bool,
    gene_resolved: bool,
    in_l2g: bool,
    in_ot_coloc: bool,
    best_h4: float,
    strong_h4: float,
    ot_coloc_complete: bool = True,
) -> tuple[str, str]:
    if not study_has_ot:
        return "NO_OT_STUDY_COVERAGE", "COVERAGE_OR_MAPPING_UNRESOLVED"
    if not locus_matched:
        if math.isfinite(best_h4) and best_h4 >= strong_h4:
            return "GWAS2M_ONLY_LOCUS_IN_OT_COVERED_STUDY", "D_POTENTIAL_NEW_LOCUS_CANDIDATE_REQUIRES_EXTERNAL_VALIDATION"
        return "GWAS2M_ONLY_LOCUS_WEAK_OR_UNTESTED", "COVERAGE_OR_MAPPING_UNRESOLVED"
    if not gene_resolved:
        return "OT_LOCUS_MATCH_EFFECTOR_UNRESOLVED", "COVERAGE_OR_MAPPING_UNRESOLVED"
    if in_ot_coloc:
        return "OT_RECAPITULATED_GENE_AND_MOLQTL", "A_ESTABLISHED_OT_GENE_AND_MOLQTL"
    # If any OT colocalisation shard is corrupt/missing, absence of a same-gene
    # colocalisation record is NOT interpretable as biological absence.
    if not ot_coloc_complete:
        if in_l2g:
            return "OT_L2G_SUPPORTED_MOLQTL_DATA_INCOMPLETE", "COVERAGE_OR_MAPPING_UNRESOLVED"
        return "OT_LOCUS_MATCH_OT_MOLQTL_DATA_INCOMPLETE", "COVERAGE_OR_MAPPING_UNRESOLVED"
    if in_l2g:
        return "OT_SUPPORTED_GENE_NO_SAME_GENE_MOLQTL", "B_OT_SUPPORTED_GENE_POTENTIAL_NEW_MECHANISM"
    return "OT_LOCUS_MATCH_GENE_NOT_IN_OT_EVIDENCE", "C_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE"


def compare_genes(
    gene_master: pd.DataFrame,
    gwas_loci: pd.DataFrame,
    matches: pd.DataFrame,
    ot_loci: pd.DataFrame,
    l2g: pd.DataFrame,
    ot_coloc: pd.DataFrame,
    l2g_top_k: int,
    ot_coloc_h4: float,
    strong_h4: float,
    ot_coloc_complete: bool = True,
) -> pd.DataFrame:
    match_map = {
        (clean_str(r["STUDY_ID"]), clean_str(r["GWAS2M_LOCUS_ID"])): clean_str(r["OT_STUDY_LOCUS_ID"])
        for _, r in matches.iterrows()
    } if not matches.empty else {}
    match_method = {
        (clean_str(r["STUDY_ID"]), clean_str(r["GWAS2M_LOCUS_ID"])): clean_str(r.get("OT_METHOD"))
        for _, r in matches.iterrows()
    } if not matches.empty else {}
    match_basis = {
        (clean_str(r["STUDY_ID"]), clean_str(r["GWAS2M_LOCUS_ID"])): clean_str(r.get("MATCH_BASIS"))
        for _, r in matches.iterrows()
    } if not matches.empty else {}

    studies_with_ot = set(ot_loci["OT_STUDY_ID"].dropna().astype(str)) if not ot_loci.empty else set()
    rows: list[dict[str, Any]] = []

    for _, r in gene_master.iterrows():
        study = clean_str(r.get("STUDY_ACCESSION")); locus = clean_str(r.get("LOCUS_ID"))
        gene = gene_base(r.get("EFFECTOR_GENE_ID"))
        otid = match_map.get((study, locus), "")
        locus_matched = bool(otid)
        study_has_ot = study in studies_with_ot

        lg = l2g[l2g["OT_STUDY_LOCUS_ID"].astype(str).eq(otid)].copy() if otid and not l2g.empty else pd.DataFrame()
        cg = ot_coloc[ot_coloc["OT_STUDY_LOCUS_ID"].astype(str).eq(otid)].copy() if otid and not ot_coloc.empty else pd.DataFrame()

        in_l2g = False; l2g_rank = np.nan; l2g_score = np.nan
        top1 = ""; top5 = ""
        if not lg.empty:
            lg = lg.sort_values("OT_L2G_RANK")
            top1 = clean_str(lg.iloc[0].get("OT_GENE_ID"))
            top5 = semicolon_unique(lg.head(l2g_top_k)["OT_GENE_ID"])
            same = lg[lg["OT_GENE_ID"].astype(str).eq(gene)] if gene else pd.DataFrame()
            if not same.empty:
                in_l2g = True
                l2g_rank = finite_float(same.iloc[0].get("OT_L2G_RANK"))
                l2g_score = finite_float(same.iloc[0].get("OT_L2G_SCORE"))

        in_coloc = False; coloc_max_h4 = np.nan; coloc_types = ""; coloc_biosamples = ""; coloc_traits = ""
        if not cg.empty and gene:
            same = cg[cg["OT_COLOC_GENE_ID"].astype(str).eq(gene)].copy()
            if not same.empty:
                h = safe_numeric(same["OT_H4"])
                coloc_max_h4 = float(h.max()) if h.notna().any() else np.nan
                in_coloc = math.isfinite(coloc_max_h4) and coloc_max_h4 >= ot_coloc_h4
                coloc_types = semicolon_unique(same.get("OT_STUDY_TYPE", pd.Series(dtype=str)), limit=50)
                coloc_biosamples = semicolon_unique(same.get("OT_BIOSAMPLE_ID", pd.Series(dtype=str)), limit=50)
                coloc_traits = semicolon_unique(same.get("OT_TRAIT", pd.Series(dtype=str)), limit=50)

        best_h4 = finite_float(r.get("BEST_H4"))
        evidence_class, novelty = classify_gene_row(
            locus_matched, study_has_ot, bool(gene), in_l2g, in_coloc,
            best_h4, strong_h4, ot_coloc_complete=ot_coloc_complete
        )

        out = r.to_dict()
        out.update({
            "OT_STUDY_COVERAGE": "YES" if study_has_ot else "NO",
            "OT_LOCUS_MATCH": "YES" if locus_matched else "NO",
            "OT_STUDY_LOCUS_ID": otid,
            "OT_METHOD": match_method.get((study, locus), ""),
            "OT_MATCH_BASIS": match_basis.get((study, locus), ""),
            "OT_L2G_GENE_PRESENT": "YES" if in_l2g else "NO",
            "OT_L2G_RANK": l2g_rank,
            "OT_L2G_SCORE": l2g_score,
            "OT_L2G_TOP1_GENE": top1,
            f"OT_L2G_TOP{l2g_top_k}_GENES": top5,
            "OT_COLOC_DATASET_COMPLETE": "YES" if ot_coloc_complete else "NO",
            "OT_COLOC_GENE_PRESENT_H4_THRESHOLD": (
                "YES" if in_coloc else ("NO" if ot_coloc_complete else "UNKNOWN_INCOMPLETE_OT_DATA")
            ),
            "OT_COLOC_MAX_H4_SAME_GENE": coloc_max_h4,
            "OT_COLOC_STUDY_TYPES": coloc_types,
            "OT_COLOC_BIOSAMPLES": coloc_biosamples,
            "OT_COLOC_TRAITS": coloc_traits,
            "OPEN_TARGETS_EVIDENCE_CLASS": evidence_class,
            "NOVELTY_TIER_STEP13": novelty,
        })
        rows.append(out)

    out = pd.DataFrame(rows)
    if not out.empty:
        tier_order = {
            "C_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE": 0,
            "D_POTENTIAL_NEW_LOCUS_CANDIDATE_REQUIRES_EXTERNAL_VALIDATION": 1,
            "B_OT_SUPPORTED_GENE_POTENTIAL_NEW_MECHANISM": 2,
            "A_ESTABLISHED_OT_GENE_AND_MOLQTL": 3,
            "COVERAGE_OR_MAPPING_UNRESOLVED": 4,
        }
        out["_TIER_ORDER"] = out["NOVELTY_TIER_STEP13"].map(tier_order).fillna(9)
        out = out.sort_values(["_TIER_ORDER", "BEST_H4", "N_SHARED_QTL_TYPES"], ascending=[True, False, False], na_position="last").drop(columns="_TIER_ORDER")
    return out.reset_index(drop=True)


# =============================================================================
# STEP12 v2 MECHANISM / CELL-CONTEXT / SIGNAL COMPARISON
# =============================================================================

def build_study_comparability(
    gwas_loci: pd.DataFrame,
    ot_loci: pd.DataFrame,
    matches: pd.DataFrame,
) -> pd.DataFrame:
    studies = sorted(
        set(gwas_loci.get("STUDY_ID", pd.Series(dtype=str)).dropna().astype(str))
        | set(ot_loci.get("OT_STUDY_ID", pd.Series(dtype=str)).dropna().astype(str))
    )
    rows: list[dict[str, Any]] = []
    for study in studies:
        ng = int(gwas_loci.get("STUDY_ID", pd.Series(dtype=str)).astype(str).eq(study).sum()) if not gwas_loci.empty else 0
        no = int(ot_loci.get("OT_STUDY_ID", pd.Series(dtype=str)).astype(str).eq(study).sum()) if not ot_loci.empty else 0
        nm = int(matches.get("STUDY_ID", pd.Series(dtype=str)).astype(str).eq(study).sum()) if not matches.empty else 0
        if ng > 0 and no > 0:
            cls = "COMPARABLE_BOTH_HAVE_LOCI"
        elif ng > 0 and no == 0:
            cls = "GWAS2M_HAS_LOCI_OT_NO_CREDIBLE_SETS"
        elif ng == 0 and no > 0:
            cls = "OT_HAS_CREDIBLE_SETS_GWAS2M_NOT_FINEMAPPED"
        else:
            cls = "NO_LOCI_BOTH"
        rows.append({
            "STUDY_ID": study,
            "N_GWAS2M_LOCI": ng,
            "N_OT_CREDIBLE_SETS": no,
            "N_PRIMARY_MATCHES": nm,
            "COMPARABILITY_CLASS": cls,
        })
    return pd.DataFrame(rows)


def build_complete_locus_inventory(
    gcomp: pd.DataFrame,
    ocomp: pd.DataFrame,
    matches: pd.DataFrame,
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    if not matches.empty:
        z = matches.copy()
        z["SOURCE"] = "MATCHED_PAIR"
        z["INVENTORY_CLASS"] = "MATCHED"
        parts.append(z)
    if not gcomp.empty:
        z = gcomp[~gcomp["COMPARISON_STATUS"].astype(str).eq("MATCHED")].copy()
        z["SOURCE"] = "GWAS2M"
        z["INVENTORY_CLASS"] = z["COMPARISON_STATUS"]
        parts.append(z)
    if not ocomp.empty:
        z = ocomp[~ocomp["COMPARISON_STATUS"].astype(str).eq("MATCHED")].copy()
        z["SOURCE"] = "OPEN_TARGETS"
        z["INVENTORY_CLASS"] = z["COMPARISON_STATUS"]
        parts.append(z)
    return pd.concat(parts, ignore_index=True, sort=False) if parts else pd.DataFrame()


def method_summary_v2(
    variant_cmp: pd.DataFrame,
    gene_cmp: pd.DataFrame,
) -> pd.DataFrame:
    if variant_cmp.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for method, v in variant_cmp.groupby("OT_METHOD", dropna=False):
        method_text = clean_str(method)
        is_pics = "pics" in method_text.lower()
        g = (
            gene_cmp[gene_cmp.get("OT_METHOD", pd.Series("", index=gene_cmp.index)).astype(str).eq(method_text)].copy()
            if not gene_cmp.empty else pd.DataFrame()
        )
        top_agree = v.get("TOP_VARIANT_AGREEMENT", pd.Series(False, index=v.index)).map(boolish)
        cs = safe_numeric(v.get("CS95_JACCARD", pd.Series(np.nan, index=v.index)))
        sp = safe_numeric(v.get("PIP_SPEARMAN", pd.Series(np.nan, index=v.index)))
        row = {
            "OT_METHOD": method_text,
            "METHOD_INTERPRETATION": (
                "PICS: lead/position/variant overlap are primary; interval/CS metrics are secondary"
                if is_pics else
                "SuSiE/other CS method: lead, interval, CS overlap and PIP agreement are interpretable"
            ),
            "N_MATCHED_PAIRS": int(len(v)),
            "TOP_VARIANT_AGREEMENT_N": int(top_agree.sum()),
            "TOP_VARIANT_AGREEMENT_RATE": float(top_agree.mean()) if len(v) else np.nan,
            "MEDIAN_VARIANT_JACCARD": float(safe_numeric(v.get("VARIANT_JACCARD", pd.Series(dtype=float))).median()),
            "MEDIAN_CS95_JACCARD": np.nan if is_pics else float(cs.median()),
            "MEDIAN_PIP_SPEARMAN": float(sp.median()),
            "N_STEP12_GENES": int(len(g)),
            "N_OT_L2G_PRESENT": int(g.get("OT_L2G_GENE_PRESENT", pd.Series("", index=g.index)).astype(str).eq("YES").sum()) if not g.empty else 0,
            "N_OT_COLOC_PRESENT": int(g.get("OT_COLOC_GENE_PRESENT_H4_THRESHOLD", pd.Series("", index=g.index)).astype(str).eq("YES").sum()) if not g.empty else 0,
        }
        rows.append(row)
    return pd.DataFrame(rows)


def _h4_flag(value: Any, threshold: float) -> bool:
    z = finite_float(value)
    return math.isfinite(z) and z >= threshold


def build_mechanism_comparison(
    gene_cmp: pd.DataFrame,
    strong_h4: float,
    suggestive_h4: float,
) -> pd.DataFrame:
    """Add explicit bulk-vs-cell and external-support interpretation.

    This is a screening layer, not a novelty claim.  In particular, absence of
    an OT same-gene molQTL record does not mean the GWAS2m context is disproved.
    """
    if gene_cmp.empty:
        return pd.DataFrame()
    out = gene_cmp.copy()
    classes: list[str] = []
    patterns: list[str] = []
    cell_screen: list[str] = []

    for _, r in out.iterrows():
        bulk_strong = _h4_flag(r.get("BULK_BEST_H4"), strong_h4)
        cell_strong = _h4_flag(r.get("CELL_BEST_H4"), strong_h4)
        bulk_shared = _h4_flag(r.get("BULK_BEST_H4"), suggestive_h4)
        cell_shared = _h4_flag(r.get("CELL_BEST_H4"), suggestive_h4)
        best_strong = _h4_flag(r.get("BEST_H4"), strong_h4)
        best_shared = _h4_flag(r.get("BEST_H4"), suggestive_h4)
        in_l2g = clean_str(r.get("OT_L2G_GENE_PRESENT")) == "YES"
        in_ot_coloc = clean_str(r.get("OT_COLOC_GENE_PRESENT_H4_THRESHOLD")) == "YES"
        locus_match = clean_str(r.get("OT_LOCUS_MATCH")) == "YES"
        coverage = clean_str(r.get("OT_STUDY_COVERAGE")) == "YES"
        ot_complete = clean_str(r.get("OT_COLOC_DATASET_COMPLETE")) == "YES"

        if bulk_strong and cell_strong:
            pattern = "BULK_AND_CELL_STRONG"
        elif cell_strong:
            pattern = "CELL_STRONG"
        elif bulk_strong:
            pattern = "BULK_STRONG"
        elif bulk_shared and cell_shared:
            pattern = "BULK_AND_CELL_SUGGESTIVE"
        elif cell_shared:
            pattern = "CELL_SUGGESTIVE"
        elif bulk_shared:
            pattern = "BULK_SUGGESTIVE"
        elif best_shared:
            pattern = "OTHER_SHARED_SIGNAL"
        else:
            pattern = "NO_SUGGESTIVE_SHARED_SIGNAL"
        patterns.append(pattern)

        if not coverage:
            cls = "UNRESOLVED_NO_OT_STUDY_COVERAGE"
        elif not locus_match:
            cls = (
                "GWAS2M_ONLY_LOCUS_STRONG_REQUIRES_EXTERNAL_VALIDATION"
                if best_strong else "GWAS2M_ONLY_LOCUS_WEAK_OR_UNTESTED"
            )
        elif in_ot_coloc:
            if bulk_strong and cell_strong:
                cls = "CONVERGENT_OT_MOLQTL_PLUS_BULK_AND_CELL"
            elif cell_strong:
                cls = "CONVERGENT_OT_MOLQTL_PLUS_CELL_SIGNAL"
            elif bulk_strong:
                cls = "CONVERGENT_OT_MOLQTL_PLUS_BULK_SIGNAL"
            else:
                cls = "CONVERGENT_OT_MOLQTL_SAME_GENE"
        elif not ot_complete:
            cls = "UNRESOLVED_OT_MOLQTL_DATA_INCOMPLETE"
        elif in_l2g and cell_strong:
            cls = "OT_GENE_SUPPORTED_NEW_CELL_CONTEXT_CANDIDATE"
        elif in_l2g and best_strong:
            cls = "OT_GENE_SUPPORTED_NEW_MOLECULAR_MECHANISM_CANDIDATE"
        elif in_l2g and best_shared:
            cls = "OT_GENE_SUPPORTED_SUGGESTIVE_GWAS2M_MECHANISM"
        elif locus_match and best_strong:
            cls = "KNOWN_OT_LOCUS_NEW_EFFECTOR_WITH_STRONG_GWAS2M_MECHANISM"
        else:
            cls = "OT_LOCUS_MATCH_NO_STRONG_EXTERNAL_GENE_SUPPORT"
        classes.append(cls)

        if finite_float(r.get("CELL_BEST_H4")) != finite_float(r.get("CELL_BEST_H4")):
            cell_screen.append("NO_FINITE_CELL_H4")
        elif in_ot_coloc:
            cell_screen.append("SAME_GENE_OT_MOLQTL_PRESENT_CONTEXT_NEEDS_MANUAL_ALIGNMENT")
        elif not ot_complete:
            cell_screen.append("OT_MOLQTL_DATA_INCOMPLETE")
        elif in_l2g and cell_strong:
            cell_screen.append("POTENTIAL_NEW_CELL_CONTEXT_FOR_OT_SUPPORTED_GENE")
        elif cell_strong:
            cell_screen.append("STRONG_CELL_CONTEXT_REQUIRES_EXTERNAL_VALIDATION")
        else:
            cell_screen.append("CELL_CONTEXT_NOT_STRONG")

    out["GWAS2M_MECHANISM_PATTERN"] = patterns
    out["EXTERNAL_MECHANISM_CLASS"] = classes
    out["CELL_CONTEXT_SCREEN"] = cell_screen
    return out


def _context_text_status(context: Any, biosamples: Any, traits: Any) -> str:
    c = norm_text(context)
    if not c:
        return "NO_GWAS2M_CELL_CONTEXT"
    external = norm_text(f"{clean_str(biosamples)} {clean_str(traits)}")
    if not external:
        return "OT_CONTEXT_METADATA_UNAVAILABLE"
    # Conservative textual screening only; it is not ontology alignment.
    if len(c) >= 4 and (c in external or external in c):
        return "TEXTUAL_CONTEXT_MATCH"
    c_tokens = {t for t in c.split() if len(t) >= 4}
    e_tokens = {t for t in external.split() if len(t) >= 4}
    if c_tokens and len(c_tokens & e_tokens) >= max(1, min(2, len(c_tokens))):
        return "PARTIAL_TEXTUAL_CONTEXT_MATCH"
    return "NO_TEXTUAL_MATCH_NOT_NEGATIVE"


def build_cell_context_comparison(mechanism: pd.DataFrame) -> pd.DataFrame:
    if mechanism.empty:
        return pd.DataFrame()
    ncell = safe_numeric(mechanism.get("CELL_N_CANDIDATES", pd.Series(0, index=mechanism.index))).fillna(0)
    ch4 = safe_numeric(mechanism.get("CELL_BEST_H4", pd.Series(np.nan, index=mechanism.index)))
    x = mechanism[(ncell > 0) | ch4.notna()].copy()
    if x.empty:
        return x
    x["OT_CONTEXT_TEXT_SCREEN"] = [
        _context_text_status(c, b, t)
        for c, b, t in zip(
            x.get("CELL_BEST_CONTEXT", pd.Series("", index=x.index)),
            x.get("OT_COLOC_BIOSAMPLES", pd.Series("", index=x.index)),
            x.get("OT_COLOC_TRAITS", pd.Series("", index=x.index)),
        )
    ]
    return x.sort_values("CELL_BEST_H4", ascending=False, na_position="last").reset_index(drop=True)


def build_signal_external_comparison(
    signal_master: pd.DataFrame,
    gene_cmp: pd.DataFrame,
) -> pd.DataFrame:
    if signal_master.empty:
        return pd.DataFrame()
    sig = signal_master.copy()
    sig["__GENE_KEY"] = sig.get("EFFECTOR_GENE_ID", pd.Series("", index=sig.index)).map(gene_base)

    if gene_cmp.empty:
        return sig.drop(columns="__GENE_KEY")

    ext_cols = [c for c in [
        "STUDY_ACCESSION", "LOCUS_ID", "EFFECTOR_GENE_ID",
        "OT_STUDY_COVERAGE", "OT_LOCUS_MATCH", "OT_STUDY_LOCUS_ID", "OT_METHOD",
        "OT_MATCH_BASIS", "OT_L2G_GENE_PRESENT", "OT_L2G_RANK", "OT_L2G_SCORE",
        "OT_L2G_TOP1_GENE", "OT_COLOC_DATASET_COMPLETE",
        "OT_COLOC_GENE_PRESENT_H4_THRESHOLD", "OT_COLOC_MAX_H4_SAME_GENE",
        "OT_COLOC_STUDY_TYPES", "OT_COLOC_BIOSAMPLES", "OT_COLOC_TRAITS",
        "OPEN_TARGETS_EVIDENCE_CLASS", "NOVELTY_TIER_STEP13",
        "EXTERNAL_MECHANISM_CLASS", "CELL_CONTEXT_SCREEN",
    ] if c in gene_cmp.columns]
    gm = gene_cmp[ext_cols].copy()
    gm["__GENE_KEY"] = gm.get("EFFECTOR_GENE_ID", pd.Series("", index=gm.index)).map(gene_base)
    gm = gm.sort_values(
        "BEST_H4" if "BEST_H4" in gm.columns else "STUDY_ACCESSION",
        ascending=False if "BEST_H4" in gm.columns else True,
        na_position="last",
    ).drop_duplicates(["STUDY_ACCESSION", "LOCUS_ID", "__GENE_KEY"], keep="first")
    if "EFFECTOR_GENE_ID" in gm.columns:
        gm = gm.drop(columns="EFFECTOR_GENE_ID")

    out = sig.merge(
        gm,
        on=["STUDY_ACCESSION", "LOCUS_ID", "__GENE_KEY"],
        how="left",
        validate="many_to_one",
    )
    return out.drop(columns="__GENE_KEY")



# =============================================================================
# SUMMARY
# =============================================================================

def build_locus_comparison(gwas_loci: pd.DataFrame, ot_loci: pd.DataFrame, matches: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mm = matches.copy() if not matches.empty else pd.DataFrame(columns=["STUDY_ID", "GWAS2M_LOCUS_ID", "OT_STUDY_LOCUS_ID"])
    gkey = set(zip(mm.get("STUDY_ID", pd.Series(dtype=str)).astype(str), mm.get("GWAS2M_LOCUS_ID", pd.Series(dtype=str)).astype(str)))
    oids = set(mm.get("OT_STUDY_LOCUS_ID", pd.Series(dtype=str)).astype(str))
    studies_with_ot = set(ot_loci["OT_STUDY_ID"].astype(str)) if not ot_loci.empty else set()

    grow: list[dict[str, Any]] = []
    for _, r in gwas_loci.iterrows():
        study = clean_str(r.get("STUDY_ID")); locus = clean_str(r.get("GWAS2M_LOCUS_ID"))
        matched = (study, locus) in gkey
        z = r.to_dict()
        z["OT_STUDY_COVERAGE"] = "YES" if study in studies_with_ot else "NO"
        z["COMPARISON_STATUS"] = "MATCHED" if matched else ("GWAS2M_ONLY_COMPARABLE" if study in studies_with_ot else "GWAS2M_ONLY_NO_OT_STUDY_COVERAGE")
        grow.append(z)
    gcomp = pd.DataFrame(grow)

    orow: list[dict[str, Any]] = []
    studies_with_g = set(gwas_loci["STUDY_ID"].astype(str)) if not gwas_loci.empty else set()
    for _, r in ot_loci.iterrows():
        otid = clean_str(r.get("OT_STUDY_LOCUS_ID")); study = clean_str(r.get("OT_STUDY_ID"))
        z = r.to_dict()
        z["COMPARISON_STATUS"] = "MATCHED" if otid in oids else ("OPENTARGETS_ONLY_COMPARABLE" if study in studies_with_g else "OPENTARGETS_ONLY_NO_GWAS2M_STUDY")
        orow.append(z)
    ocomp = pd.DataFrame(orow)
    return gcomp, ocomp, mm


def summary_table(gene_cmp: pd.DataFrame) -> pd.DataFrame:
    if gene_cmp.empty:
        return pd.DataFrame()
    return (
        gene_cmp.groupby("NOVELTY_TIER_STEP13", dropna=False)
        .agg(
            N_GENES=("EFFECTOR_KEY", "size"),
            N_LOCI=("LOCUS_ID", "nunique"),
            N_STUDIES=("STUDY_ACCESSION", "nunique"),
            MAX_H4=("BEST_H4", "max"),
            MEDIAN_H4=("BEST_H4", "median"),
        )
        .reset_index()
        .sort_values("N_GENES", ascending=False)
    )


def locus_benchmark_summary(
    gcomp: pd.DataFrame,
    ocomp: pd.DataFrame,
    matches: pd.DataFrame,
    variant_cmp: pd.DataFrame,
) -> pd.DataFrame:
    """Human-readable locus benchmark counts.

    IMPORTANT: GWAS2m loci and OT credible sets are distinct objects, so the
    unmatched totals are reported separately rather than forced into a single
    denominator. Primary matching is one-to-one within each study.
    """
    def nstatus(df: pd.DataFrame, value: str) -> int:
        if df.empty or "COMPARISON_STATUS" not in df.columns:
            return 0
        return int(df["COMPARISON_STATUS"].astype(str).eq(value).sum())

    n_g = int(len(gcomp))
    n_o = int(len(ocomp))
    n_m = int(len(matches))
    exact = int(matches.get("EXACT_LEAD", pd.Series(False, index=matches.index)).map(boolish).sum()) if not matches.empty else 0
    shared = int((safe_numeric(matches.get("N_SHARED_VARIANTS", pd.Series(0, index=matches.index))).fillna(0) > 0).sum()) if not matches.empty else 0

    top_agree = 0
    cs_overlap = 0
    detail = 0
    if not variant_cmp.empty:
        detail = int(variant_cmp.get("DETAIL_AVAILABLE", pd.Series(False, index=variant_cmp.index)).map(boolish).sum())
        top_agree = int(variant_cmp.get("TOP_VARIANT_AGREEMENT", pd.Series(False, index=variant_cmp.index)).map(boolish).sum())
        cs_overlap = int((safe_numeric(variant_cmp.get("N_SHARED_CS95", pd.Series(0, index=variant_cmp.index))).fillna(0) > 0).sum())

    rows = [
        {"LEVEL": "LOCUS", "CATEGORY": "GWAS2M_TOTAL_LOCI", "COUNT": n_g, "INTERPRETATION": "Independent GWAS2m locus rows"},
        {"LEVEL": "LOCUS", "CATEGORY": "OPEN_TARGETS_TOTAL_CREDIBLE_SETS", "COUNT": n_o, "INTERPRETATION": "OT credible sets for the same queried GWAS study IDs"},
        {"LEVEL": "LOCUS", "CATEGORY": "COMMON_PRIMARY_MATCHED", "COUNT": n_m, "INTERPRETATION": "One-to-one primary GWAS2m <-> OT matches"},
        {"LEVEL": "LOCUS", "CATEGORY": "GWAS2M_ONLY_TOTAL", "COUNT": max(0, n_g - n_m), "INTERPRETATION": "GWAS2m loci without a primary OT match"},
        {"LEVEL": "LOCUS", "CATEGORY": "GWAS2M_ONLY_COMPARABLE", "COUNT": nstatus(gcomp, "GWAS2M_ONLY_COMPARABLE"), "INTERPRETATION": "OT covers the study, but no OT credible set matched this GWAS2m locus"},
        {"LEVEL": "LOCUS", "CATEGORY": "GWAS2M_ONLY_NO_OT_STUDY_COVERAGE", "COUNT": nstatus(gcomp, "GWAS2M_ONLY_NO_OT_STUDY_COVERAGE"), "INTERPRETATION": "Cannot call an OT miss; the study has no OT credible-set coverage"},
        {"LEVEL": "LOCUS", "CATEGORY": "OPEN_TARGETS_ONLY_TOTAL", "COUNT": max(0, n_o - n_m), "INTERPRETATION": "OT credible sets without a primary GWAS2m match"},
        {"LEVEL": "LOCUS", "CATEGORY": "OPEN_TARGETS_ONLY_COMPARABLE", "COUNT": nstatus(ocomp, "OPENTARGETS_ONLY_COMPARABLE"), "INTERPRETATION": "OT locus present in a GWAS2m-represented study but not recovered/matched by GWAS2m"},
        {"LEVEL": "LOCUS", "CATEGORY": "OPEN_TARGETS_ONLY_NO_GWAS2M_STUDY", "COUNT": nstatus(ocomp, "OPENTARGETS_ONLY_NO_GWAS2M_STUDY"), "INTERPRETATION": "OT study was not represented in the GWAS2m locus set"},
        {"LEVEL": "MATCH_DETAIL", "CATEGORY": "EXACT_LEAD_VARIANT_MATCH", "COUNT": exact, "INTERPRETATION": "Matched loci with identical lead variant"},
        {"LEVEL": "MATCH_DETAIL", "CATEGORY": "ANY_SHARED_VARIANT", "COUNT": shared, "INTERPRETATION": "Primary matches sharing at least one variant"},
        {"LEVEL": "MATCH_DETAIL", "CATEGORY": "VARIANT_DETAIL_AVAILABLE", "COUNT": detail, "INTERPRETATION": "Primary matches with variant-level data available in both pipelines"},
        {"LEVEL": "MATCH_DETAIL", "CATEGORY": "TOP_PIP_VARIANT_AGREEMENT", "COUNT": top_agree, "INTERPRETATION": "Matched loci whose highest-PIP variant is identical"},
        {"LEVEL": "MATCH_DETAIL", "CATEGORY": "CS95_HAS_SHARED_VARIANT", "COUNT": cs_overlap, "INTERPRETATION": "Matched loci whose 95% credible sets overlap by at least one variant"},
    ]
    return pd.DataFrame(rows)


def match_basis_summary(matches: pd.DataFrame) -> pd.DataFrame:
    if matches.empty or "MATCH_BASIS" not in matches.columns:
        return pd.DataFrame(columns=["MATCH_BASIS", "N_MATCHED_LOCI"])
    return (
        matches["MATCH_BASIS"].fillna("UNKNOWN").astype(str)
        .value_counts(dropna=False)
        .rename_axis("MATCH_BASIS")
        .reset_index(name="N_MATCHED_LOCI")
    )


def gene_benchmark_summary(gene_cmp: pd.DataFrame, strong_h4: float) -> pd.DataFrame:
    if gene_cmp.empty:
        return pd.DataFrame(columns=["CATEGORY", "COUNT", "INTERPRETATION"])

    matched = gene_cmp.get("OT_LOCUS_MATCH", pd.Series("NO", index=gene_cmp.index)).astype(str).eq("YES")
    gene_resolved = gene_cmp.get("EFFECTOR_GENE_ID", pd.Series("", index=gene_cmp.index)).map(gene_base).astype(str).str.len() > 0
    l2g_yes = gene_cmp.get("OT_L2G_GENE_PRESENT", pd.Series("NO", index=gene_cmp.index)).astype(str).eq("YES")
    coloc_state = gene_cmp.get("OT_COLOC_GENE_PRESENT_H4_THRESHOLD", pd.Series("", index=gene_cmp.index)).astype(str)
    h4 = safe_numeric(gene_cmp.get("BEST_H4", pd.Series(np.nan, index=gene_cmp.index)))
    strong = h4 >= strong_h4

    def n(mask):
        return int(pd.Series(mask, index=gene_cmp.index).fillna(False).sum())

    rows = [
        {"CATEGORY": "STEP12_LOCUS_GENE_ROWS", "COUNT": len(gene_cmp), "INTERPRETATION": "All Step12 locus x candidate-effector rows"},
        {"CATEGORY": "ROWS_AT_MATCHED_OT_LOCI", "COUNT": n(matched), "INTERPRETATION": "Candidate-gene rows whose GWAS2m locus has an OT primary match"},
        {"CATEGORY": "RESOLVED_EFFECTOR_GENE_ROWS", "COUNT": n(gene_resolved), "INTERPRETATION": "Rows with a resolvable effector gene ID"},
        {"CATEGORY": "SAME_GENE_FOUND_IN_OT_L2G", "COUNT": n(matched & gene_resolved & l2g_yes), "INTERPRETATION": "GWAS2m effector also appears in OT L2G for the matched locus"},
        {"CATEGORY": "SAME_GENE_NOT_IN_OT_L2G", "COUNT": n(matched & gene_resolved & ~l2g_yes), "INTERPRETATION": "Potentially different effector prioritisation; review strong-H4 rows first"},
        {"CATEGORY": "SAME_GENE_OT_COLOC_H4_SUPPORTED", "COUNT": n(matched & gene_resolved & coloc_state.eq("YES")), "INTERPRETATION": "OT also has same-gene molecular-QTL colocalisation above threshold"},
        {"CATEGORY": "OT_COLOC_COMPARISON_UNKNOWN_INCOMPLETE_DATA", "COUNT": n(coloc_state.str.startswith("UNKNOWN")), "INTERPRETATION": "Cannot interpret OT molQTL absence because OT colocalisation data are incomplete"},
        {"CATEGORY": "STRONG_STEP12_H4_ROWS", "COUNT": n(strong), "INTERPRETATION": f"Step12 BEST_H4 >= {strong_h4:g}"},
        {"CATEGORY": "STRONG_H4_AT_MATCHED_OT_LOCUS", "COUNT": n(strong & matched), "INTERPRETATION": "Strong GWAS2m molecular evidence at a locus also present in OT"},
        {"CATEGORY": "STRONG_H4_AT_GWAS2M_ONLY_LOCUS", "COUNT": n(strong & ~matched), "INTERPRETATION": "High-priority candidates for external novelty validation"},
    ]
    return pd.DataFrame(rows)


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    args = arguments()
    if PIPELINE["study"] and not args.study:
        args.study = PIPELINE["study"]
    if not 0 <= args.strong_h4 <= 1 or not 0 <= args.suggestive_h4 <= 1 or not 0 <= args.ot_coloc_h4 <= 1:
        raise SystemExit("H4 thresholds must be between 0 and 1")
    if args.suggestive_h4 > args.strong_h4:
        raise SystemExit("--suggestive-h4 must be <= --strong-h4")
    if args.l2g_top_k < 1:
        raise SystemExit("--l2g-top-k must be >= 1")

    root = Path(args.root).resolve()
    phenotype = args.phenotype.strip()
    ancestry_code, ancestry_label = canonical_ancestry(args.ancestry)
    paths = pipeline_paths(root, phenotype, ancestry_label, args.ot_cache or None, args.ot_release)
    paths["OUT"].mkdir(parents=True, exist_ok=True)

    banner("GWAS2m STEP13 - OPEN TARGETS COMPARISON")
    print(f"Version       : {VERSION}")
    print(f"Root          : {root}")
    print(f"Phenotype     : {phenotype}")
    print(f"Ancestry      : {ancestry_label} ({ancestry_code})")
    print(f"Step12        : {paths['S12']}")
    print(f"Open Targets  : {paths['OT_CACHE']}")
    print(f"Output        : {paths['OUT']}")

    (
        genes, locus_master, step12_variants, step12_signals,
        step12_molecular, step12_strong, step12_pathway,
    ) = load_step12(paths, args.study)
    studies = sorted(set(genes["STUDY_ACCESSION"].dropna().astype(str)))
    gwas_variants = load_gwas_variants(paths, step12_variants, studies)
    gwas_loci = build_gwas_loci(locus_master, gwas_variants)

    # Resolve OT local datasets; download missing release datasets when enabled.
    resolved = ensure_ot_datasets(
        paths["OT_CACHE"], args.ot_release, args.download_ot, args.rsync_base
    )
    if resolved["credible_set"] is None or resolved["l2g_prediction"] is None:
        missing = [k for k in ["credible_set", "l2g_prediction"] if resolved[k] is None]
        raise SystemExit(
            "Missing required Open Targets dataset(s): " + ", ".join(missing) +
            f"\nExpected local Parquet folders under:\n  {paths['OT_CACHE']}\n"
            "Pass --ot-cache if your downloaded files are elsewhere, or download them once with\n"
            "  envs/pipeline/bin/python Step00_Download_Resources.py --groups opentargets "
            "--output Step00_Download_OpenTargets.sh && sbatch Step00_Download_OpenTargets.sh"
        )
    if resolved["colocalisation"] is None and not args.allow_no_coloc:
        raise SystemExit(
            "Open Targets colocalisation dataset was not found.\n"
            "Either place it in the OT cache or rerun with --allow-no-coloc."
        )

    require_duckdb()
    banner("OPEN TARGETS DATASETS")
    for k, v in resolved.items():
        print(f"{k:18s}: {v if v else 'NOT FOUND'}")

    # Validate every resolved Parquet dataset shard-by-shard before analysis.
    # When enabled, automatically retry corrupt/truncated shards from the
    # official Open Targets release before benchmarking.
    dataset_complete: dict[str, bool] = {}
    dataset_has_healthy: dict[str, bool] = {}
    for k, v in resolved.items():
        if v is not None:
            good = healthy_parquet_files(v, allow_all_bad=True)
            dataset_has_healthy[k] = bool(good)
            dataset_complete[k] = parquet_dataset_complete(v)
        else:
            dataset_has_healthy[k] = False
            dataset_complete[k] = False

    incomplete_before = [k for k, ok in dataset_complete.items() if resolved.get(k) is not None and not ok]
    if args.repair_ot and incomplete_before:
        banner("AUTOMATIC OPEN TARGETS REPAIR / RETRY")
        if shutil.which("rsync") is None:
            print("[WARN] rsync is not available; automatic repair cannot run.")
        else:
            for k in incomplete_before:
                v = resolved[k]
                if v is None:
                    continue
                result = repair_bad_parquet_shards(
                    v, args.ot_release, args.repair_retries, args.rsync_base
                )
                print(
                    f"{k:18s}: "
                    f"{'REPAIRED' if result.get('SUCCESS') else 'STILL INCOMPLETE'} - "
                    f"{result.get('MESSAGE', '')}"
                )

        # Recompute health after repair attempts.
        for k, v in resolved.items():
            if v is not None:
                _clear_parquet_health(v)
                good = healthy_parquet_files(v, allow_all_bad=True)
                dataset_has_healthy[k] = bool(good)
                dataset_complete[k] = parquet_dataset_complete(v)

    # Required core datasets must retain at least one readable shard.
    for k in ["credible_set", "l2g_prediction"]:
        if not dataset_has_healthy.get(k, False):
            raise RuntimeError(
                f"Open Targets required dataset {k!r} has no readable Parquet shards after repair attempts."
            )

    print("\nOpen Targets Parquet health:")
    for k in sorted(dataset_complete):
        if resolved.get(k) is not None:
            print(f"  {k:18s}: {'COMPLETE' if dataset_complete[k] else 'INCOMPLETE - see QC audit'}")

    # Credible sets for exactly the GWAS studies represented by Step12.
    credible_subset, study_col = subset_parquet(
        resolved["credible_set"], ["studyId", "study_id"], studies
    )
    if credible_subset.empty:
        raise SystemExit(
            "No Open Targets credible sets were found for the Step12 GWAS study IDs.\n"
            "This may indicate study-ID mismatch or no OT coverage for these studies."
        )
    ot_loci, ot_variants = flatten_ot_credible_sets(credible_subset)
    if args.study:
        ot_loci = ot_loci[ot_loci["OT_STUDY_ID"].astype(str).eq(args.study)].copy()
        if not ot_variants.empty:
            keep_ids = set(ot_loci["OT_STUDY_LOCUS_ID"].astype(str))
            ot_variants = ot_variants[ot_variants["OT_STUDY_LOCUS_ID"].astype(str).isin(keep_ids)].copy()

    banner("LOCUS MATCHING")
    pairs_all = candidate_locus_pairs(gwas_loci, ot_loci, gwas_variants, ot_variants, args.locus_padding)
    matches = select_primary_matches(pairs_all)
    print(f"GWAS2m loci          : {len(gwas_loci):,}")
    print(f"OT credible sets     : {len(ot_loci):,}")
    print(f"Candidate pairs      : {len(pairs_all):,}")
    print(f"Primary matched pairs: {len(matches):,}")

    g_locus_cmp, ot_locus_cmp, _ = build_locus_comparison(gwas_loci, ot_loci, matches)
    study_comparability = build_study_comparability(gwas_loci, ot_loci, matches)
    complete_inventory = build_complete_locus_inventory(g_locus_cmp, ot_locus_cmp, matches)
    variant_cmp = compare_variants(matches, gwas_variants, ot_variants)
    variant_detail = build_variant_level_detail(matches, gwas_variants, ot_variants)
    missed_loci = build_missed_loci_diagnostic(
        g_locus_cmp, ot_locus_cmp, pairs_all, gwas_loci, ot_loci, args.locus_padding
    )

    # L2G only needs matched OT locus IDs.
    matched_otids = set(matches.get("OT_STUDY_LOCUS_ID", pd.Series(dtype=str)).dropna().astype(str))
    l2g_subset, _ = subset_parquet(
        resolved["l2g_prediction"], ["studyLocusId", "study_locus_id"], matched_otids
    )
    l2g = flatten_l2g(l2g_subset)

    # OT colocalisation + resolution of other QTL loci.
    coloc_resolved = pd.DataFrame()
    coloc_subset = pd.DataFrame()
    other_cred_subset = pd.DataFrame()
    study_subset = pd.DataFrame()
    if resolved["colocalisation"] is not None and matched_otids:
        coloc_subset, _ = subset_parquet(
            resolved["colocalisation"], ["studyLocusId", "study_locus_id"], matched_otids
        )
        coloc = flatten_colocalisation(coloc_subset)
        if not coloc.empty:
            other_ids = set(coloc["OTHER_STUDY_LOCUS_ID"].dropna().astype(str))
            if other_ids:
                other_cred_subset, _ = subset_parquet(
                    resolved["credible_set"], ["studyLocusId", "study_locus_id"], other_ids
                )
                other_meta, _ = flatten_ot_credible_sets(other_cred_subset)
                other_studies = set(other_meta.get("OT_STUDY_ID", pd.Series(dtype=str)).dropna().astype(str))
                smeta = pd.DataFrame()
                if (
                    resolved["study"] is not None
                    and dataset_has_healthy.get("study", False)
                    and other_studies
                ):
                    study_subset, _ = subset_parquet(
                        resolved["study"], ["studyId", "study_id", "id"], other_studies
                    )
                    smeta = study_metadata(study_subset)
                elif other_studies:
                    print(
                        "[WARN] Open Targets study metadata unavailable/unreadable; "
                        "continuing with qtlGeneId from credible_set where available."
                    )
                coloc_resolved = resolve_ot_coloc(coloc, other_meta, smeta)

    gene_cmp = compare_genes(
        genes, gwas_loci, matches, ot_loci, l2g, coloc_resolved,
        args.l2g_top_k, args.ot_coloc_h4, args.strong_h4,
        ot_coloc_complete=(
            dataset_complete.get("colocalisation", False)
            and dataset_has_healthy.get("study", False)
        ),
    )
    mechanism_cmp = build_mechanism_comparison(gene_cmp, args.strong_h4, args.suggestive_h4)
    # Use the richer mechanism table as the canonical gene comparison from here onward.
    gene_cmp = mechanism_cmp
    cell_context_cmp = build_cell_context_comparison(mechanism_cmp)
    signal_external_cmp = build_signal_external_comparison(step12_signals, mechanism_cmp)
    method_summary = method_summary_v2(variant_cmp, mechanism_cmp)

    tiers = summary_table(gene_cmp)
    locus_summary = locus_benchmark_summary(g_locus_cmp, ot_locus_cmp, matches, variant_cmp)
    basis_summary = match_basis_summary(matches)
    gene_summary = gene_benchmark_summary(gene_cmp, args.strong_h4)
    ot_coloc_fully_interpretable = (
        dataset_complete.get("colocalisation", False)
        and dataset_has_healthy.get("study", False)
    )
    performance = build_performance_metrics(
        g_locus_cmp, ot_locus_cmp, matches, variant_cmp, gene_cmp,
        args.strong_h4, ot_coloc_fully_interpretable,
    )
    strong_gene_benchmark = gene_cmp[
        safe_numeric(gene_cmp.get("BEST_H4", pd.Series(np.nan, index=gene_cmp.index))) >= args.strong_h4
    ].copy() if not gene_cmp.empty else pd.DataFrame()

    # Save raw OT subsets and comparisons for audit/reproducibility.
    tables = {
        "Step13_Study_Comparability.tsv": study_comparability,
        "Step13_Complete_Locus_Inventory.tsv": complete_inventory,
        "Step13_Locus_Pair_Candidates.tsv": pairs_all,
        "Step13_Locus_Matches.tsv": matches,
        "Step13_GWAS2m_Locus_Comparison.tsv": g_locus_cmp,
        "Step13_OpenTargets_Locus_Comparison.tsv": ot_locus_cmp,
        "Step13_Variant_Comparison.tsv": variant_cmp,
        "Step13_Variant_Level_Detail.tsv": variant_detail,
        "Step13_Missed_Loci_Why.tsv": missed_loci,
        "Step13_Performance_Metrics.tsv": performance,
        "Step13_Strong_Gene_Benchmark.tsv": strong_gene_benchmark,
        "Step13_Gene_Comparison.tsv": gene_cmp,
        "Step13_Molecular_Mechanism_Comparison.tsv": mechanism_cmp,
        "Step13_Cell_Context_Comparison.tsv": cell_context_cmp,
        "Step13_Signal_External_Comparison.tsv": signal_external_cmp,
        "Step13_Method_Summary.tsv": method_summary,
        "Step13_Novelty_Tier_Summary.tsv": tiers,
        "Step13_Locus_Benchmark_Summary.tsv": locus_summary,
        "Step13_Match_Basis_Summary.tsv": basis_summary,
        "Step13_Gene_Benchmark_Summary.tsv": gene_summary,
        "OpenTargets_credible_set_subset.tsv": credible_subset,
        "OpenTargets_l2g_subset.tsv": l2g,
        "OpenTargets_colocalisation_subset.tsv": coloc_subset,
        "OpenTargets_colocalisation_resolved.tsv": coloc_resolved,
        "Step13_OpenTargets_Parquet_Health.tsv": pd.DataFrame(_PARQUET_HEALTH_ROWS),
        "Step13_OpenTargets_Repair_Audit.tsv": pd.DataFrame(_OT_REPAIR_ROWS),
    }
    for name, df in tables.items():
        df.to_csv(paths["OUT"] / name, sep="\t", index=False)

    # Human-friendly CSV copies of the benchmark tables most likely to be
    # inspected outside Python/R.
    save_tsv_and_csv(performance, paths["OUT"], "Step13_Performance_Metrics")
    save_tsv_and_csv(variant_cmp, paths["OUT"], "Step13_Matched_Locus_Detail")
    save_tsv_and_csv(variant_detail, paths["OUT"], "Step13_Variant_Level_Detail")
    save_tsv_and_csv(missed_loci, paths["OUT"], "Step13_Missed_Loci_Why")
    save_tsv_and_csv(strong_gene_benchmark, paths["OUT"], "Step13_Strong_Gene_Benchmark")

    # High-value interpretation subsets.
    gene_cmp[gene_cmp["NOVELTY_TIER_STEP13"].eq("A_ESTABLISHED_OT_GENE_AND_MOLQTL")].to_csv(
        paths["OUT"] / "Step13_A_Recapitulated_Known_Biology.tsv", sep="\t", index=False
    )
    gene_cmp[gene_cmp["NOVELTY_TIER_STEP13"].eq("B_OT_SUPPORTED_GENE_POTENTIAL_NEW_MECHANISM")].to_csv(
        paths["OUT"] / "Step13_B_Potential_New_Mechanisms.tsv", sep="\t", index=False
    )
    gene_cmp[gene_cmp["NOVELTY_TIER_STEP13"].eq("C_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE")].to_csv(
        paths["OUT"] / "Step13_C_New_Effector_Candidates.tsv", sep="\t", index=False
    )
    gene_cmp[gene_cmp["NOVELTY_TIER_STEP13"].eq("D_POTENTIAL_NEW_LOCUS_CANDIDATE_REQUIRES_EXTERNAL_VALIDATION")].to_csv(
        paths["OUT"] / "Step13_D_Potential_New_Locus_Candidates.tsv", sep="\t", index=False
    )
    gene_cmp[gene_cmp["NOVELTY_TIER_STEP13"].eq("COVERAGE_OR_MAPPING_UNRESOLVED")].to_csv(
        paths["OUT"] / "Step13_Unresolved_Or_Coverage_Gaps.tsv", sep="\t", index=False
    )

    # Mechanism-aware screening subsets. These are candidates for follow-up,
    # not definitive novelty claims.
    mechanism_cmp[mechanism_cmp["EXTERNAL_MECHANISM_CLASS"].astype(str).str.contains("NEW_CELL_CONTEXT", na=False)].to_csv(
        paths["OUT"] / "Step13_Potential_New_Cell_Contexts.tsv", sep="\t", index=False
    )
    mechanism_cmp[mechanism_cmp["EXTERNAL_MECHANISM_CLASS"].astype(str).str.contains("NEW_MOLECULAR_MECHANISM", na=False)].to_csv(
        paths["OUT"] / "Step13_Potential_New_Molecular_Mechanisms.tsv", sep="\t", index=False
    )
    mechanism_cmp[mechanism_cmp["EXTERNAL_MECHANISM_CLASS"].astype(str).str.startswith("CONVERGENT", na=False)].to_csv(
        paths["OUT"] / "Step13_Externally_Convergent_Mechanisms.tsv", sep="\t", index=False
    )

    metadata = {
        "step": 13,
        "version": VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "phenotype": phenotype,
        "ancestry_code": ancestry_code,
        "ancestry_label": ancestry_label,
        "study_filter": args.study or None,
        "open_targets_release_label": args.ot_release,
        "open_targets_cache": str(paths["OT_CACHE"]),
        "resolved_datasets": {k: str(v) if v else None for k, v in resolved.items()},
        "n_step12_gene_rows": len(genes),
        "n_step12_signal_rows": len(step12_signals),
        "n_step12_molecular_rows": len(step12_molecular),
        "n_step12_strong_rows": len(step12_strong),
        "n_step12_pathway_genes": len(step12_pathway),
        "n_gwas2m_loci": len(gwas_loci),
        "n_ot_credible_sets": len(ot_loci),
        "n_matched_loci": len(matches),
        "n_gwas2m_only_loci": max(0, len(gwas_loci) - len(matches)),
        "n_opentargets_only_credible_sets": max(0, len(ot_loci) - len(matches)),
        "gwas2m_locus_status_counts": g_locus_cmp["COMPARISON_STATUS"].value_counts(dropna=False).to_dict() if not g_locus_cmp.empty else {},
        "opentargets_locus_status_counts": ot_locus_cmp["COMPARISON_STATUS"].value_counts(dropna=False).to_dict() if not ot_locus_cmp.empty else {},
        "match_basis_counts": matches["MATCH_BASIS"].value_counts(dropna=False).to_dict() if not matches.empty and "MATCH_BASIS" in matches.columns else {},
        "n_l2g_rows": len(l2g),
        "n_ot_coloc_rows": len(coloc_resolved),
        "novelty_tier_counts": gene_cmp["NOVELTY_TIER_STEP13"].value_counts(dropna=False).to_dict() if not gene_cmp.empty else {},
        "mechanism_class_counts": gene_cmp["EXTERNAL_MECHANISM_CLASS"].value_counts(dropna=False).to_dict() if not gene_cmp.empty and "EXTERNAL_MECHANISM_CLASS" in gene_cmp.columns else {},
        "cell_context_screen_counts": cell_context_cmp["CELL_CONTEXT_SCREEN"].value_counts(dropna=False).to_dict() if not cell_context_cmp.empty and "CELL_CONTEXT_SCREEN" in cell_context_cmp.columns else {},
        "open_targets_dataset_complete": dataset_complete,
        "open_targets_dataset_has_healthy_shards": dataset_has_healthy,
        "n_corrupt_or_truncated_parquet_shards": int(sum(1 for r in _PARQUET_HEALTH_ROWS if r.get("STATUS") == "CORRUPT_OR_TRUNCATED")),
        "automatic_repair_enabled": bool(args.repair_ot),
        "n_repair_audit_rows": len(_OT_REPAIR_ROWS),
        "important_note": "Open Targets is an external benchmark, not ground truth. Tier D requires GWAS Catalog/literature/independent validation.",
    }
    (paths["OUT"] / "Step13_summary.json").write_text(
        json.dumps(metadata, indent=2, default=json_safe) + "\n", encoding="utf-8"
    )

    banner("STEP13 COMPLETE - LOCUS BENCHMARK")
    n_g_only = max(0, len(gwas_loci) - len(matches))
    n_o_only = max(0, len(ot_loci) - len(matches))
    print(f"GWAS2m total loci                    : {len(gwas_loci):,}")
    print(f"Open Targets total credible sets     : {len(ot_loci):,}")
    print(f"COMMON / primary matched loci        : {len(matches):,}")
    print(f"GWAS2m ONLY / no primary OT match    : {n_g_only:,}")
    print(f"Open Targets ONLY / missed by GWAS2m : {n_o_only:,}")
    print()
    if not g_locus_cmp.empty:
        print("GWAS2m LOCUS STATUS")
        print(g_locus_cmp["COMPARISON_STATUS"].value_counts(dropna=False).rename_axis("STATUS").reset_index(name="N_LOCI").to_string(index=False))
    print()
    if not ot_locus_cmp.empty:
        print("OPEN TARGETS LOCUS STATUS")
        print(ot_locus_cmp["COMPARISON_STATUS"].value_counts(dropna=False).rename_axis("STATUS").reset_index(name="N_CREDIBLE_SETS").to_string(index=False))
    print()
    if not basis_summary.empty:
        print("HOW THE COMMON LOCI MATCHED")
        print(basis_summary.to_string(index=False))
    print()
    if not variant_cmp.empty:
        exact_top = int(variant_cmp.get("TOP_VARIANT_AGREEMENT", pd.Series(False, index=variant_cmp.index)).map(boolish).sum())
        cs_shared = int((safe_numeric(variant_cmp.get("N_SHARED_CS95", pd.Series(0, index=variant_cmp.index))).fillna(0) > 0).sum())
        print(f"Same highest-PIP/top variant         : {exact_top:,} / {len(matches):,}")
        print(f"95% credible sets share variant(s)   : {cs_shared:,} / {len(matches):,}")

    banner("STEP13 BENCHMARK PERFORMANCE - OPEN TARGETS AS EXTERNAL REFERENCE")
    print("NOTE: these are concordance metrics against Open Targets, NOT absolute biological accuracy.")
    if not performance.empty:
        display_cols = ["LEVEL", "METRIC", "VALUE", "NUMERATOR", "DENOMINATOR", "PERCENT", "INTERPRETATION"]
        with pd.option_context("display.max_columns", None, "display.width", 220, "display.max_colwidth", 70):
            print(performance[display_cols].to_string(index=False))

    if not variant_cmp.empty:
        banner("MATCHED LOCUS DETAIL - VARIANT / CREDIBLE-SET / PIP AGREEMENT")
        cols = [c for c in [
            "STUDY_ID", "GWAS2M_LOCUS_ID", "OT_STUDY_LOCUS_ID", "MATCH_BASIS",
            "GWAS2M_TOP_VARIANT", "OT_TOP_VARIANT", "TOP_VARIANT_AGREEMENT",
            "N_SHARED_VARIANTS", "VARIANT_JACCARD", "N_SHARED_CS95", "CS95_JACCARD",
            "PIP_SPEARMAN", "PIP_PEARSON", "PIP_MAE",
        ] if c in variant_cmp.columns]
        with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 260):
            print(variant_cmp[cols].to_string(index=False))

    if not missed_loci.empty:
        banner("UNMATCHED LOCI - WHAT WAS MISSED AND WHY")
        cols = [c for c in [
            "SOURCE", "COMPARISON_STATUS", "STUDY_ID", "GWAS2M_LOCUS_ID",
            "OT_STUDY_LOCUS_ID", "CHR", "LEAD_VARIANT",
            "HAS_NONPRIMARY_CANDIDATE_PAIR", "N_NONPRIMARY_CANDIDATE_PAIRS",
            "NEAREST_OTHER_LOCUS_ID", "NEAREST_GAP_BP", "WHY_NOT_PRIMARY_MATCH",
        ] if c in missed_loci.columns]
        with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 260, "display.max_colwidth", 80):
            print(missed_loci[cols].to_string(index=False))

    banner("STEP13 COMPLETE - GENE / BIOLOGICAL BENCHMARK")
    print(f"Step12 locus-gene rows               : {len(gene_cmp):,}")
    if not gene_summary.empty:
        print(gene_summary.to_string(index=False))
    if not tiers.empty:
        print("\nNOVELTY / BENCHMARK TIERS")
        print(tiers.to_string(index=False))

    if not strong_gene_benchmark.empty:
        banner(f"STRONG STEP12 GENE BENCHMARK - H4 >= {args.strong_h4:g}")
        cols = [c for c in [
            "STUDY_ACCESSION", "LOCUS_ID", "EFFECTOR_GENE_ID", "EFFECTOR_GENE_SYMBOL",
            "BEST_H4", "BEST_EVIDENCE_SOURCE", "BEST_QTL_TYPE", "BEST_CONTEXT",
            "BEST_CONDITION", "BEST_GWAS_VARIANT", "BEST_QTL_VARIANT",
            "BULK_BEST_H4", "BULK_BEST_CONTEXT", "CELL_BEST_H4",
            "CELL_BEST_CONTEXT", "CELL_BEST_CONDITION", "OT_LOCUS_MATCH",
            "OT_STUDY_LOCUS_ID", "OT_L2G_GENE_PRESENT", "OT_L2G_RANK",
            "OT_L2G_SCORE", "OT_L2G_TOP1_GENE",
            "OT_COLOC_GENE_PRESENT_H4_THRESHOLD", "OT_COLOC_MAX_H4_SAME_GENE",
            "OPEN_TARGETS_EVIDENCE_CLASS", "EXTERNAL_MECHANISM_CLASS",
            "CELL_CONTEXT_SCREEN",
        ] if c in strong_gene_benchmark.columns]
        with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 280, "display.max_colwidth", 80):
            print(strong_gene_benchmark[cols].to_string(index=False))

    if not args.compact_screen:
        if not method_summary.empty:
            banner("METHOD-AWARE OPEN TARGETS BENCHMARK SUMMARY")
            with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 320, "display.max_colwidth", 100):
                print(method_summary.to_string(index=False))

        if not mechanism_cmp.empty:
            banner("TOP 50 INTEGRATED GWAS2m <-> OPEN TARGETS MECHANISMS")
            topm = mechanism_cmp.sort_values("BEST_H4", ascending=False, na_position="last").head(50)
            cols = [c for c in [
                "STUDY_ACCESSION", "LOCUS_ID", "EFFECTOR_GENE_ID", "EFFECTOR_GENE_SYMBOL",
                "BEST_H4", "BEST_EVIDENCE_SOURCE", "BEST_CONTEXT", "BEST_CONDITION",
                "BEST_GWAS_VARIANT", "BEST_QTL_VARIANT",
                "BULK_BEST_H4", "BULK_BEST_CONTEXT",
                "CELL_BEST_H4", "CELL_BEST_CONTEXT", "CELL_BEST_CONDITION",
                "OT_LOCUS_MATCH", "OT_METHOD", "OT_L2G_GENE_PRESENT", "OT_L2G_RANK",
                "OT_L2G_SCORE", "OT_COLOC_GENE_PRESENT_H4_THRESHOLD",
                "OT_COLOC_MAX_H4_SAME_GENE", "EXTERNAL_MECHANISM_CLASS",
                "NOVELTY_TIER_STEP13",
            ] if c in topm.columns]
            with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 420, "display.max_colwidth", 90):
                print(topm[cols].to_string(index=False))

        if not cell_context_cmp.empty:
            banner("TOP 40 CELL-TYPE MECHANISMS VS OPEN TARGETS")
            z = cell_context_cmp.head(40)
            cols = [c for c in [
                "STUDY_ACCESSION", "LOCUS_ID", "EFFECTOR_GENE_ID", "EFFECTOR_GENE_SYMBOL",
                "CELL_BEST_H4", "CELL_BEST_CONTEXT", "CELL_BEST_CONDITION",
                "CELL_BEST_GWAS_VARIANT", "CELL_BEST_QTL_VARIANT",
                "OT_L2G_GENE_PRESENT", "OT_COLOC_GENE_PRESENT_H4_THRESHOLD",
                "OT_COLOC_MAX_H4_SAME_GENE", "OT_COLOC_BIOSAMPLES", "OT_COLOC_TRAITS",
                "OT_CONTEXT_TEXT_SCREEN", "CELL_CONTEXT_SCREEN", "EXTERNAL_MECHANISM_CLASS",
            ] if c in z.columns]
            with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 420, "display.max_colwidth", 100):
                print(z[cols].to_string(index=False))

    if not dataset_complete.get("colocalisation", False) or not dataset_has_healthy.get("study", False):
        print("\nIMPORTANT: Open Targets colocalisation/study data are incomplete.")
        print("Tier A/B/C absence-based interpretation is intentionally suppressed.")
        print("Repair the corrupt OT shards, rerun Step13, and only then interpret gene/mechanism novelty.")

    print(f"\nPerformance metrics:\n  {paths['OUT'] / 'Step13_Performance_Metrics.csv'}")
    print(f"Matched locus detail:\n  {paths['OUT'] / 'Step13_Matched_Locus_Detail.csv'}")
    print(f"Variant-level detail:\n  {paths['OUT'] / 'Step13_Variant_Level_Detail.csv'}")
    print(f"Missed loci + reasons:\n  {paths['OUT'] / 'Step13_Missed_Loci_Why.csv'}")
    print(f"Strong gene benchmark:\n  {paths['OUT'] / 'Step13_Strong_Gene_Benchmark.csv'}")
    print(f"OT Parquet QC audit:\n  {paths['OUT'] / 'Step13_OpenTargets_Parquet_Health.tsv'}")
    print(f"OT repair audit:\n  {paths['OUT'] / 'Step13_OpenTargets_Repair_Audit.tsv'}")
    print(f"Locus summary:\n  {paths['OUT'] / 'Step13_Locus_Benchmark_Summary.tsv'}")
    print(f"Match basis summary:\n  {paths['OUT'] / 'Step13_Match_Basis_Summary.tsv'}")
    print(f"Gene summary:\n  {paths['OUT'] / 'Step13_Gene_Benchmark_Summary.tsv'}")
    print(f"Main gene output:\n  {paths['OUT'] / 'Step13_Gene_Comparison.tsv'}")
    print(f"Mechanism comparison:\n  {paths['OUT'] / 'Step13_Molecular_Mechanism_Comparison.tsv'}")
    print(f"Cell-context comparison:\n  {paths['OUT'] / 'Step13_Cell_Context_Comparison.tsv'}")
    print(f"Signal-level external comparison:\n  {paths['OUT'] / 'Step13_Signal_External_Comparison.tsv'}")
    print(f"Study comparability:\n  {paths['OUT'] / 'Step13_Study_Comparability.tsv'}")
    print(f"Complete locus inventory:\n  {paths['OUT'] / 'Step13_Complete_Locus_Inventory.tsv'}")
    print(f"New-effector candidates:\n  {paths['OUT'] / 'Step13_C_New_Effector_Candidates.tsv'}")
    print(f"Potential new-locus candidates:\n  {paths['OUT'] / 'Step13_D_Potential_New_Locus_Candidates.tsv'}")
    print("\nNext: review OT coverage/QC, then inspect convergent mechanisms, potential new cell contexts, new effectors, and only after that proceed to pathway analysis.")


if __name__ == "__main__":
    pipeline_main(main)
