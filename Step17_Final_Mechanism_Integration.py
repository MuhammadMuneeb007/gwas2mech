#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================================
GWAS2Mech STEP17 v2.0
PATHWAY-FIRST MECHANISM INTEGRATION / EVIDENCE TIERS / AUDITABLE TSV OUTPUT
================================================================================

GWAS2Mech PIPELINE USAGE (same flags in Steps 12-17)
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR                       # run on this machine
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR --slurm               # or write a one-job SLURM script
    sbatch Step17_Mechanism_parkinson_s_disease_european.sh
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR --stats
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR --index 2             # study 2 only -> <ancestry>/<GCST>/
    # 1000 Genomes 30x LD panel: every input / output folder is <ancestry>_1000G30x
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --slurm
    sbatch Step17_Mechanism_parkinson_s_disease_european_1000G30x.sh
    envs/pipeline/bin/python Step17_Final_Mechanism_Integration.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

PURPOSE
-------
Step17 is the final integration layer of GWAS2m.

It combines the evidence generated in Steps 12-16:

    Step12  molecular-QTL / colocalisation master evidence
    Step13  Open Targets comparison
    Step14  deduplicated genes + pathway enrichment
    Step15  PPI / network evidence
    Step16  directional experimental perturbation evidence (nonconfirmatory unless replicated)

The goal is NOT to invent a new causal probability.

The goal is to produce a transparent final table answering:

    * Which genes have the strongest genetic/molecular support?
    * Which are independently supported by Open Targets?
    * Which are GWAS2m alternative effectors at known loci?
    * Which occur where Open Targets has a coverage gap?
    * Which genes drive significant disease-relevant pathways?
    * Which genes physically connect to known/OT-supported disease genes?
    * Which targets are tractable or clinically precedented?
    * Is therapeutic direction known, conflicting, or unavailable?
    * What evidence tier should each candidate receive?
    * What is the evidence basis for the final prioritisation?

SCIENTIFIC GUARDRAILS
---------------------
1. H4 is evidence for a shared signal, NOT proof that a gene is causal.
2. Open Targets is an external comparator, NOT a biological gold standard.
3. Absence from Open Targets does NOT automatically mean novelty.
4. Pathway enrichment is supportive mechanism context, NOT causal proof.
5. STRING network proximity is supportive mechanism context, NOT causal proof.
6. Added STRING connector proteins are NOT promoted to GWAS2m genetic targets.
7. ChEMBL drug precedence for another disease is NOT proof of efficacy for the
   phenotype being analysed.
8. Therapeutic direction remains UNKNOWN unless signed, allele-harmonised
   molecular evidence supports it.
9. The final Step17 score is a transparent heuristic evidence-ranking score,
   NOT a probability of causality, clinical success, or treatment efficacy.
10. The evidence tier is rule-based and should be interpreted together with the
    underlying columns, never in isolation.

PRIMARY OUTPUT
--------------
17_mechanism/<phenotype>/<ancestry>/

    Step17_Final_Target_Ranking.tsv
    Step17_Final_Evidence_Matrix.tsv
    Step17_Evidence_Tiers.tsv

    Step17_Pathway_Mechanism_Map.tsv
    Step17_Network_Mechanism_Map.tsv
    Step17_Drug_Mechanism_Map.tsv
    Step17_External_Concordance.tsv

    Step17_High_Priority_Candidates.tsv
    Step17_GWAS2m_Alternative_Effectors.tsv
    Step17_OpenTargets_Convergent_Targets.tsv
    Step17_OT_Coverage_Gap_Targets.tsv

    Step17_QC_Audit.tsv
    Step17_summary.json

RUN
---
python Step17_Final_Mechanism_Integration.py \
  --phenotype "parkinson's disease" \
  --ancestry EUR

Include suggestive candidates (H4 >= 0.50) in addition to CORE_STRONG:
python Step17_Final_Mechanism_Integration.py \
  --phenotype "parkinson's disease" \
  --ancestry EUR \
  --include-suggestive

Compact terminal output:
python Step17_Final_Mechanism_Integration.py \
  --phenotype "parkinson's disease" \
  --ancestry EUR \
  --compact-screen

DEPENDENCIES
------------
python -m pip install pandas numpy
================================================================================
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd



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
    "NAME": "Step17", "TITLE": "final disease-mechanism integration", "FOLDER": "17_mechanism",
    "SCRIPT_PREFIX": "Step17_Mechanism", "SUMMARY": "Step17_summary.json",
    "PREVIOUS": "Steps 12-16",
    "INPUTS": [
        ("Step12", "12_master_table", ["Step12_Locus_Gene_Master.tsv"], True, True),
        ("Step13", "13_opentargets_comparison", ["Step13_Gene_Comparison.tsv"], False, True),
        ("Step14", "14_pathways", ["Step14_Gene_Ranking.tsv", "Step14_Pathway_Significant.tsv"], True, True),
        ("Step15", "15_network", ["Step15_Gene_Pathway_Handoff.tsv", "Step15_Gene_Mapping_Audit.tsv"], True, True),
        ("Step16", "16_perturbation", ["Step16_summary.json", "Step16_Directional_Pathway_Responses_All.tsv"], True, True),
    ],
    "MODULES": [("pandas", True, "tables")],
    "STATS_TABLES": [
        ("Step17_Pathway_Master.tsv", ["PATHWAY_ID", "PATHWAY_NAME", "SOURCE", "STEP14_BEST_FDR", "N_PRIORITISED_DRIVER_GENES", "N_PERTURBED_DRIVER_GENES", "HIGHEST_MECHANISM_TIER"], "", "PATHWAY-FIRST EVIDENCE"),
        ("Step17_Mechanism_Master.tsv", ["MECHANISM_ID", "PATHWAY_NAME", "GENE_SYMBOL", "DATASET", "CELL_TYPE", "PATHWAY_NES", "GENE_SET_PERMUTATION_FDR", "MECHANISM_TIER", "DIRECTION_ALIGNMENT", "MISSING_EVIDENCE"], "", "GENE-PATHWAY-CONTEXT EVIDENCE"),
        ("Step17_Missing_Evidence.tsv", ["MECHANISM_ID", "PATHWAY_NAME", "GENE_SYMBOL", "MISSING_EVIDENCE"], "", "EVIDENCE GAPS"),
    ],
}


def frame_extra_checks(ns, P) -> list:
    checks=[]
    for name, (folder, filename) in M17_FILES.items():
        if name not in M17_CORE_INPUTS:
            continue
        f=P["root"]/folder/P["p"]/results_slug(P["a"])/filename
        checks.append((name+" "+filename, f.is_file() and f.stat().st_size>10 if f.is_file() else False, str(f), True))
    summary=P["root"]/"16_perturbation"/P["p"]/results_slug(P["a"])/"Step16_summary.json"
    good=False
    if summary.is_file():
        try: good=str(json.loads(summary.read_text()).get("version","")).startswith("3.")
        except Exception: pass
    checks.append(("Completed Step16 v3 summary",good,str(summary),True))
    return checks


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


VERSION = "3.0.0-step16v3-source-id-mechanisms"

ANCESTRY = {
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

def norm(value) -> str:
    x = re.sub(r"[^a-z0-9]+", " ", str(value).lower())
    return re.sub(r"\s+", " ", x).strip()


def slug(value) -> str:
    return norm(value).replace(" ", "_")


def ancestry_info(value):
    key = norm(value)

    if key not in ANCESTRY:
        raise SystemExit(
            f"Unsupported ancestry {value!r}; use EUR/AFR/EAS/SAS/AMR."
        )

    return ANCESTRY[key]


def section(title, char="="):
    print()
    print(char * 158)
    print(title)
    print(char * 158)


def show(df, title=None, max_rows=None, columns=None):
    if title:
        section(title)

    if df is None or df.empty:
        print("(no rows)")
        return

    x = df.copy()

    if columns:
        cols = [c for c in columns if c in x.columns]
        x = x[cols]

    total = len(x)

    if max_rows is not None:
        x = x.head(max_rows)

    with pd.option_context(
        "display.max_rows", None,
        "display.max_columns", None,
        "display.width", 1400,
        "display.max_colwidth", 180,
        "display.expand_frame_repr", False,
        "display.float_format", lambda v: f"{v:.6g}",
    ):
        print(x.to_string(index=False))

    if max_rows is not None and total > len(x):
        print(f"\n[showing {len(x):,}/{total:,} rows]")


def clean_gene_id(value):
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass

    x = str(value).strip()

    if x.lower() in {"", "nan", "none", "null", "na", ".", "-"}:
        return ""

    m = re.search(r"(ENSG\d+)", x, flags=re.I)

    return m.group(1).upper() if m else x


def clean_symbol(value):
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass

    x = str(value).strip()

    if x.lower() in {"", "nan", "none", "null", "na", ".", "-"}:
        return ""

    return x


def safe_numeric(values):
    return pd.to_numeric(values, errors="coerce")


def safe_float(value, default=np.nan):
    try:
        if pd.isna(value):
            return default
    except Exception:
        pass

    try:
        return float(value)
    except Exception:
        return default


def boolish(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)

    if value is None:
        return False

    try:
        if pd.isna(value):
            return False
    except Exception:
        pass

    return str(value).strip().lower() in {
        "1", "true", "t", "yes", "y",
    }


def first_existing(columns: Iterable[str], candidates: Iterable[str]):
    columns = list(columns)

    lower = {
        str(c).lower(): c
        for c in columns
    }

    for candidate in candidates:
        if candidate in columns:
            return candidate

        found = lower.get(
            str(candidate).lower()
        )

        if found is not None:
            return found

    return None


def split_tokens(value):
    if value is None:
        return []

    try:
        if pd.isna(value):
            return []
    except Exception:
        pass

    return [
        x.strip()
        for x in re.split(r"[;,|]", str(value))
        if x.strip()
        and x.strip().lower() not in {
            "nan", "none", "null", "na",
        }
    ]


def join_unique(values):
    out = []

    for value in values:
        for token in split_tokens(value):
            if token not in out:
                out.append(token)

    return ";".join(out)


def dedupe_keep_order(values):
    seen = set()
    out = []

    for value in values:
        if value in seen:
            continue

        seen.add(value)
        out.append(value)

    return out


def percentile(values):
    x = safe_numeric(pd.Series(values))

    if x.notna().sum() <= 1:
        return pd.Series(
            0.0,
            index=x.index,
        )

    return (
        x.rank(
            pct=True,
            method="average",
        )
        .fillna(0.0)
    )


def read_tsv(path, required=False):
    path = Path(path)

    if not path.exists():
        if required:
            raise SystemExit(
                f"Required input missing:\n  {path}"
            )

        return pd.DataFrame()

    try:
        return pd.read_csv(
            path,
            sep="\t",
            low_memory=False,
        )

    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def ensure_gene_id(df, candidates):
    if df.empty:
        return df, None

    col = first_existing(
        df.columns,
        candidates,
    )

    if col is None:
        return df, None

    df = df.copy()

    df[col] = df[col].map(
        clean_gene_id
    )

    return df, col


# =============================================================================
# INPUT PATHS / LOADERS
# =============================================================================

def load_all_inputs(
    root,
    phenotype_slug,
    ancestry_slug,
):
    step12_dir = (
        root
        /
        "12_master_table"
        /
        phenotype_slug
        /
        ancestry_slug
    )

    step13_dir = (
        root
        /
        "13_opentargets_comparison"
        /
        phenotype_slug
        /
        ancestry_slug
    )

    step14_dir = (
        root
        /
        "14_pathways"
        /
        phenotype_slug
        /
        ancestry_slug
    )

    step15_dir = (
        root
        /
        "15_network"
        /
        phenotype_slug
        /
        ancestry_slug
    )

    step16_dir = (
        root
        /
        "16_drug_targets"
        /
        phenotype_slug
        /
        ancestry_slug
    )

    data = {
        "step12_gene_master":
            read_tsv(
                step12_dir
                /
                "Step12_Locus_Gene_Master.tsv",
                required=True,
            ),

        "step13_gene":
            read_tsv(
                step13_dir
                /
                "Step13_Gene_Comparison.tsv",
                required=False,
            ),

        "step13_strong":
            read_tsv(
                step13_dir
                /
                "Step13_Strong_Gene_Benchmark.tsv",
                required=False,
            ),

        "step14_ranking":
            read_tsv(
                step14_dir
                /
                "Step14_Gene_Ranking.tsv",
                required=True,
            ),

        "step14_membership":
            read_tsv(
                step14_dir
                /
                "Step14_Gene_Set_Membership.tsv",
                required=True,
            ),

        "step13_pathways":
            read_tsv(
                step14_dir
                /
                "Step14_Pathway_Significant.tsv",
                required=False,
            ),

        "step14_drivers":
            read_tsv(
                step14_dir
                /
                "Step14_Pathway_Gene_Drivers.tsv",
                required=False,
            ),

        "step15_centrality":
            read_tsv(
                step15_dir
                /
                "Step15_Network_Centrality.tsv",
                required=False,
            ),

        "step15_direct":
            read_tsv(
                step15_dir
                /
                "Step15_Direct_Physical_Edges.tsv",
                required=False,
            ),

        "step15_bridges":
            read_tsv(
                step15_dir
                /
                "Step15_OT_Bridge_Summary.tsv",
                required=False,
            ),

        "step15_connectors":
            read_tsv(
                step15_dir
                /
                "Step15_Connector_Genes.tsv",
                required=False,
            ),

        "step16_targets":
            read_tsv(
                step16_dir
                /
                "Step16_Target_Prioritisation.tsv",
                required=True,
            ),

        "step16_mechanisms":
            read_tsv(
                step16_dir
                /
                "Step16_Drug_Mechanisms.tsv",
                required=False,
            ),

        "step16_direction":
            read_tsv(
                step16_dir
                /
                "Step16_Direction_Summary.tsv",
                required=False,
            ),

        "step16_phenotype_drugs":
            read_tsv(
                step16_dir
                /
                "Step16_Phenotype_Matched_Drugs.tsv",
                required=False,
            ),
    }

    data["dirs"] = {
        "step12": step12_dir,
        "step13": step13_dir,
        "step14": step14_dir,
        "step15": step15_dir,
        "step16": step16_dir,
    }

    return data


def gene_sets_from_membership(membership):
    result = {}

    if membership.empty:
        return result

    if not {
        "GENE_SET",
        "GENE_ID",
    }.issubset(
        membership.columns
    ):
        return result

    for name, group in membership.groupby(
        "GENE_SET",
        sort=False,
    ):
        result[str(name)] = dedupe_keep_order([
            clean_gene_id(value)
            for value
            in group[
                "GENE_ID"
            ]
            if clean_gene_id(value)
        ])

    return result


# =============================================================================
# PATHWAY EVIDENCE
# =============================================================================

def make_pathway_evidence(
    pathways,
    drivers,
):
    if pathways.empty:
        return pd.DataFrame(), pd.DataFrame()

    significant = pathways.copy()

    if "SIGNIFICANT" in significant.columns:
        significant = significant[
            significant[
                "SIGNIFICANT"
            ]
            .map(boolish)
        ].copy()

    if significant.empty:
        return pd.DataFrame(), pd.DataFrame()

    primary_sets = {
        "CORE_STRONG",
        "CORE_STRONG_PRIOR_ROBUST",
        "CORE_STRONG_MULTI_QTL",
        "LOCUS_TOP_STRONG",
    }

    significant["IS_PRIMARY_GWAS2M_PATHWAY_SET"] = (
        significant[
            "GENE_SET"
        ]
        .astype(str)
        .isin(primary_sets)
    )

    if drivers.empty:
        rows = []

        for row in significant.itertuples(
            index=False
        ):
            for gid in split_tokens(
                getattr(
                    row,
                    "INTERSECTION_GENES",
                    "",
                )
            ):
                rows.append({
                    "GENE_SET":
                        getattr(
                            row,
                            "GENE_SET",
                            "",
                        ),

                    "SOURCE":
                        getattr(
                            row,
                            "SOURCE",
                            "",
                        ),

                    "TERM_ID":
                        getattr(
                            row,
                            "TERM_ID",
                            "",
                        ),

                    "TERM_NAME":
                        getattr(
                            row,
                            "TERM_NAME",
                            "",
                        ),

                    "FDR":
                        safe_float(
                            getattr(
                                row,
                                "FDR",
                                np.nan,
                            )
                        ),

                    "GENE_ID":
                        clean_gene_id(
                            gid
                        ),

                    "IS_PRIMARY_GWAS2M_PATHWAY_SET":
                        getattr(
                            row,
                            "IS_PRIMARY_GWAS2M_PATHWAY_SET",
                            False,
                        ),
                })

        driver_table = pd.DataFrame(
            rows
        )

    else:
        driver_table = drivers.copy()

        driver_table, gid_col = ensure_gene_id(
            driver_table,
            [
                "GENE_ID",
                "EFFECTOR_GENE_ID",
            ],
        )

        if gid_col and gid_col != "GENE_ID":
            driver_table = driver_table.rename(
                columns={
                    gid_col:
                        "GENE_ID",
                }
            )

        if "GENE_SET" in driver_table.columns:
            driver_table[
                "IS_PRIMARY_GWAS2M_PATHWAY_SET"
            ] = (
                driver_table[
                    "GENE_SET"
                ]
                .astype(str)
                .isin(
                    primary_sets
                )
            )

    if driver_table.empty:
        return significant, pd.DataFrame()

    rows = []

    for gid, group in driver_table.groupby(
        "GENE_ID",
        sort=False,
    ):
        primary = group[
            group[
                "IS_PRIMARY_GWAS2M_PATHWAY_SET"
            ]
            .map(boolish)
        ].copy()

        rows.append({
            "GENE_ID":
                clean_gene_id(
                    gid
                ),

            "N_SIGNIFICANT_PATHWAYS_ALL_SETS":
                group[
                    [
                        "SOURCE",
                        "TERM_NAME",
                    ]
                ]
                .drop_duplicates()
                .shape[0],

            "N_SIGNIFICANT_PATHWAYS_PRIMARY":
                primary[
                    [
                        "SOURCE",
                        "TERM_NAME",
                    ]
                ]
                .drop_duplicates()
                .shape[0]
                if not primary.empty
                else
                0,

            "N_PATHWAY_DATABASES_PRIMARY":
                primary[
                    "SOURCE"
                ]
                .astype(str)
                .nunique()
                if (
                    not primary.empty
                    and
                    "SOURCE"
                    in primary.columns
                )
                else
                0,

            "BEST_PATHWAY_FDR":
                safe_numeric(
                    primary[
                        "FDR"
                    ]
                ).min()
                if (
                    not primary.empty
                    and
                    "FDR"
                    in primary.columns
                    and
                    safe_numeric(
                        primary[
                            "FDR"
                        ]
                    )
                    .notna()
                    .any()
                )
                else
                np.nan,

            "PRIMARY_PATHWAYS":
                join_unique(
                    primary[
                        "TERM_NAME"
                    ]
                )
                if (
                    not primary.empty
                    and
                    "TERM_NAME"
                    in primary.columns
                )
                else
                "",

            "PRIMARY_PATHWAY_SOURCES":
                join_unique(
                    primary[
                        "SOURCE"
                    ]
                )
                if (
                    not primary.empty
                    and
                    "SOURCE"
                    in primary.columns
                )
                else
                "",

            "IS_PRIMARY_PATHWAY_DRIVER":
                (
                    "YES"
                    if len(primary)
                    else
                    "NO"
                ),
        })

    return (
        significant,
        pd.DataFrame(
            rows
        ),
    )


# =============================================================================
# NETWORK EVIDENCE
# =============================================================================

def make_network_evidence(
    centrality,
    direct_edges,
    bridges,
):
    centrality, gid_col = ensure_gene_id(
        centrality,
        [
            "GENE_ID",
            "EFFECTOR_GENE_ID",
        ],
    )

    if gid_col and gid_col != "GENE_ID":
        centrality = centrality.rename(
            columns={
                gid_col:
                    "GENE_ID",
            }
        )

    direct_genes = set()

    if not direct_edges.empty:
        for col in [
            "GENE_A_ID",
            "GENE_B_ID",
        ]:
            if col in direct_edges.columns:
                direct_genes.update(
                    clean_gene_id(value)
                    for value
                    in direct_edges[col]
                    if clean_gene_id(value)
                )

    bridge_genes = set()

    if not bridges.empty:
        if (
            "GWAS2M_DISCORDANT_GENE_ID"
            in
            bridges.columns
        ):
            yes = bridges[
                bridges.get(
                    "CONNECTED",
                    ""
                )
                .astype(str)
                .eq("YES")
            ].copy()

            bridge_genes.update(
                clean_gene_id(v)
                for v
                in yes[
                    "GWAS2M_DISCORDANT_GENE_ID"
                ]
                if clean_gene_id(v)
            )

    rows = []

    all_genes = set(
        centrality[
            "GENE_ID"
        ].dropna().astype(str)
    ) if (
        not centrality.empty
        and
        "GENE_ID"
        in centrality.columns
    ) else set()

    all_genes |= direct_genes
    all_genes |= bridge_genes

    for gid in sorted(
        all_genes
    ):
        z = (
            centrality[
                centrality[
                    "GENE_ID"
                ]
                ==
                gid
            ].copy()
            if not centrality.empty
            else pd.DataFrame()
        )

        row = (
            z.iloc[0]
            if not z.empty
            else
            pd.Series(
                dtype=object
            )
        )

        rows.append({
            "GENE_ID":
                gid,

            "NETWORK_DEGREE":
                safe_float(
                    row.get(
                        "DEGREE",
                        np.nan,
                    )
                ),

            "NETWORK_WEIGHTED_DEGREE":
                safe_float(
                    row.get(
                        "WEIGHTED_DEGREE",
                        np.nan,
                    )
                ),

            "NETWORK_BETWEENNESS":
                safe_float(
                    row.get(
                        "BETWEENNESS",
                        np.nan,
                    )
                ),

            "NETWORK_PAGERANK":
                safe_float(
                    row.get(
                        "PAGERANK",
                        np.nan,
                    )
                ),

            "NETWORK_COMPONENT_ID":
                row.get(
                    "COMPONENT_ID",
                    np.nan,
                ),

            "NETWORK_COMMUNITY_ID":
                row.get(
                    "COMMUNITY_ID",
                    np.nan,
                ),

            "IN_DIRECT_PHYSICAL_NETWORK":
                (
                    "YES"
                    if gid
                    in direct_genes
                    else
                    "NO"
                ),

            "HAS_PATH_TO_OT_SUPPORTED_GENE":
                (
                    "YES"
                    if gid
                    in bridge_genes
                    else
                    "NO"
                ),
        })

    return pd.DataFrame(
        rows
    )


# =============================================================================
# DRUG / THERAPEUTIC MAP
# =============================================================================

def make_drug_map(
    mechanisms,
    phenotype_drugs,
):
    if mechanisms.empty:
        return pd.DataFrame()

    mechanisms, gid_col = ensure_gene_id(
        mechanisms,
        [
            "GENE_ID",
            "EFFECTOR_GENE_ID",
        ],
    )

    if gid_col and gid_col != "GENE_ID":
        mechanisms = mechanisms.rename(
            columns={
                gid_col:
                    "GENE_ID",
            }
        )

    rows = []

    for gid, group in mechanisms.groupby(
        "GENE_ID",
        sort=False,
    ):
        aligned = (
            group[
                group[
                    "DIRECTION_ALIGNMENT"
                ]
                .astype(str)
                .eq("ALIGNED")
            ]
            if "DIRECTION_ALIGNMENT"
            in group.columns
            else
            pd.DataFrame()
        )

        pheno = (
            group[
                group[
                    "PHENOTYPE_INDICATION_MATCH"
                ]
                .astype(str)
                .eq("YES")
            ]
            if "PHENOTYPE_INDICATION_MATCH"
            in group.columns
            else
            pd.DataFrame()
        )

        rows.append({
            "GENE_ID":
                clean_gene_id(
                    gid
                ),

            "N_DRUG_MECHANISM_ROWS":
                len(
                    group
                ),

            "N_UNIQUE_DRUGS":
                group[
                    "MOLECULE_CHEMBL_ID"
                ]
                .astype(str)
                .nunique()
                if "MOLECULE_CHEMBL_ID"
                in group.columns
                else
                0,

            "MAX_CLINICAL_PHASE":
                safe_numeric(
                    group[
                        "MAX_PHASE_OVERALL"
                    ]
                ).max()
                if (
                    "MAX_PHASE_OVERALL"
                    in group.columns
                    and
                    safe_numeric(
                        group[
                            "MAX_PHASE_OVERALL"
                        ]
                    )
                    .notna()
                    .any()
                )
                else
                np.nan,

            "DRUG_NAMES":
                join_unique(
                    group[
                        "MOLECULE_NAME"
                    ]
                )
                if "MOLECULE_NAME"
                in group.columns
                else
                "",

            "MECHANISMS_OF_ACTION":
                join_unique(
                    group[
                        "MECHANISM_OF_ACTION"
                    ]
                )
                if "MECHANISM_OF_ACTION"
                in group.columns
                else
                "",

            "ACTION_TYPES":
                join_unique(
                    group[
                        "ACTION_TYPE"
                    ]
                )
                if "ACTION_TYPE"
                in group.columns
                else
                "",

            "N_PHENOTYPE_MATCHED_DRUGS":
                pheno[
                    "MOLECULE_CHEMBL_ID"
                ]
                .astype(str)
                .nunique()
                if (
                    not pheno.empty
                    and
                    "MOLECULE_CHEMBL_ID"
                    in pheno.columns
                )
                else
                0,

            "PHENOTYPE_MATCHED_DRUGS":
                join_unique(
                    pheno[
                        "MOLECULE_NAME"
                    ]
                )
                if (
                    not pheno.empty
                    and
                    "MOLECULE_NAME"
                    in pheno.columns
                )
                else
                "",

            "N_DIRECTIONALLY_ALIGNED_DRUGS":
                aligned[
                    "MOLECULE_CHEMBL_ID"
                ]
                .astype(str)
                .nunique()
                if (
                    not aligned.empty
                    and
                    "MOLECULE_CHEMBL_ID"
                    in aligned.columns
                )
                else
                0,

            "DIRECTIONALLY_ALIGNED_DRUGS":
                join_unique(
                    aligned[
                        "MOLECULE_NAME"
                    ]
                )
                if (
                    not aligned.empty
                    and
                    "MOLECULE_NAME"
                    in aligned.columns
                )
                else
                "",
        })

    return pd.DataFrame(
        rows
    )


# =============================================================================
# EXTERNAL CONCORDANCE / CLASSIFICATION
# =============================================================================

def external_concordance_class(
    row,
    strong_h4,
):
    h4 = safe_float(
        row.get(
            "MAX_H4",
            np.nan,
        )
    )

    ot_l2g = boolish(
        row.get(
            "OT_L2G_SUPPORTED_ANY",
            False,
        )
    )

    ot_match = boolish(
        row.get(
            "OT_ANY_LOCUS_MATCH",
            False,
        )
    )

    if pd.isna(h4):
        return "UNRESOLVED"

    if h4 < strong_h4:
        return "BELOW_STRONG_THRESHOLD"

    if ot_l2g:
        return "EXTERNALLY_CONVERGENT_OT_L2G"

    if ot_match:
        return "GWAS2M_ALTERNATIVE_EFFECTOR_AT_OT_MATCHED_LOCUS"

    return "OT_COVERAGE_GAP_OR_UNMATCHED_LOCUS"


# =============================================================================
# FINAL SCORE
# =============================================================================

def compute_final_score(
    table,
):
    x = table.copy()

    # -------------------------------------------------------------------------
    # 45 points: genetic / molecular evidence
    # -------------------------------------------------------------------------

    x["SCORE_GENETIC_MOLECULAR"] = (
        safe_numeric(
            x[
                "MAX_H4"
            ]
        )
        .fillna(0.0)
        .clip(
            0,
            1,
        )
        *
        30.0
    )

    x["SCORE_GENETIC_MOLECULAR"] += (
        x[
            "PRIOR_ROBUST"
        ]
        .map(
            boolish
        )
        .astype(int)
        *
        5.0
    )

    x["SCORE_GENETIC_MOLECULAR"] += (
        (
            safe_numeric(
                x[
                    "N_QTL_TYPES"
                ]
            )
            .fillna(0)
            >=
            2
        )
        .astype(int)
        *
        5.0
    )

    x["SCORE_GENETIC_MOLECULAR"] += (
        x[
            "IS_LOCUS_TOP_STRONG"
        ]
        .map(
            boolish
        )
        .astype(int)
        *
        5.0
    )

    # -------------------------------------------------------------------------
    # 15 points: external concordance
    # -------------------------------------------------------------------------

    x["SCORE_EXTERNAL"] = (
        x[
            "OT_L2G_SUPPORTED_ANY"
        ]
        .map(
            boolish
        )
        .astype(int)
        *
        10.0
    )

    x["SCORE_EXTERNAL"] += (
        x[
            "OT_ANY_LOCUS_MATCH"
        ]
        .map(
            boolish
        )
        .astype(int)
        *
        5.0
    )

    # -------------------------------------------------------------------------
    # 15 points: pathway mechanism evidence
    # -------------------------------------------------------------------------

    n_pathways = (
        safe_numeric(
            x[
                "N_SIGNIFICANT_PATHWAYS_PRIMARY"
            ]
        )
        .fillna(0)
    )

    pathway_component = np.minimum(
        n_pathways,
        5,
    ) / 5.0

    x["SCORE_PATHWAY"] = (
        pathway_component
        *
        10.0
    )

    x["SCORE_PATHWAY"] += (
        x[
            "IS_PRIMARY_PATHWAY_DRIVER"
        ]
        .map(
            boolish
        )
        .astype(int)
        *
        5.0
    )

    # -------------------------------------------------------------------------
    # 10 points: physical/network context
    # -------------------------------------------------------------------------

    pagerank_pct = percentile(
        x[
            "NETWORK_PAGERANK"
        ]
    )

    x["SCORE_NETWORK"] = (
        x[
            "IN_DIRECT_PHYSICAL_NETWORK"
        ]
        .map(
            boolish
        )
        .astype(int)
        *
        4.0
    )

    x["SCORE_NETWORK"] += (
        pagerank_pct
        *
        4.0
    )

    x["SCORE_NETWORK"] += (
        x[
            "HAS_PATH_TO_OT_SUPPORTED_GENE"
        ]
        .map(
            boolish
        )
        .astype(int)
        *
        2.0
    )

    # -------------------------------------------------------------------------
    # 10 points: tractability / clinical precedence
    # -------------------------------------------------------------------------

    any_tractability = (
        x[
            [
                "OT_SM_TRACTABLE",
                "OT_AB_TRACTABLE",
                "OT_PROTAC_TRACTABLE",
                "OT_OTHER_CLINICAL_TRACTABLE",
            ]
        ]
        .apply(
            lambda col: col.map(
                boolish
            )
        )
        .any(
            axis=1
        )
    )

    x["SCORE_THERAPEUTIC"] = (
        any_tractability
        .astype(int)
        *
        3.0
    )

    x["SCORE_THERAPEUTIC"] += (
        (
            safe_numeric(
                x[
                    "N_UNIQUE_DRUGS"
                ]
            )
            .fillna(0)
            >
            0
        )
        .astype(int)
        *
        3.0
    )

    x["SCORE_THERAPEUTIC"] += (
        (
            safe_numeric(
                x[
                    "N_PHENOTYPE_MATCHED_DRUGS"
                ]
            )
            .fillna(0)
            >
            0
        )
        .astype(int)
        *
        4.0
    )

    # -------------------------------------------------------------------------
    # 5 points: direction
    # -------------------------------------------------------------------------

    x["SCORE_DIRECTION"] = 0.0

    x.loc[
        x[
            "DIRECTION_STATUS"
        ]
        .astype(str)
        .eq(
            "CONSENSUS"
        ),
        "SCORE_DIRECTION",
    ] += 2.0

    x.loc[
        safe_numeric(
            x[
                "N_DIRECTIONALLY_ALIGNED_DRUGS"
            ]
        )
        .fillna(0)
        >
        0,
        "SCORE_DIRECTION",
    ] += 3.0

    x["FINAL_EVIDENCE_SCORE"] = (
        x[
            "SCORE_GENETIC_MOLECULAR"
        ]
        +
        x[
            "SCORE_EXTERNAL"
        ]
        +
        x[
            "SCORE_PATHWAY"
        ]
        +
        x[
            "SCORE_NETWORK"
        ]
        +
        x[
            "SCORE_THERAPEUTIC"
        ]
        +
        x[
            "SCORE_DIRECTION"
        ]
    )

    x[
        "FINAL_EVIDENCE_SCORE"
    ] = (
        x[
            "FINAL_EVIDENCE_SCORE"
        ]
        .clip(
            0,
            100,
        )
    )

    return x


# =============================================================================
# RULE-BASED EVIDENCE TIER
# =============================================================================

def evidence_tier(
    row,
    strong_h4,
    suggestive_h4,
):
    h4 = safe_float(
        row.get(
            "MAX_H4",
            np.nan,
        )
    )

    prior = boolish(
        row.get(
            "PRIOR_ROBUST",
            False,
        )
    )

    multi_qtl = (
        safe_float(
            row.get(
                "N_QTL_TYPES",
                0,
            ),
            default=0,
        )
        >=
        2
    )

    ot = boolish(
        row.get(
            "OT_L2G_SUPPORTED_ANY",
            False,
        )
    )

    pathway = boolish(
        row.get(
            "IS_PRIMARY_PATHWAY_DRIVER",
            False,
        )
    )

    direct_network = boolish(
        row.get(
            "IN_DIRECT_PHYSICAL_NETWORK",
            False,
        )
    )

    tractable = any([
        boolish(
            row.get(
                "OT_SM_TRACTABLE",
                False,
            )
        ),
        boolish(
            row.get(
                "OT_AB_TRACTABLE",
                False,
            )
        ),
        boolish(
            row.get(
                "OT_PROTAC_TRACTABLE",
                False,
            )
        ),
        boolish(
            row.get(
                "OT_OTHER_CLINICAL_TRACTABLE",
                False,
            )
        ),
    ])

    drug = (
        safe_float(
            row.get(
                "N_UNIQUE_DRUGS",
                0,
            ),
            default=0,
        )
        >
        0
    )

    support_layers = sum([
        int(pathway),
        int(direct_network),
        int(tractable),
        int(drug),
    ])

    if pd.isna(h4):
        return (
            "U_UNRESOLVED",
            "No interpretable H4."
        )

    if h4 >= strong_h4:
        if ot and support_layers >= 1:
            return (
                "A1_EXTERNAL_CONVERGENT",
                "Strong molecular evidence + OT L2G concordance + at least one additional mechanism/therapeutic layer."
            )

        if (
            prior
            and
            multi_qtl
            and
            support_layers >= 2
        ):
            return (
                "A2_GWAS2M_MULTI_LAYER",
                "Strong prior-robust multi-QTL evidence + at least two pathway/network/therapeutic layers."
            )

        if (
            prior
            and
            (
                multi_qtl
                or
                support_layers >= 1
            )
        ):
            return (
                "B_STRONG_ROBUST",
                "Strong prior-robust molecular evidence with additional support."
            )

        return (
            "C_STRONG_MOLECULAR",
            "Strong H4 molecular evidence but limited independent/multi-layer support."
        )

    if h4 >= suggestive_h4:
        return (
            "D_SUGGESTIVE",
            "Suggestive H4 evidence; retain for sensitivity analysis."
        )

    return (
        "E_BACKGROUND",
        "Below Step17 prioritisation threshold."
    )


# =============================================================================
# MECHANISM SUMMARY
# =============================================================================

def mechanism_summary(
    row,
):
    symbol = clean_symbol(
        row.get(
            "EFFECTOR_GENE_SYMBOL",
            "",
        )
    )

    gid = clean_gene_id(
        row.get(
            "EFFECTOR_GENE_ID",
            "",
        )
    )

    label = (
        f"{symbol} ({gid})"
        if symbol
        else
        gid
    )

    h4 = safe_float(
        row.get(
            "MAX_H4",
            np.nan,
        )
    )

    parts = []

    if pd.notna(h4):
        parts.append(
            f"max H4={h4:.4f}"
        )

    qtl = str(
        row.get(
            "QTL_TYPES",
            "",
        )
    ).strip()

    if qtl:
        parts.append(
            f"QTL={qtl}"
        )

    if boolish(
        row.get(
            "PRIOR_ROBUST",
            False,
        )
    ):
        parts.append(
            "prior-robust"
        )

    if boolish(
        row.get(
            "OT_L2G_SUPPORTED_ANY",
            False,
        )
    ):
        parts.append(
            "OT-L2G concordant"
        )
    elif boolish(
        row.get(
            "OT_ANY_LOCUS_MATCH",
            False,
        )
    ):
        parts.append(
            "alternative effector at OT-matched locus"
        )
    else:
        parts.append(
            "OT coverage-gap/unmatched context"
        )

    pathways = split_tokens(
        row.get(
            "PRIMARY_PATHWAYS",
            "",
        )
    )

    if pathways:
        parts.append(
            "pathways="
            +
            " | ".join(
                pathways[:3]
            )
        )

    if boolish(
        row.get(
            "IN_DIRECT_PHYSICAL_NETWORK",
            False,
        )
    ):
        parts.append(
            "direct physical-network seed"
        )

    if boolish(
        row.get(
            "HAS_PATH_TO_OT_SUPPORTED_GENE",
            False,
        )
    ):
        parts.append(
            "network path to OT-supported gene"
        )

    n_drugs = safe_float(
        row.get(
            "N_UNIQUE_DRUGS",
            0,
        ),
        default=0,
    )

    if n_drugs > 0:
        parts.append(
            f"ChEMBL drugs={int(n_drugs)}"
        )

    n_pheno = safe_float(
        row.get(
            "N_PHENOTYPE_MATCHED_DRUGS",
            0,
        ),
        default=0,
    )

    if n_pheno > 0:
        parts.append(
            f"phenotype-matched drugs={int(n_pheno)}"
        )

    direction = str(
        row.get(
            "THERAPEUTIC_DIRECTION",
            "UNKNOWN",
        )
    )

    if direction and direction != "UNKNOWN":
        parts.append(
            f"direction={direction}"
        )
    else:
        parts.append(
            "therapeutic direction unresolved"
        )

    return (
        label
        +
        ": "
        +
        "; ".join(
            parts
        )
    )


# =============================================================================
# TERMINAL REPORT
# =============================================================================

def print_candidate_detail(
    table,
    compact=False,
):
    section(
        "FINAL CANDIDATE-BY-CANDIDATE EVIDENCE"
    )

    max_rows = (
        20
        if compact
        else
        None
    )

    z = table.head(
        max_rows
    ) if max_rows else table

    for row in z.to_dict(
        "records"
    ):
        symbol = clean_symbol(
            row.get(
                "EFFECTOR_GENE_SYMBOL",
                "",
            )
        )

        gid = clean_gene_id(
            row.get(
                "EFFECTOR_GENE_ID",
                "",
            )
        )

        print()

        print(
            f"#{int(row.get('STEP17_RANK', 0))} "
            f"{symbol or gid} ({gid})"
        )

        print(
            "-" * 125
        )

        print(
            f"  Evidence tier                    : "
            f"{row.get('EVIDENCE_TIER', '')}"
        )

        print(
            f"  Final evidence score             : "
            f"{safe_float(row.get('FINAL_EVIDENCE_SCORE', np.nan)):.3f}/100"
        )

        print(
            f"  Max H4                           : "
            f"{safe_float(row.get('MAX_H4', np.nan)):.6f}"
        )

        print(
            f"  Prior robust                     : "
            f"{row.get('PRIOR_ROBUST', '')}"
        )

        print(
            f"  QTL modalities                   : "
            f"{row.get('QTL_TYPES', '')}"
        )

        print(
            f"  Tissues                          : "
            f"{row.get('TISSUES', '')}"
        )

        print(
            f"  OT external class                : "
            f"{row.get('EXTERNAL_CONCORDANCE_CLASS', '')}"
        )

        print(
            f"  OT L2G same-gene support         : "
            f"{row.get('OT_L2G_SUPPORTED_ANY', '')}"
        )

        print(
            f"  Significant primary pathways     : "
            f"{int(safe_float(row.get('N_SIGNIFICANT_PATHWAYS_PRIMARY', 0), 0))}"
        )

        print(
            f"  Pathways                         : "
            f"{row.get('PRIMARY_PATHWAYS', '')}"
        )

        print(
            f"  Direct physical network          : "
            f"{row.get('IN_DIRECT_PHYSICAL_NETWORK', '')}"
        )

        print(
            f"  Path to OT-supported target      : "
            f"{row.get('HAS_PATH_TO_OT_SUPPORTED_GENE', '')}"
        )

        print(
            f"  Network PageRank                 : "
            f"{safe_float(row.get('NETWORK_PAGERANK', np.nan))}"
        )

        print(
            f"  Small-molecule tractable         : "
            f"{row.get('OT_SM_TRACTABLE', '')}"
        )

        print(
            f"  Antibody tractable               : "
            f"{row.get('OT_AB_TRACTABLE', '')}"
        )

        print(
            f"  PROTAC tractable                 : "
            f"{row.get('OT_PROTAC_TRACTABLE', '')}"
        )

        print(
            f"  ChEMBL drugs                     : "
            f"{int(safe_float(row.get('N_UNIQUE_DRUGS', 0), 0))}"
        )

        print(
            f"  Phenotype-matched drugs          : "
            f"{int(safe_float(row.get('N_PHENOTYPE_MATCHED_DRUGS', 0), 0))}"
        )

        print(
            f"  Drug names                       : "
            f"{row.get('DRUG_NAMES', '')}"
        )

        print(
            f"  Therapeutic direction            : "
            f"{row.get('THERAPEUTIC_DIRECTION', 'UNKNOWN')}"
        )

        print(
            f"  Mechanism summary                : "
            f"{row.get('MECHANISM_SUMMARY', '')}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--phenotype",
        required=True,
    )

    parser.add_argument(
        "--ancestry",
        required=True,
    )

    parser.add_argument(
        "--strong-h4",
        type=float,
        default=0.80,
    )

    parser.add_argument(
        "--suggestive-h4",
        type=float,
        default=0.50,
    )

    parser.add_argument(
        "--include-suggestive",
        action="store_true",
        help="Include H4 >= suggestive threshold in addition to CORE_STRONG.",
    )

    parser.add_argument(
        "--output-dir",
        default="",
    )

    parser.add_argument(
        "--compact-screen",
        action="store_true",
    )

    args = parser.parse_args()

    root = Path.cwd().resolve()

    ancestry_code, ancestry_label = ancestry_info(
        args.ancestry
    )

    phenotype_slug = slug(
        args.phenotype
    )

    ancestry_slug = results_slug(slug(ancestry_label))

    output = (
        Path(
            args.output_dir
        ).resolve()
        if args.output_dir
        else
        root
        /
        "17_mechanism"
        /
        phenotype_slug
        /
        ancestry_slug
    )

    output.mkdir(
        parents=True,
        exist_ok=True,
    )

    section(
        "GWAS2m STEP17 v3.0 (pathway-first, Step16 v3 source-ID integration) - FINAL DISEASE-MECHANISM INTEGRATION"
    )

    print(
        f"Root                     : {root}"
    )

    print(
        f"Phenotype                : {args.phenotype}"
    )

    print(
        f"Ancestry                 : {ancestry_label} ({ancestry_code})"
    )

    print(
        f"Strong H4               : >= {args.strong_h4}"
    )

    print(
        f"Suggestive H4           : >= {args.suggestive_h4}"
    )

    print(
        f"Include suggestive      : {args.include_suggestive}"
    )

    print(
        f"Output                   : {output}"
    )

    data = load_all_inputs(
        root,
        phenotype_slug,
        ancestry_slug,
    )

    gene_sets = gene_sets_from_membership(
        data[
            "step14_membership"
        ]
    )

    ranking = data[
        "step14_ranking"
    ].copy()

    ranking[
        "EFFECTOR_GENE_ID"
    ] = ranking[
        "EFFECTOR_GENE_ID"
    ].map(
        clean_gene_id
    )

    if args.include_suggestive:
        selected_ids = set(
            gene_sets.get(
                "EXPANDED_SUGGESTIVE",
                [],
            )
        )
    else:
        selected_ids = set(
            gene_sets.get(
                "CORE_STRONG",
                [],
            )
        )

    if not selected_ids:
        raise SystemExit(
            "No target genes found in the selected Step14 gene set."
        )

    final = ranking[
        ranking[
            "EFFECTOR_GENE_ID"
        ]
        .isin(
            selected_ids
        )
    ].copy()

    # -------------------------------------------------------------------------
    # Add set membership flags
    # -------------------------------------------------------------------------

    membership_flags = [
        "CORE_STRONG",
        "CORE_STRONG_PRIOR_ROBUST",
        "CORE_STRONG_MULTI_QTL",
        "LOCUS_TOP_STRONG",
        "CORE_STRONG_OT_SUPPORTED",
        "CORE_STRONG_OT_DISCORDANT",
        "CORE_STRONG_OT_COVERAGE_GAP",
    ]

    for name in membership_flags:
        final[
            "IS_"
            +
            name
        ] = (
            final[
                "EFFECTOR_GENE_ID"
            ]
            .isin(
                set(
                    gene_sets.get(
                        name,
                        [],
                    )
                )
            )
            .map({
                True: "YES",
                False: "NO",
            })
        )

    final[
        "IS_LOCUS_TOP_STRONG"
    ] = final[
        "IS_LOCUS_TOP_STRONG"
    ]

    # -------------------------------------------------------------------------
    # Pathway layer
    # -------------------------------------------------------------------------

    (
        pathway_significant,
        pathway_gene,
    ) = make_pathway_evidence(
        data[
            "step13_pathways"
        ],
        data[
            "step14_drivers"
        ],
    )

    if not pathway_gene.empty:
        final = final.merge(
            pathway_gene,
            left_on="EFFECTOR_GENE_ID",
            right_on="GENE_ID",
            how="left",
        ).drop(
            columns=[
                "GENE_ID",
            ],
            errors="ignore",
        )

    for col, default in {
        "N_SIGNIFICANT_PATHWAYS_ALL_SETS": 0,
        "N_SIGNIFICANT_PATHWAYS_PRIMARY": 0,
        "N_PATHWAY_DATABASES_PRIMARY": 0,
        "BEST_PATHWAY_FDR": np.nan,
        "PRIMARY_PATHWAYS": "",
        "PRIMARY_PATHWAY_SOURCES": "",
        "IS_PRIMARY_PATHWAY_DRIVER": "NO",
    }.items():
        if col not in final.columns:
            final[
                col
            ] = default
        else:
            final[
                col
            ] = final[
                col
            ].fillna(
                default
            )

    # -------------------------------------------------------------------------
    # Network layer
    # -------------------------------------------------------------------------

    network = make_network_evidence(
        data[
            "step15_centrality"
        ],
        data[
            "step15_direct"
        ],
        data[
            "step15_bridges"
        ],
    )

    if not network.empty:
        final = final.merge(
            network,
            left_on="EFFECTOR_GENE_ID",
            right_on="GENE_ID",
            how="left",
        ).drop(
            columns=[
                "GENE_ID",
            ],
            errors="ignore",
        )

    for col, default in {
        "NETWORK_DEGREE": np.nan,
        "NETWORK_WEIGHTED_DEGREE": np.nan,
        "NETWORK_BETWEENNESS": np.nan,
        "NETWORK_PAGERANK": np.nan,
        "NETWORK_COMPONENT_ID": np.nan,
        "NETWORK_COMMUNITY_ID": np.nan,
        "IN_DIRECT_PHYSICAL_NETWORK": "NO",
        "HAS_PATH_TO_OT_SUPPORTED_GENE": "NO",
    }.items():
        if col not in final.columns:
            final[
                col
            ] = default
        else:
            final[
                col
            ] = final[
                col
            ].fillna(
                default
            )

    # -------------------------------------------------------------------------
    # Step16 therapeutic layer
    # -------------------------------------------------------------------------

    step16 = data[
        "step16_targets"
    ].copy()

    step16, s16_gid = ensure_gene_id(
        step16,
        [
            "EFFECTOR_GENE_ID",
            "GENE_ID",
        ],
    )

    if s16_gid is not None:
        if s16_gid != "EFFECTOR_GENE_ID":
            step16 = step16.rename(
                columns={
                    s16_gid:
                        "EFFECTOR_GENE_ID",
                }
            )

        keep = [
            "EFFECTOR_GENE_ID",
            "OT_SM_TRACTABLE",
            "OT_AB_TRACTABLE",
            "OT_PROTAC_TRACTABLE",
            "OT_OTHER_CLINICAL_TRACTABLE",
            "CHEMBL_TARGET_FOUND",
            "N_UNIQUE_DRUGS",
            "N_APPROVED_OR_PHASE4_DRUGS",
            "CHEMBL_MAX_PHASE",
            "N_PHENOTYPE_MATCHED_DRUGS",
            "PHENOTYPE_MAX_PHASE",
            "N_DIRECTIONALLY_ALIGNED_DRUGS",
            "DIRECTIONALLY_ALIGNED_DRUGS",
            "PHENOTYPE_MATCHED_DRUGS",
            "DIRECTION_STATUS",
            "THERAPEUTIC_DIRECTION",
            "THERAPEUTIC_PRIORITY_SCORE",
            "THERAPEUTIC_PRIORITY_TIER",
        ]

        keep = [
            c
            for c in keep
            if c in step16.columns
        ]

        s16_small = (
            step16[
                keep
            ]
            .drop_duplicates(
                "EFFECTOR_GENE_ID"
            )
        )

        final = final.merge(
            s16_small,
            on="EFFECTOR_GENE_ID",
            how="left",
        )

    # Drug map
    drug_map = make_drug_map(
        data[
            "step16_mechanisms"
        ],
        data[
            "step16_phenotype_drugs"
        ],
    )

    if not drug_map.empty:
        # Add detailed drug strings only; numeric columns may already come from Step16.
        detailed_cols = [
            "GENE_ID",
            "DRUG_NAMES",
            "MECHANISMS_OF_ACTION",
            "ACTION_TYPES",
        ]

        detail = drug_map[
            [
                c
                for c in detailed_cols
                if c in drug_map.columns
            ]
        ].copy()

        final = final.merge(
            detail,
            left_on="EFFECTOR_GENE_ID",
            right_on="GENE_ID",
            how="left",
        ).drop(
            columns=[
                "GENE_ID",
            ],
            errors="ignore",
        )

    # Defaults
    defaults = {
        "OT_SM_TRACTABLE": "NO",
        "OT_AB_TRACTABLE": "NO",
        "OT_PROTAC_TRACTABLE": "NO",
        "OT_OTHER_CLINICAL_TRACTABLE": "NO",
        "CHEMBL_TARGET_FOUND": "NO",
        "N_UNIQUE_DRUGS": 0,
        "N_APPROVED_OR_PHASE4_DRUGS": 0,
        "CHEMBL_MAX_PHASE": np.nan,
        "N_PHENOTYPE_MATCHED_DRUGS": 0,
        "PHENOTYPE_MAX_PHASE": np.nan,
        "N_DIRECTIONALLY_ALIGNED_DRUGS": 0,
        "DIRECTIONALLY_ALIGNED_DRUGS": "",
        "PHENOTYPE_MATCHED_DRUGS": "",
        "DIRECTION_STATUS": "UNAVAILABLE",
        "THERAPEUTIC_DIRECTION": "UNKNOWN",
        "DRUG_NAMES": "",
        "MECHANISMS_OF_ACTION": "",
        "ACTION_TYPES": "",
    }

    for col, default in defaults.items():
        if col not in final.columns:
            final[
                col
            ] = default

        final[
            col
        ] = final[
            col
        ].fillna(
            default
        )

    # -------------------------------------------------------------------------
    # External concordance
    # -------------------------------------------------------------------------

    if "OT_ANY_LOCUS_MATCH" not in final.columns:
        final[
            "OT_ANY_LOCUS_MATCH"
        ] = "UNKNOWN"

    if "OT_L2G_SUPPORTED_ANY" not in final.columns:
        final[
            "OT_L2G_SUPPORTED_ANY"
        ] = "UNKNOWN"

    final[
        "EXTERNAL_CONCORDANCE_CLASS"
    ] = final.apply(
        lambda row:
        external_concordance_class(
            row,
            args.strong_h4,
        ),
        axis=1,
    )

    # -------------------------------------------------------------------------
    # Score + tier
    # -------------------------------------------------------------------------

    final = compute_final_score(
        final
    )

    tier_rows = final.apply(
        lambda row:
        evidence_tier(
            row,
            args.strong_h4,
            args.suggestive_h4,
        ),
        axis=1,
    )

    final[
        "EVIDENCE_TIER"
    ] = [
        x[0]
        for x
        in tier_rows
    ]

    final[
        "EVIDENCE_TIER_REASON"
    ] = [
        x[1]
        for x
        in tier_rows
    ]

    final[
        "MECHANISM_SUMMARY"
    ] = final.apply(
        mechanism_summary,
        axis=1,
    )

    tier_priority = {
        "A1_EXTERNAL_CONVERGENT": 0,
        "A2_GWAS2M_MULTI_LAYER": 1,
        "B_STRONG_ROBUST": 2,
        "C_STRONG_MOLECULAR": 3,
        "D_SUGGESTIVE": 4,
        "E_BACKGROUND": 5,
        "U_UNRESOLVED": 9,
    }

    final[
        "_TIER_SORT"
    ] = final[
        "EVIDENCE_TIER"
    ].map(
        tier_priority
    ).fillna(
        8
    )

    final = final.sort_values(
        [
            "_TIER_SORT",
            "FINAL_EVIDENCE_SCORE",
            "MAX_H4",
        ],
        ascending=[
            True,
            False,
            False,
        ],
        kind="stable",
    ).drop(
        columns=[
            "_TIER_SORT",
        ]
    ).reset_index(
        drop=True
    )

    final[
        "STEP17_RANK"
    ] = np.arange(
        1,
        len(final) + 1,
    )

    # Put rank first
    cols = [
        "STEP17_RANK",
    ] + [
        c
        for c in final.columns
        if c != "STEP17_RANK"
    ]

    final = final[
        cols
    ]

    # -------------------------------------------------------------------------
    # Derived output subsets
    # -------------------------------------------------------------------------

    external = final[
        [
            c
            for c in [
                "STEP17_RANK",
                "EFFECTOR_GENE_ID",
                "EFFECTOR_GENE_SYMBOL",
                "MAX_H4",
                "OT_ANY_LOCUS_MATCH",
                "OT_L2G_SUPPORTED_ANY",
                "OT_L2G_BEST_RANK",
                "EXTERNAL_CONCORDANCE_CLASS",
                "EVIDENCE_TIER",
                "FINAL_EVIDENCE_SCORE",
            ]
            if c in final.columns
        ]
    ].copy()

    tiers = final[
        [
            "STEP17_RANK",
            "EFFECTOR_GENE_ID",
            "EFFECTOR_GENE_SYMBOL",
            "MAX_H4",
            "EVIDENCE_TIER",
            "EVIDENCE_TIER_REASON",
            "FINAL_EVIDENCE_SCORE",
            "EXTERNAL_CONCORDANCE_CLASS",
        ]
    ].copy()

    high = final[
        final[
            "EVIDENCE_TIER"
        ]
        .isin(
            {
                "A1_EXTERNAL_CONVERGENT",
                "A2_GWAS2M_MULTI_LAYER",
                "B_STRONG_ROBUST",
            }
        )
    ].copy()

    alternatives = final[
        final[
            "EXTERNAL_CONCORDANCE_CLASS"
        ]
        ==
        "GWAS2M_ALTERNATIVE_EFFECTOR_AT_OT_MATCHED_LOCUS"
    ].copy()

    convergent = final[
        final[
            "EXTERNAL_CONCORDANCE_CLASS"
        ]
        ==
        "EXTERNALLY_CONVERGENT_OT_L2G"
    ].copy()

    coverage_gap = final[
        final[
            "EXTERNAL_CONCORDANCE_CLASS"
        ]
        ==
        "OT_COVERAGE_GAP_OR_UNMATCHED_LOCUS"
    ].copy()

    # Mechanism maps
    pathway_map = final[
        [
            c
            for c in [
                "STEP17_RANK",
                "EFFECTOR_GENE_ID",
                "EFFECTOR_GENE_SYMBOL",
                "MAX_H4",
                "N_SIGNIFICANT_PATHWAYS_PRIMARY",
                "N_PATHWAY_DATABASES_PRIMARY",
                "BEST_PATHWAY_FDR",
                "PRIMARY_PATHWAY_SOURCES",
                "PRIMARY_PATHWAYS",
                "IS_PRIMARY_PATHWAY_DRIVER",
            ]
            if c in final.columns
        ]
    ].copy()

    network_map = final[
        [
            c
            for c in [
                "STEP17_RANK",
                "EFFECTOR_GENE_ID",
                "EFFECTOR_GENE_SYMBOL",
                "MAX_H4",
                "IN_DIRECT_PHYSICAL_NETWORK",
                "HAS_PATH_TO_OT_SUPPORTED_GENE",
                "NETWORK_DEGREE",
                "NETWORK_WEIGHTED_DEGREE",
                "NETWORK_BETWEENNESS",
                "NETWORK_PAGERANK",
                "NETWORK_COMPONENT_ID",
                "NETWORK_COMMUNITY_ID",
            ]
            if c in final.columns
        ]
    ].copy()

    drug_out = final[
        [
            c
            for c in [
                "STEP17_RANK",
                "EFFECTOR_GENE_ID",
                "EFFECTOR_GENE_SYMBOL",
                "MAX_H4",
                "OT_SM_TRACTABLE",
                "OT_AB_TRACTABLE",
                "OT_PROTAC_TRACTABLE",
                "CHEMBL_TARGET_FOUND",
                "N_UNIQUE_DRUGS",
                "CHEMBL_MAX_PHASE",
                "N_PHENOTYPE_MATCHED_DRUGS",
                "PHENOTYPE_MATCHED_DRUGS",
                "DIRECTION_STATUS",
                "THERAPEUTIC_DIRECTION",
                "N_DIRECTIONALLY_ALIGNED_DRUGS",
                "DRUG_NAMES",
                "MECHANISMS_OF_ACTION",
            ]
            if c in final.columns
        ]
    ].copy()

    # -------------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------------

    final.to_csv(
        output
        /
        "Step17_Final_Target_Ranking.tsv",
        sep="\t",
        index=False,
    )

    final.to_csv(
        output
        /
        "Step17_Final_Evidence_Matrix.tsv",
        sep="\t",
        index=False,
    )

    tiers.to_csv(
        output
        /
        "Step17_Evidence_Tiers.tsv",
        sep="\t",
        index=False,
    )

    pathway_map.to_csv(
        output
        /
        "Step17_Pathway_Mechanism_Map.tsv",
        sep="\t",
        index=False,
    )

    network_map.to_csv(
        output
        /
        "Step17_Network_Mechanism_Map.tsv",
        sep="\t",
        index=False,
    )

    drug_out.to_csv(
        output
        /
        "Step17_Drug_Mechanism_Map.tsv",
        sep="\t",
        index=False,
    )

    external.to_csv(
        output
        /
        "Step17_External_Concordance.tsv",
        sep="\t",
        index=False,
    )

    high.to_csv(
        output
        /
        "Step17_High_Priority_Candidates.tsv",
        sep="\t",
        index=False,
    )

    alternatives.to_csv(
        output
        /
        "Step17_GWAS2m_Alternative_Effectors.tsv",
        sep="\t",
        index=False,
    )

    convergent.to_csv(
        output
        /
        "Step17_OpenTargets_Convergent_Targets.tsv",
        sep="\t",
        index=False,
    )

    coverage_gap.to_csv(
        output
        /
        "Step17_OT_Coverage_Gap_Targets.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Terminal report
    # -------------------------------------------------------------------------

    show(
        final,
        title="FINAL STEP17 TARGET RANKING",
        max_rows=None,
        columns=[
            "STEP17_RANK",
            "EFFECTOR_GENE_SYMBOL",
            "EFFECTOR_GENE_ID",
            "MAX_H4",
            "PRIOR_ROBUST",
            "N_QTL_TYPES",
            "OT_L2G_SUPPORTED_ANY",
            "EXTERNAL_CONCORDANCE_CLASS",
            "N_SIGNIFICANT_PATHWAYS_PRIMARY",
            "IN_DIRECT_PHYSICAL_NETWORK",
            "HAS_PATH_TO_OT_SUPPORTED_GENE",
            "OT_SM_TRACTABLE",
            "OT_AB_TRACTABLE",
            "OT_PROTAC_TRACTABLE",
            "N_UNIQUE_DRUGS",
            "N_PHENOTYPE_MATCHED_DRUGS",
            "THERAPEUTIC_DIRECTION",
            "EVIDENCE_TIER",
            "FINAL_EVIDENCE_SCORE",
        ],
    )

    show(
        tiers,
        title="STEP17 EVIDENCE TIERS",
    )

    show(
        external,
        title="OPEN TARGETS EXTERNAL CONCORDANCE",
    )

    show(
        pathway_map[
            safe_numeric(
                pathway_map.get(
                    "N_SIGNIFICANT_PATHWAYS_PRIMARY",
                    pd.Series(
                        0,
                        index=pathway_map.index,
                    ),
                )
            )
            >
            0
        ],
        title="GENES DRIVING PRIMARY SIGNIFICANT PATHWAYS",
    )

    show(
        network_map[
            (
                network_map.get(
                    "IN_DIRECT_PHYSICAL_NETWORK",
                    ""
                )
                .astype(str)
                .eq("YES")
            )
            |
            (
                network_map.get(
                    "HAS_PATH_TO_OT_SUPPORTED_GENE",
                    ""
                )
                .astype(str)
                .eq("YES")
            )
        ],
        title="NETWORK-SUPPORTED CANDIDATES",
    )

    show(
        drug_out[
            (
                safe_numeric(
                    drug_out.get(
                        "N_UNIQUE_DRUGS",
                        pd.Series(
                            0,
                            index=drug_out.index,
                        ),
                    )
                )
                >
                0
            )
            |
            (
                safe_numeric(
                    drug_out.get(
                        "N_PHENOTYPE_MATCHED_DRUGS",
                        pd.Series(
                            0,
                            index=drug_out.index,
                        ),
                    )
                )
                >
                0
            )
        ],
        title="DRUG / CLINICAL-PRECEDENCE CANDIDATES",
    )

    print_candidate_detail(
        final,
        compact=args.compact_screen,
    )

    # -------------------------------------------------------------------------
    # Plain-English summary
    # -------------------------------------------------------------------------

    section(
        "PLAIN-ENGLISH STEP17 INTERPRETATION"
    )

    tier_counts = (
        final[
            "EVIDENCE_TIER"
        ]
        .value_counts()
    )

    print(
        f"1. Step17 integrated {len(final)} final candidate genes."
    )

    print(
        f"2. Externally convergent OT-L2G targets: {len(convergent)}."
    )

    print(
        f"3. GWAS2m alternative effectors at OT-matched loci: {len(alternatives)}."
    )

    print(
        f"4. Strong candidates in OT coverage-gap/unmatched contexts: {len(coverage_gap)}."
    )

    print(
        f"5. Candidates driving >=1 primary significant pathway: "
        f"{int((safe_numeric(final['N_SIGNIFICANT_PATHWAYS_PRIMARY']) > 0).sum())}."
    )

    print(
        f"6. Candidates in the direct physical network: "
        f"{int(final['IN_DIRECT_PHYSICAL_NETWORK'].map(boolish).sum())}."
    )

    print(
        f"7. Candidates with a network path to an OT-supported strong gene: "
        f"{int(final['HAS_PATH_TO_OT_SUPPORTED_GENE'].map(boolish).sum())}."
    )

    any_tract = (
        final[
            [
                "OT_SM_TRACTABLE",
                "OT_AB_TRACTABLE",
                "OT_PROTAC_TRACTABLE",
                "OT_OTHER_CLINICAL_TRACTABLE",
            ]
        ]
        .apply(
            lambda col:
            col.map(
                boolish
            )
        )
        .any(
            axis=1
        )
    )

    print(
        f"8. Candidates with at least one positive tractability modality: "
        f"{int(any_tract.sum())}."
    )

    print(
        f"9. Candidates with ChEMBL drug mechanisms: "
        f"{int((safe_numeric(final['N_UNIQUE_DRUGS']) > 0).sum())}."
    )

    print(
        f"10. Candidates with phenotype-matched ChEMBL drugs: "
        f"{int((safe_numeric(final['N_PHENOTYPE_MATCHED_DRUGS']) > 0).sum())}."
    )

    print(
        f"11. Candidates with resolved therapeutic direction: "
        f"{int(final['DIRECTION_STATUS'].astype(str).eq('CONSENSUS').sum())}."
    )

    print()

    print(
        "Evidence tier counts:"
    )

    for tier, count in tier_counts.items():
        print(
            f"  {tier}: {count}"
        )

    if not final.empty:
        top = final.iloc[0]

        print()

        print(
            f"Top integrated candidate: "
            f"{clean_symbol(top.get('EFFECTOR_GENE_SYMBOL', '')) or top.get('EFFECTOR_GENE_ID', '')}"
        )

        print(
            f"  Tier : {top.get('EVIDENCE_TIER', '')}"
        )

        print(
            f"  Score: {safe_float(top.get('FINAL_EVIDENCE_SCORE', np.nan)):.3f}/100"
        )

        print(
            f"  Basis: {top.get('MECHANISM_SUMMARY', '')}"
        )

    print()

    print(
        "IMPORTANT: this is an integrated evidence ranking. "
        "It is not a causal probability and it does not establish clinical efficacy."
    )

    # -------------------------------------------------------------------------
    # QC
    # -------------------------------------------------------------------------

    qc_rows = [
        {
            "CHECK":
                "FINAL_CANDIDATE_GENES",
            "STATUS":
                "INFO",
            "VALUE":
                len(final),
            "DETAIL":
                (
                    "EXPANDED_SUGGESTIVE"
                    if args.include_suggestive
                    else
                    "CORE_STRONG"
                ),
        },
        {
            "CHECK":
                "STEP14_PATHWAY_LAYER",
            "STATUS":
                (
                    "OK"
                    if not pathway_gene.empty
                    else
                    "WARN"
                ),
            "VALUE":
                len(pathway_gene),
            "DETAIL":
                "Gene-level significant pathway driver summaries.",
        },
        {
            "CHECK":
                "STEP15_NETWORK_LAYER",
            "STATUS":
                (
                    "OK"
                    if not network.empty
                    else
                    "WARN"
                ),
            "VALUE":
                len(network),
            "DETAIL":
                "Gene-level network summaries.",
        },
        {
            "CHECK":
                "STEP16_DRUG_LAYER",
            "STATUS":
                (
                    "OK"
                    if not data[
                        "step16_targets"
                    ].empty
                    else
                    "WARN"
                ),
            "VALUE":
                len(
                    data[
                        "step16_targets"
                    ]
                ),
            "DETAIL":
                "Step16 target-prioritisation rows.",
        },
        {
            "CHECK":
                "THERAPEUTIC_DIRECTION",
            "STATUS":
                (
                    "OK"
                    if final[
                        "DIRECTION_STATUS"
                    ]
                    .astype(str)
                    .eq(
                        "CONSENSUS"
                    )
                    .any()
                    else
                    "WARN"
                ),
            "VALUE":
                int(
                    final[
                        "DIRECTION_STATUS"
                    ]
                    .astype(str)
                    .eq(
                        "CONSENSUS"
                    )
                    .sum()
                ),
            "DETAIL":
                (
                    "Direction is expected to remain unavailable when "
                    "signed allele-harmonised GWAS/QTL effects are absent."
                ),
        },
        {
            "CHECK":
                "FINAL_SCORE_IS_HEURISTIC",
            "STATUS":
                "INFO",
            "VALUE":
                "YES",
            "DETAIL":
                (
                    "FINAL_EVIDENCE_SCORE is not a probability; "
                    "all component scores are retained for auditability."
                ),
        },
    ]

    qc = pd.DataFrame(
        qc_rows
    )

    qc.to_csv(
        output
        /
        "Step17_QC_Audit.tsv",
        sep="\t",
        index=False,
    )

    show(
        qc,
        title="STEP17 QC AUDIT",
    )

    # -------------------------------------------------------------------------
    # Summary JSON
    # -------------------------------------------------------------------------

    summary = {
        "version":
            VERSION,

        "created_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),

        "phenotype":
            args.phenotype,

        "ancestry_code":
            ancestry_code,

        "ancestry_label":
            ancestry_label,

        "include_suggestive":
            args.include_suggestive,

        "strong_h4":
            args.strong_h4,

        "suggestive_h4":
            args.suggestive_h4,

        "n_final_candidates":
            len(
                final
            ),

        "evidence_tier_counts":
            {
                str(k):
                    int(v)
                for k, v
                in tier_counts.items()
            },

        "n_ot_convergent":
            len(
                convergent
            ),

        "n_gwas2m_alternative_effectors":
            len(
                alternatives
            ),

        "n_ot_coverage_gap":
            len(
                coverage_gap
            ),

        "n_pathway_supported":
            int(
                (
                    safe_numeric(
                        final[
                            "N_SIGNIFICANT_PATHWAYS_PRIMARY"
                        ]
                    )
                    >
                    0
                )
                .sum()
            ),

        "n_direct_physical_network":
            int(
                final[
                    "IN_DIRECT_PHYSICAL_NETWORK"
                ]
                .map(
                    boolish
                )
                .sum()
            ),

        "n_tractable":
            int(
                any_tract.sum()
            ),

        "n_targets_with_drugs":
            int(
                (
                    safe_numeric(
                        final[
                            "N_UNIQUE_DRUGS"
                        ]
                    )
                    >
                    0
                )
                .sum()
            ),

        "n_phenotype_matched_drugs":
            int(
                (
                    safe_numeric(
                        final[
                            "N_PHENOTYPE_MATCHED_DRUGS"
                        ]
                    )
                    >
                    0
                )
                .sum()
            ),

        "n_direction_resolved":
            int(
                final[
                    "DIRECTION_STATUS"
                ]
                .astype(str)
                .eq(
                    "CONSENSUS"
                )
                .sum()
            ),

        "scientific_notes": [
            "Step17 evidence tiers are rule-based summaries, not causal probabilities.",
            "Open Targets is an external comparator rather than a biological gold standard.",
            "Network and pathway support provide mechanism context, not causal proof.",
            "Drug precedence is target tractability/clinical context, not evidence of efficacy for the phenotype.",
            "Therapeutic direction remains unknown without signed allele-harmonised molecular evidence.",
        ],
    }

    (
        output
        /
        "Step17_summary.json"
    ).write_text(
        json.dumps(
            summary,
            indent=2,
            default=str,
        )
        +
        "\n",
        encoding="utf-8",
    )

    # -------------------------------------------------------------------------
    # Finish
    # -------------------------------------------------------------------------

    section(
        "STEP17 COMPLETE"
    )

    print(
        f"Output directory:"
    )

    print(
        f"  {output}"
    )

    print()

    print(
        "Main outputs:"
    )

    for name in [
        "Step17_Final_Target_Ranking.tsv",
        "Step17_Final_Evidence_Matrix.tsv",
        "Step17_Evidence_Tiers.tsv",
        "Step17_Pathway_Mechanism_Map.tsv",
        "Step17_Network_Mechanism_Map.tsv",
        "Step17_Drug_Mechanism_Map.tsv",
        "Step17_External_Concordance.tsv",
        "Step17_High_Priority_Candidates.tsv",
        "Step17_GWAS2m_Alternative_Effectors.tsv",
        "Step17_OpenTargets_Convergent_Targets.tsv",
        "Step17_OT_Coverage_Gap_Targets.tsv",
        "Step17_QC_Audit.tsv",
        "Step17_summary.json",
    ]:
        print(
            f"  {output / name}"
        )

    print()

    print(
        "PIPELINE STATUS:"
    )

    print(
        "  GWAS -> fine mapping -> annotation -> molecular QTL coloc (bulk + cell type)"
    )

    print(
        "  -> Open Targets benchmark -> pathways -> PPI/network -> druggability"
    )

    print(
        "  -> FINAL INTEGRATED DISEASE-MECHANISM PRIORITISATION COMPLETE."
    )


# The complete historical Step17 implementation above is retained as legacy_main.
legacy_main = main


# =============================================================================
# PATHWAY-FIRST MECHANISTIC INTEGRATION (Step17 v2)
# All files under 17_mechanism/<phenotype>/<ancestry>[_<panel>]. No HTML.
# Original Step17 gene/drug ranking helpers above are retained for provenance.
# =============================================================================

import hashlib as _mhash
import os as _mos
import time as _mtime

M17_VERSION = "3.0.0-step16v3-source-id-mechanisms"

M17_FILES = {
    "S12": ("12_master_table", "Step12_Locus_Gene_Master.tsv"),
    "S12LONG": ("12_master_table", "Step12_Molecular_Evidence_Long.tsv"),
    "S12LEAD": ("12_master_table", "Step12_Lead_Variant_Locus_Map.tsv"),
    "S12CS": ("12_master_table", "Step12_Credible_Set_Master.tsv"),
    "S13": ("13_opentargets_comparison", "Step13_Gene_Comparison.tsv"),
    "S14": ("14_pathways", "Step14_Pathway_Significant.tsv"),
    "S14GENES": ("14_pathways", "Step14_Gene_Ranking.tsv"),
    "S14DRIVERS": ("14_pathways", "Step14_Pathway_Gene_Drivers.tsv"),
    "S14SENS": ("14_pathways", "Step14_Pathway_Sensitivity_Summary.tsv"),
    "S15": ("15_network", "Step15_Gene_Pathway_Handoff.tsv"),
    "S15MAP": ("15_network", "Step15_Gene_Mapping_Audit.tsv"),
    "S15NET": ("15_network", "Step15_Pathway_Network_Evidence.tsv"),
    "S15CENTRAL": ("15_network", "Step15_Primary_Network_Centrality.tsv"),
    "S16": ("16_perturbation", "Step16_Directional_Pathway_Responses_All.tsv"),
    "S16MATCH": ("16_perturbation", "Step16_Step14_Pathway_Responses.tsv"),
    "S16MAP": ("16_perturbation", "Step16_Step14_Pathway_Mapping_Audit.tsv"),
    "S16LIB": ("16_perturbation", "Step16_Pathway_Library_Inventory.tsv"),
    "S16RUN": ("16_perturbation", "Step16_Run_Audit.tsv"),
    "S16COVER": ("16_perturbation", "Step16_Gene_Coverage_Summary.tsv"),
    "S16LINCS": ("16_perturbation", "Step16_LINCS_Gene_Coverage.tsv"),
}

M17_CORE_INPUTS = ("S12", "S14", "S14GENES", "S15", "S15MAP", "S16")
M17_NULLS = {"", "NAN", "NA", "NONE", "NULL", "UNKNOWN", "NOT_AVAILABLE", "NOT_MEASURED", "-"}


def m17_s(v):
    if v is None:
        return ""
    try:
        if pd.isna(v):
            return ""
    except (TypeError, ValueError):
        pass
    z = str(v).strip()
    return "" if z.upper() in M17_NULLS else z


def m17_float(v):
    try:
        f = float(v)
        return f if np.isfinite(f) else np.nan
    except (TypeError, ValueError):
        return np.nan


def m17_col(df, *options):
    names = {str(k).lower(): k for k in df.columns}
    return next((names[z.lower()] for z in options if z.lower() in names), None)


def m17_val(row, *options):
    for key in options:
        if key in row:
            v = m17_s(row[key])
            if v:
                return v
    return ""


def m17_norm(v):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", m17_s(v).lower())).strip()


def m17_token(v):
    return sorted(set(t.strip() for t in re.split(r"[;|,]", m17_s(v)) if t.strip() and t.strip().upper() not in M17_NULLS))


def m17_pathkey(row):
    src = m17_s(row.get("SOURCE", ""))
    ident = m17_s(row.get("TERM_ID", ""))
    term = m17_s(row.get("TERM_NAME", ""))
    return (m17_norm(src), m17_norm(ident) if ident else m17_norm(term))


def m17_sign(v):
    if m17_s(v).upper() in ("INCREASE", "UP", "POSITIVE", "+1", "RISK_UP", "INCREASED"):
        return 1
    if m17_s(v).upper() in ("DECREASE", "DOWN", "NEGATIVE", "-1", "RISK_DOWN", "DECREASED"):
        return -1
    try:
        x = float(v)
        return 1 if x > 0 else -1 if x < 0 else 0
    except (ValueError, TypeError):
        return None


def m17_read(path, mandatory=False):
    if not path.is_file():
        if mandatory:
            raise FileNotFoundError(f"Required upstream file missing: {path}")
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", low_memory=False, dtype={"TERM_ID": str, "GENE_ID": str})
    except (ValueError, pd.errors.EmptyDataError):
        return pd.read_csv(path, sep="\t", low_memory=False)


def m17_write(out, filename, rows, columns=None):
    out.mkdir(parents=True, exist_ok=True)
    df = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    if columns is not None:
        df = df.reindex(columns=columns)
    path = out / filename
    tmp = out / f".{filename}.tmp.{_mos.getpid()}"
    df.to_csv(tmp, sep="\t", index=False)
    _mos.replace(tmp, path)
    print(f"[SAVE] {filename}: {len(df):,} rows", flush=True)
    return df


def m17_paths(args):
    # Use the existing Step17 pipeline's ancestry and panel conventions, incl. --index.
    P = frame_paths(args)
    d = {}
    for key, (folder, filename) in M17_FILES.items():
        d[key] = P["root"] / folder / P["p"] / results_slug(P["a"]) / filename
    return P, d


def m17_gene_map(rank, mapping):
    gm = {}
    for df in (rank, mapping):
        if df.empty:
            continue
        idc = m17_col(df, "EFFECTOR_GENE_ID", "GENE_ID")
        symc = m17_col(df, "EFFECTOR_GENE_SYMBOL", "GENE_SYMBOL")
        if not idc:
            continue
        for d in df.to_dict("records"):
            gid = m17_s(d.get(idc, "")).split(".")[0]
            if gid:
                sym = m17_s(d.get(symc, "")) if symc else ""
                gm[gid] = sym or gm.get(gid, gid)
    # IDs and unique symbols only; ambiguous symbols are never silently guessed.
    lookup = {}
    for gid, sym in gm.items():
        lookup[gid.upper()] = gid
    by_symbol = {}
    for gid, sym in gm.items():
        if sym:
            by_symbol.setdefault(sym.upper(), set()).add(gid)
    for sym, ids in by_symbol.items():
        if len(ids) == 1:
            lookup[sym] = next(iter(ids))
    return gm, lookup


def m17_identifier(value, lookup):
    key = m17_s(value).split(".")[0].upper()
    return lookup.get(key, key if key.startswith("ENSG") else "")


def m17_derived_genes(d, lookup):
    # Preserve gene membership from the reported Step14 significant enrichment.
    vals = []
    for fld in ("INTERSECTION_GENES", "INTERSECTION_GENE_IDS"):
        for raw in m17_token(d.get(fld, "")):
            gid = m17_identifier(raw, lookup)
            if gid:
                vals.append(gid)
    return set(vals)


def m17_risk_effects(long_df, gm, lookup):
    # Do not infer QTL effect direction from H4 (it contains no effect sign).
    # Preserve only explicitly recorded, allele-harmonised risk directions.
    effects = {}
    if long_df.empty:
        return effects
    gidc = m17_col(long_df, "EFFECTOR_GENE_ID", "GENE_ID")
    dirc = m17_col(long_df, "RISK_ASSOCIATED_MOLECULAR_DIRECTION", "RISK_ALLELE_MOLECULAR_DIRECTION")
    harmc = m17_col(long_df, "ALLELE_HARMONIZATION_STATUS", "HARMONIZATION_STATUS")
    if not gidc or not dirc or not harmc:
        return effects
    for d in long_df.to_dict("records"):
        gid = m17_identifier(d.get(gidc), lookup)
        state = m17_s(d.get(harmc)).upper()
        sign = m17_sign(d.get(dirc))
        if gid and sign in (-1, 1) and state in {"MATCHED", "FLIPPED", "HARMONIZED", "HARMONISED", "PASS"}:
            effects.setdefault(gid, set()).add(sign)
    return effects


def m17_canonical_source(raw):
    """Normalize source names without crossing ontology boundaries."""
    v = re.sub(r"[\s_\-]+", "", m17_s(raw).upper())
    if v in {"GO:BP", "GOBP", "GOBIOLOGICALPROCESS"}: return "GO:BP"
    if v in {"WP", "WIKIPATHWAYS", "WIKIPATHWAY"}: return "WP"
    if v in {"KEGG", "KEGGLEGACY", "KEGGMEDICUS"}: return "KEGG"
    if v in {"REAC", "REACTOME", "REACTOMEPATHWAYS"}: return "REAC"
    return v


def m17_native_term_id(raw, source):
    """Extract source-native IDs using the same rules as Step16 v3."""
    z = m17_s(raw)
    source = m17_canonical_source(source)
    if source == "GO:BP":
        match = re.search(r"GO\s*:\s*(\d{7})", z, re.I)
        return f"GO:{match.group(1)}" if match else ""
    if source == "WP":
        match = re.search(r"\bWP\s*[-_: ]?\s*(\d{1,6})\b", z, re.I)
        return f"WP{match.group(1)}" if match else ""
    if source == "KEGG":
        match = re.search(r"(?:\bhsa|\bKEGG\s*:\s*)(\d{5})\b", z, re.I)
        return f"hsa{match.group(1)}" if match else ""
    if source == "REAC":
        match = re.search(r"R-HSA-\d+", z, re.I)
        return match.group(0).upper() if match else ""
    return ""


def m17_matching(exp, path, keys=None):
    """Accept only same-source, same-native-pathway-ID evidence.

    An explicit Step16 matching flag is cross-checked; string-name similarity
    NEVER establishes equivalence between GO / WikiPathways / KEGG / Reactome.
    Different gene-set releases may still have different memberships.
    """
    p_source = m17_canonical_source(path.get("SOURCE", ""))
    p_id = m17_native_term_id(path.get("TERM_ID", ""), p_source)
    e_source = m17_canonical_source(exp.get("PATHWAY_SOURCE", ""))
    e_id = m17_native_term_id(exp.get("PATHWAY_ID", ""), e_source)
    if p_id and p_source == e_source and p_id == e_id:
        return "EXACT_SOURCE_NATIVE_ID"
    if m17_s(exp.get("STEP14_PATHWAY_MATCH", "")).upper() == "YES" and m17_s(exp.get("PATHWAY_MATCH_METHOD", "")).upper() == "EXACT_SOURCE_NATIVE_ID":
        match_source = m17_canonical_source(exp.get("STEP14_SOURCE", ""))
        match_id = m17_native_term_id(exp.get("STEP14_TERM_ID", ""), match_source)
        if p_id and match_source == p_source and match_id == p_id and e_source == p_source and e_id == p_id:
            return "STEP16_V3_PREVERIFIED_EXACT_SOURCE_NATIVE_ID"
    return ""


def m17_tier(exp, alignment):
    if exp is None:
        return "MT1_UNRESOLVED_NO_MATCHING_EXPERIMENT"
    # FDR from permutations of gene labels is NOT experiment-level significance.
    independent = m17_s(exp.get("INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED", "")).upper() == "YES"
    q = m17_float(exp.get("EXPERIMENT_LEVEL_FDR"))
    replicated = independent and np.isfinite(q) and q <= 0.05
    if replicated and alignment == "DIRECTIONALLY_CONSISTENT_VALIDATED":
        return "MT4_GENETIC_EXPERIMENTAL_DIRECTIONAL_CONCORDANCE"
    if replicated:
        return "MT3_REPLICATE_SUPPORTED_PERTURBATION"
    return "MT2_EXPLORATORY_PERTURBATION_RESPONSE"


def m17_evidence_scores(row, exp):
    # Transparent *heuristic component* score, never a probability of causality.
    h4 = m17_float(row.get("MAX_H4"))
    gene = round(30 * max(0, min(1, h4)), 2) if np.isfinite(h4) else 0.0
    fdr = m17_float(row.get("STEP14_PATHWAY_FDR"))
    pathway = 15.0 if np.isfinite(fdr) and fdr <= .05 else 0.0
    ppi = 5.0 if m17_s(row.get("DIRECT_PHYSICAL_PPI")).upper() == "YES" else 0.0
    # Mere discovery is very weak; candidate experiment effect cannot be confirmed here.
    perturb = 5.0 if exp is not None else 0.0
    valid = 20.0 if m17_s(row.get("MECHANISM_TIER")).startswith(("MT3_", "MT4_")) else 0.0
    aligned = 25.0 if m17_s(row.get("MECHANISM_TIER")).startswith("MT4_") else 0.0
    return {"SCORE_GENETIC_H4_0_30": gene, "SCORE_PATHWAY_FDR_0_15": pathway,
            "SCORE_PHYSICAL_PPI_0_5": ppi, "SCORE_PERTURBATION_COVERAGE_0_5": perturb,
            "SCORE_EXPERIMENT_VALIDATION_0_20": valid, "SCORE_DIRECTIONAL_CONCORDANCE_0_25": aligned,
            "MECHANISM_EVIDENCE_SCORE_HEURISTIC": round(gene+pathway+ppi+perturb+valid+aligned, 2)}


def mechanism_main():
    parser = argparse.ArgumentParser(description="GWAS2Mech Step17 pathway-first mechanistic evidence tables")
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--include-suggestive", action="store_true",
                        help="Include suggestive gene candidates in the unassigned-gene audit")
    parser.add_argument("--strong-h4", type=float, default=.8)
    parser.add_argument("--suggestive-h4", type=float, default=.5)
    parser.add_argument("--compact-screen", action="store_true")
    parser.add_argument("--version", action="version", version="GWAS2Mech Step17 "+M17_VERSION)
    args = parser.parse_args()
    # pipeline_main already canonicalised --reference-panel and --index to PIPELINE.
    args.index = None   # --index was already resolved to PIPELINE["study"] by the pipeline frame
    args.reference_panel = PIPELINE["panel"]
    P, paths = m17_paths(args)
    out = Path(args.output_dir).resolve() if args.output_dir else P["OUT"]
    missing = [str(paths[k]) for k in M17_CORE_INPUTS if not paths[k].is_file()]
    if missing:
        raise SystemExit("Step17 required inputs missing (no placeholder mechanism results will be written):\n  " + "\n  ".join(missing))
    all_d = {k: m17_read(v, mandatory=k in M17_CORE_INPUTS) for k,v in paths.items()}
    step16_summary = paths["S16"].parent / "Step16_summary.json"
    if not step16_summary.is_file():
        raise SystemExit(f"Step16 completed summary missing: {step16_summary}")
    summary16 = json.loads(step16_summary.read_text())
    if not str(summary16.get("version", "")).startswith("3."):
        raise SystemExit(f"Step17 v3 needs COMPLETED Step16 v3 source-ID outputs; observed Step16 summary version: {summary16.get('version')}. Run Step16 v3 first.")
    if not (paths["S16"].parent / "Step16_Step14_Pathway_Mapping_Audit.tsv").is_file():
        raise SystemExit("Missing Step16 v3 pathway ID mapping audit. Step16 v3 must complete before Step17.")
    print(f"GWAS2Mech Step17 {M17_VERSION}; panel={PIPELINE['panel']}; output={out}", flush=True)
    genes = all_d["S14GENES"]; mapping = all_d["S15MAP"]
    gm, lookup = m17_gene_map(genes, mapping)
    if not gm:
        raise SystemExit("No usable gene identifiers from Step14/15")
    primary_gene_ids=set()
    for d in mapping.to_dict("records"):
        if m17_s(d.get("IS_GWAS2M_PRIMARY")).upper()=="YES":
            gid=m17_identifier(m17_val(d,"GENE_ID","EFFECTOR_GENE_ID"),lookup)
            if gid:primary_gene_ids.add(gid)
    gene_records = {}
    for d in genes.to_dict("records"):
        gid = m17_identifier(m17_val(d,"EFFECTOR_GENE_ID", "GENE_ID"), lookup)
        if gid:
            gene_records[gid] = d
    if all_d["S14"].empty:
        raise SystemExit("Step14_Pathway_Significant.tsv has zero rows; cannot make a pathway-first mechanism table")
    # Source + term ID + gene set yields genuine Step14 rows, then collapse by unique pathway.
    pathways = {}
    for d in all_d["S14"].to_dict("records"):
        key = m17_pathkey(d)
        if not all(key):
            continue
        if key not in pathways:
            pathways[key] = {"SOURCE":m17_s(d.get("SOURCE")),
                "TERM_ID":m17_s(d.get("TERM_ID")), "TERM_NAME":m17_s(d.get("TERM_NAME")),
                "SET_NAMES":set(), "DRIVER_IDS":set(), "STEP14_FDRS":[], "RAW_ROWS":[]}
        a = pathways[key]
        a["SET_NAMES"].add(m17_s(d.get("GENE_SET")))
        a["DRIVER_IDS"].update(m17_derived_genes(d,lookup))
        fdr = m17_float(d.get("FDR"))
        if np.isfinite(fdr): a["STEP14_FDRS"].append(fdr)
        a["RAW_ROWS"].append(d)
    # The Step14 driver table may enumerate genes even when INTERSECTION_GENES is not populated.
    driver_df=all_d["S14DRIVERS"]
    if not driver_df.empty:
        for d in driver_df.to_dict("records"):
            entry=pathways.get(m17_pathkey(d))
            if entry:
                for token in (m17_val(d,"GENE_ID","EFFECTOR_GENE_ID"),
                              m17_val(d,"GENE_SYMBOL","EFFECTOR_GENE_SYMBOL")):
                    gid=m17_identifier(token,lookup)
                    if gid:entry["DRIVER_IDS"].add(gid)
    # Step15 retains the curated gene -> significant pathway handoff.
    handoff = all_d["S15"]
    for d in handoff.to_dict("records"):
        entry=pathways.get(m17_pathkey(d))
        if entry:
            gid=m17_identifier(m17_val(d,"GENE_ID","EFFECTOR_GENE_ID","GENE_SYMBOL"),lookup)
            if gid:entry["DRIVER_IDS"].add(gid)
    # Keep source IDs separate; no fuzzy cross-ontology matching.
    native_lookup={}
    pathway_names={}
    for key, pt in pathways.items():
        src=m17_canonical_source(pt["SOURCE"])
        pid=m17_native_term_id(pt["TERM_ID"],src)
        if pid:
            native_lookup.setdefault((src,pid),key)
        pathway_names.setdefault(m17_norm(pt["TERM_NAME"]),set()).add(src)
    exp_map={}; mapping_audit=[]; experiment_rows=[]; raw_unmapped_responses=[]
    for i,d in enumerate(all_d["S16"].to_dict("records")):
        gene=m17_identifier(m17_val(d,"GENE_ID","GENE_SYMBOL"),lookup)
        exp_src=m17_canonical_source(d.get("PATHWAY_SOURCE"))
        exp_pid=m17_native_term_id(d.get("PATHWAY_ID"),exp_src)
        key=native_lookup.get((exp_src,exp_pid)) if gene and exp_pid else None
        reason=m17_matching(d,pathways[key]) if key is not None else ""
        matches=[(key,reason)] if key is not None and reason else []
        if not matches:
            raw_unmapped_responses.append(d)
            same_name_other_ontology = any(
                source!=exp_src for source in pathway_names.get(m17_norm(d.get("PATHWAY_NAME")),set()))
            same_source_name_only = exp_src in pathway_names.get(m17_norm(d.get("PATHWAY_NAME")),set())
            unmapped_status = (
                "UNRESOLVED_GENE_MAPPING" if not gene else
                "CROSS_ONTOLOGY_NAME_ONLY_NOT_ACCEPTED" if same_name_other_ontology else
                "SAME_SOURCE_NAME_ONLY_NOT_ACCEPTED" if same_source_name_only else
                "STEP16_PATHWAY_NATIVE_ID_MISSING" if not exp_pid else
                "NO_STEP14_SIGNIFICANT_NATIVE_ID_MATCH")
            mapping_audit.append({"GENE_SYMBOL":m17_val(d,"GENE_SYMBOL"),
              "GENE_ID":m17_val(d,"GENE_ID"),
              "DATASET":m17_val(d,"DATASET"),"STEP16_SOURCE":exp_src,"STEP16_TERM_ID":exp_pid,
              "STEP16_PATHWAY_NAME":m17_val(d,"PATHWAY_NAME"),
              "STEP16_PATHWAY_DESCRIPTION":m17_val(d,"PATHWAY_DESCRIPTION"),
              "MATCH_STATUS":unmapped_status,"MATCH_METHOD":"",
              "NOTE":"Unmatched does not mean biological non-response"})
        for key,how in matches:
            exp_map.setdefault((key,gene),[]).append(d)
            mapping_audit.append({"GENE_SYMBOL":m17_val(d,"GENE_SYMBOL"),
              "GENE_ID":gene,"DATASET":m17_val(d,"DATASET"),
              "STEP16_SOURCE":exp_src,"STEP16_TERM_ID":exp_pid,
              "STEP16_PATHWAY_NAME":m17_val(d,"PATHWAY_NAME"),
              "STEP16_PATHWAY_DESCRIPTION":m17_val(d,"PATHWAY_DESCRIPTION"),
              "MATCH_STATUS":"MATCHED_EXACT", "MATCH_METHOD":how,
              "NOTE":"Same pathway ID; gene-set release definitions may differ; no causal inference"})
            experiment_rows.append({"PATHWAY_SOURCE":pathways[key]["SOURCE"],"PATHWAY_ID":pathways[key]["TERM_ID"],
              "PATHWAY_NAME":pathways[key]["TERM_NAME"],"GENE_ID":gene,
              "GENE_SYMBOL":gm.get(gene,m17_val(d,"GENE_SYMBOL")),"DATASET":d.get("DATASET", ""),
              "CELL_TYPE":d.get("CELL_TYPE", ""),"NES":m17_float(d.get("NES")),
              "GENE_SET_PERMUTATION_FDR":m17_float(d.get("GENE_SET_PERMUTATION_FDR")),
              "GENE_SET_NOMINAL_P":m17_float(d.get("GENE_SET_NOMINAL_P")),
              "EXPERIMENT_LEVEL_FDR":m17_float(d.get("EXPERIMENT_LEVEL_FDR")),
              "RESPONSE_DIRECTION":d.get("RESPONSE_DIRECTION", ""),"EVIDENCE_CLASS":d.get("EVIDENCE_CLASS", ""),
              "TARGET_ACTIVITY_DIRECTION":d.get("PERTURBATION_ACTIVITY_DIRECTION", ""),
              "TARGET_EXPRESSION_DIFF":m17_float(d.get("TARGET_EXPRESSION_DIFF")),
              "TARGET_EXPRESSION_MEASURED":d.get("TARGET_EXPRESSION_MEASURED", ""),
              "INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED":d.get("INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED", ""),
              "LEADING_EDGE_GENES":d.get("LEADING_EDGE_GENES", ""),
              "LIBRARY_NAME":d.get("LIBRARY_NAME", ""),"LIBRARY_SHA256":d.get("LIBRARY_SHA256", ""),
              "PATHWAY_GENE_SET_SIZE":d.get("PATHWAY_GENE_SET_SIZE", ""),
              "FDR_SCOPE":d.get("FDR_SCOPE", ""),
              "INTERSECTION_WITH_STEP14_DRIVER": "YES" if gene in pathways[key]["DRIVER_IDS"] else "NO",
              "MATCH_METHOD":how})
    # Gene-level genetic and Open Targets evidence; do not count many rows as
    # multiple independent signals. Preserve observed sources without making effects up.
    master_by_gene = {}
    for d in all_d["S12"].to_dict("records"):
        gid=m17_identifier(m17_val(d,"EFFECTOR_GENE_ID","GENE_ID"),lookup)
        if gid:master_by_gene.setdefault(gid,[]).append(d)
    ot_by_gene={}
    for d in all_d["S13"].to_dict("records"):
        gid=m17_identifier(m17_val(d,"EFFECTOR_GENE_ID","GENE_ID"),lookup)
        if gid:ot_by_gene.setdefault(gid,[]).append(d)
    pathway_net={}
    for d in all_d["S15NET"].to_dict("records"):
        key=m17_pathkey(d)
        if key in pathways:pathway_net[key]=d
    pathway_sensitivity={}
    for d in all_d["S14SENS"].to_dict("records"):
        key=m17_pathkey(d)
        if key in pathways:
            atype=m17_val(d,"ANALYSIS_TYPE")
            pathway_sensitivity.setdefault(key,{})[atype]=d

    # Network support and context metadata, strictly gene-specific.
    network_by_gene = {}
    central = all_d["S15CENTRAL"]
    for d in central.to_dict("records"):
        gid=m17_identifier(m17_val(d,"GENE_ID","EFFECTOR_GENE_ID","GENE_SYMBOL"),lookup)
        if gid:network_by_gene[gid]=d
    risk=m17_risk_effects(all_d["S12LONG"],gm,lookup)
    exp_coverage={}; lin_coverage={}
    for d in all_d["S16COVER"].to_dict("records"):
        gid=m17_identifier(m17_val(d,"GENE_ID","GENE_SYMBOL"),lookup)
        if gid:exp_coverage[gid]=d
    for d in all_d["S16LINCS"].to_dict("records"):
        gid=m17_identifier(m17_val(d,"GENE_ID","GENE_SYMBOL"),lookup)
        if gid:lin_coverage[gid]=d
    gene_evidence=[]; mechanisms=[]; pathway_master=[]; missing_rows=[]
    for key, p in sorted(pathways.items(),key=lambda kv:(kv[1]["SOURCE"],kv[1]["TERM_NAME"])):
        drivers=p["DRIVER_IDS"]
        fdr=min(p["STEP14_FDRS"]) if p["STEP14_FDRS"] else np.nan
        # Also preserve perturbation-only links, but distinguish them from Step14 membership.
        perturbed_for_path={gene for (k,gene) in exp_map if k==key}
        union_genes=drivers | perturbed_for_path
        perpath=[]
        if not union_genes:
            union_genes={""}
        for gid in sorted(union_genes):
            gene_detail=gene_records.get(gid,{})
            sym=gm.get(gid, gid) if gid else ""
            exp_list=exp_map.get((key,gid),[]) if gid else []
            if not exp_list: exp_list=[None]
            for exp in exp_list:
                risk_dirs=risk.get(gid,set())
                risk_state=("CONFLICTING_HARMONISED_QTL_SIGNS" if len(risk_dirs)>1 else
                            "RISK_ASSOCIATED_MOLECULAR_INCREASE" if risk_dirs=={1} else
                            "RISK_ASSOCIATED_MOLECULAR_DECREASE" if risk_dirs=={-1} else
                            "UNKNOWN_NO_HARMONISED_SIGNED_EFFECT")
                direction="UNRESOLVED_NO_VERIFIED_TARGET_ENGAGEMENT_AND_GENETIC_EFFECT"
                # Experimental target direction is deliberately not inferred from CRISPRi label.
                if exp is not None and len(risk_dirs)==1:
                    direction="UNRESOLVED_TARGET_ENGAGEMENT_OR_COMPARABILITY"
                elif exp is None:
                    direction="NOT_TESTED"
                tier=m17_tier(exp,direction)
                network=network_by_gene.get(gid,{})
                gm_sources=master_by_gene.get(gid,[])
                ot_sources=ot_by_gene.get(gid,[])
                loci={m17_val(e,"LOCUS_ID") for e in gm_sources if m17_val(e,"LOCUS_ID")}
                studies={m17_val(e,"STUDY_ACCESSION") for e in gm_sources if m17_val(e,"STUDY_ACCESSION")}
                h4_s12=[m17_float(m17_val(e,"PP_H4","MAX_H4")) for e in gm_sources]
                h4_s12=[v for v in h4_s12 if np.isfinite(v)]
                ot_l2g="YES" if any(m17_s(e.get("OT_L2G_SUPPORTED_ANY")).upper()=="YES" for e in ot_sources) else "NOT_OBSERVED_OR_NOT_QUERIED"
                net_path=pathway_net.get(key,{})
                sens=pathway_sensitivity.get(key,{})
                network_degree=m17_float(m17_val(network,"DEGREE","NETWORK_DEGREE"))
                miss=[]
                if not gid:miss.append("NO_PRIORITISED_DRIVER_GENE_IDENTIFIED")
                if gid and gid not in drivers: miss.append("PERTURBATION_ONLY_NOT_A_STEP14_PATHWAY_DRIVER")
                if exp is None:
                    if m17_float(exp_coverage.get(gid,{}).get("N_DATASETS")) > 0:
                        miss.append("GENE_PERTURBED_BUT_NO_IDENTICAL_SOURCE_NATIVE_PATHWAY_ID_RESPONSE")
                    else:
                        miss.append("NO_EXACT_PATHWAY_MATCHED_PERTURBATION_RESULT")
                if len(risk_dirs)!=1:miss.append("MISSING_OR_CONFLICTING_HARMONISED_GWAS_QTL_DIRECTION")
                if exp is not None:
                    if m17_s(exp.get("PERTURBATION_ACTIVITY_DIRECTION", "")).upper() in {"", "NOT_MEASURED_DO_NOT_INFER_FROM_MODE", "NOT_VERIFIED_FROM_TARGET_TRANSCRIPT_DIFFERENCE"}:
                        miss.append("PERTURBATION_TARGET_ENGAGEMENT_NOT_MEASURED")
                    if not (m17_s(exp.get("INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED", "")).upper()=="YES" and
                        np.isfinite(m17_float(exp.get("EXPERIMENT_LEVEL_FDR")))):
                        miss.append("NO_VERIFIED_INDEPENDENT_REPLICATE_LEVEL_FDR")
                if np.isnan(network_degree):miss.append("DIRECT_PPI_NOT_QUANTIFIED_IN_MAPPED_NETWORK")
                no_inference=";".join(sorted(set(miss)))
                dataset=m17_s(exp.get("DATASET")) if exp else ""
                ctype=m17_s(exp.get("CELL_TYPE")) if exp else ""
                mech_key="|".join([p["SOURCE"],p["TERM_ID"] or p["TERM_NAME"],gid,dataset,ctype])
                mid="MECH_"+_mhash.sha1(mech_key.encode()).hexdigest()[:16]
                row={"MECHANISM_ID":mid,"PHENOTYPE":args.phenotype,"ANCESTRY":args.ancestry,"REFERENCE_PANEL":PIPELINE["panel"],
                    "PATHWAY_SOURCE":p["SOURCE"],"PATHWAY_ID":p["TERM_ID"],"PATHWAY_NAME":p["TERM_NAME"],
                    "STEP14_GENE_SETS":";".join(sorted(x for x in p["SET_NAMES"] if x)),
                    "STEP14_PATHWAY_FDR":fdr,"GENE_ID":gid,"GENE_SYMBOL":sym,
                    "IS_GWAS2M_PRIMARY_GENE":"YES" if gid in primary_gene_ids else "NO_OR_UNVERIFIED",
                    "PATHWAY_RELATION":"STEP14_GENETIC_DRIVER" if gid in drivers else
                      "PERTURBATION_RESPONSE_ONLY_NO_STEP14_MEMBERSHIP" if gid else "NO_GENE_DRIVER",
                    "MAX_H4":m17_float(m17_val(gene_detail,"MAX_H4")),
                    "STEP14_EVIDENCE_SCORE":m17_float(m17_val(gene_detail,"EVIDENCE_SCORE_STEP14")),
                    "PRIOR_ROBUST":m17_val(gene_detail,"PRIOR_ROBUST"),
                    "QTL_MODALITIES":m17_val(gene_detail,"QTL_TYPES","BEST_QTL_TYPE"),
                    "N_STEP12_LOCUS_IDS":len(loci),"N_STEP12_STUDIES":len(studies),
                    "STEP12_LOCUS_IDS":";".join(sorted(loci)),"STEP12_STUDY_ACCESSIONS":";".join(sorted(studies)),
                    "STEP12_MAX_PP_H4":max(h4_s12) if h4_s12 else np.nan,
                    "OT_L2G_SUPPORT":ot_l2g,
                    "STEP15_N_DIRECT_PRIMARY_PPI_EDGES_IN_PATHWAY":m17_float(net_path.get("N_DIRECT_PRIMARY_STRING_EDGES")),
                    "STEP15_PRIMARY_STRING_GENE_PAIRS":m17_s(net_path.get("PRIMARY_STRING_GENE_PAIRS")),
                    "STEP14_REGION_LOO_FRAC_SIGNIFICANT":m17_float(sens.get("REGION_LOO",{}).get("FRAC_SIGNIFICANT")),
                    "STEP14_GENE_LOO_FRAC_SIGNIFICANT":m17_float(sens.get("GENE_LOO",{}).get("FRAC_SIGNIFICANT")),
                    "RISK_ASSOCIATED_MOLECULAR_DIRECTION":risk_state,
                    "DIRECT_PHYSICAL_PPI":"YES" if np.isfinite(network_degree) and network_degree>0 else "UNRESOLVED",
                    "NETWORK_DEGREE":network_degree,
                    "N_GENE_SCPERTURB_DATASETS":m17_float(exp_coverage.get(gid,{}).get("N_DATASETS")),
                    "N_GENE_SCPERTURB_CONTRASTS":m17_float(exp_coverage.get(gid,{}).get("N_CONTRASTS")),
                    "N_GENE_LINCS_SIGNATURES_METADATA_ONLY":m17_float(lin_coverage.get(gid,{}).get("N_GENETIC_L1000_SIGNATURES")),
                    "DATASET":dataset,"CELL_TYPE":ctype,"EXPERIMENT_MODE":m17_s(exp.get("EXPERIMENT_MODE")) if exp else "",
                    "N_TARGET_CELLS":m17_float(exp.get("TARGET_CELLS")) if exp else np.nan,
                    "N_CONTROL_CELLS":m17_float(exp.get("CONTROL_CELLS")) if exp else np.nan,
                    "PATHWAY_NES":m17_float(exp.get("NES")) if exp else np.nan,
                    "GENE_SET_PERMUTATION_FDR":m17_float(exp.get("GENE_SET_PERMUTATION_FDR")) if exp else np.nan,
                    "EXPERIMENT_LEVEL_FDR":m17_float(exp.get("EXPERIMENT_LEVEL_FDR")) if exp else np.nan,
                    "GENE_SET_NOMINAL_P":m17_float(exp.get("GENE_SET_NOMINAL_P")) if exp else np.nan,
                    "STEP16_LEADING_EDGE_GENES":m17_s(exp.get("LEADING_EDGE_GENES")) if exp else "",
                    "STEP16_EXPRESSION_SCALE":m17_s(exp.get("EXPRESSION_SCALE")) if exp else "",
                    "STEP16_EXPRESSION_TRANSFORM":m17_s(exp.get("EXPRESSION_TRANSFORM")) if exp else "",
                    "STEP16_FDR_SCOPE":m17_s(exp.get("FDR_SCOPE")) if exp else "",
                    "STEP16_CONTROL_SAMPLING":m17_s(exp.get("CONTROL_SAMPLING")) if exp else "",
                    "STEP16_LIBRARY_NAME":m17_s(exp.get("LIBRARY_NAME")) if exp else "",
                    "STEP16_LIBRARY_SHA256":m17_s(exp.get("LIBRARY_SHA256")) if exp else "",
                    "STEP16_SOURCE_NATIVE_ID_MATCH_METHOD":m17_s(exp.get("PATHWAY_MATCH_METHOD")) if exp else "",
                    "STEP16_PATHWAY_GENE_SET_SIZE":m17_float(exp.get("PATHWAY_GENE_SET_SIZE")) if exp else np.nan,
                    "STEP16_TARGET_EXPRESSION_DIFF_DESCRIPTIVE":m17_float(exp.get("TARGET_EXPRESSION_DIFF")) if exp else np.nan,
                    "STEP16_TARGET_EXPRESSION_MEASURED":m17_s(exp.get("TARGET_EXPRESSION_MEASURED")) if exp else "",
                    "STEP16_INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED":m17_s(exp.get("INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED")) if exp else "",
                    "STEP16_INFERENCE_SCOPE":m17_s(exp.get("INFERENCE_SCOPE")) if exp else "",
                    "PERTURBATION_RESPONSE_DIRECTION":m17_s(exp.get("RESPONSE_DIRECTION")) if exp else "NOT_TESTED",
                    "EXPERIMENTAL_TARGET_ACTIVITY_DIRECTION":m17_s(exp.get("PERTURBATION_ACTIVITY_DIRECTION")) if exp else "NOT_TESTED",
                    "DIRECTION_ALIGNMENT":direction,"MECHANISM_TIER":tier,"MISSING_EVIDENCE":no_inference,
                    "EVIDENCE_SCOPE":"EXPLORATORY_PERMUTATION_ONLY_NOT_EXPERIMENT_LEVEL" if exp else "NO_MATCHED_EXPERIMENT",
                    "INTERPRETATION":"Context-specific transcriptional pathway response; not GWAS allele causality" if exp else
                        "Genetic candidate/pathway association without matched experimental response"}
                row.update(m17_evidence_scores(row,exp))
                mechanisms.append(row); perpath.append(row)
                if miss:
                    missing_rows.append({"MECHANISM_ID":mid,"PATHWAY_ID":p["TERM_ID"],"PATHWAY_NAME":p["TERM_NAME"],
                        "GENE_ID":gid,"GENE_SYMBOL":sym,"DATASET":dataset,"MECHANISM_TIER":tier,"MISSING_EVIDENCE":no_inference})
        unique_path_exp={(r["GENE_ID"],r["DATASET"]) for r in perpath if r["DATASET"]}
        covered={r["GENE_ID"] for r in perpath if r["DATASET"] and r["GENE_ID"] in drivers}
        pathway_master.append({"PATHWAY_SOURCE":p["SOURCE"],"PATHWAY_ID":p["TERM_ID"],
          "PATHWAY_NAME":p["TERM_NAME"],"SOURCE":p["SOURCE"],"STEP14_BEST_FDR":fdr,
          "STEP14_GENE_SETS":";".join(sorted(x for x in p["SET_NAMES"] if x)),
          "N_PRIORITISED_DRIVER_GENES":len(drivers),"DRIVER_GENES":";".join(gm.get(x,x) for x in sorted(drivers)),
          "N_PERTURBED_DRIVER_GENES":len(covered),"N_EXPERIMENT_CONTRASTS":len(unique_path_exp),
          "N_MECHANISM_ROWS":len(perpath),
          "STEP15_N_DIRECT_PRIMARY_PPI_EDGES":m17_float(pathway_net.get(key,{}).get("N_DIRECT_PRIMARY_STRING_EDGES")),
          "STEP15_PRIMARY_STRING_GENE_PAIRS":m17_s(pathway_net.get(key,{}).get("PRIMARY_STRING_GENE_PAIRS")),
          "STEP14_REGION_LOO_FRAC_SIGNIFICANT":m17_float(pathway_sensitivity.get(key,{}).get("REGION_LOO",{}).get("FRAC_SIGNIFICANT")),
          "STEP14_GENE_LOO_FRAC_SIGNIFICANT":m17_float(pathway_sensitivity.get(key,{}).get("GENE_LOO",{}).get("FRAC_SIGNIFICANT")),
          "HIGHEST_MECHANISM_TIER":max((r["MECHANISM_TIER"] for r in perpath),key=lambda t:int(t[2]) if t.startswith("MT") else 0),
          "BEST_EXPERIMENT_PERMUTATION_FDR":min((m17_float(r["GENE_SET_PERMUTATION_FDR"]) for r in perpath if np.isfinite(m17_float(r["GENE_SET_PERMUTATION_FDR"]))), default=np.nan),
          "PATHWAY_MAPPING_STATUS":"EXACT_SOURCE_NATIVE_ID_EXPERIMENT_MATCH" if unique_path_exp else
             "NO_MATCHED_EXPERIMENT_OR_INSUFFICIENT_SOURCE_NATIVE_ID_MAPPING",
          "MISSING_EVIDENCE_COUNT":sum(bool(r["MISSING_EVIDENCE"]) for r in perpath)})
    # All primary genes lacking any pathway mapping must remain visible in their own audit.
    genes_on_paths={g for p in pathways.values() for g in p["DRIVER_IDS"]}
    unassigned=[]
    for gid,gd in sorted(gene_records.items()):
        if gid not in genes_on_paths:
            unassigned.append({"GENE_ID":gid,"GENE_SYMBOL":gm.get(gid,gid),
                 "MAX_H4":m17_float(gd.get("MAX_H4")),"REASON":"NO_STEP14_SIGNIFICANT_PATHWAY_DRIVER_MEMBERSHIP",
                 "SC_PERTURB_DATASETS":m17_float(exp_coverage.get(gid,{}).get("N_DATASETS")),
                 "LINCS_SIGNATURES_METADATA_ONLY":m17_float(lin_coverage.get(gid,{}).get("N_GENETIC_L1000_SIGNATURES"))})
    out.mkdir(parents=True,exist_ok=True)
    m17_write(out,"Step17_Pathway_Master.tsv",pathway_master)
    # Source-faithful genetic / locus / annotation evidence. Never equate the
    # highest-PIP variant, representative coloc variant and GWAS lead allele.
    lead_index={}
    for z in all_d["S12LEAD"].to_dict("records"):
        lead_index.setdefault((m17_val(z,"STUDY_ACCESSION"),m17_val(z,"LOCUS_ID")),[]).append(z)
    cs_index={}
    for z in all_d["S12CS"].to_dict("records"):
        cs_index.setdefault((m17_val(z,"STUDY_ACCESSION"),m17_val(z,"LOCUS_ID")),[]).append(z)
    locus_rows=[]; best_locus_by_gene={}
    for gid, entries in master_by_gene.items():
        for z in entries:
            study=m17_val(z,"STUDY_ACCESSION"); locus=m17_val(z,"LOCUS_ID")
            lead=lead_index.get((study,locus),[])
            credible=cs_index.get((study,locus),[])
            row={"GENE_ID":gid,"GENE_SYMBOL":gm.get(gid,gid),
                 "STUDY_ACCESSION":study,"LOCUS_ID":locus,
                 "BEST_COLOC_PP_H4":m17_float(m17_val(z,"BEST_PP_H4","BEST_H4")),
                 "COLOC_CLASS":m17_val(z,"BEST_EVIDENCE_CLASS"),
                 "GWAS_SIGNAL_REPRESENTATIVE_VARIANT":m17_val(z,"BEST_GWAS_VARIANT"),
                 "QTL_SIGNAL_REPRESENTATIVE_VARIANT":m17_val(z,"BEST_QTL_VARIANT"),
                 "GWAS_LEAD_VARIANT_IDS":";".join(sorted({m17_val(q,"LEAD_IDS","LEAD_VARIANTS","TOP_VARIANT") for q in lead if m17_val(q,"LEAD_IDS","LEAD_VARIANTS","TOP_VARIANT")})),
                 "MOLECULAR_QTL_TYPE":m17_val(z,"BEST_QTL_TYPE"),
                 "QTL_TISSUE_CONTEXT":m17_val(z,"BEST_CONTEXT"),
                 "QTL_CONDITION":m17_val(z,"BEST_CONDITION"),
                 "MOLECULAR_TRAIT":m17_val(z,"BEST_MOLECULAR_TRAIT_ID"),
                 "BEST_QTL_DATASET_ID":m17_val(z,"BEST_DATASET_ID"),
                 "PRIOR_ROBUST_STRONG":m17_val(z,"PRIOR_ROBUST_STRONG"),
                 "N_CREDIBLE_SETS_IN_LOCUS":len(credible),
                 "CREDIBLE_SET_IDS":";".join(sorted({m17_val(c,"CREDIBLE_SET") for c in credible if m17_val(c,"CREDIBLE_SET")})),
                 "CREDIBLE_SET_TOP_VARIANTS":";".join(sorted({m17_val(c,"CS_TOP_VARIANT") for c in credible if m17_val(c,"CS_TOP_VARIANT")})),
                 "GENE_VARIANT_TOP_PIP_VARIANT":m17_val(z,"GENE_VARIANT_TOP_PIP_VARIANT"),
                 "GENE_VARIANT_MAX_PIP":m17_float(z.get("GENE_VARIANT_MAX_PIP")),
                 "LOCUS_TOP_PIP_VARIANT":m17_val(z,"TOP_PIP_VARIANT"),
                 "LOCUS_TOP_PIP_CONSEQUENCE":m17_val(z,"TOP_PIP_CONSEQUENCE"),
                 "LOCUS_TOP_PIP_GENE":m17_val(z,"TOP_PIP_SYMBOL"),
                 "GENE_VARIANT_CONSEQUENCES":m17_val(z,"GENE_VARIANT_CONSEQUENCES"),
                 "GENE_VARIANT_FUNCTIONAL_CATEGORIES":m17_val(z,"GENE_VARIANT_FUNCTIONAL_CATEGORIES"),
                 "LOCUS_VARIANTS_ANNOTATED":m17_float(z.get("ANNOT_N_VARIANTS")),
                 "LOCUS_CREDIBLE_VARIANTS_ANNOTATED":m17_float(z.get("ANNOT_N_IN_95_CREDIBLE_SET")),
                 "LOCUS_SPLICEAI_GE_020":m17_float(z.get("ANNOT_N_SPLICEAI_GE_0_20")),
                 "SOURCE_FILE":str(paths["S12"])}
            locus_rows.append(row)
            previous=best_locus_by_gene.get(gid)
            if previous is None or (np.isfinite(row["BEST_COLOC_PP_H4"]) and (
                    not np.isfinite(previous["BEST_COLOC_PP_H4"]) or
                    row["BEST_COLOC_PP_H4"]>previous["BEST_COLOC_PP_H4"])):
                best_locus_by_gene[gid]=row
    m17_write(out,"Step17_Genetic_Locus_Evidence.tsv",locus_rows)
    for row in mechanisms:
        strongest=best_locus_by_gene.get(row["GENE_ID"],{})
        for key in ("STUDY_ACCESSION","LOCUS_ID","BEST_COLOC_PP_H4","COLOC_CLASS",
            "GWAS_SIGNAL_REPRESENTATIVE_VARIANT","QTL_SIGNAL_REPRESENTATIVE_VARIANT",
            "GWAS_LEAD_VARIANT_IDS","MOLECULAR_QTL_TYPE","QTL_TISSUE_CONTEXT",
            "QTL_CONDITION","MOLECULAR_TRAIT","BEST_QTL_DATASET_ID",
            "N_CREDIBLE_SETS_IN_LOCUS","CREDIBLE_SET_TOP_VARIANTS",
            "GENE_VARIANT_TOP_PIP_VARIANT","GENE_VARIANT_MAX_PIP",
            "LOCUS_TOP_PIP_VARIANT","LOCUS_TOP_PIP_CONSEQUENCE","LOCUS_TOP_PIP_GENE",
            "GENE_VARIANT_CONSEQUENCES","GENE_VARIANT_FUNCTIONAL_CATEGORIES"):
            row["REPRESENTATIVE_"+key]=strongest.get(key,"")
        fdr=m17_float(row.get("GENE_SET_PERMUTATION_FDR"))
        p_fdr=m17_float(row.get("STEP14_PATHWAY_FDR"))
        gene=row.get("GENE_SYMBOL") or row.get("GENE_ID") or "unidentified candidate"
        narrative=(f"The {row.get('PATHWAY_SOURCE')} pathway '{row.get('PATHWAY_NAME')}' "
            f"({row.get('PATHWAY_ID')}) had Step14 enrichment FDR "
            f"{p_fdr:.4g} " if np.isfinite(p_fdr) else
            f"The {row.get('PATHWAY_SOURCE')} pathway '{row.get('PATHWAY_NAME')}' "
            f"({row.get('PATHWAY_ID')}) lacked a reported Step14 FDR ")
        narrative+=f"and included gene {gene} as {row.get('PATHWAY_RELATION')} "
        if strongest:
            narrative+=(f"in a genetically informed evidence chain. Example GWAS study "
              f"{strongest.get('STUDY_ACCESSION')}, locus {strongest.get('LOCUS_ID')}, "
              f"representative GWAS signal variant {strongest.get('GWAS_SIGNAL_REPRESENTATIVE_VARIANT') or 'unavailable'}, "
              f"QTL signal variant {strongest.get('QTL_SIGNAL_REPRESENTATIVE_VARIANT') or 'unavailable'}, "
              f"molecular modality {strongest.get('MOLECULAR_QTL_TYPE') or 'unavailable'} in "
              f"{strongest.get('QTL_TISSUE_CONTEXT') or 'unspecified context'}, "
              f"and coloc PP.H4 {strongest.get('BEST_COLOC_PP_H4')}.")
        else: narrative+=". No gene-specific Step12 locus record was resolved."
        if row.get("DATASET"):
            narrative+=(f" Step16 perturbation in {row.get('DATASET')} gave signed pathway NES "
              f"{m17_float(row.get('PATHWAY_NES')):.3g} and exploratory gene-set permutation "
              f"FDR {fdr:.3g}." if np.isfinite(fdr) else
              f" Step16 perturbation in {row.get('DATASET')} had no finite pathway permutation FDR.")
        else:
            narrative+=" No identical-source, identical-native-ID Step16 response for this gene-pathway pair was observed."
        narrative+=(f" Classification: {row.get('MECHANISM_TIER')}. Missing evidence: "
                   f"{row.get('MISSING_EVIDENCE') or 'none explicitly recorded'}. "
                   "A shared QTL signal and enrichment do not prove causal pathway mediation; "
                   "a descriptive NES is neither risk-allele direction nor clinical efficacy.")
        row["MECHANISM_NARRATIVE"] = narrative
    mech_df=m17_write(out,"Step17_Mechanism_Master.tsv",mechanisms)
    m17_write(out,"Step17_Mechanism_Narratives.tsv",[
        {"MECHANISM_ID":z["MECHANISM_ID"],"PATHWAY_SOURCE":z["PATHWAY_SOURCE"],
         "PATHWAY_ID":z["PATHWAY_ID"],"GENE_SYMBOL":z["GENE_SYMBOL"],
         "MECHANISM_TIER":z["MECHANISM_TIER"],"NARRATIVE":z["MECHANISM_NARRATIVE"]}
        for z in mechanisms])
    m17_write(out,"Step17_Experiment_Evidence.tsv",experiment_rows,
        columns=["PATHWAY_SOURCE","PATHWAY_ID","PATHWAY_NAME","GENE_ID","GENE_SYMBOL","DATASET",
                 "CELL_TYPE","NES","GENE_SET_NOMINAL_P","GENE_SET_PERMUTATION_FDR","EXPERIMENT_LEVEL_FDR",
                 "RESPONSE_DIRECTION","EVIDENCE_CLASS","TARGET_ACTIVITY_DIRECTION","TARGET_EXPRESSION_DIFF",
                 "TARGET_EXPRESSION_MEASURED","INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED",
                 "LEADING_EDGE_GENES","LIBRARY_NAME","LIBRARY_SHA256","PATHWAY_GENE_SET_SIZE",
                 "FDR_SCOPE","INTERSECTION_WITH_STEP14_DRIVER","MATCH_METHOD"])
    m17_write(out,"Step17_Library_Mapping_Audit.tsv",all_d["S16MAP"])
    m17_write(out,"Step17_Perturbation_Libraries.tsv",all_d["S16LIB"])
    m17_write(out,"Step17_All_Step16_Responses.tsv",all_d["S16"])
    # Retain perturbation-only GO/WP/KEGG gene sets separately, never promote
    # them to genetically implicated Step14 mechanisms on NES alone.
    perturb_only={}
    for row in raw_unmapped_responses:
        key=(m17_val(row,"PATHWAY_SOURCE"),m17_val(row,"PATHWAY_ID"),
             m17_val(row,"PATHWAY_KEY","PATHWAY_NAME"))
        if not key[0]:continue
        item=perturb_only.setdefault(key,{"PATHWAY_SOURCE":key[0],"PATHWAY_ID":key[1],"PATHWAY_KEY":key[2],
              "N_RESPONSE_ROWS":0,"PERTURBED_GENES":set(),"CELL_MODELS":set(),"FDR_VALUES":[]})
        item["N_RESPONSE_ROWS"]+=1
        item["PERTURBED_GENES"].add(m17_val(row,"GENE_SYMBOL"))
        item["CELL_MODELS"].add(m17_val(row,"CELL_TYPE"))
        val=m17_float(row.get("GENE_SET_PERMUTATION_FDR"))
        if np.isfinite(val):item["FDR_VALUES"].append(val)
    m17_write(out,"Step17_Perturbation_Only_Pathways.tsv",[
        {"PATHWAY_SOURCE":q["PATHWAY_SOURCE"],"PATHWAY_ID":q["PATHWAY_ID"],"PATHWAY_KEY":q["PATHWAY_KEY"],
         "N_RESPONSE_ROWS":q["N_RESPONSE_ROWS"],"PERTURBED_GENES":";".join(sorted(x for x in q["PERTURBED_GENES"] if x)),
         "CELL_MODELS":";".join(sorted(x for x in q["CELL_MODELS"] if x)),
         "BEST_GENE_SET_PERMUTATION_FDR":min(q["FDR_VALUES"]) if q["FDR_VALUES"] else np.nan,
         "GENETICALLY_ENRICHED_STEP14":"NOT_ESTABLISHED_OR_UNMAPPED",
         "EVIDENCE_SCOPE":"EXPLORATORY_EXPERIMENTAL_NOT_GWAS_PATHWAY_MECHANISM"}
        for q in perturb_only.values()])
    m17_write(out,"Step17_Pathway_Mapping_Audit.tsv",mapping_audit)
    m17_write(out,"Step17_Missing_Evidence.tsv",missing_rows)
    m17_write(out,"Step17_Final_Evidence_Matrix.tsv",mech_df)
    m17_write(out,"Step17_Pathway_Mechanism_Map.tsv",pathway_master)
    m17_write(out,"Step17_Evidence_Tiers.tsv",[
        {"MECHANISM_ID":r["MECHANISM_ID"],"PATHWAY_NAME":r["PATHWAY_NAME"],
         "GENE_ID":r["GENE_ID"],"GENE_SYMBOL":r["GENE_SYMBOL"],"MECHANISM_TIER":r["MECHANISM_TIER"],
         "DIRECTION_ALIGNMENT":r["DIRECTION_ALIGNMENT"],"MISSING_EVIDENCE":r["MISSING_EVIDENCE"]}
        for r in mechanisms])
    m17_write(out,"Step17_Unassigned_Genes.tsv",unassigned)
    gene_summary = m17_write(out,"Step17_Gene_Evidence_Summary.tsv",pd.DataFrame([
      {"GENE_ID":gid,"GENE_SYMBOL":gm.get(gid,gid),"MAX_H4":m17_float(gd.get("MAX_H4")),
       "STEP14_RANK":m17_float(gd.get("RANK_STEP14")),
       "N_PATHWAYS_WITH_DRIVER_MEMBERSHIP":sum(gid in p["DRIVER_IDS"] for p in pathways.values()),
       "IS_GWAS2M_PRIMARY_GENE":"YES" if gid in primary_gene_ids else "NO_OR_UNVERIFIED",
       "N_SCPERTURB_DATASETS":m17_float(exp_coverage.get(gid,{}).get("N_DATASETS")),
       "N_LINCS_METADATA_SIGNATURES":m17_float(lin_coverage.get(gid,{}).get("N_GENETIC_L1000_SIGNATURES"))}
      for gid,gd in sorted(gene_records.items())]))
    # Backward-compatible gene-facing ranking: derived from the Step14 order and
    # the strongest *heuristic* mechanism row, not a causal probability.
    best_mech = (mech_df.sort_values("MECHANISM_EVIDENCE_SCORE_HEURISTIC",ascending=False)
                 .drop_duplicates("GENE_ID") if not mech_df.empty else pd.DataFrame())
    if not best_mech.empty:
        usable_best=best_mech[["GENE_ID","MECHANISM_EVIDENCE_SCORE_HEURISTIC","MECHANISM_TIER"]]
        gene_rank=gene_summary.merge(usable_best,on="GENE_ID",how="left")
    else:
        gene_rank=gene_summary.copy()
        gene_rank["MECHANISM_EVIDENCE_SCORE_HEURISTIC"]=np.nan
        gene_rank["MECHANISM_TIER"]="MT1_UNRESOLVED_NO_PATHWAY"
    gene_rank=gene_rank.sort_values(["STEP14_RANK","GENE_SYMBOL"],na_position="last").reset_index(drop=True)
    gene_rank.insert(0,"STEP17_RANK",np.arange(1,len(gene_rank)+1))
    gene_rank["RANK_INTERPRETATION"]="STEP14_GENE_RANK_RETAINED_NOT_EXPERIMENTAL_CAUSALITY"
    m17_write(out,"Step17_Final_Target_Ranking.tsv",gene_rank)
    # Preserve the source long evidence with no invented effect alleles or beta signs.
    if not all_d["S12LONG"].empty:
        m17_write(out,"Step17_Source_Molecular_Evidence_Long.tsv",all_d["S12LONG"])
    # Machine-readable QC: exact mapping coverage and known evidence deficiencies.
    qc_rows=[
      {"CHECK":"STEP14_SIGNIFICANT_PATHWAYS_PRESERVED", "STATUS":"PASS" if len(pathways)==len(pathway_master) else "FAIL",
       "VALUE":len(pathway_master),"DETAIL":"One row per source-native pathway identity"},
      {"CHECK":"STEP16_SOURCE_VERSION", "STATUS":"PASS","VALUE":summary16.get("version"),
       "DETAIL":"Recorded in Step16_summary.json"},
      {"CHECK":"STEP16_DESCRIPTIVE_CONTRASTS", "STATUS":"EXPLORATORY","VALUE":summary16.get("n_successful_descriptive_contrasts",0),
       "DETAIL":"Not replicated causal experiments"},
      {"CHECK":"SOURCE_ONTOLOGY_CONSISTENCY", "STATUS":"STRICT","VALUE":len(experiment_rows),
       "DETAIL":"Strict same-source native GO:BP / WikiPathways / KEGG ID; no name matching"},
      {"CHECK":"STEP16_SOURCE_NATIVE_ID_MAPPABILITY","STATUS":"AUDIT",
       "VALUE":summary16.get("n_step14_pathways_exactly_mappable",0),
       "DETAIL":"Count of Step14 terms with source-native IDs identified in Step16 v3 libraries"},
      {"CHECK":"STEP16_UNMAPPED_RESPONSE_ROWS", "STATUS":"INFORMATIONAL","VALUE":len(raw_unmapped_responses),
       "DETAIL":"All raw experiment rows preserved separately"},
      {"CHECK":"GWAS_QTL_ALLELE_HARMONIZATION", "STATUS":"UNRESOLVED" if not risk else "PARTIALLY_AVAILABLE",
       "VALUE":len(risk),"DETAIL":"H4 and NES do not indicate GWAS risk-allele direction"},
      {"CHECK":"EXPERIMENT_REPLICATE_AWARE_FDR", "STATUS":"NOT_AVAILABLE" if not any(x["MECHANISM_TIER"].startswith(("MT3", "MT4")) for x in mechanisms) else "PARTIALLY_AVAILABLE",
       "VALUE":sum(x["MECHANISM_TIER"].startswith(("MT3", "MT4")) for x in mechanisms),
       "DETAIL":"Gene-set permutation FDR is distinct from replicate-aware experiment FDR"},
    ]
    m17_write(out,"Step17_QC_Audit.tsv",qc_rows)
    source_manifest=[]
    for k,v in paths.items():
        source_manifest.append({"INPUT":k,"PATH":str(v),"EXISTS":v.is_file(),
          "SIZE_BYTES":v.stat().st_size if v.is_file() else 0,
          "MODIFIED_EPOCH":v.stat().st_mtime if v.is_file() else np.nan})
    m17_write(out,"Step17_Source_Manifest.tsv",source_manifest)
    for name in ("Step17_Pathway_Master", "Step17_Mechanism_Master"):
        source=out/(name+".tsv")
        pd.read_csv(source,sep="\t",low_memory=False).to_csv(out/(name+".csv"),index=False)
    summary={"version":M17_VERSION,"status":"COMPLETE_EVIDENCE_INTEGRATION_NOT_CAUSAL_VALIDATION",
      "n_unique_pathways":len(pathways),"n_ranked_genes":len(gene_records),
      "n_primary_gene_ids":len(primary_gene_ids),
      "n_genes_with_pathway_membership":len(genes_on_paths),
      "n_mechanism_rows":len(mech_df),"n_step16_exploratory_responses_mapped":len(experiment_rows),
      "n_step16_descriptive_contrasts":summary16.get("n_successful_descriptive_contrasts"),
      "n_mechanism_rows_MT1":sum(x["MECHANISM_TIER"].startswith("MT1") for x in mechanisms),
      "n_mechanism_rows_MT2":sum(x["MECHANISM_TIER"].startswith("MT2") for x in mechanisms),
      "n_mechanism_rows_MT3":sum(x["MECHANISM_TIER"].startswith("MT3") for x in mechanisms),
      "n_mechanism_rows_MT4":sum(x["MECHANISM_TIER"].startswith("MT4") for x in mechanisms),
      "n_unmapped_step16_pathway_rows":len(raw_unmapped_responses),
      "n_cross_ontology_name_only_unaccepted":sum(x["MATCH_STATUS"]=="CROSS_ONTOLOGY_NAME_ONLY_NOT_ACCEPTED" for x in mapping_audit),
      "n_step16_native_id_missing":sum(x["MATCH_STATUS"]=="STEP16_PATHWAY_NATIVE_ID_MISSING" for x in mapping_audit),
      "n_step16_same_source_name_only_unaccepted":sum(x["MATCH_STATUS"]=="SAME_SOURCE_NAME_ONLY_NOT_ACCEPTED" for x in mapping_audit),
      "n_step16_matched_unique_pathway_ids":len({(r["PATHWAY_SOURCE"],r["PATHWAY_ID"]) for r in experiment_rows}),
      "n_step16_matched_unique_genes":len({r["GENE_ID"] for r in experiment_rows}),
      "n_perturbation_only_unique_gene_sets":len(perturb_only),
      "n_missing_evidence_rows":len(missing_rows),
      "score_interpretation":"Heuristic evidence display, not a calibrated probability, odds ratio, statistical p value, or causal certainty",
      "direction_interpretation":"Allele-harmonised QTL effects and verified target engagement required for mechanistic alignment",
      "experiment_interpretation":"GSEA gene-set permutation FDR is not independent biological-replicate-level FDR",
      "no_html":True}
    tmp=out/(".Step17_summary.json.tmp."+str(_mos.getpid()))
    tmp.write_text(json.dumps(summary,indent=2,allow_nan=False,default=str)+"\n")
    _mos.replace(tmp,out/"Step17_summary.json")
    print("[COMPLETE] "+json.dumps({k:v for k,v in summary.items() if isinstance(v,(int,float))},indent=2))
    if not args.compact_screen:
        print("[NOTE] Same source-native pathway ID does not establish identical gene-set membership across releases.")
        print("[NOTE] No genetic risk-to-pathway direction can be inferred from NES alone.")
        print("[NOTE] MT2 is exploratory, even if gene-set permutation FDR < 0.05.")
    print("[OUTPUT]",out,flush=True)


if __name__ == "__main__":
    pipeline_main(mechanism_main)
