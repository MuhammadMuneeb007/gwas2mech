#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================================
GWAS2m STEP 14 v3.0
PUBLICATION-GRADE GENE PRIORITISATION + PATHWAY ENRICHMENT + OT COMPARISON
================================================================================

GWAS2Mech PIPELINE USAGE (same flags in Steps 12-17)
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR                       # run on this machine
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR --slurm               # or write a one-job SLURM script
    sbatch Step14_Pathways_parkinson_s_disease_european.sh
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR --stats
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR --index 2             # study 2 only -> <ancestry>/<GCST>/
    # 1000 Genomes 30x LD panel: every input / output folder is <ancestry>_1000G30x
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --slurm
    sbatch Step14_Pathways_parkinson_s_disease_european_1000G30x.sh
    envs/pipeline/bin/python Step14_Pathway_Analysis.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

WHY THIS STEP EXISTS
--------------------
Step12/Step13 contain many rows per biological gene because the same gene may
appear in:
    * multiple GWAS studies
    * multiple loci/signals
    * multiple tissues
    * multiple QTL modalities
    * multiple molecular traits

Those repeated rows MUST NOT be used as repeated pathway votes.

Step14 therefore performs:

    Step12 locus-gene evidence
        -> deduplicated gene evidence
        -> transparent evidence ranking
        -> biologically defined gene sets
        -> pathway enrichment
        -> pathway redundancy/driver-gene audit
        -> Open Targets pathway comparison
        -> complete terminal report

SCIENTIFIC PRINCIPLES
---------------------
1. One gene = one pathway-enrichment vote.
2. Tissue/QTL recurrence is evidence about a gene, not extra gene counts.
3. The main set is CORE_STRONG: MAX_H4 >= 0.80.
4. EXPANDED_SUGGESTIVE (MAX_H4 >= 0.50) is a sensitivity analysis.
5. Open Targets evidence is NOT used to rank GWAS2m genes. It remains an
   independent external comparator.
6. Pathways are considered significant after FDR correction.
7. A pathway needs at least --min-overlap distinct input genes (default 2).
8. Very tiny / extremely broad terms can be filtered using term-size limits.
9. An enrichment result is not proof of pathway causality.
10. A GWAS2m-only pathway is not automatically a novel pathway.

DEFAULT INPUTS
--------------
12_master_table/<phenotype>/<ancestry>/Step12_Locus_Gene_Master.tsv
13_opentargets_comparison/<phenotype>/<ancestry>/Step13_Gene_Comparison.tsv

DEFAULT OUTPUTS
---------------
14_pathways/<phenotype>/<ancestry>/

    Step14_Gene_Ranking.tsv
    Step14_Gene_Tissue_Evidence.tsv
    Step14_Gene_CellType_Evidence.tsv
    Step14_Gene_Context_Evidence.tsv
    Step14_Gene_Modality_Evidence.tsv
    Step14_Locus_Top_Genes.tsv
    Step14_Gene_Set_Membership.tsv
    Step14_Publication_Priority_Genes.tsv

    Step14_Pathway_All_Results.tsv
    Step14_Pathway_Significant.tsv
    Step14_Pathway_Set_Summary.tsv
    Step14_Pathway_Gene_Drivers.tsv
    Step14_Pathway_Redundancy.tsv

    Step14_OpenTargets_Pathway_Comparison.tsv
    Step14_OpenTargets_Pathway_Comparison_Summary.tsv

    Step14_QC_Audit.tsv
    Step14_summary.json

INSTALL
-------
Online enrichment:
    python -m pip install gprofiler-official pandas numpy scipy statsmodels

Offline GMT enrichment:
    python -m pip install pandas numpy scipy statsmodels

RUN
---
Normal:
    python Step14_Pathway_Analysis.py \
      --phenotype "parkinson's disease" \
      --ancestry EUR

Build/rank gene sets only:
    python Step14_Pathway_Analysis.py \
      --phenotype "parkinson's disease" \
      --ancestry EUR \
      --skip-enrichment

Offline GMT:
    python Step14_Pathway_Analysis.py \
      --phenotype "parkinson's disease" \
      --ancestry EUR \
      --gmt resources/pathways/Reactome.gmt \
      --gmt resources/pathways/GO_BP.gmt

Less terminal output:
    ... --compact-screen

================================================================================
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from scipy.stats import hypergeom
except Exception:
    hypergeom = None

try:
    from statsmodels.stats.multitest import multipletests
except Exception:
    multipletests = None



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
    "NAME": "Step14", "TITLE": "gene prioritisation and pathway enrichment", "FOLDER": "14_pathways",
    "SCRIPT_PREFIX": "Step14_Pathways", "SUMMARY": "Step14_summary.json",
    "PREVIOUS": "Step12_Master_Table.py (and Step13_Compare_OpenTargets.py for the Open Targets comparison)",
    "INPUTS": [
        ("Step12", "12_master_table", ["Step12_Locus_Gene_Master.tsv"], True, True),
        ("Step13", "13_opentargets_comparison", ["Step13_Gene_Comparison.tsv"], False, True),
    ],
    "MODULES": [("scipy", True, "hypergeometric test (GMT mode)"), ("statsmodels", True, "FDR"),
                ("gprofiler", False, "g:Profiler online enrichment (default backend)")],
    "STATS_TABLES": [
        ("Step14_Gene_Ranking.tsv", ["EFFECTOR_GENE_SYMBOL", "EFFECTOR_GENE_ID", "MAX_H4", "PRIOR_ROBUST",
                                     "BEST_QTL_TYPE", "BEST_CONTEXT", "N_STUDIES", "PUBLICATION_PRIORITY_CLASS",
                                     "EVIDENCE_SCORE_STEP14"], "EVIDENCE_SCORE_STEP14", "GENE RANKING"),
        ("Step14_Pathway_Significant.tsv", ["GENE_SET", "SOURCE", "TERM_ID", "TERM_NAME", "FDR", "INTERSECTION_SIZE",
                                            "INTERSECTION_GENE_LABELS"], "", "SIGNIFICANT PATHWAYS"),
    ],
}


def frame_extra_checks(ns, P) -> list:
    if "--skip-enrichment" in _fsys.argv:
        return [("pathway enrichment", True, "skipped (--skip-enrichment)", False)]
    gmts = [_fsys.argv[i + 1] for i, a in enumerate(_fsys.argv[:-1]) if a == "--gmt"]
    if gmts:
        return [(f"GMT {g}", _FPath(g).is_file(), "offline gene-set library", True) for g in gmts]
    import importlib.util as _fimp
    ok_module = _fimp.find_spec("gprofiler") is not None
    reachable, detail = False, ""
    try:
        import urllib.request as _fu
        with _fu.urlopen("https://biit.cs.ut.ee/gprofiler/api/util/organisms_list/", timeout=20) as r:
            reachable, detail = r.status == 200, f"HTTP {r.status}"
    except Exception as exc:
        detail = f"{type(exc).__name__}: {str(exc)[:80]}"
    return [("g:Profiler (online)", ok_module and reachable,
             (detail if ok_module else "pip install gprofiler-official") +
             " - or run offline with --gmt <library.gmt>", True)]


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


VERSION = "3.0.0"

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

QTL_TYPES = ["eQTL", "sQTL", "isoQTL", "exonQTL", "pQTL"]

DEFAULT_GPROFILER_SOURCES = [
    "GO:BP",
    "REAC",
    "KEGG",
    "WP",
]


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
            f"Unsupported ancestry {value!r}. Use EUR, AFR, EAS, SAS or AMR."
        )
    return ANCESTRY[key]


def section(title: str, char: str = "="):
    print()
    print(char * 150)
    print(title)
    print(char * 150)


def subsection(title: str):
    section(title, "-")


def show(
    df: pd.DataFrame | None,
    title: str | None = None,
    *,
    max_rows: int | None = None,
    columns: list[str] | None = None,
):
    if title:
        section(title)

    if df is None or df.empty:
        print("(no rows)")
        return

    x = df.copy()

    if columns:
        columns = [c for c in columns if c in x.columns]
        x = x[columns]

    if max_rows is not None:
        x = x.head(max_rows)

    with pd.option_context(
        "display.max_rows", None,
        "display.max_columns", None,
        "display.width", 1000,
        "display.max_colwidth", 140,
        "display.expand_frame_repr", False,
        "display.float_format", lambda v: f"{v:.6g}",
    ):
        print(x.to_string(index=False))

    if max_rows is not None and len(x) < len(df):
        print(f"\n[showing {len(x):,}/{len(df):,} rows]")


def first_existing(columns: Iterable[str], candidates: Iterable[str]):
    columns = list(columns)
    lower = {str(c).lower(): c for c in columns}

    for candidate in candidates:
        if candidate in columns:
            return candidate

        found = lower.get(str(candidate).lower())
        if found is not None:
            return found

    return None


def safe_numeric(values):
    return pd.to_numeric(values, errors="coerce")


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


def clean_gene_id(value) -> str:
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass

    x = str(value).strip()

    if x.lower() in {"", "nan", "none", "null", "na", "."}:
        return ""

    m = re.search(r"(ENSG\d+)", x, flags=re.I)

    return m.group(1).upper() if m else x


def clean_symbol(value) -> str:
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
        and x.strip().lower() not in {"nan", "none", "null", "na"}
    ]


def join_unique(values) -> str:
    output = []

    for value in values:
        for token in split_tokens(value):
            if token not in output:
                output.append(token)

    return ";".join(output)


def safe_max(values):
    x = safe_numeric(values)
    return x.max() if x.notna().any() else np.nan


def safe_min(values):
    x = safe_numeric(values)
    return x.min() if x.notna().any() else np.nan


def safe_median(values):
    x = safe_numeric(values)
    return x.median() if x.notna().any() else np.nan


def bh_fdr(pvalues):
    p = safe_numeric(pd.Series(pvalues))

    out = pd.Series(np.nan, index=p.index, dtype=float)

    good = p.notna()

    if not good.any():
        return out

    if multipletests is None:
        # Manual Benjamini-Hochberg fallback.
        vals = p[good].to_numpy(dtype=float)
        order = np.argsort(vals)
        ranked = vals[order]
        n = len(ranked)
        adjusted = ranked * n / np.arange(1, n + 1)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        adjusted = np.clip(adjusted, 0, 1)
        inverse = np.empty_like(order)
        inverse[order] = np.arange(n)
        out.loc[good] = adjusted[inverse]
        return out

    out.loc[good] = multipletests(
        p[good].to_numpy(dtype=float),
        method="fdr_bh",
    )[1]

    return out


def jaccard(a, b):
    a = set(a)
    b = set(b)

    if not a and not b:
        return np.nan

    union = a | b

    return len(a & b) / len(union) if union else np.nan


def pair_key(study, locus):
    return f"{study}|{locus}"


# =============================================================================
# INPUT LOADING
# =============================================================================

def load_step12(path: Path):
    if not path.exists():
        raise SystemExit(
            "Step12 locus-gene master is missing:\n"
            f"  {path}"
        )

    x = pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
    )

    required = [
        "STUDY_ACCESSION",
        "LOCUS_ID",
        "EFFECTOR_GENE_ID",
        "BEST_H4",
    ]

    missing = [
        col
        for col in required
        if col not in x.columns
    ]

    if missing:
        raise SystemExit(
            "Step12 locus-gene master is missing required columns:\n"
            f"  {missing}"
        )

    x["EFFECTOR_GENE_ID"] = (
        x["EFFECTOR_GENE_ID"]
        .map(clean_gene_id)
    )

    if "EFFECTOR_GENE_SYMBOL" not in x.columns:
        x["EFFECTOR_GENE_SYMBOL"] = ""

    x["EFFECTOR_GENE_SYMBOL"] = (
        x["EFFECTOR_GENE_SYMBOL"]
        .map(clean_symbol)
    )

    x["BEST_H4"] = safe_numeric(
        x["BEST_H4"]
    )

    x["_LOCUS_STUDY_KEY"] = [
        pair_key(study, locus)
        for study, locus
        in zip(
            x["STUDY_ACCESSION"],
            x["LOCUS_ID"],
        )
    ]

    x = x[
        x["EFFECTOR_GENE_ID"]
        .astype(str)
        .str.len()
        >
        0
    ].copy()

    return x


def load_step13(path: Path):
    if not path.exists():
        return pd.DataFrame()

    try:
        x = pd.read_csv(
            path,
            sep="\t",
            low_memory=False,
        )
    except Exception:
        return pd.DataFrame()

    if "EFFECTOR_GENE_ID" in x.columns:
        x["EFFECTOR_GENE_ID"] = (
            x["EFFECTOR_GENE_ID"]
            .map(clean_gene_id)
        )

    return x


# =============================================================================
# STEP12 -> DEDUPLICATED GENE EVIDENCE
# =============================================================================

def qtl_type_column(df):
    return first_existing(
        df.columns,
        [
            "BEST_QTL_TYPE",
            "QTL_TYPE",
        ],
    )


def context_column(df):
    """Best available biological context (bulk tissue or cell type)."""
    return first_existing(
        df.columns,
        [
            "BEST_CONTEXT",
            "CONTEXT_NAME",
            "BEST_TISSUE",
            "TISSUE",
        ],
    )


def tissue_column(df):
    """Bulk-tissue context only when the Step12 v2 schema is available."""
    return first_existing(
        df.columns,
        [
            "BULK_BEST_CONTEXT",
            "BEST_TISSUE",
            "TISSUE",
        ],
    )


def cell_column(df):
    return first_existing(
        df.columns,
        [
            "CELL_BEST_CONTEXT",
            "CELL_TYPE",
        ],
    )


def prior_robust_column(df):
    return first_existing(
        df.columns,
        [
            "PRIOR_ROBUST_STRONG",
            "PRIOR_ROBUST",
        ],
    )


def supported_modality_columns(df):
    result = {}

    for qtl in QTL_TYPES:
        candidates = [
            f"{qtl.upper()}_SUPPORTED",
            f"{qtl}_SUPPORTED",
        ]

        found = first_existing(
            df.columns,
            candidates,
        )

        if found:
            result[qtl] = found

    return result


def make_symbol_lookup(step12):
    counts = defaultdict(lambda: defaultdict(int))

    for row in step12.itertuples(index=False):
        gid = clean_gene_id(
            getattr(
                row,
                "EFFECTOR_GENE_ID",
                "",
            )
        )

        symbol = clean_symbol(
            getattr(
                row,
                "EFFECTOR_GENE_SYMBOL",
                "",
            )
        )

        if gid and symbol:
            counts[gid][symbol] += 1

    result = {}

    for gid, options in counts.items():
        result[gid] = sorted(
            options.items(),
            key=lambda kv: (-kv[1], kv[0]),
        )[0][0]

    return result


def build_context_evidence(step12):
    """Collapse Step12 v2 bulk and cell-type evidence into explicit contexts.

    One output row = gene x evidence source x context x condition.  This keeps
    cell-type evidence distinct from bulk tissue evidence while still ensuring
    that a gene contributes only once to any pathway gene set.
    """
    rows = []

    source_specs = [
        (
            "STEP09_BULK",
            "BULK_BEST_H4",
            "BULK_BEST_CONTEXT",
            "BULK_BEST_CONDITION",
            "BULK_BEST_QTL_TYPE",
        ),
        (
            "STEP10_CELLTYPE",
            "CELL_BEST_H4",
            "CELL_BEST_CONTEXT",
            "CELL_BEST_CONDITION",
            "CELL_BEST_QTL_TYPE",
        ),
    ]

    for source, h4_col, context_col, condition_col, qtl_col in source_specs:
        if h4_col not in step12.columns or context_col not in step12.columns:
            continue

        for row in step12.itertuples(index=False):
            gid = clean_gene_id(getattr(row, "EFFECTOR_GENE_ID", ""))
            context = str(getattr(row, context_col, "") or "").strip()
            condition = str(getattr(row, condition_col, "") or "").strip() if condition_col in step12.columns else ""
            qtl = str(getattr(row, qtl_col, "") or "").strip() if qtl_col in step12.columns else ""
            h4 = pd.to_numeric(pd.Series([getattr(row, h4_col, np.nan)]), errors="coerce").iloc[0]

            if not gid or not context or context.lower() in {"nan", "none"} or pd.isna(h4):
                continue

            rows.append(
                {
                    "EFFECTOR_GENE_ID": gid,
                    "EFFECTOR_GENE_SYMBOL": clean_symbol(getattr(row, "EFFECTOR_GENE_SYMBOL", "")),
                    "EVIDENCE_SOURCE": source,
                    "CONTEXT": context,
                    "CONDITION": "" if condition.lower() in {"nan", "none"} else condition,
                    "QTL_TYPE": "" if qtl.lower() in {"nan", "none"} else qtl,
                    "H4": h4,
                    "STUDY_ACCESSION": str(getattr(row, "STUDY_ACCESSION", "")),
                    "LOCUS_ID": str(getattr(row, "LOCUS_ID", "")),
                }
            )

    if not rows:
        # Legacy fallback: retain the original best-context representation.
        ccol = context_column(step12)
        qcol = qtl_type_column(step12)
        if not ccol:
            return pd.DataFrame()
        for row in step12.itertuples(index=False):
            gid = clean_gene_id(getattr(row, "EFFECTOR_GENE_ID", ""))
            context = str(getattr(row, ccol, "") or "").strip()
            h4 = pd.to_numeric(pd.Series([getattr(row, "BEST_H4", np.nan)]), errors="coerce").iloc[0]
            if not gid or not context or context.lower() in {"nan", "none"} or pd.isna(h4):
                continue
            rows.append(
                {
                    "EFFECTOR_GENE_ID": gid,
                    "EFFECTOR_GENE_SYMBOL": clean_symbol(getattr(row, "EFFECTOR_GENE_SYMBOL", "")),
                    "EVIDENCE_SOURCE": str(getattr(row, "BEST_EVIDENCE_SOURCE", "LEGACY")),
                    "CONTEXT": context,
                    "CONDITION": str(getattr(row, "BEST_CONDITION", "")),
                    "QTL_TYPE": str(getattr(row, qcol, "")) if qcol else "",
                    "H4": h4,
                    "STUDY_ACCESSION": str(getattr(row, "STUDY_ACCESSION", "")),
                    "LOCUS_ID": str(getattr(row, "LOCUS_ID", "")),
                }
            )

    raw = pd.DataFrame(rows)
    if raw.empty:
        return raw

    out_rows = []
    group_cols = ["EFFECTOR_GENE_ID", "EVIDENCE_SOURCE", "CONTEXT", "CONDITION"]
    for key, group in raw.groupby(group_cols, dropna=False, sort=False):
        gid, source, context, condition = key
        out_rows.append(
            {
                "EFFECTOR_GENE_ID": clean_gene_id(gid),
                "EFFECTOR_GENE_SYMBOL": join_unique(group["EFFECTOR_GENE_SYMBOL"]),
                "EVIDENCE_SOURCE": source,
                "CONTEXT": context,
                "CONDITION": condition,
                "MAX_H4": safe_max(group["H4"]),
                "MEDIAN_H4": safe_median(group["H4"]),
                "N_LOCUS_GENE_ROWS": len(group),
                "N_STUDIES": group["STUDY_ACCESSION"].astype(str).nunique(),
                "N_LOCUS_STUDY_PAIRS": pd.Series(
                    [pair_key(s, l) for s, l in zip(group["STUDY_ACCESSION"], group["LOCUS_ID"])]
                ).nunique(),
                "QTL_TYPES": join_unique(group["QTL_TYPE"]),
            }
        )

    out = pd.DataFrame(out_rows)
    return out.sort_values(
        ["EFFECTOR_GENE_ID", "EVIDENCE_SOURCE", "MAX_H4"],
        ascending=[True, True, False],
        kind="stable",
    )


def build_tissue_evidence(step12):
    contexts = build_context_evidence(step12)
    if contexts.empty:
        return contexts
    out = contexts[contexts["EVIDENCE_SOURCE"].eq("STEP09_BULK")].copy()
    if out.empty:
        return out
    return out.rename(columns={"CONTEXT": "TISSUE"})


def build_cell_evidence(step12):
    contexts = build_context_evidence(step12)
    if contexts.empty:
        return contexts
    out = contexts[contexts["EVIDENCE_SOURCE"].eq("STEP10_CELLTYPE")].copy()
    if out.empty:
        return out
    return out.rename(columns={"CONTEXT": "CELL_TYPE"})


def build_modality_evidence(step12):
    qcol = qtl_type_column(step12)
    ccol = context_column(step12)

    if not qcol:
        return pd.DataFrame()

    rows = []

    for (
        gid,
        qtl,
    ), group in step12.groupby(
        [
            "EFFECTOR_GENE_ID",
            qcol,
        ],
        dropna=False,
        sort=False,
    ):
        contexts = (
            sorted(
                {
                    str(v).strip()
                    for v in group[ccol]
                    if str(v).strip()
                    and str(v).lower() not in {"nan", "none"}
                }
            )
            if ccol else []
        )

        rows.append(
            {
                "EFFECTOR_GENE_ID": clean_gene_id(gid),
                "EFFECTOR_GENE_SYMBOL": join_unique(group["EFFECTOR_GENE_SYMBOL"]),
                "QTL_TYPE": str(qtl),
                "MAX_H4": safe_max(group["BEST_H4"]),
                "MEDIAN_H4": safe_median(group["BEST_H4"]),
                "N_LOCUS_GENE_ROWS": len(group),
                "N_STUDIES": group["STUDY_ACCESSION"].astype(str).nunique(),
                "N_LOCUS_STUDY_PAIRS": group["_LOCUS_STUDY_KEY"].astype(str).nunique(),
                "N_CONTEXTS": len(contexts),
                "CONTEXTS": ";".join(contexts),
                # Legacy aliases retained so older downstream notebooks do not break.
                "N_TISSUES": len(contexts),
                "TISSUES": ";".join(contexts),
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(
            ["EFFECTOR_GENE_ID", "MAX_H4"],
            ascending=[True, False],
        )
    return out


def build_gene_ranking(step12):
    qcol = qtl_type_column(step12)
    prior_col = prior_robust_column(step12)
    supported_cols = supported_modality_columns(step12)
    symbol_lookup = make_symbol_lookup(step12)

    rows = []

    def numeric_max(group, col):
        if col not in group.columns:
            return np.nan
        return safe_max(group[col])

    def best_row(group, col):
        if col not in group.columns:
            return None
        vals = safe_numeric(group[col])
        if not vals.notna().any():
            return None
        return group.loc[vals.idxmax()]

    for gid, group in step12.groupby("EFFECTOR_GENE_ID", sort=False):
        gid = clean_gene_id(gid)
        qtls = set()

        if qcol:
            for value in group[qcol]:
                for token in split_tokens(value):
                    qtls.add(token)

        for qtl, col in supported_cols.items():
            if group[col].map(boolish).any():
                qtls.add(qtl)

        max_h4 = safe_max(group["BEST_H4"])
        median_h4 = safe_median(group["BEST_H4"])
        overall_best = best_row(group, "BEST_H4")

        prior_robust = (
            bool(group[prior_col].map(boolish).any())
            if prior_col else False
        )

        n_locus_study = group["_LOCUS_STUDY_KEY"].astype(str).nunique()
        n_studies = group["STUDY_ACCESSION"].astype(str).nunique()
        n_qtl = len(qtls)

        bulk_h4 = numeric_max(group, "BULK_BEST_H4")
        cell_h4 = numeric_max(group, "CELL_BEST_H4")
        bulk_best = best_row(group, "BULK_BEST_H4")
        cell_best = best_row(group, "CELL_BEST_H4")

        bulk_contexts = sorted({
            str(v).strip() for v in group.get("BULK_BEST_CONTEXT", pd.Series(dtype=str))
            if str(v).strip() and str(v).lower() not in {"nan", "none"}
        })
        cell_contexts = sorted({
            str(v).strip() for v in group.get("CELL_BEST_CONTEXT", pd.Series(dtype=str))
            if str(v).strip() and str(v).lower() not in {"nan", "none"}
        })

        max_pip = numeric_max(group, "GENE_VARIANT_MAX_PIP")
        max_spliceai = numeric_max(group, "GENE_VARIANT_MAX_SPLICEAI")
        max_pangolin = numeric_max(group, "GENE_VARIANT_MAX_PANGOLIN")
        pip_best = best_row(group, "GENE_VARIANT_MAX_PIP")
        splice_best = best_row(group, "GENE_VARIANT_MAX_SPLICEAI")
        pangolin_best = best_row(group, "GENE_VARIANT_MAX_PANGOLIN")

        # Core GWAS2m score remains independent of Open Targets.
        score = (
            (0.0 if pd.isna(max_h4) else float(max_h4))
            + 0.15 * int(prior_robust)
            + 0.08 * max(0, min(n_qtl - 1, 4))
            + 0.06 * math.log1p(max(0, n_locus_study - 1))
            + 0.04 * math.log1p(max(0, n_studies - 1))
        )

        bulk_strong = pd.notna(bulk_h4) and float(bulk_h4) >= 0.80
        cell_strong = pd.notna(cell_h4) and float(cell_h4) >= 0.80

        rows.append(
            {
                "EFFECTOR_GENE_ID": gid,
                "EFFECTOR_GENE_SYMBOL": symbol_lookup.get(gid, ""),
                "MAX_H4": max_h4,
                "MEDIAN_H4": median_h4,
                "PRIOR_ROBUST": "YES" if prior_robust else "NO",
                "BEST_EVIDENCE_SOURCE": clean_symbol(overall_best.get("BEST_EVIDENCE_SOURCE", "")) if overall_best is not None else "",
                "BEST_QTL_TYPE": clean_symbol(overall_best.get("BEST_QTL_TYPE", "")) if overall_best is not None else "",
                "BEST_CONTEXT": str(overall_best.get("BEST_CONTEXT", "")) if overall_best is not None else "",
                "BEST_CONDITION": str(overall_best.get("BEST_CONDITION", "")) if overall_best is not None else "",
                "BEST_GWAS_VARIANT": str(overall_best.get("BEST_GWAS_VARIANT", "")) if overall_best is not None else "",
                "BEST_QTL_VARIANT": str(overall_best.get("BEST_QTL_VARIANT", "")) if overall_best is not None else "",
                "BULK_MAX_H4": bulk_h4,
                "BULK_BEST_CONTEXT": str(bulk_best.get("BULK_BEST_CONTEXT", "")) if bulk_best is not None else "",
                "BULK_BEST_QTL_TYPE": str(bulk_best.get("BULK_BEST_QTL_TYPE", "")) if bulk_best is not None else "",
                "N_BULK_CONTEXTS": len(bulk_contexts),
                "BULK_CONTEXTS": ";".join(bulk_contexts),
                "CELL_MAX_H4": cell_h4,
                "CELL_BEST_CONTEXT": str(cell_best.get("CELL_BEST_CONTEXT", "")) if cell_best is not None else "",
                "CELL_BEST_CONDITION": str(cell_best.get("CELL_BEST_CONDITION", "")) if cell_best is not None else "",
                "CELL_BEST_QTL_TYPE": str(cell_best.get("CELL_BEST_QTL_TYPE", "")) if cell_best is not None else "",
                "N_CELL_CONTEXTS": len(cell_contexts),
                "CELL_CONTEXTS": ";".join(cell_contexts),
                "BULK_STRONG": "YES" if bulk_strong else "NO",
                "CELL_STRONG": "YES" if cell_strong else "NO",
                "BULK_CELL_STRONG": "YES" if bulk_strong and cell_strong else "NO",
                "MAX_PIP": max_pip,
                "TOP_PIP_VARIANT": str(pip_best.get("GENE_VARIANT_TOP_PIP_VARIANT", "")) if pip_best is not None else "",
                "MAX_SPLICEAI": max_spliceai,
                "TOP_SPLICEAI_VARIANT": str(splice_best.get("GENE_VARIANT_TOP_SPLICEAI_VARIANT", "")) if splice_best is not None else "",
                "MAX_PANGOLIN": max_pangolin,
                "TOP_PANGOLIN_VARIANT": str(pangolin_best.get("GENE_VARIANT_TOP_PANGOLIN_VARIANT", "")) if pangolin_best is not None else "",
                "N_QTL_TYPES": n_qtl,
                "QTL_TYPES": ";".join(sorted(qtls)),
                "N_TISSUES": len(bulk_contexts),
                "TISSUES": ";".join(bulk_contexts),
                "N_CELL_TYPES": len(cell_contexts),
                "CELL_TYPES": ";".join(cell_contexts),
                "N_STUDIES": n_studies,
                "STUDIES": join_unique(group["STUDY_ACCESSION"]),
                "N_INDEPENDENT_LOCUS_STUDY_PAIRS": n_locus_study,
                "LOCUS_STUDY_PAIRS": join_unique(group["_LOCUS_STUDY_KEY"]),
                "N_LOCUS_GENE_ROWS": len(group),
                "N_H4_GE_0_80_ROWS": int((safe_numeric(group["BEST_H4"]) >= 0.80).sum()),
                "N_H4_GE_0_50_ROWS": int((safe_numeric(group["BEST_H4"]) >= 0.50).sum()),
                "EVIDENCE_SCORE_STEP14": score,
            }
        )

    genes = pd.DataFrame(rows)
    if genes.empty:
        return genes

    genes = genes.sort_values(
        [
            "MAX_H4",
            "PRIOR_ROBUST",
            "BULK_CELL_STRONG",
            "N_QTL_TYPES",
            "N_INDEPENDENT_LOCUS_STUDY_PAIRS",
            "EVIDENCE_SCORE_STEP14",
        ],
        ascending=[False, False, False, False, False, False],
        kind="stable",
    ).reset_index(drop=True)

    genes["RANK_STEP14"] = np.arange(1, len(genes) + 1)
    return genes


def build_locus_top_genes(step12):
    rows = []

    for (
        study,
        locus,
    ), group in step12.groupby(
        [
            "STUDY_ACCESSION",
            "LOCUS_ID",
        ],
        sort=True,
    ):
        z = group.copy()

        z["_H4"] = safe_numeric(
            z[
                "BEST_H4"
            ]
        )

        z = z.sort_values(
            [
                "_H4",
            ],
            ascending=[
                False,
            ],
            kind="stable",
        )

        best = z.iloc[0]

        rows.append(
            {
                "STUDY_ACCESSION":
                    study,

                "LOCUS_ID":
                    locus,

                "TOP_GENE_ID":
                    clean_gene_id(
                        best[
                            "EFFECTOR_GENE_ID"
                        ]
                    ),

                "TOP_GENE_SYMBOL":
                    clean_symbol(
                        best.get(
                            "EFFECTOR_GENE_SYMBOL",
                            "",
                        )
                    ),

                "TOP_GENE_H4":
                    best[
                        "_H4"
                    ],

                "N_CANDIDATE_GENES":
                    group[
                        "EFFECTOR_GENE_ID"
                    ]
                    .astype(str)
                    .nunique(),

                "ALL_CANDIDATE_GENES":
                    join_unique(
                        group[
                            "EFFECTOR_GENE_ID"
                        ]
                    ),
            }
        )

    return pd.DataFrame(rows)


# =============================================================================
# ATTACH OPEN TARGETS EVIDENCE WITHOUT USING IT IN GWAS2m SCORE
# =============================================================================

def attach_step13_evidence(genes, step13):
    """Attach external evidence without changing the GWAS2m ranking score."""
    out = genes.copy()

    defaults = {
        "OT_ANY_LOCUS_MATCH": "UNKNOWN",
        "OT_L2G_SUPPORTED_ANY": "UNKNOWN",
        "OT_L2G_BEST_RANK": np.nan,
        "OT_L2G_MAX_SCORE": np.nan,
        "OT_MOLQTL_SUPPORTED_ANY": "UNKNOWN",
        "OT_MOLQTL_MAX_H4": np.nan,
        "OT_EVIDENCE_CLASSES": "",
        "OT_EXTERNAL_MECHANISM_CLASSES": "",
        "OT_CELL_CONTEXT_SCREENS": "",
        "OT_NOVELTY_TIERS": "",
    }

    for col, default in defaults.items():
        out[col] = default

    if step13.empty:
        return out

    gid_col = first_existing(step13.columns, ["EFFECTOR_GENE_ID"])
    if not gid_col:
        return out

    locus_match_col = first_existing(step13.columns, ["OT_LOCUS_MATCH"])
    present_col = first_existing(step13.columns, ["OT_L2G_GENE_PRESENT"])
    rank_col = first_existing(step13.columns, ["OT_L2G_RANK"])
    score_col = first_existing(step13.columns, ["OT_L2G_SCORE"])
    molqtl_col = first_existing(step13.columns, ["OT_COLOC_GENE_PRESENT_H4_THRESHOLD"])
    molqtl_h4_col = first_existing(step13.columns, ["OT_COLOC_MAX_H4_SAME_GENE"])
    class_col = first_existing(step13.columns, ["OPEN_TARGETS_EVIDENCE_CLASS"])
    mechanism_col = first_existing(step13.columns, ["EXTERNAL_MECHANISM_CLASS"])
    cell_screen_col = first_existing(step13.columns, ["CELL_CONTEXT_SCREEN"])
    novelty_col = first_existing(step13.columns, ["NOVELTY_TIER_STEP13"])

    rows = []
    for gid, group in step13.groupby(gid_col, sort=False):
        gid = clean_gene_id(gid)
        if not gid:
            continue
        rows.append(
            {
                "EFFECTOR_GENE_ID": gid,
                "OT_ANY_LOCUS_MATCH": "YES" if locus_match_col and group[locus_match_col].map(boolish).any() else "NO",
                "OT_L2G_SUPPORTED_ANY": "YES" if present_col and group[present_col].map(boolish).any() else "NO",
                "OT_L2G_BEST_RANK": safe_numeric(group[rank_col]).min() if rank_col and safe_numeric(group[rank_col]).notna().any() else np.nan,
                "OT_L2G_MAX_SCORE": safe_numeric(group[score_col]).max() if score_col and safe_numeric(group[score_col]).notna().any() else np.nan,
                "OT_MOLQTL_SUPPORTED_ANY": "YES" if molqtl_col and group[molqtl_col].map(boolish).any() else "NO",
                "OT_MOLQTL_MAX_H4": safe_numeric(group[molqtl_h4_col]).max() if molqtl_h4_col and safe_numeric(group[molqtl_h4_col]).notna().any() else np.nan,
                "OT_EVIDENCE_CLASSES": join_unique(group[class_col]) if class_col else "",
                "OT_EXTERNAL_MECHANISM_CLASSES": join_unique(group[mechanism_col]) if mechanism_col else "",
                "OT_CELL_CONTEXT_SCREENS": join_unique(group[cell_screen_col]) if cell_screen_col else "",
                "OT_NOVELTY_TIERS": join_unique(group[novelty_col]) if novelty_col else "",
            }
        )

    aggregate = pd.DataFrame(rows)
    if aggregate.empty:
        return out

    out = out.drop(columns=list(defaults), errors="ignore").merge(
        aggregate,
        on="EFFECTOR_GENE_ID",
        how="left",
    )

    for col, default in defaults.items():
        if col not in out.columns:
            out[col] = default
        elif isinstance(default, str):
            out[col] = out[col].fillna(default)
    return out


def add_publication_priority(genes, strong_h4=0.80, suggestive_h4=0.50):
    out = genes.copy()
    tags_all = []
    classes = []

    for _, row in out.iterrows():
        max_h4 = pd.to_numeric(pd.Series([row.get("MAX_H4")]), errors="coerce").iloc[0]
        strong = pd.notna(max_h4) and max_h4 >= strong_h4
        suggestive = pd.notna(max_h4) and max_h4 >= suggestive_h4
        external = str(row.get("OT_EXTERNAL_MECHANISM_CLASSES", ""))
        novelty = str(row.get("OT_NOVELTY_TIERS", ""))
        tags = []

        if strong:
            tags.append("STRONG_H4")
        elif suggestive:
            tags.append("SUGGESTIVE_H4")
        if row.get("PRIOR_ROBUST") == "YES":
            tags.append("PRIOR_ROBUST")
        if row.get("BULK_STRONG") == "YES":
            tags.append("STRONG_BULK")
        if row.get("CELL_STRONG") == "YES":
            tags.append("STRONG_CELL")
        if row.get("BULK_CELL_STRONG") == "YES":
            tags.append("BULK_CELL_CONVERGENCE")
        if row.get("OT_L2G_SUPPORTED_ANY") == "YES":
            tags.append("OT_L2G_SUPPORTED")
        if row.get("OT_MOLQTL_SUPPORTED_ANY") == "YES":
            tags.append("OT_MOLQTL_SUPPORTED")
        if "NEW_CELL_CONTEXT" in external:
            tags.append("NEW_CELL_CONTEXT_CANDIDATE")
        if "NEW_MOLECULAR_MECHANISM" in external:
            tags.append("NEW_MOLECULAR_MECHANISM_CANDIDATE")
        if "C_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE" in novelty:
            tags.append("NEW_EFFECTOR_CANDIDATE")
        if "D_POTENTIAL_NEW_LOCUS_CANDIDATE_REQUIRES_EXTERNAL_VALIDATION" in novelty:
            tags.append("POTENTIAL_NEW_LOCUS_REQUIRES_VALIDATION")

        if not suggestive:
            cls = "BACKGROUND_OR_WEAK_EVIDENCE"
        elif strong and "NEW_CELL_CONTEXT" in external and row.get("OT_L2G_SUPPORTED_ANY") == "YES":
            cls = "P1_OT_SUPPORTED_NEW_CELL_CONTEXT_CANDIDATE"
        elif strong and "NEW_MOLECULAR_MECHANISM" in external and row.get("OT_L2G_SUPPORTED_ANY") == "YES":
            cls = "P1_OT_SUPPORTED_NEW_MOLECULAR_MECHANISM_CANDIDATE"
        elif strong and row.get("BULK_CELL_STRONG") == "YES":
            cls = "P1_STRONG_BULK_CELL_CONVERGENCE"
        elif strong and row.get("OT_L2G_SUPPORTED_ANY") == "YES":
            cls = "P1_EXTERNAL_GENE_SUPPORT"
        elif strong and "D_POTENTIAL_NEW_LOCUS_CANDIDATE_REQUIRES_EXTERNAL_VALIDATION" in novelty:
            cls = "P2_POTENTIAL_NEW_LOCUS_REQUIRES_VALIDATION"
        elif strong and "C_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE" in novelty:
            cls = "P2_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE"
        elif strong and row.get("OT_ANY_LOCUS_MATCH") == "NO":
            cls = "P3_STRONG_OT_COVERAGE_GAP"
        elif strong:
            cls = "P3_STRONG_INTERNAL_EVIDENCE"
        else:
            cls = "P4_SUGGESTIVE_SENSITIVITY"

        tags_all.append(";".join(tags))
        classes.append(cls)

    out["PUBLICATION_EVIDENCE_TAGS"] = tags_all
    out["PUBLICATION_PRIORITY_CLASS"] = classes
    return out


def build_publication_priority_table(genes, suggestive_h4=0.50):
    if genes.empty:
        return genes.copy()
    out = genes[safe_numeric(genes["MAX_H4"]) >= suggestive_h4].copy()
    priority_order = {
        "P1_OT_SUPPORTED_NEW_CELL_CONTEXT_CANDIDATE": 1,
        "P1_OT_SUPPORTED_NEW_MOLECULAR_MECHANISM_CANDIDATE": 2,
        "P1_STRONG_BULK_CELL_CONVERGENCE": 3,
        "P1_EXTERNAL_GENE_SUPPORT": 4,
        "P2_POTENTIAL_NEW_LOCUS_REQUIRES_VALIDATION": 5,
        "P2_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE": 6,
        "P3_STRONG_OT_COVERAGE_GAP": 7,
        "P3_STRONG_INTERNAL_EVIDENCE": 8,
        "P4_SUGGESTIVE_SENSITIVITY": 9,
    }
    out["_PRIORITY_ORDER"] = out["PUBLICATION_PRIORITY_CLASS"].map(priority_order).fillna(99)
    out = out.sort_values(
        ["_PRIORITY_ORDER", "MAX_H4", "PRIOR_ROBUST", "EVIDENCE_SCORE_STEP14"],
        ascending=[True, False, False, False],
        kind="stable",
    ).drop(columns="_PRIORITY_ORDER")
    return out.reset_index(drop=True)


def extract_ot_top1(step13):
    if step13.empty:
        return []

    col = first_existing(
        step13.columns,
        [
            "OT_L2G_TOP1_GENE",
        ],
    )

    if not col:
        return []

    return sorted(
        {
            clean_gene_id(value)
            for value
            in step13[col]
            if clean_gene_id(value)
        }
    )


# =============================================================================
# GENE SETS
# =============================================================================

def modality_strong_genes(
    modality,
    threshold,
):
    result = {}

    if modality.empty:
        for qtl in QTL_TYPES:
            result[qtl] = []
        return result

    for qtl in QTL_TYPES:
        z = modality[
            (
                modality[
                    "QTL_TYPE"
                ]
                .astype(str)
                .str.lower()
                ==
                qtl.lower()
            )
            &
            (
                safe_numeric(
                    modality[
                        "MAX_H4"
                    ]
                )
                >=
                threshold
            )
        ]

        result[qtl] = sorted(
            {
                clean_gene_id(value)
                for value
                in z[
                    "EFFECTOR_GENE_ID"
                ]
                if clean_gene_id(value)
            }
        )

    return result


def build_gene_sets(
    genes,
    modality,
    locus_top,
    step13,
    strong_h4,
    suggestive_h4,
):
    result = {}

    strong = safe_numeric(genes["MAX_H4"]) >= strong_h4
    suggestive = safe_numeric(genes["MAX_H4"]) >= suggestive_h4

    result["CORE_STRONG"] = genes.loc[strong, "EFFECTOR_GENE_ID"].tolist()
    result["CORE_STRONG_PRIOR_ROBUST"] = genes.loc[
        strong & genes["PRIOR_ROBUST"].eq("YES"), "EFFECTOR_GENE_ID"
    ].tolist()
    result["CORE_STRONG_MULTI_QTL"] = genes.loc[
        strong & (safe_numeric(genes["N_QTL_TYPES"]) >= 2), "EFFECTOR_GENE_ID"
    ].tolist()
    result["EXPANDED_SUGGESTIVE"] = genes.loc[suggestive, "EFFECTOR_GENE_ID"].tolist()

    if "BULK_MAX_H4" in genes.columns:
        result["CORE_STRONG_BULK"] = genes.loc[
            safe_numeric(genes["BULK_MAX_H4"]) >= strong_h4, "EFFECTOR_GENE_ID"
        ].tolist()
    if "CELL_MAX_H4" in genes.columns:
        result["CORE_STRONG_CELL"] = genes.loc[
            safe_numeric(genes["CELL_MAX_H4"]) >= strong_h4, "EFFECTOR_GENE_ID"
        ].tolist()
    if "BULK_CELL_STRONG" in genes.columns:
        result["CORE_STRONG_BULK_AND_CELL"] = genes.loc[
            strong & genes["BULK_CELL_STRONG"].eq("YES"), "EFFECTOR_GENE_ID"
        ].tolist()

    if "OT_L2G_SUPPORTED_ANY" in genes.columns:
        result["CORE_STRONG_OT_SUPPORTED"] = genes.loc[
            strong & genes["OT_L2G_SUPPORTED_ANY"].eq("YES"), "EFFECTOR_GENE_ID"
        ].tolist()
        result["CORE_STRONG_OT_DISCORDANT"] = genes.loc[
            strong & genes["OT_ANY_LOCUS_MATCH"].eq("YES") & genes["OT_L2G_SUPPORTED_ANY"].eq("NO"),
            "EFFECTOR_GENE_ID",
        ].tolist()
        result["CORE_STRONG_OT_COVERAGE_GAP"] = genes.loc[
            strong & genes["OT_ANY_LOCUS_MATCH"].eq("NO"), "EFFECTOR_GENE_ID"
        ].tolist()

    if "OT_MOLQTL_SUPPORTED_ANY" in genes.columns:
        result["CORE_STRONG_OT_MOLQTL_SUPPORTED"] = genes.loc[
            strong & genes["OT_MOLQTL_SUPPORTED_ANY"].eq("YES"), "EFFECTOR_GENE_ID"
        ].tolist()

    if "OT_EXTERNAL_MECHANISM_CLASSES" in genes.columns:
        external = genes["OT_EXTERNAL_MECHANISM_CLASSES"].fillna("").astype(str)
        result["CORE_STRONG_NEW_CELL_CONTEXT_CANDIDATE"] = genes.loc[
            strong & external.str.contains("NEW_CELL_CONTEXT", regex=False), "EFFECTOR_GENE_ID"
        ].tolist()
        result["CORE_STRONG_NEW_MOLECULAR_MECHANISM_CANDIDATE"] = genes.loc[
            strong & external.str.contains("NEW_MOLECULAR_MECHANISM", regex=False), "EFFECTOR_GENE_ID"
        ].tolist()

    if "OT_NOVELTY_TIERS" in genes.columns:
        novelty = genes["OT_NOVELTY_TIERS"].fillna("").astype(str)
        result["CORE_STRONG_NEW_EFFECTOR_CANDIDATE"] = genes.loc[
            strong & novelty.str.contains("C_KNOWN_OT_LOCUS_NEW_EFFECTOR_CANDIDATE", regex=False),
            "EFFECTOR_GENE_ID",
        ].tolist()
        result["CORE_STRONG_POTENTIAL_NEW_LOCUS"] = genes.loc[
            strong & novelty.str.contains("D_POTENTIAL_NEW_LOCUS_CANDIDATE_REQUIRES_EXTERNAL_VALIDATION", regex=False),
            "EFFECTOR_GENE_ID",
        ].tolist()

    if not locus_top.empty:
        result["LOCUS_TOP_STRONG"] = locus_top.loc[
            safe_numeric(locus_top["TOP_GENE_H4"]) >= strong_h4,
            "TOP_GENE_ID",
        ].tolist()

    modality_sets = modality_strong_genes(modality, strong_h4)
    for qtl, values in modality_sets.items():
        result[f"STRONG_{qtl}"] = values

    ot = extract_ot_top1(step13)
    if ot:
        result["OPEN_TARGETS_TOP1"] = ot

    clean = {}
    for name, values in result.items():
        clean[name] = sorted({clean_gene_id(v) for v in values if clean_gene_id(v)})
    return clean


def gene_set_membership_table(
    gene_sets,
    genes,
):
    symbol_lookup = dict(
        zip(
            genes[
                "EFFECTOR_GENE_ID"
            ],
            genes[
                "EFFECTOR_GENE_SYMBOL"
            ],
        )
    )

    h4_lookup = dict(
        zip(
            genes[
                "EFFECTOR_GENE_ID"
            ],
            genes[
                "MAX_H4"
            ],
        )
    )

    rows = []

    for gene_set, members in gene_sets.items():
        for rank, gid in enumerate(
            members,
            start=1,
        ):
            rows.append(
                {
                    "GENE_SET":
                        gene_set,

                    "SET_MEMBER_INDEX":
                        rank,

                    "GENE_ID":
                        gid,

                    "GENE_SYMBOL":
                        symbol_lookup.get(
                            gid,
                            "",
                        ),

                    "GWAS2M_MAX_H4":
                        h4_lookup.get(
                            gid,
                            np.nan,
                        ),
                }
            )

    return pd.DataFrame(rows)


# =============================================================================
# BACKGROUND
# =============================================================================

def load_background(
    genes,
    background_file,
):
    if not background_file:
        return sorted(
            set(
                genes[
                    "EFFECTOR_GENE_ID"
                ]
                .dropna()
                .astype(str)
            )
        ), "STEP12_UNIQUE_GENES"

    path = Path(
        background_file
    ).resolve()

    if not path.exists():
        raise SystemExit(
            f"Background file does not exist: {path}"
        )

    if path.suffix.lower() in {
        ".tsv",
        ".csv",
        ".txt",
    }:
        sep = (
            "\t"
            if path.suffix.lower()
            in {
                ".tsv",
                ".txt",
            }
            else
            ","
        )

        x = pd.read_csv(
            path,
            sep=sep,
            low_memory=False,
        )

        col = first_existing(
            x.columns,
            [
                "GENE_ID",
                "EFFECTOR_GENE_ID",
                "gene",
                "gene_id",
                "ensembl_gene_id",
            ],
        )

        if col is None:
            # One-column plain text read may have been interpreted as header.
            values = [
                clean_gene_id(v)
                for v
                in path.read_text(
                    encoding="utf-8",
                    errors="replace",
                ).splitlines()
            ]
        else:
            values = [
                clean_gene_id(v)
                for v
                in x[col]
            ]

    else:
        values = [
            clean_gene_id(v)
            for v
            in path.read_text(
                encoding="utf-8",
                errors="replace",
            ).splitlines()
        ]

    values = sorted(
        {
            v
            for v
            in values
            if v
        }
    )

    return values, str(path)


# =============================================================================
# GMT OFFLINE ENRICHMENT
# =============================================================================

def parse_gmt(path: Path):
    pathways = {}

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            parts = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if len(parts) < 3:
                continue

            name = parts[0].strip()

            genes = {
                clean_gene_id(x)
                for x
                in parts[2:]
                if clean_gene_id(x)
            }

            if name and genes:
                pathways[
                    name
                ] = genes

    return pathways


def offline_enrichment(
    query,
    background,
    gmt_files,
    gene_set_name,
):
    if hypergeom is None:
        raise RuntimeError(
            "Offline enrichment requires scipy."
        )

    query = (
        set(query)
        &
        set(background)
    )

    universe = set(
        background
    )

    M = len(
        universe
    )

    N = len(
        query
    )

    if not M or not N:
        return pd.DataFrame()

    rows = []

    for gmt_path in gmt_files:
        gmt_path = Path(
            gmt_path
        ).resolve()

        library_name = gmt_path.stem

        pathways = parse_gmt(
            gmt_path
        )

        for term_name, members in pathways.items():
            members = (
                members
                &
                universe
            )

            K = len(
                members
            )

            overlap = (
                query
                &
                members
            )

            k = len(
                overlap
            )

            if not k or not K:
                continue

            pvalue = float(
                hypergeom.sf(
                    k - 1,
                    M,
                    K,
                    N,
                )
            )

            rows.append(
                {
                    "GENE_SET":
                        gene_set_name,

                    "SOURCE":
                        library_name,

                    "TERM_ID":
                        "",

                    "TERM_NAME":
                        term_name,

                    "P_VALUE":
                        pvalue,

                    "QUERY_SIZE":
                        N,

                    "TERM_SIZE":
                        K,

                    "INTERSECTION_SIZE":
                        k,

                    "INTERSECTION_GENES":
                        ";".join(
                            sorted(
                                overlap
                            )
                        ),

                    "BACKEND":
                        "OFFLINE_GMT",
                }
            )

    result = pd.DataFrame(
        rows
    )

    if result.empty:
        return result

    result[
        "FDR"
    ] = bh_fdr(
        result[
            "P_VALUE"
        ]
    )

    return result


# =============================================================================
# G:PROFILER
# =============================================================================

def gprofiler_enrichment(
    query,
    background,
    gene_set_name,
    sources,
    organism,
):
    try:
        from gprofiler import GProfiler
    except ImportError as exc:
        raise RuntimeError(
            "g:Profiler backend requires gprofiler-official.\n"
            "Install with:\n"
            "  python -m pip install gprofiler-official"
        ) from exc

    query = sorted(
        set(query)
    )

    background = sorted(
        set(background)
    )

    if not query:
        return pd.DataFrame()

    gp = GProfiler(
        return_dataframe=True
    )

    result = gp.profile(
        organism=organism,
        query=query,
        background=background,
        sources=sources,
        user_threshold=1.0,
        significance_threshold_method="fdr",
        no_evidences=False,
    )

    if result is None or len(result) == 0:
        return pd.DataFrame()

    result = pd.DataFrame(
        result
    )

    rows = []

    for _, row in result.iterrows():
        intersections = row.get(
            "intersections",
            [],
        )

        if not isinstance(
            intersections,
            (
                list,
                tuple,
                set,
                np.ndarray,
            ),
        ):
            intersections = []

        rows.append(
            {
                "GENE_SET":
                    gene_set_name,

                "SOURCE":
                    row.get(
                        "source",
                        "",
                    ),

                "TERM_ID":
                    row.get(
                        "native",
                        "",
                    ),

                "TERM_NAME":
                    row.get(
                        "name",
                        "",
                    ),

                "P_VALUE":
                    row.get(
                        "p_value",
                        np.nan,
                    ),

                # g:Profiler p_value already reflects selected correction.
                "FDR":
                    row.get(
                        "p_value",
                        np.nan,
                    ),

                "QUERY_SIZE":
                    row.get(
                        "query_size",
                        len(query),
                    ),

                "TERM_SIZE":
                    row.get(
                        "term_size",
                        np.nan,
                    ),

                "INTERSECTION_SIZE":
                    row.get(
                        "intersection_size",
                        len(intersections),
                    ),

                "INTERSECTION_GENES":
                    ";".join(
                        map(
                            str,
                            intersections,
                        )
                    ),

                "PRECISION":
                    row.get(
                        "precision",
                        np.nan,
                    ),

                "RECALL":
                    row.get(
                        "recall",
                        np.nan,
                    ),

                "BACKEND":
                    "GPROFILER",
            }
        )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# PATHWAY FILTER / QC
# =============================================================================

def filter_pathway_results(
    results,
    fdr,
    min_overlap,
    min_term_size,
    max_term_size,
):
    if results.empty:
        return results

    out = results.copy()

    out[
        "FDR"
    ] = safe_numeric(
        out[
            "FDR"
        ]
    )

    out[
        "INTERSECTION_SIZE"
    ] = safe_numeric(
        out[
            "INTERSECTION_SIZE"
        ]
    )

    out[
        "TERM_SIZE"
    ] = safe_numeric(
        out[
            "TERM_SIZE"
        ]
    )

    out[
        "PASS_FDR"
    ] = (
        out[
            "FDR"
        ]
        <=
        fdr
    )

    out[
        "PASS_MIN_OVERLAP"
    ] = (
        out[
            "INTERSECTION_SIZE"
        ]
        >=
        min_overlap
    )

    out[
        "PASS_TERM_SIZE"
    ] = (
        out[
            "TERM_SIZE"
        ]
        .between(
            min_term_size,
            max_term_size,
            inclusive="both",
        )
        |
        out[
            "TERM_SIZE"
        ]
        .isna()
    )

    out[
        "SIGNIFICANT"
    ] = (
        out[
            "PASS_FDR"
        ]
        &
        out[
            "PASS_MIN_OVERLAP"
        ]
        &
        out[
            "PASS_TERM_SIZE"
        ]
    )

    return out


def add_gene_labels(
    pathways,
    genes,
):
    if pathways.empty:
        return pathways

    symbol_map = dict(
        zip(
            genes[
                "EFFECTOR_GENE_ID"
            ],
            genes[
                "EFFECTOR_GENE_SYMBOL"
            ],
        )
    )

    def labels(value):
        ids = split_tokens(
            value
        )

        out = []

        for gid in ids:
            gid = clean_gene_id(
                gid
            )

            symbol = clean_symbol(
                symbol_map.get(
                    gid,
                    "",
                )
            )

            out.append(
                (
                    f"{symbol} ({gid})"
                    if symbol
                    else
                    gid
                )
            )

        return "; ".join(
            out
        )

    pathways = pathways.copy()

    pathways[
        "INTERSECTION_GENE_LABELS"
    ] = pathways[
        "INTERSECTION_GENES"
    ].map(
        labels
    )

    return pathways


def pathway_set_summary(
    pathway_results,
    gene_sets,
):
    rows = []

    for gene_set, members in gene_sets.items():
        z = (
            pathway_results[
                pathway_results[
                    "GENE_SET"
                ]
                ==
                gene_set
            ].copy()
            if not pathway_results.empty
            else
            pd.DataFrame()
        )

        sig = (
            z[
                z[
                    "SIGNIFICANT"
                ]
                .map(boolish)
            ].copy()
            if (
                not z.empty
                and
                "SIGNIFICANT"
                in z.columns
            )
            else
            pd.DataFrame()
        )

        rows.append(
            {
                "GENE_SET":
                    gene_set,

                "N_UNIQUE_INPUT_GENES":
                    len(members),

                "N_TERMS_RETURNED":
                    len(z),

                "N_SIGNIFICANT_TERMS":
                    len(sig),

                "MIN_FDR":
                    (
                        safe_min(
                            z[
                                "FDR"
                            ]
                        )
                        if not z.empty
                        else
                        np.nan
                    ),

                "BEST_PATHWAY":
                    (
                        z.sort_values(
                            "FDR",
                            ascending=True,
                        )
                        .iloc[0][
                            "TERM_NAME"
                        ]
                        if (
                            not z.empty
                            and
                            safe_numeric(
                                z[
                                    "FDR"
                                ]
                            )
                            .notna()
                            .any()
                        )
                        else
                        ""
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# PATHWAY DRIVER GENES
# =============================================================================

def explode_pathway_drivers(
    significant,
    genes,
):
    if significant.empty:
        return pd.DataFrame()

    gene_lookup = genes.set_index(
        "EFFECTOR_GENE_ID"
    )

    rows = []

    for _, row in significant.iterrows():
        for gid in split_tokens(
            row.get(
                "INTERSECTION_GENES",
                "",
            )
        ):
            gid = clean_gene_id(
                gid
            )

            gene_row = (
                gene_lookup.loc[
                    gid
                ]
                if gid
                in
                gene_lookup.index
                else
                None
            )

            rows.append(
                {
                    "GENE_SET":
                        row.get(
                            "GENE_SET",
                            "",
                        ),

                    "SOURCE":
                        row.get(
                            "SOURCE",
                            "",
                        ),

                    "TERM_ID":
                        row.get(
                            "TERM_ID",
                            "",
                        ),

                    "TERM_NAME":
                        row.get(
                            "TERM_NAME",
                            "",
                        ),

                    "FDR":
                        row.get(
                            "FDR",
                            np.nan,
                        ),

                    "GENE_ID":
                        gid,

                    "GENE_SYMBOL":
                        (
                            clean_symbol(
                                gene_row[
                                    "EFFECTOR_GENE_SYMBOL"
                                ]
                            )
                            if gene_row is not None
                            else
                            ""
                        ),

                    "GENE_MAX_H4":
                        (
                            gene_row[
                                "MAX_H4"
                            ]
                            if gene_row is not None
                            else
                            np.nan
                        ),

                    "GENE_QTL_TYPES":
                        (
                            gene_row[
                                "QTL_TYPES"
                            ]
                            if gene_row is not None
                            else
                            ""
                        ),

                    "GENE_N_TISSUES":
                        (
                            gene_row[
                                "N_TISSUES"
                            ]
                            if gene_row is not None
                            else
                            np.nan
                        ),

                    "GENE_OT_L2G_SUPPORTED":
                        (
                            gene_row.get(
                                "OT_L2G_SUPPORTED_ANY",
                                "",
                            )
                            if gene_row is not None
                            else
                            ""
                        ),
                }
            )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# SIMPLE PATHWAY REDUNDANCY AUDIT
# =============================================================================

def pathway_redundancy(
    significant,
    overlap_threshold,
):
    if significant.empty:
        return pd.DataFrame()

    rows = []

    for gene_set, group in significant.groupby(
        "GENE_SET",
        sort=False,
    ):
        z = group.sort_values(
            [
                "FDR",
                "INTERSECTION_SIZE",
            ],
            ascending=[
                True,
                False,
            ],
        ).reset_index(
            drop=True
        )

        clusters = []
        cluster_id = 0

        for i, row in z.iterrows():
            members = set(
                split_tokens(
                    row[
                        "INTERSECTION_GENES"
                    ]
                )
            )

            assigned = None

            best_j = 0.0
            best_rep = None

            for cluster in clusters:
                rep = cluster[
                    "rep_members"
                ]

                score = jaccard(
                    members,
                    rep,
                )

                if (
                    pd.notna(score)
                    and
                    score
                    >=
                    overlap_threshold
                    and
                    score
                    >
                    best_j
                ):
                    best_j = score
                    assigned = cluster
                    best_rep = cluster[
                        "rep_name"
                    ]

            if assigned is None:
                cluster_id += 1

                assigned = {
                    "cluster_id":
                        cluster_id,

                    "rep_name":
                        row[
                            "TERM_NAME"
                        ],

                    "rep_members":
                        members,
                }

                clusters.append(
                    assigned
                )

                best_rep = row[
                    "TERM_NAME"
                ]

            rows.append(
                {
                    "GENE_SET":
                        gene_set,

                    "REDUNDANCY_CLUSTER":
                        assigned[
                            "cluster_id"
                        ],

                    "IS_CLUSTER_REPRESENTATIVE":
                        (
                            row[
                                "TERM_NAME"
                            ]
                            ==
                            assigned[
                                "rep_name"
                            ]
                        ),

                    "REPRESENTATIVE_PATHWAY":
                        assigned[
                            "rep_name"
                        ],

                    "SOURCE":
                        row[
                            "SOURCE"
                        ],

                    "TERM_NAME":
                        row[
                            "TERM_NAME"
                        ],

                    "FDR":
                        row[
                            "FDR"
                        ],

                    "INTERSECTION_SIZE":
                        row[
                            "INTERSECTION_SIZE"
                        ],

                    "INTERSECTION_GENES":
                        row[
                            "INTERSECTION_GENES"
                        ],

                    "JACCARD_TO_REPRESENTATIVE":
                        (
                            1.0
                            if row[
                                "TERM_NAME"
                            ]
                            ==
                            assigned[
                                "rep_name"
                            ]
                            else
                            best_j
                        ),
                }
            )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# GWAS2m vs OPEN TARGETS PATHWAYS
# =============================================================================

def compare_pathway_sets(
    significant,
    left_name,
    right_name,
):
    if significant.empty:
        return (
            pd.DataFrame(),
            pd.DataFrame(),
        )

    left = significant[
        significant[
            "GENE_SET"
        ]
        ==
        left_name
    ]

    right = significant[
        significant[
            "GENE_SET"
        ]
        ==
        right_name
    ]

    def keys(df):
        return {
            (
                str(row.SOURCE),
                str(row.TERM_NAME),
            )
            for row
            in df.itertuples()
        }

    a = keys(
        left
    )

    b = keys(
        right
    )

    union = (
        a
        |
        b
    )

    common = (
        a
        &
        b
    )

    rows = []

    for source, name in sorted(
        union
    ):
        rows.append(
            {
                "SOURCE":
                    source,

                "TERM_NAME":
                    name,

                "GWAS2M_SIGNIFICANT":
                    (
                        source,
                        name,
                    )
                    in
                    a,

                "OPEN_TARGETS_SIGNIFICANT":
                    (
                        source,
                        name,
                    )
                    in
                    b,

                "PATHWAY_CLASS":
                    (
                        "COMMON"
                        if (
                            source,
                            name,
                        )
                        in
                        common
                        else
                        "GWAS2M_ONLY"
                        if (
                            source,
                            name,
                        )
                        in
                        a
                        else
                        "OPEN_TARGETS_ONLY"
                    ),
            }
        )

    summary = pd.DataFrame(
        [
            {
                "GWAS2M_GENE_SET":
                    left_name,

                "OPEN_TARGETS_GENE_SET":
                    right_name,

                "GWAS2M_SIGNIFICANT_PATHWAYS":
                    len(a),

                "OPEN_TARGETS_SIGNIFICANT_PATHWAYS":
                    len(b),

                "COMMON_SIGNIFICANT_PATHWAYS":
                    len(common),

                "UNION_SIGNIFICANT_PATHWAYS":
                    len(union),

                "PATHWAY_JACCARD":
                    (
                        len(common)
                        /
                        len(union)
                        if union
                        else
                        np.nan
                    ),
            }
        ]
    )

    return (
        pd.DataFrame(rows),
        summary,
    )


# =============================================================================
# SCREEN REPORT
# =============================================================================

def gene_label(
    gid,
    genes,
):
    z = genes[
        genes[
            "EFFECTOR_GENE_ID"
        ]
        ==
        gid
    ]

    if z.empty:
        return gid

    symbol = clean_symbol(
        z.iloc[0][
            "EFFECTOR_GENE_SYMBOL"
        ]
    )

    h4 = z.iloc[0][
        "MAX_H4"
    ]

    label = (
        f"{symbol} ({gid})"
        if symbol
        else
        gid
    )

    if pd.notna(h4):
        label += (
            f" [H4={float(h4):.6f}]"
        )

    return label


def print_pipeline_collapse(step12, genes, gene_sets):
    section("WHAT HAPPENED TO THE LARGE MASTER TABLE?")
    print(f"Step12 locus x gene rows                       : {len(step12):,}")
    print(f"Unique genes after deduplication               : {len(genes):,}")
    print(f"Unique strong genes (MAX_H4 >= 0.80)           : {len(gene_sets.get('CORE_STRONG', [])):,}")
    print(f"Unique suggestive+ genes (MAX_H4 >= 0.50)      : {len(gene_sets.get('EXPANDED_SUGGESTIVE', [])):,}")
    print(f"Strong + prior robust genes                    : {len(gene_sets.get('CORE_STRONG_PRIOR_ROBUST', [])):,}")
    print(f"Strong bulk genes                              : {len(gene_sets.get('CORE_STRONG_BULK', [])):,}")
    print(f"Strong cell-type genes                         : {len(gene_sets.get('CORE_STRONG_CELL', [])):,}")
    print(f"Strong in both bulk + cell                     : {len(gene_sets.get('CORE_STRONG_BULK_AND_CELL', [])):,}")
    print(f"OT-supported strong genes                      : {len(gene_sets.get('CORE_STRONG_OT_SUPPORTED', [])):,}")
    print(f"OT-supported new-cell-context candidates       : {len(gene_sets.get('CORE_STRONG_NEW_CELL_CONTEXT_CANDIDATE', [])):,}")
    print(f"Strong potential new-locus candidates          : {len(gene_sets.get('CORE_STRONG_POTENTIAL_NEW_LOCUS', [])):,}")
    print()
    print("RULE: repeated tissues / cells / QTL datasets / studies do NOT duplicate a gene in enrichment.")
    print("They are retained as evidence columns supporting that one gene.")


def print_all_gene_ranking(genes, compact):
    section("DEDUPLICATED GENE RANKING - BULK + CELL + VARIANT + OPEN TARGETS CONTEXT")
    columns = [
        "RANK_STEP14", "EFFECTOR_GENE_ID", "EFFECTOR_GENE_SYMBOL",
        "MAX_H4", "PRIOR_ROBUST", "BEST_EVIDENCE_SOURCE", "BEST_QTL_TYPE",
        "BEST_CONTEXT", "BEST_CONDITION", "BEST_GWAS_VARIANT", "BEST_QTL_VARIANT",
        "BULK_MAX_H4", "BULK_BEST_CONTEXT", "CELL_MAX_H4", "CELL_BEST_CONTEXT",
        "CELL_BEST_CONDITION", "BULK_CELL_STRONG",
        "MAX_PIP", "TOP_PIP_VARIANT", "MAX_SPLICEAI", "MAX_PANGOLIN",
        "N_QTL_TYPES", "N_STUDIES", "N_INDEPENDENT_LOCUS_STUDY_PAIRS",
        "OT_L2G_SUPPORTED_ANY", "OT_L2G_BEST_RANK", "OT_MOLQTL_SUPPORTED_ANY",
        "OT_EXTERNAL_MECHANISM_CLASSES", "OT_NOVELTY_TIERS",
        "PUBLICATION_PRIORITY_CLASS", "EVIDENCE_SCORE_STEP14",
    ]
    show(genes, max_rows=(50 if compact else None), columns=columns)


def print_gene_sets(gene_sets, genes):
    section("EXACT GENE SETS USED FOR PATHWAY ANALYSIS")
    print("Each list below is DEDUPLICATED. A gene appears once per set.\n")
    for name, members in gene_sets.items():
        print()
        print(f"{name}  [{len(members)} unique gene(s)]")
        print("-" * 100)
        if not members:
            print("  (empty)")
            continue
        for i, gid in enumerate(members, start=1):
            print(f"  {i:>3}. {gene_label(gid, genes)}")


def print_strong_support(genes, tissue, cell, modality, strong_h4):
    section("WHY EACH CORE-STRONG GENE WAS SELECTED")
    strong = genes[safe_numeric(genes["MAX_H4"]) >= strong_h4].sort_values("RANK_STEP14")

    for row in strong.itertuples(index=False):
        gid = row.EFFECTOR_GENE_ID
        symbol = clean_symbol(row.EFFECTOR_GENE_SYMBOL)
        print()
        print(f"{row.RANK_STEP14}. {symbol + ' ' if symbol else ''}({gid})")
        print("-" * 110)
        print(f"  MAX H4                          : {row.MAX_H4:.6f}")
        print(f"  Prior robust                    : {row.PRIOR_ROBUST}")
        print(f"  Best evidence source            : {getattr(row, 'BEST_EVIDENCE_SOURCE', '')}")
        print(f"  Best context / condition        : {getattr(row, 'BEST_CONTEXT', '')} / {getattr(row, 'BEST_CONDITION', '')}")
        print(f"  Best GWAS / QTL variants        : {getattr(row, 'BEST_GWAS_VARIANT', '')} / {getattr(row, 'BEST_QTL_VARIANT', '')}")
        print(f"  Bulk max H4 / context           : {getattr(row, 'BULK_MAX_H4', np.nan)} / {getattr(row, 'BULK_BEST_CONTEXT', '')}")
        print(f"  Cell max H4 / context           : {getattr(row, 'CELL_MAX_H4', np.nan)} / {getattr(row, 'CELL_BEST_CONTEXT', '')}")
        print(f"  Cell condition                  : {getattr(row, 'CELL_BEST_CONDITION', '')}")
        print(f"  Fine-map max PIP                : {getattr(row, 'MAX_PIP', np.nan)}")
        print(f"  Max SpliceAI                    : {getattr(row, 'MAX_SPLICEAI', np.nan)}")
        print(f"  Max Pangolin                    : {getattr(row, 'MAX_PANGOLIN', np.nan)}")
        print(f"  QTL modalities                  : {row.QTL_TYPES}")
        print(f"  Number studies                  : {row.N_STUDIES}")
        print(f"  Independent locus-study pairs   : {row.N_INDEPENDENT_LOCUS_STUDY_PAIRS}")
        print(f"  OT L2G supports same gene        : {getattr(row, 'OT_L2G_SUPPORTED_ANY', '')}")
        print(f"  OT best L2G rank                 : {getattr(row, 'OT_L2G_BEST_RANK', np.nan)}")
        print(f"  OT same-gene molQTL support      : {getattr(row, 'OT_MOLQTL_SUPPORTED_ANY', '')}")
        print(f"  External mechanism class         : {getattr(row, 'OT_EXTERNAL_MECHANISM_CLASSES', '')}")
        print(f"  Publication priority             : {getattr(row, 'PUBLICATION_PRIORITY_CLASS', '')}")
        print(f"  Evidence score (not probability) : {row.EVIDENCE_SCORE_STEP14:.6f}")

        m = modality[modality["EFFECTOR_GENE_ID"].eq(gid)].sort_values("MAX_H4", ascending=False)
        if not m.empty:
            print("  QTL-modality evidence:")
            for mr in m.itertuples(index=False):
                print(f"    - {mr.QTL_TYPE}: max H4={mr.MAX_H4:.6f}; contexts={getattr(mr, 'CONTEXTS', '')}; locus-study pairs={mr.N_LOCUS_STUDY_PAIRS}")

        t = tissue[tissue["EFFECTOR_GENE_ID"].eq(gid)].sort_values("MAX_H4", ascending=False) if not tissue.empty else pd.DataFrame()
        if not t.empty:
            print("  Bulk tissue evidence:")
            for tr in t.itertuples(index=False):
                print(f"    - {tr.TISSUE}: max H4={tr.MAX_H4:.6f}; QTL={tr.QTL_TYPES}; studies={tr.N_STUDIES}")

        c = cell[cell["EFFECTOR_GENE_ID"].eq(gid)].sort_values("MAX_H4", ascending=False) if not cell.empty else pd.DataFrame()
        if not c.empty:
            print("  Cell-type evidence:")
            for cr in c.itertuples(index=False):
                print(f"    - {cr.CELL_TYPE} [{cr.CONDITION}]: max H4={cr.MAX_H4:.6f}; QTL={cr.QTL_TYPES}; studies={cr.N_STUDIES}")


def print_pathway_results(
    all_results,
    significant,
    gene_sets,
    compact,
):
    section(
        "PATHWAY ENRICHMENT RESULTS BY GENE SET"
    )

    if (
        all_results is None
        or all_results.empty
        or "GENE_SET" not in all_results.columns
    ):
        print("(no pathway results available; enrichment may have been skipped or returned no terms)")
        for name, members in gene_sets.items():
            print(f"  {name}: {len(members)} unique input gene(s)")
        return

    for name, members in gene_sets.items():
        print()

        print(
            f"{name}: {len(members)} unique input gene(s)"
        )

        print(
            "-" * 130
        )

        if len(members) < 2:
            print(
                "  NOT TESTED: fewer than 2 genes."
            )
            continue

        z = all_results[
            all_results[
                "GENE_SET"
            ]
            ==
            name
        ].copy()

        s = significant[
            significant[
                "GENE_SET"
            ]
            ==
            name
        ].copy()

        print(
            f"  Returned terms       : {len(z):,}"
        )

        print(
            f"  Significant terms    : {len(s):,}"
        )

        if s.empty:
            if z.empty:
                print(
                    "  No terms returned."
                )
                continue

            print(
                "  No term passed all FDR / overlap / term-size filters."
            )

            best = z.sort_values(
                "FDR"
            ).head(
                10
            )

            show(
                best,
                columns=[
                    "SOURCE",
                    "TERM_NAME",
                    "FDR",
                    "INTERSECTION_SIZE",
                    "TERM_SIZE",
                    "INTERSECTION_GENE_LABELS",
                    "PASS_FDR",
                    "PASS_MIN_OVERLAP",
                    "PASS_TERM_SIZE",
                ],
            )

            continue

        s = s.sort_values(
            [
                "FDR",
                "INTERSECTION_SIZE",
            ],
            ascending=[
                True,
                False,
            ],
        )

        show(
            s,
            max_rows=(
                30
                if compact
                else
                None
            ),
            columns=[
                "SOURCE",
                "TERM_ID",
                "TERM_NAME",
                "FDR",
                "INTERSECTION_SIZE",
                "TERM_SIZE",
                "INTERSECTION_GENE_LABELS",
            ],
        )


def print_pathway_comparison(
    comparison,
    summary,
    compact,
):
    section(
        "GWAS2m CORE_STRONG vs OPEN TARGETS TOP1 PATHWAYS"
    )

    if summary.empty:
        print(
            "(comparison unavailable)"
        )
        return

    print(
        summary.to_string(
            index=False
        )
    )

    print()

    print(
        "INTERPRETATION:"
    )

    print(
        "  COMMON          = significant in both gene-set analyses."
    )

    print(
        "  GWAS2m_ONLY     = significant only for GWAS2m CORE_STRONG."
    )

    print(
        "  OPEN_TARGETS_ONLY = significant only for OT TOP1 genes."
    )

    print(
        "  These are pathway-set differences, not proof of biological novelty."
    )

    if comparison.empty:
        return

    for klass in [
        "COMMON",
        "GWAS2m_ONLY",
        "OPEN_TARGETS_ONLY",
    ]:
        z = comparison[
            comparison[
                "PATHWAY_CLASS"
            ]
            ==
            klass
        ]

        print()

        print(
            f"{klass}: {len(z)} pathway(s)"
        )

        print(
            "-" * 100
        )

        if z.empty:
            print(
                "  (none)"
            )
        else:
            show(
                z,
                max_rows=(
                    50
                    if compact
                    else
                    None
                ),
                columns=[
                    "SOURCE",
                    "TERM_NAME",
                ],
            )


def print_plain_english(
    step12,
    genes,
    gene_sets,
    significant,
    comparison_summary,
    strong_h4,
    suggestive_h4,
):
    section(
        "PLAIN-ENGLISH BIOLOGICAL INTERPRETATION"
    )

    print(
        f"1. Step12 contained {len(step12):,} locus x gene evidence rows."
    )

    print(
        f"2. Those collapse to {len(genes):,} unique genes."
    )

    print(
        f"3. {len(gene_sets.get('CORE_STRONG', []))} unique genes have MAX_H4 >= {strong_h4}."
    )

    print(
        f"4. {len(gene_sets.get('CORE_STRONG_PRIOR_ROBUST', []))} strong genes are prior-robust."
    )

    print(
        f"5. {len(gene_sets.get('CORE_STRONG_MULTI_QTL', []))} strong genes have >=2 QTL modalities."
    )

    print(
        f"6. {len(gene_sets.get('EXPANDED_SUGGESTIVE', []))} unique genes have MAX_H4 >= {suggestive_h4}."
    )

    if "OT_L2G_SUPPORTED_ANY" in genes.columns:
        strong = genes[
            genes[
                "MAX_H4"
            ]
            >=
            strong_h4
        ]

        print(
            f"7. Among strong unique genes, "
            f"{int((strong['OT_L2G_SUPPORTED_ANY'] == 'YES').sum())} "
            f"have Open Targets L2G support."
        )

    if not significant.empty:
        core = significant[
            significant[
                "GENE_SET"
            ]
            ==
            "CORE_STRONG"
        ]

        print(
            f"8. CORE_STRONG produced {len(core)} significant pathways after filters."
        )

        if not core.empty:
            print(
                "   Main CORE_STRONG pathways:"
            )

            for row in core.sort_values(
                "FDR"
            ).itertuples(
                index=False
            ):
                print(
                    f"   - {row.SOURCE}: {row.TERM_NAME}; "
                    f"FDR={row.FDR:.5g}; "
                    f"genes={row.INTERSECTION_GENE_LABELS}"
                )

    if not comparison_summary.empty:
        row = comparison_summary.iloc[0]

        print(
            f"9. GWAS2m-vs-OT significant pathway Jaccard = "
            f"{row['PATHWAY_JACCARD']:.4f}."
        )

        print(
            f"   Common={int(row['COMMON_SIGNIFICANT_PATHWAYS'])}; "
            f"GWAS2m-only={int(row['GWAS2M_SIGNIFICANT_PATHWAYS'] - row['COMMON_SIGNIFICANT_PATHWAYS'])}; "
            f"OT-only={int(row['OPEN_TARGETS_SIGNIFICANT_PATHWAYS'] - row['COMMON_SIGNIFICANT_PATHWAYS'])}."
        )

    print()

    print(
        "IMPORTANT: pathway enrichment is supportive interpretation. "
        "It does not establish gene/pathway causality."
    )


# =============================================================================
# QC
# =============================================================================

def make_qc(
    step12,
    genes,
    gene_sets,
    background,
    background_source,
    pathway_results,
    enrichment_errors,
):
    rows = []

    rows.extend(
        [
            {
                "CHECK":
                    "STEP12_LOCUS_GENE_ROWS",

                "STATUS":
                    "INFO",

                "VALUE":
                    len(step12),

                "DETAIL":
                    "Input rows before gene deduplication.",
            },
            {
                "CHECK":
                    "UNIQUE_GENES",

                "STATUS":
                    "INFO",

                "VALUE":
                    len(genes),

                "DETAIL":
                    "Unique genes after collapsing repeated tissues/QTL/studies.",
            },
            {
                "CHECK":
                    "BACKGROUND_GENES",

                "STATUS":
                    (
                        "WARN"
                        if background_source
                        ==
                        "STEP12_UNIQUE_GENES"
                        else
                        "OK"
                    ),

                "VALUE":
                    len(background),

                "DETAIL":
                    (
                        "Default background is Step12 unique genes. "
                        "For final publication, use a broader tested/eligible gene universe if available."
                        if background_source
                        ==
                        "STEP12_UNIQUE_GENES"
                        else
                        f"Custom background: {background_source}"
                    ),
            },
            {
                "CHECK":
                    "CORE_STRONG_GENES",

                "STATUS":
                    "INFO",

                "VALUE":
                    len(
                        gene_sets.get(
                            "CORE_STRONG",
                            [],
                        )
                    ),

                "DETAIL":
                    "Primary pathway input.",
            },
            {
                "CHECK":
                    "PATHWAY_RESULTS",

                "STATUS":
                    (
                        "OK"
                        if len(pathway_results)
                        else
                        "WARN"
                    ),

                "VALUE":
                    len(pathway_results),

                "DETAIL":
                    "Raw pathway rows returned.",
            },
        ]
    )

    for name, err in enrichment_errors.items():
        rows.append(
            {
                "CHECK":
                    f"ENRICHMENT_{name}",

                "STATUS":
                    "ERROR",

                "VALUE":
                    "",

                "DETAIL":
                    err,
            }
        )

    missing_symbols = int(
        (
            genes[
                "EFFECTOR_GENE_SYMBOL"
            ]
            .fillna("")
            .astype(str)
            .str.len()
            ==
            0
        )
        .sum()
    )

    rows.append(
        {
            "CHECK":
                "GENES_MISSING_SYMBOL",

            "STATUS":
                (
                    "WARN"
                    if missing_symbols
                    else
                    "OK"
                ),

            "VALUE":
                missing_symbols,

            "DETAIL":
                "Ensembl IDs remain usable even when symbol is missing.",
        }
    )

    return pd.DataFrame(
        rows
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
        "--fdr",
        type=float,
        default=0.05,
    )

    parser.add_argument(
        "--min-overlap",
        type=int,
        default=2,
        help="Minimum number of distinct query genes required to call a pathway significant.",
    )

    parser.add_argument(
        "--min-term-size",
        type=int,
        default=5,
    )

    parser.add_argument(
        "--max-term-size",
        type=int,
        default=500,
    )

    parser.add_argument(
        "--redundancy-jaccard",
        type=float,
        default=0.70,
        help="Gene-overlap Jaccard used to group highly redundant significant pathways.",
    )

    parser.add_argument(
        "--organism",
        default="hsapiens",
    )

    parser.add_argument(
        "--sources",
        default="GO:BP,REAC,KEGG,WP",
        help="Comma-separated g:Profiler pathway sources.",
    )

    parser.add_argument(
        "--gmt",
        action="append",
        default=[],
        help="Offline GMT file. Repeat for multiple GMT libraries.",
    )

    parser.add_argument(
        "--background-file",
        default="",
        help=(
            "Optional custom tested/eligible-gene universe. "
            "Default is all unique Step12 genes."
        ),
    )

    parser.add_argument(
        "--step12-file",
        default="",
    )

    parser.add_argument(
        "--step13-file",
        default="",
    )

    parser.add_argument(
        "--skip-enrichment",
        action="store_true",
    )

    parser.add_argument(
        "--compact-screen",
        action="store_true",
        help="Print fewer rows. Default prints the complete interpretation report.",
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

    step12_path = (
        Path(
            args.step12_file
        ).resolve()
        if args.step12_file
        else
        root
        /
        "12_master_table"
        /
        phenotype_slug
        /
        ancestry_slug
        /
        "Step12_Locus_Gene_Master.tsv"
    )

    step13_path = (
        Path(
            args.step13_file
        ).resolve()
        if args.step13_file
        else
        root
        /
        "13_opentargets_comparison"
        /
        phenotype_slug
        /
        ancestry_slug
        /
        "Step13_Gene_Comparison.tsv"
    )

    output = (
        root
        /
        "14_pathways"
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
        "GWAS2m STEP14 v3.0 - INTEGRATED GENE PRIORITISATION + PATHWAY ANALYSIS"
    )

    print(
        f"Root                   : {root}"
    )

    print(
        f"Phenotype              : {args.phenotype}"
    )

    print(
        f"Ancestry               : {ancestry_label} ({ancestry_code})"
    )

    print(
        f"Step12                 : {step12_path}"
    )

    print(
        f"Step13                 : {step13_path if step13_path.exists() else 'NOT AVAILABLE'}"
    )

    print(
        f"Output                 : {output}"
    )

    print(
        f"Strong H4              : >= {args.strong_h4}"
    )

    print(
        f"Suggestive H4          : >= {args.suggestive_h4}"
    )

    print(
        f"FDR                    : <= {args.fdr}"
    )

    print(
        f"Minimum pathway overlap: >= {args.min_overlap} genes"
    )

    print(
        f"Term-size filter       : {args.min_term_size} - {args.max_term_size}"
    )

    # -------------------------------------------------------------------------
    # Load and collapse evidence
    # -------------------------------------------------------------------------

    step12 = load_step12(
        step12_path
    )

    step13 = load_step13(
        step13_path
    )

    contexts = build_context_evidence(
        step12
    )

    tissues = build_tissue_evidence(
        step12
    )

    cells = build_cell_evidence(
        step12
    )

    modalities = build_modality_evidence(
        step12
    )

    genes = build_gene_ranking(
        step12
    )

    genes = attach_step13_evidence(
        genes,
        step13,
    )

    genes = add_publication_priority(
        genes,
        args.strong_h4,
        args.suggestive_h4,
    )

    publication_priority = build_publication_priority_table(
        genes,
        args.suggestive_h4,
    )

    locus_top = build_locus_top_genes(
        step12
    )

    gene_sets = build_gene_sets(
        genes,
        modalities,
        locus_top,
        step13,
        args.strong_h4,
        args.suggestive_h4,
    )

    membership = gene_set_membership_table(
        gene_sets,
        genes,
    )

    background, background_source = load_background(
        genes,
        args.background_file,
    )

    # -------------------------------------------------------------------------
    # Save gene evidence
    # -------------------------------------------------------------------------

    genes.to_csv(
        output
        /
        "Step14_Gene_Ranking.tsv",
        sep="\t",
        index=False,
    )

    tissues.to_csv(
        output
        /
        "Step14_Gene_Tissue_Evidence.tsv",
        sep="\t",
        index=False,
    )

    cells.to_csv(
        output
        /
        "Step14_Gene_CellType_Evidence.tsv",
        sep="\t",
        index=False,
    )

    contexts.to_csv(
        output
        /
        "Step14_Gene_Context_Evidence.tsv",
        sep="\t",
        index=False,
    )

    publication_priority.to_csv(
        output
        /
        "Step14_Publication_Priority_Genes.tsv",
        sep="\t",
        index=False,
    )

    modalities.to_csv(
        output
        /
        "Step14_Gene_Modality_Evidence.tsv",
        sep="\t",
        index=False,
    )

    locus_top.to_csv(
        output
        /
        "Step14_Locus_Top_Genes.tsv",
        sep="\t",
        index=False,
    )

    membership.to_csv(
        output
        /
        "Step14_Gene_Set_Membership.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Print gene evidence
    # -------------------------------------------------------------------------

    print_pipeline_collapse(
        step12,
        genes,
        gene_sets,
    )

    print_all_gene_ranking(
        genes,
        args.compact_screen,
    )

    print_gene_sets(
        gene_sets,
        genes,
    )

    print_strong_support(
        genes,
        tissues,
        cells,
        modalities,
        args.strong_h4,
    )

    # -------------------------------------------------------------------------
    # Enrichment
    # -------------------------------------------------------------------------

    pathway_parts = []
    enrichment_errors = {}

    if not args.skip_enrichment:
        section(
            "RUNNING PATHWAY ENRICHMENT"
        )

        backend = (
            "OFFLINE_GMT"
            if args.gmt
            else
            "GPROFILER"
        )

        print(
            f"Backend                : {backend}"
        )

        print(
            f"Background              : {len(background):,} genes"
        )

        print(
            f"Background source       : {background_source}"
        )

        print(
            f"g:Profiler sources      : {args.sources}"
        )

        print()

        preferred_order = [
            "CORE_STRONG",
            "CORE_STRONG_PRIOR_ROBUST",
            "CORE_STRONG_MULTI_QTL",
            "CORE_STRONG_BULK",
            "CORE_STRONG_CELL",
            "CORE_STRONG_BULK_AND_CELL",
            "LOCUS_TOP_STRONG",
            "EXPANDED_SUGGESTIVE",
            "CORE_STRONG_OT_SUPPORTED",
            "CORE_STRONG_OT_MOLQTL_SUPPORTED",
            "CORE_STRONG_NEW_CELL_CONTEXT_CANDIDATE",
            "CORE_STRONG_NEW_MOLECULAR_MECHANISM_CANDIDATE",
            "CORE_STRONG_NEW_EFFECTOR_CANDIDATE",
            "CORE_STRONG_POTENTIAL_NEW_LOCUS",
            "CORE_STRONG_OT_DISCORDANT",
            "CORE_STRONG_OT_COVERAGE_GAP",
            "OPEN_TARGETS_TOP1",
        ]

        preferred_order.extend(
            [
                f"STRONG_{qtl}"
                for qtl
                in QTL_TYPES
            ]
        )

        ordered_sets = []

        for name in preferred_order:
            if name in gene_sets:
                ordered_sets.append(name)

        for name in gene_sets:
            if name not in ordered_sets:
                ordered_sets.append(name)

        sources = [
            x.strip()
            for x
            in args.sources.split(",")
            if x.strip()
        ]

        for name in ordered_sets:
            query = gene_sets[
                name
            ]

            if len(query) < 2:
                print(
                    f"[SKIP] {name}: only {len(query)} unique gene(s)"
                )
                continue

            print(
                f"[RUN ] {name}: {len(query)} unique gene(s)"
            )

            try:
                if args.gmt:
                    result = offline_enrichment(
                        query,
                        background,
                        args.gmt,
                        name,
                    )

                else:
                    result = gprofiler_enrichment(
                        query,
                        background,
                        name,
                        sources,
                        args.organism,
                    )

                if not result.empty:
                    pathway_parts.append(
                        result
                    )

            except Exception as exc:
                message = (
                    f"{type(exc).__name__}: {exc}"
                )

                enrichment_errors[
                    name
                ] = message

                print(
                    f"[FAIL] {name}: {message}"
                )

    raw_pathways = (
        pd.concat(
            pathway_parts,
            ignore_index=True,
            sort=False,
        )
        if pathway_parts
        else
        pd.DataFrame()
    )

    if not raw_pathways.empty:
        raw_pathways = filter_pathway_results(
            raw_pathways,
            args.fdr,
            args.min_overlap,
            args.min_term_size,
            args.max_term_size,
        )

        raw_pathways = add_gene_labels(
            raw_pathways,
            genes,
        )

        raw_pathways = raw_pathways.sort_values(
            [
                "GENE_SET",
                "FDR",
                "INTERSECTION_SIZE",
            ],
            ascending=[
                True,
                True,
                False,
            ],
        )

    significant = (
        raw_pathways[
            raw_pathways[
                "SIGNIFICANT"
            ]
            .map(boolish)
        ].copy()
        if (
            not raw_pathways.empty
            and
            "SIGNIFICANT"
            in raw_pathways.columns
        )
        else
        pd.DataFrame()
    )

    set_summary = pathway_set_summary(
        raw_pathways,
        gene_sets,
    )

    drivers = explode_pathway_drivers(
        significant,
        genes,
    )

    redundancy = pathway_redundancy(
        significant,
        args.redundancy_jaccard,
    )

    pathway_comparison = pd.DataFrame()
    pathway_comparison_summary = pd.DataFrame()

    if (
        "OPEN_TARGETS_TOP1"
        in
        gene_sets
    ):
        (
            pathway_comparison,
            pathway_comparison_summary,
        ) = compare_pathway_sets(
            significant,
            "CORE_STRONG",
            "OPEN_TARGETS_TOP1",
        )

    # -------------------------------------------------------------------------
    # Save pathway outputs
    # -------------------------------------------------------------------------

    raw_pathways.to_csv(
        output
        /
        "Step14_Pathway_All_Results.tsv",
        sep="\t",
        index=False,
    )

    significant.to_csv(
        output
        /
        "Step14_Pathway_Significant.tsv",
        sep="\t",
        index=False,
    )

    set_summary.to_csv(
        output
        /
        "Step14_Pathway_Set_Summary.tsv",
        sep="\t",
        index=False,
    )

    drivers.to_csv(
        output
        /
        "Step14_Pathway_Gene_Drivers.tsv",
        sep="\t",
        index=False,
    )

    redundancy.to_csv(
        output
        /
        "Step14_Pathway_Redundancy.tsv",
        sep="\t",
        index=False,
    )

    pathway_comparison.to_csv(
        output
        /
        "Step14_OpenTargets_Pathway_Comparison.tsv",
        sep="\t",
        index=False,
    )

    pathway_comparison_summary.to_csv(
        output
        /
        "Step14_OpenTargets_Pathway_Comparison_Summary.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Screen report
    # -------------------------------------------------------------------------

    show(
        set_summary,
        title="PATHWAY ENRICHMENT SUMMARY",
    )

    print_pathway_results(
        raw_pathways,
        significant,
        gene_sets,
        args.compact_screen,
    )

    print_pathway_comparison(
        pathway_comparison,
        pathway_comparison_summary,
        args.compact_screen,
    )

    print_plain_english(
        step12,
        genes,
        gene_sets,
        significant,
        pathway_comparison_summary,
        args.strong_h4,
        args.suggestive_h4,
    )

    # -------------------------------------------------------------------------
    # QC
    # -------------------------------------------------------------------------

    qc = make_qc(
        step12,
        genes,
        gene_sets,
        background,
        background_source,
        raw_pathways,
        enrichment_errors,
    )

    qc.to_csv(
        output
        /
        "Step14_QC_Audit.tsv",
        sep="\t",
        index=False,
    )

    show(
        qc,
        title="STEP14 QC AUDIT",
    )

    metadata = {
        "version":
            VERSION,

        "created_at_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),

        "phenotype":
            args.phenotype,

        "ancestry_code":
            ancestry_code,

        "ancestry_label":
            ancestry_label,

        "step12":
            str(
                step12_path
            ),

        "step13":
            (
                str(
                    step13_path
                )
                if step13_path.exists()
                else
                None
            ),

        "n_step12_rows":
            int(
                len(
                    step12
                )
            ),

        "n_unique_genes":
            int(
                len(
                    genes
                )
            ),

        "n_publication_priority_genes":
            int(
                len(
                    publication_priority
                )
            ),

        "n_background_genes":
            int(
                len(
                    background
                )
            ),

        "background_source":
            background_source,

        "strong_h4":
            args.strong_h4,

        "suggestive_h4":
            args.suggestive_h4,

        "fdr":
            args.fdr,

        "min_overlap":
            args.min_overlap,

        "term_size_range":
            [
                args.min_term_size,
                args.max_term_size,
            ],

        "gene_sets":
            {
                name:
                    len(
                        values
                    )
                for name, values
                in gene_sets.items()
            },

        "n_pathway_rows":
            int(
                len(
                    raw_pathways
                )
            ),

        "n_significant_pathways":
            int(
                len(
                    significant
                )
            ),

        "enrichment_errors":
            enrichment_errors,

        "scientific_notes":
            [
                "Each gene contributes once to pathway enrichment regardless of repeated tissue/QTL evidence.",
                "Open Targets does not contribute to the GWAS2m evidence ranking score.",
                "Bulk-tissue and cell-type evidence are retained separately before gene-level deduplication.",
                "PIP, SpliceAI and Pangolin are carried as orthogonal evidence fields and do not automatically inflate H4.",
                "Pathway enrichment is associative/supportive, not causal proof.",
                "Default background is the set of unique genes reaching Step12; a broader tested/eligible-gene universe is preferable for final publication if available.",
            ],
    }

    (
        output
        /
        "Step14_summary.json"
    ).write_text(
        json.dumps(
            metadata,
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
        "STEP14 COMPLETE"
    )

    print(
        f"Unique Step12 genes                      : {len(genes):,}"
    )

    print(
        f"CORE_STRONG genes                        : {len(gene_sets.get('CORE_STRONG', [])):,}"
    )

    print(
        f"CORE_STRONG_PRIOR_ROBUST genes           : {len(gene_sets.get('CORE_STRONG_PRIOR_ROBUST', [])):,}"
    )

    print(
        f"CORE_STRONG_MULTI_QTL genes              : {len(gene_sets.get('CORE_STRONG_MULTI_QTL', [])):,}"
    )

    print(
        f"CORE_STRONG_BULK genes                   : {len(gene_sets.get('CORE_STRONG_BULK', [])):,}"
    )

    print(
        f"CORE_STRONG_CELL genes                   : {len(gene_sets.get('CORE_STRONG_CELL', [])):,}"
    )

    print(
        f"CORE_STRONG_BULK_AND_CELL genes          : {len(gene_sets.get('CORE_STRONG_BULK_AND_CELL', [])):,}"
    )

    print(
        f"EXPANDED_SUGGESTIVE genes                : {len(gene_sets.get('EXPANDED_SUGGESTIVE', [])):,}"
    )

    print(
        f"Significant pathway rows                 : {len(significant):,}"
    )

    print()

    print(
        "Outputs:"
    )

    outputs = [
        "Step14_Gene_Ranking.tsv",
        "Step14_Gene_Tissue_Evidence.tsv",
        "Step14_Gene_CellType_Evidence.tsv",
        "Step14_Gene_Context_Evidence.tsv",
        "Step14_Gene_Modality_Evidence.tsv",
        "Step14_Publication_Priority_Genes.tsv",
        "Step14_Locus_Top_Genes.tsv",
        "Step14_Gene_Set_Membership.tsv",
        "Step14_Pathway_All_Results.tsv",
        "Step14_Pathway_Significant.tsv",
        "Step14_Pathway_Set_Summary.tsv",
        "Step14_Pathway_Gene_Drivers.tsv",
        "Step14_Pathway_Redundancy.tsv",
        "Step14_OpenTargets_Pathway_Comparison.tsv",
        "Step14_OpenTargets_Pathway_Comparison_Summary.tsv",
        "Step14_QC_Audit.tsv",
        "Step14_summary.json",
    ]

    for name in outputs:
        print(
            f"  {output / name}"
        )

    print()

    print(
        "NEXT PIPELINE STEP:"
    )

    print(
        "  Step15 = PPI / network analysis using CORE_STRONG + PRIOR_ROBUST"
    )

    print(
        "  genes, with known disease genes as external anchors."
    )


if __name__ == "__main__":
    pipeline_main(main)
