#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GWAS2m SETUP 11B -- SINGLE-CELL eQTL RESOURCE DOWNLOADER
=========================================================

SCHEDULING MODEL
----------------
ONE SLURM ARRAY JOB for <=1000 missing physical files.
ONE ARRAY INDEX = ONE PHYSICAL FILE.
NO %N throttle.

Example:
    167 biological contexts
    x 4 physical artifacts/context
    = 668 expected physical files

If 591 files are still missing:
    #SBATCH --array=1-591

and:
    task 1   -> exactly one file
    task 2   -> exactly one file
    ...
    task 591 -> exactly one file

Artifacts per provider-compatible single-cell context:
    1. dense eQTL summary statistics
    2. dense tabix index
    3. provider SuSiE LBF
    4. provider SuSiE credible sets

Downloads use aria2c, resume <file>.part, and retry until completion
(or until the SLURM walltime ends). Re-run safely to resume anything left.

Typical use:
    python Setup_11B_Download_ContextSpecific_eQTL.py \
      --profile singlecell \
      --qtl-types eQTL \
      --partition general
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd


try:
    import Step11_Colocalize_GTEx_SuSiE as core
except Exception as exc:
    raise SystemExit(
        "Could not import Step11_Colocalize_GTEx_SuSiE.py\n"
        f"{type(exc).__name__}: {exc}"
    )

try:
    import Step11_Provider_SuSiE_Formal as provider
except Exception as exc:
    raise SystemExit(
        "Could not import Step11_Provider_SuSiE_Formal.py\n"
        f"{type(exc).__name__}: {exc}"
    )


VERSION = "2.1.0-singlecell-one-array-index-per-file"
DEFAULT_ARRAY_LIMIT = 1000


SINGLECELL_FAMILIES = [
    ["Jerber_2021"],
    [
        "Randolph_2021_reannotated",
        "Randoplh_2021_reannotated",
        "Randolph_2021",
        "Randoplh_2021",
    ],
    [
        "Perez_2022_reannotated",
        "Perez_2022",
    ],
    [
        "OneK1K_reannotated",
        "OneK1K",
    ],
    [
        "Nathan_2022_reannotated",
        "Nathan_2022",
    ],
    [
        "Cytoimmgen_reannotated",
        "Cytoimmgen",
    ],
    ["IBDverse"],
]


def banner(text: str) -> None:
    print()
    print("=" * 120)
    print(text)
    print("=" * 120, flush=True)


def clean(value) -> str:
    value = str(value or "").strip()
    if value.lower() in {"", "nan", "na", "none", "."}:
        return ""
    return value


def url_to_https(url: str) -> str:
    url = clean(url)
    prefix = "ftp://ftp.ebi.ac.uk/"
    if url.startswith(prefix):
        return "https://ftp.ebi.ac.uk/" + url[len(prefix):]
    return url


def basename_from_url(url: str) -> str:
    value = clean(url).split("?", 1)[0].rstrip("/")
    return Path(urlparse(value).path).name


def nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def artifact_complete(path: Path) -> bool:
    if not nonempty(path):
        return False

    for suffix in [".part", ".aria2", ".tmp"]:
        if Path(str(path) + suffix).exists():
            return False

    return True


def gzip_ok(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(16 * 1024 * 1024):
                pass
        return True
    except Exception:
        return False


def active_job(job_id: str) -> bool:
    job_id = str(job_id or "").strip()
    if not job_id:
        return False

    proc = subprocess.run(
        ["squeue", "-h", "-j", job_id, "-o", "%T"],
        text=True,
        capture_output=True,
    )

    if proc.returncode != 0:
        return False

    return any(
        x.strip() in {
            "PENDING",
            "RUNNING",
            "CONFIGURING",
            "COMPLETING",
            "SUSPENDED",
        }
        for x in proc.stdout.splitlines()
    )


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(
        json.dumps(payload, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    os.replace(tmp, path)


def choose_alias(
    available: list[str],
    aliases: list[str],
) -> str | None:
    lookup = {
        x.casefold(): x
        for x in available
    }

    for alias in aliases:
        hit = lookup.get(alias.casefold())
        if hit:
            return hit

    return None


def select_studies(catalogue: pd.DataFrame):
    available = sorted(
        {
            clean(x)
            for x in catalogue["study_label"].tolist()
            if clean(x)
        }
    )

    selected = []
    unavailable = []

    for family in SINGLECELL_FAMILIES:
        hit = choose_alias(
            available,
            family,
        )

        if hit:
            selected.append(hit)
        else:
            unavailable.append(family)

    return (
        list(dict.fromkeys(selected)),
        unavailable,
    )


def provider_supported(row: pd.Series) -> bool:
    return bool(
        clean(row.get("ftp_path", ""))
        and clean(row.get("ftp_lbf_path", ""))
        and clean(row.get("ftp_cs_path", ""))
    )


def build_contexts(
    root: Path,
    profile: str,
    qtl_types: set[str],
):
    if profile != "singlecell":
        raise RuntimeError(
            "This Step11B setup is intentionally single-cell only. "
            "Use --profile singlecell."
        )

    catalogue = core.load_eqtl_catalogue(
        root,
        False,
    ).copy()

    if catalogue.empty:
        raise RuntimeError(
            "eQTL Catalogue metadata is empty."
        )

    if "QTL_TYPE" not in catalogue.columns:
        if "quant_method" not in catalogue.columns:
            raise RuntimeError(
                "Metadata lacks QTL_TYPE and quant_method. "
                f"Columns={list(catalogue.columns)}"
            )

        catalogue["QTL_TYPE"] = (
            catalogue["quant_method"]
            .map(core.quant_method_to_qtl_type)
        )

    required = [
        "study_label",
        "dataset_id",
        "QTL_TYPE",
        "ftp_path",
        "ftp_lbf_path",
        "ftp_cs_path",
    ]

    missing = [
        c
        for c in required
        if c not in catalogue.columns
    ]

    if missing:
        raise RuntimeError(
            f"eQTL Catalogue metadata missing columns: {missing}"
        )

    selected, unavailable = (
        select_studies(catalogue)
    )

    if not selected:
        raise RuntimeError(
            "No configured single-cell studies found in local metadata."
        )

    x = catalogue[
        catalogue["study_label"]
        .astype(str)
        .isin(selected)
    ].copy()

    allowed = {
        str(v).lower()
        for v in qtl_types
    }

    x = x[
        x["QTL_TYPE"]
        .astype(str)
        .str.lower()
        .isin(allowed)
    ].copy()

    if x.empty:
        raise RuntimeError(
            f"No single-cell contexts matched QTL types {qtl_types}"
        )

    x["PROVIDER_COMPATIBLE"] = x.apply(
        provider_supported,
        axis=1,
    )

    eligible = (
        x[x["PROVIDER_COMPATIBLE"]]
        .drop_duplicates("dataset_id")
        .reset_index(drop=True)
    )

    excluded = (
        x[~x["PROVIDER_COMPATIBLE"]]
        .drop_duplicates("dataset_id")
        .reset_index(drop=True)
    )

    return (
        eligible,
        excluded,
        selected,
        unavailable,
    )


def build_artifacts(
    root: Path,
    contexts: pd.DataFrame,
):
    dense_root = (
        root
        / "resources"
        / "coloc"
        / "eqtl_catalogue"
        / "dense"
    )

    dense_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    artifact_rows = []
    registry_rows = []

    for _, row in contexts.iterrows():
        dataset_id = clean(row["dataset_id"])
        study_label = clean(row["study_label"])
        qtl_type = clean(row["QTL_TYPE"])

        tissue_label = clean(
            row.get(
                "tissue_label",
                row.get("tissue", ""),
            )
        )

        condition_label = clean(
            row.get(
                "condition_label",
                row.get("condition", ""),
            )
        )

        dense_url = url_to_https(
            row["ftp_path"]
        )

        lbf_url = url_to_https(
            row["ftp_lbf_path"]
        )

        cs_url = url_to_https(
            row["ftp_cs_path"]
        )

        dense_name = basename_from_url(
            dense_url
        )

        dense_path = (
            dense_root
            / dataset_id
            / dense_name
        )

        index_url = dense_url + ".tbi"

        index_path = Path(
            str(dense_path) + ".tbi"
        )

        lbf_path = (
            provider.provider_local_path(
                root,
                dataset_id,
                lbf_url,
            )
        )

        cs_path = (
            provider.provider_local_path(
                root,
                dataset_id,
                cs_url,
            )
        )

        common = {
            "DATASET_ID": dataset_id,
            "STUDY_LABEL": study_label,
            "TISSUE_LABEL": tissue_label,
            "CONDITION_LABEL": condition_label,
            "QTL_TYPE": qtl_type,
        }

        for (
            kind,
            url,
            local_path,
        ) in [
            (
                "DENSE_SUMSTATS",
                dense_url,
                dense_path,
            ),
            (
                "DENSE_INDEX",
                index_url,
                index_path,
            ),
            (
                "PROVIDER_LBF",
                lbf_url,
                lbf_path,
            ),
            (
                "PROVIDER_CS",
                cs_url,
                cs_path,
            ),
        ]:
            artifact_rows.append(
                {
                    **common,
                    "ARTIFACT_TYPE": kind,
                    "URL": url,
                    "LOCAL_PATH": str(local_path),
                }
            )

        registry_rows.append(
            {
                **common,
                "DENSE_PATH": str(dense_path),
                "DENSE_INDEX_PATH": str(index_path),
                "PROVIDER_LBF_PATH": str(lbf_path),
                "PROVIDER_CS_PATH": str(cs_path),
            }
        )

    artifacts = (
        pd.DataFrame(artifact_rows)
        .drop_duplicates("LOCAL_PATH")
        .reset_index(drop=True)
    )

    registry = (
        pd.DataFrame(registry_rows)
        .drop_duplicates("DATASET_ID")
        .reset_index(drop=True)
    )

    artifacts["COMPLETE_NOW"] = [
        artifact_complete(
            Path(str(p))
        )
        for p in artifacts["LOCAL_PATH"]
    ]

    pending = (
        artifacts[
            ~artifacts["COMPLETE_NOW"]
        ]
        .copy()
        .reset_index(drop=True)
    )

    if not pending.empty:
        pending.insert(
            0,
            "TASK_ID",
            range(1, len(pending) + 1),
        )

    return (
        pending,
        registry,
        artifacts,
    )


def read_task(
    manifest: Path,
    task_id: int,
):
    x = pd.read_csv(
        manifest,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    ids = pd.to_numeric(
        x["TASK_ID"],
        errors="coerce",
    )

    hit = x[
        ids.eq(int(task_id))
    ]

    if len(hit) != 1:
        raise RuntimeError(
            f"TASK_ID={task_id}: expected 1 row, found {len(hit)}"
        )

    return hit.iloc[0].to_dict()


def aria2_attempt(
    url: str,
    destination: Path,
) -> None:
    if not shutil.which("aria2c"):
        raise RuntimeError(
            "aria2c is required."
        )

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    partial = Path(
        str(destination) + ".part"
    )

    cmd = [
        "aria2c",
        "--continue=true",
        "--allow-overwrite=true",
        "--auto-file-renaming=false",
        "--file-allocation=none",
        "--max-connection-per-server=1",
        "--split=1",
        "--retry-wait=10",
        "--max-tries=20",
        "--connect-timeout=60",
        "--timeout=60",
        "--summary-interval=60",
        "--console-log-level=notice",
        "--dir",
        str(destination.parent),
        "--out",
        partial.name,
        url,
    ]

    print(
        "$ "
        + " ".join(
            shlex.quote(str(x))
            for x in cmd
        ),
        flush=True,
    )

    subprocess.run(
        cmd,
        check=True,
    )

    if not nonempty(partial):
        raise RuntimeError(
            f"aria2c returned without a non-empty file: {partial}"
        )

    if str(destination).endswith(".gz"):
        if not gzip_ok(partial):
            raise RuntimeError(
                f"gzip integrity check failed: {partial}"
            )

    os.replace(
        partial,
        destination,
    )

    for control in [
        Path(str(partial) + ".aria2"),
        Path(str(destination) + ".aria2"),
    ]:
        try:
            if control.exists():
                control.unlink()
        except OSError:
            pass


def worker(
    manifest: Path,
    task_id: int,
    retry_sleep: int,
) -> None:
    row = read_task(
        manifest,
        task_id,
    )

    destination = Path(
        row["LOCAL_PATH"]
    )

    url = str(
        row["URL"]
    )

    kind = str(
        row["ARTIFACT_TYPE"]
    )

    banner(
        f"STEP11B DOWNLOAD ARRAY TASK {task_id}"
    )

    for key in [
        "STUDY_LABEL",
        "DATASET_ID",
        "TISSUE_LABEL",
        "CONDITION_LABEL",
        "QTL_TYPE",
        "ARTIFACT_TYPE",
    ]:
        print(
            f"{key:18s}: {row.get(key, '')}"
        )

    print(
        f"{'URL':18s}: {url}"
    )

    print(
        f"{'LOCAL_PATH':18s}: {destination}"
    )

    if artifact_complete(
        destination
    ):
        print(
            "[CACHE] Already complete."
        )
        return

    attempt = 0

    while True:
        attempt += 1

        try:
            print(
                f"[ATTEMPT {attempt}] "
                f"{kind} -> {destination}",
                flush=True,
            )

            aria2_attempt(
                url,
                destination,
            )

            if not artifact_complete(
                destination
            ):
                raise RuntimeError(
                    "Final artifact failed completion check."
                )

            print()
            print(
                "[COMPLETE]"
            )

            print(
                f"Size: {destination.stat().st_size / 1024**3:.3f} GB"
            )

            return

        except KeyboardInterrupt:
            raise

        except Exception as exc:
            print(
                f"[RETRY] {type(exc).__name__}: {exc}",
                flush=True,
            )

            print(
                f"[RETRY] sleeping {retry_sleep}s; "
                "the next attempt resumes the same .part file",
                flush=True,
            )

            time.sleep(
                max(1, int(retry_sleep))
            )


def verify(
    registry_file: Path,
) -> bool:
    registry = pd.read_csv(
        registry_file,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    tabix = shutil.which(
        "tabix"
    )

    rows = []
    all_paths = []

    for _, row in registry.iterrows():
        dense = Path(
            row["DENSE_PATH"]
        )

        index = Path(
            row["DENSE_INDEX_PATH"]
        )

        lbf = Path(
            row["PROVIDER_LBF_PATH"]
        )

        cs = Path(
            row["PROVIDER_CS_PATH"]
        )

        all_paths.extend(
            [
                dense,
                index,
                lbf,
                cs,
            ]
        )

        dense_ok = artifact_complete(
            dense
        )

        index_ok = artifact_complete(
            index
        )

        lbf_ok = artifact_complete(
            lbf
        )

        cs_ok = artifact_complete(
            cs
        )

        tabix_ok = False

        if dense_ok and index_ok:
            if tabix:
                proc = subprocess.run(
                    [
                        tabix,
                        "-l",
                        str(dense),
                    ],
                    text=True,
                    capture_output=True,
                )

                tabix_ok = (
                    proc.returncode == 0
                    and bool(proc.stdout.strip())
                )
            else:
                tabix_ok = True

        ready = (
            dense_ok
            and index_ok
            and tabix_ok
            and lbf_ok
            and cs_ok
        )

        rows.append(
            {
                "DATASET_ID": row["DATASET_ID"],
                "STUDY_LABEL": row["STUDY_LABEL"],
                "TISSUE_LABEL": row.get("TISSUE_LABEL", ""),
                "CONDITION_LABEL": row.get("CONDITION_LABEL", ""),
                "QTL_TYPE": row["QTL_TYPE"],
                "DENSE_OK": dense_ok,
                "INDEX_OK": index_ok,
                "TABIX_OK": tabix_ok,
                "LBF_OK": lbf_ok,
                "CS_OK": cs_ok,
                "READY_FOR_STEP11B": ready,
            }
        )

    status = pd.DataFrame(
        rows
    )

    status_file = (
        registry_file.parent
        / "resource_status.tsv"
    )

    status.to_csv(
        status_file,
        sep="\t",
        index=False,
    )

    unique_paths = list(
        dict.fromkeys(
            map(str, all_paths)
        )
    )

    physical_total = len(
        unique_paths
    )

    physical_complete = sum(
        artifact_complete(
            Path(p)
        )
        for p in unique_paths
    )

    n_ready = int(
        status[
            "READY_FOR_STEP11B"
        ].sum()
    )

    banner(
        "STEP11B SINGLE-CELL RESOURCE STATUS"
    )

    print(
        status.to_string(
            index=False
        )
    )

    print()
    print(
        "=" * 100
    )
    print(
        "BIOLOGICAL CONTEXTS"
    )
    print(
        "=" * 100
    )

    print(
        f"Contexts total       : {len(status)}"
    )
    print(
        f"Contexts ready       : {n_ready}"
    )
    print(
        f"Contexts not ready   : {len(status) - n_ready}"
    )

    print()
    print(
        "=" * 100
    )
    print(
        "PHYSICAL FILES"
    )
    print(
        "=" * 100
    )

    print(
        f"Files expected       : {physical_total}"
    )
    print(
        f"Files complete       : {physical_complete}"
    )
    print(
        f"Files remaining      : {physical_total - physical_complete}"
    )

    print()
    print(
        f"Status table         : {status_file}"
    )

    all_ready = (
        len(status) > 0
        and bool(
            status[
                "READY_FOR_STEP11B"
            ].all()
        )
    )

    marker = (
        registry_file.parent
        / ".complete"
    )

    if all_ready:
        write_json(
            marker,
            {
                "STATUS": "COMPLETE",
                "N_CONTEXTS": len(status),
                "N_PHYSICAL_FILES": physical_total,
                "VERIFIED_UNIX": time.time(),
            },
        )

    elif marker.exists():
        marker.unlink()

    return all_ready


def sbatch(
    script: Path,
    dependency: str = "",
) -> str:
    cmd = [
        "sbatch",
        "--parsable",
    ]

    if dependency:
        cmd.append(
            f"--dependency={dependency}"
        )

    cmd.append(
        str(script)
    )

    out = subprocess.check_output(
        cmd,
        text=True,
    ).strip()

    return (
        out
        .split(";", 1)[0]
        .strip()
    )


def submit_arrays(
    root: Path,
    pending: pd.DataFrame,
    registry_file: Path,
    partition: str,
    walltime: str,
    mem: str,
    cpus: int,
    array_limit: int,
    retry_sleep: int,
) -> None:
    if pending.empty:
        print(
            "[CACHE] No missing physical files."
        )

        verify(
            registry_file
        )

        return

    work = (
        registry_file.parent
    )

    manifest = (
        work
        / "download_manifest.tsv"
    )

    pending.to_csv(
        manifest,
        sep="\t",
        index=False,
    )

    logs = (
        root
        / "logs"
        / "resource_setup"
        / "step11b_singlecell"
    )

    logs.mkdir(
        parents=True,
        exist_ok=True,
    )

    scripts = (
        work
        / "slurm"
    )

    scripts.mkdir(
        parents=True,
        exist_ok=True,
    )

    state_file = (
        work
        / "active_download_submission.json"
    )

    old = {}

    if state_file.exists():
        try:
            old = json.loads(
                state_file.read_text()
            )
        except Exception:
            old = {}

    active = [
        str(j)
        for j in old.get(
            "ARRAY_JOB_IDS",
            [],
        )
        if active_job(str(j))
    ]

    if active:
        banner(
            "DOWNLOAD ARRAY ALREADY ACTIVE"
        )

        print(
            "Active job(s): "
            + ", ".join(active)
        )

        print(
            "No duplicate array submitted."
        )

        return

    py = Path(
        sys.executable
    ).resolve()

    program = Path(
        __file__
    ).resolve()

    array_limit = min(
        1000,
        max(
            1,
            int(array_limit),
        ),
    )

    total = len(
        pending
    )

    job_ids = []

    for (
        batch_no,
        global_start,
    ) in enumerate(
        range(
            1,
            total + 1,
            array_limit,
        ),
        start=1,
    ):
        global_end = min(
            total,
            global_start
            + array_limit
            - 1,
        )

        n_tasks = (
            global_end
            - global_start
            + 1
        )

        offset = (
            global_start
            - 1
        )

        script = (
            scripts
            / f"download_batch_{batch_no:03d}.sh"
        )

        script.write_text(
            f"""#!/bin/bash
#SBATCH --job-name=ctQTLdl_{batch_no:03d}
#SBATCH --partition={partition}
#SBATCH --time={walltime}
#SBATCH --mem={mem}
#SBATCH --cpus-per-task={cpus}
#SBATCH --ntasks=1
#SBATCH --array=1-{n_tasks}
#SBATCH --output={logs}/download.B{batch_no:03d}.%A_%a.out
#SBATCH --error={logs}/download.B{batch_no:03d}.%A_%a.err

set -euo pipefail

cd {shlex.quote(str(root))}

GLOBAL_TASK=$(({offset} + SLURM_ARRAY_TASK_ID))

exec {shlex.quote(str(py))} {shlex.quote(str(program))} \\
    --worker \\
    --manifest {shlex.quote(str(manifest))} \\
    --task-id "$GLOBAL_TASK" \\
    --worker-retry-sleep {int(retry_sleep)}
""",
            encoding="utf-8",
        )

        script.chmod(
            0o755
        )

        job_id = sbatch(
            script
        )

        job_ids.append(
            job_id
        )

        print(
            f"[SUBMIT] files {global_start}-{global_end} "
            f"-> array job {job_id}"
        )

    verify_script = (
        scripts
        / "verify_after_download.sh"
    )

    verify_script.write_text(
        f"""#!/bin/bash
#SBATCH --job-name=ctQTLverify
#SBATCH --partition={partition}
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output={logs}/verify.%j.out
#SBATCH --error={logs}/verify.%j.err

set -euo pipefail

cd {shlex.quote(str(root))}

{shlex.quote(str(py))} {shlex.quote(str(program))} \\
    --verify-only \\
    --registry {shlex.quote(str(registry_file))}
""",
        encoding="utf-8",
    )

    verify_script.chmod(
        0o755
    )

    dependency = (
        "afterany:"
        + ":".join(job_ids)
    )

    verify_job = sbatch(
        verify_script,
        dependency=dependency,
    )

    write_json(
        state_file,
        {
            "VERSION": VERSION,
            "N_MISSING_PHYSICAL_FILES": total,
            "ARRAY_JOB_IDS": job_ids,
            "VERIFY_JOB_ID": verify_job,
            "SUBMITTED_UNIX": time.time(),
            "MANIFEST": str(manifest),
            "REGISTRY": str(registry_file),
        },
    )

    banner(
        "STEP11B DOWNLOAD ARRAY SUBMITTED"
    )

    print(
        f"Missing physical files : {total}"
    )

    print(
        f"Array job(s)           : {', '.join(job_ids)}"
    )

    if total <= 1000:
        print(
            f"Array indices          : 1-{total}"
        )

    print(
        "ONE ARRAY INDEX        : ONE PHYSICAL FILE"
    )

    print(
        "Array throttle         : NONE"
    )

    print(
        f"Verification job       : {verify_job}"
    )

    print(
        f"Manifest               : {manifest}"
    )


def arguments():
    p = argparse.ArgumentParser()

    p.add_argument(
        "--root",
        default=".",
    )

    p.add_argument(
        "--profile",
        choices=["singlecell"],
        default="singlecell",
    )

    p.add_argument(
        "--qtl-types",
        default="eQTL",
    )

    p.add_argument(
        "--partition",
        default="general",
    )

    p.add_argument(
        "--time",
        default="24:00:00",
    )

    p.add_argument(
        "--mem",
        default="2G",
    )

    p.add_argument(
        "--cpus",
        type=int,
        default=1,
    )

    p.add_argument(
        "--array-limit",
        type=int,
        default=DEFAULT_ARRAY_LIMIT,
    )

    p.add_argument(
        "--worker-retry-sleep",
        type=int,
        default=30,
    )

    p.add_argument(
        "--plan-only",
        action="store_true",
    )

    p.add_argument(
        "--worker",
        action="store_true",
    )

    p.add_argument(
        "--manifest",
        default="",
    )

    p.add_argument(
        "--task-id",
        type=int,
    )

    p.add_argument(
        "--verify-only",
        action="store_true",
    )

    p.add_argument(
        "--registry",
        default="",
    )

    return p.parse_args()


def main() -> None:
    args = arguments()

    root = Path(
        args.root
    ).resolve()

    if args.worker:
        if (
            not args.manifest
            or args.task_id is None
        ):
            raise SystemExit(
                "--worker requires --manifest and --task-id"
            )

        worker(
            Path(args.manifest),
            int(args.task_id),
            int(args.worker_retry_sleep),
        )

        return

    if args.verify_only:
        if not args.registry:
            raise SystemExit(
                "--verify-only requires --registry"
            )

        ok = verify(
            Path(args.registry)
        )

        raise SystemExit(
            0 if ok else 2
        )

    qtl_types = {
        x.strip()
        for x in args.qtl_types.split(",")
        if x.strip()
    }

    banner(
        "GWAS2m SETUP 11B -- SINGLE-CELL eQTL"
    )

    print(
        f"Version        : {VERSION}"
    )

    print(
        f"Root           : {root}"
    )

    print(
        f"Profile        : {args.profile}"
    )

    print(
        f"QTL types      : {sorted(qtl_types)}"
    )

    print(
        f"Partition      : {args.partition}"
    )

    print(
        "Scheduling     : ONE SLURM ARRAY JOB"
    )

    print(
        "Task unit      : ONE ARRAY INDEX = ONE PHYSICAL FILE"
    )

    print(
        "Array throttle : NONE"
    )

    (
        contexts,
        excluded,
        selected,
        unavailable,
    ) = build_contexts(
        root,
        args.profile,
        qtl_types,
    )

    banner(
        "SELECTED SINGLE-CELL STUDIES"
    )

    for study in selected:
        print(
            f"[KEEP] {study}"
        )

    for family in unavailable:
        print(
            "[NOT IN LOCAL METADATA] "
            + " | ".join(family)
        )

    banner(
        "PROVIDER-SUSIE COMPATIBILITY"
    )

    print(
        f"Biological contexts      : {len(contexts)}"
    )

    print(
        f"Excluded contexts        : {len(excluded)}"
    )

    if contexts.empty:
        raise SystemExit(
            "No provider-compatible single-cell contexts."
        )

    (
        pending,
        registry,
        artifacts,
    ) = build_artifacts(
        root,
        contexts,
    )

    work = (
        root
        / "resources"
        / "coloc"
        / "eqtl_catalogue"
        / "context_specific_sets"
        / "singlecell"
    )

    work.mkdir(
        parents=True,
        exist_ok=True,
    )

    registry_file = (
        work
        / "context_registry.tsv"
    )

    registry.to_csv(
        registry_file,
        sep="\t",
        index=False,
    )

    contexts.to_csv(
        work
        / "eligible_contexts.tsv",
        sep="\t",
        index=False,
    )

    excluded.to_csv(
        work
        / "excluded_contexts.tsv",
        sep="\t",
        index=False,
    )

    artifacts.to_csv(
        work
        / "all_artifacts.tsv",
        sep="\t",
        index=False,
    )

    n_complete = int(
        artifacts["COMPLETE_NOW"].sum()
    )

    banner(
        "DOWNLOAD PLAN"
    )

    print(
        f"Biological contexts      : {len(registry)}"
    )

    print(
        f"Physical files expected  : {len(artifacts)}"
    )

    print(
        f"Physical files complete  : {n_complete}"
    )

    print(
        f"Physical files missing   : {len(pending)}"
    )

    print(
        f"SLURM array indices      : {len(pending)}"
    )

    print(
        "ONE ARRAY INDEX         : ONE PHYSICAL FILE"
    )

    print(
        f"Registry                 : {registry_file}"
    )

    if not pending.empty:
        print()
        print(
            pending[
                "ARTIFACT_TYPE"
            ]
            .value_counts()
            .to_string()
        )

    if args.plan_only:
        print()
        print(
            "[PLAN ONLY] No jobs submitted."
        )
        return

    submit_arrays(
        root=root,
        pending=pending,
        registry_file=registry_file,
        partition=args.partition,
        walltime=args.time,
        mem=args.mem,
        cpus=args.cpus,
        array_limit=args.array_limit,
        retry_sleep=args.worker_retry_sleep,
    )


if __name__ == "__main__":
    main()
