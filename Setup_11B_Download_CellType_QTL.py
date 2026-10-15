#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GWAS2m SETUP 11B -- DOWNLOAD CELL-TYPE QTL RESOURCES
=====================================================

Purpose
-------
Prepare all local resources required by:

    Step11B_CellType_QTL_Colocalization.py

For every eQTL Catalogue context matching the requested study/QTL type,
this script downloads:

  1. dense summary statistics:
       <dataset>.all.tsv.gz

  2. tabix index:
       <dataset>.all.tsv.gz.tbi

  3. provider SuSiE LBF file
       if advertised by eQTL Catalogue metadata

  4. provider SuSiE credible-set file
       if advertised by eQTL Catalogue metadata

Parallel architecture
---------------------
ONE downloadable file = ONE SLURM ARRAY TASK.

Example:

    30 dense files
    + 30 .tbi files
    + provider LBF files
    + provider credible-set files

Each file becomes an independent SLURM task.

NO %N throttle is used.

If >1000 files are required, the downloader automatically splits
the work into multiple arrays of at most 1000 tasks.

Existing valid files are skipped.

Downloads are resumable and use temporary .part files before
atomic rename.

Default target
--------------
Jerber_2021 cell-type eQTL resources.

This intentionally downloads ALL matching Jerber cell contexts:
dopaminergic neurons, neurons, astrocytes, serotonergic neurons,
ependymal cells, floor-plate progenitors, neuroblasts, etc.

Typical use
-----------

    python Setup_11B_Download_CellType_QTL.py

Equivalent explicit form:

    python Setup_11B_Download_CellType_QTL.py \
        --study-regex '^Jerber_2021$' \
        --qtl-types eQTL \
        --partition general

The default controller mode automatically submits the SLURM arrays.
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
from types import SimpleNamespace
from urllib.parse import urlparse

import pandas as pd


# =============================================================================
# GWAS2m imports
# =============================================================================

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


VERSION = "1.0.0-celltype-qtl-array-download"

DEFAULT_ARRAY_LIMIT = 1000


# =============================================================================
# Generic helpers
# =============================================================================

def banner(text: str) -> None:

    print()
    print("=" * 120)
    print(text)
    print("=" * 120, flush=True)


def slug(value: str) -> str:

    return core.slugify(
        str(value)
    )


def nonempty(
    path: Path,
) -> bool:

    return (
        path.exists()
        and path.is_file()
        and path.stat().st_size > 0
    )


def clean_url(
    value,
) -> str:

    value = str(
        value or ""
    ).strip()

    if (
        not value
        or value.lower()
        in {
            "nan",
            "na",
            "none",
            ".",
        }
    ):
        return ""

    return value


def https_if_ebi_ftp(
    url: str,
) -> str:

    """
    Convert EBI ftp:// URLs to HTTPS.

    This is generally more reliable from compute nodes.
    """

    url = str(
        url
    ).strip()

    prefix = (
        "ftp://ftp.ebi.ac.uk/"
    )

    if url.startswith(
        prefix
    ):

        return (
            "https://ftp.ebi.ac.uk/"
            + url[len(prefix):]
        )

    return url


def url_basename(
    url: str,
) -> str:

    clean = url.split(
        "?",
        1,
    )[0]

    path = urlparse(
        clean
    ).path.rstrip("/")

    return Path(
        path
    ).name


# =============================================================================
# Validation
# =============================================================================

def gzip_test(
    path: Path,
) -> bool:

    if not nonempty(
        path
    ):
        return False

    if not str(
        path
    ).endswith(
        ".gz"
    ):
        return True

    try:

        with gzip.open(
            path,
            "rb",
        ) as fh:

            while fh.read(
                8 * 1024 * 1024
            ):
                pass

        return True

    except Exception:

        return False


def artifact_valid(
    kind: str,
    path: Path,
) -> bool:

    if not nonempty(
        path
    ):
        return False

    # aria2 incomplete marker
    if Path(
        str(path) + ".aria2"
    ).exists():

        return False

    if kind in {
        "DENSE_SUMSTATS",
        "PROVIDER_LBF",
        "PROVIDER_CS",
    }:

        if str(
            path
        ).endswith(
            ".gz"
        ):

            return gzip_test(
                path
            )

        return True

    if kind == "DENSE_INDEX":

        return (
            path.stat().st_size
            > 16
        )

    return True


# =============================================================================
# Command runner
# =============================================================================

def run(
    cmd: list[str],
) -> subprocess.CompletedProcess:

    print(
        "$ "
        + " ".join(
            shlex.quote(
                str(x)
            )
            for x in cmd
        ),
        flush=True,
    )

    return subprocess.run(
        [
            str(x)
            for x in cmd
        ],
        check=True,
        text=True,
    )


# =============================================================================
# Downloader
# =============================================================================

def download_one(
    url: str,
    target: Path,
) -> None:

    target.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    url = https_if_ebi_ftp(
        url
    )

    part = Path(
        str(target)
        + ".part"
    )

    # -------------------------------------------------------------------------
    # Prefer aria2c
    # -------------------------------------------------------------------------

    if shutil.which(
        "aria2c"
    ):

        cmd = [

            "aria2c",

            "--continue=true",

            "--allow-overwrite=true",

            "--auto-file-renaming=false",

            "--file-allocation=none",

            "--max-connection-per-server=4",

            "--split=4",

            "--min-split-size=20M",

            "--retry-wait=5",

            "--max-tries=10",

            "--timeout=60",

            "--dir",
            str(
                target.parent
            ),

            "--out",
            part.name,

            url,
        ]

        run(
            cmd
        )

    # -------------------------------------------------------------------------
    # curl fallback
    # -------------------------------------------------------------------------

    elif shutil.which(
        "curl"
    ):

        cmd = [

            "curl",

            "-fL",

            "--retry",
            "10",

            "--retry-delay",
            "5",

            "--connect-timeout",
            "60",

            "-C",
            "-",

            "-o",
            str(
                part
            ),

            url,
        ]

        run(
            cmd
        )

    # -------------------------------------------------------------------------
    # wget fallback
    # -------------------------------------------------------------------------

    elif shutil.which(
        "wget"
    ):

        cmd = [

            "wget",

            "--continue",

            "--tries=10",

            "--timeout=60",

            "-O",
            str(
                part
            ),

            url,
        ]

        run(
            cmd
        )

    else:

        raise RuntimeError(
            "Need aria2c, curl, or wget."
        )

    if not nonempty(
        part
    ):

        raise RuntimeError(
            f"Download produced empty file: {part}"
        )

    # Atomic completion.
    os.replace(
        part,
        target,
    )


# =============================================================================
# Manifest helpers
# =============================================================================

def read_manifest_row(
    manifest: Path,
    task_id: int,
) -> dict[str, str]:

    df = pd.read_csv(
        manifest,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    if (
        "TASK_ID"
        not in df.columns
    ):

        raise RuntimeError(
            f"Manifest lacks TASK_ID: {manifest}"
        )

    ids = pd.to_numeric(
        df["TASK_ID"],
        errors="coerce",
    )

    hit = df[
        ids
        == int(
            task_id
        )
    ]

    if len(
        hit
    ) != 1:

        raise RuntimeError(
            "Expected one manifest row for "
            f"TASK_ID={task_id}; "
            f"found {len(hit)}"
        )

    return (
        hit.iloc[0]
        .to_dict()
    )


# =============================================================================
# Worker
# =============================================================================

def worker(
    manifest: Path,
    task_id: int,
) -> None:

    row = read_manifest_row(
        manifest,
        task_id,
    )

    kind = row[
        "ARTIFACT_TYPE"
    ]

    url = row[
        "URL"
    ]

    target = Path(
        row[
            "LOCAL_PATH"
        ]
    )

    banner(
        f"CELL-TYPE QTL DOWNLOAD TASK {task_id}"
    )

    print(
        f"Dataset       : "
        f"{row['DATASET_ID']}"
    )

    print(
        f"Study         : "
        f"{row['STUDY_LABEL']}"
    )

    print(
        f"Cell/context  : "
        f"{row['TISSUE_LABEL']}"
    )

    print(
        f"QTL type      : "
        f"{row['QTL_TYPE']}"
    )

    print(
        f"Artifact      : "
        f"{kind}"
    )

    print(
        f"URL           : "
        f"{url}"
    )

    print(
        f"Local         : "
        f"{target}"
    )

    # -------------------------------------------------------------------------
    # Resume
    # -------------------------------------------------------------------------

    if artifact_valid(
        kind,
        target,
    ):

        print()
        print(
            "[CACHE] Existing valid file."
        )

        return

    # Existing invalid file is preserved for inspection.
    if target.exists():

        bad = target.with_name(
            target.name
            + f".invalid.{int(time.time())}"
        )

        print(
            "[WARN] Existing file failed "
            f"validation; moving to {bad}"
        )

        target.rename(
            bad
        )

    # -------------------------------------------------------------------------
    # Download
    # -------------------------------------------------------------------------

    download_one(
        url,
        target,
    )

    # -------------------------------------------------------------------------
    # Validate
    # -------------------------------------------------------------------------

    if not artifact_valid(
        kind,
        target,
    ):

        raise RuntimeError(
            "Downloaded file failed "
            f"validation: {target}"
        )

    print()
    print(
        f"[OK] {kind}"
    )

    print(
        f"     {target}"
    )

    print(
        f"     "
        f"{target.stat().st_size / 1024**2:.2f} MB"
    )


# =============================================================================
# Build download plan
# =============================================================================

def build_plans(
    root: Path,
    study_regex: str,
    qtl_types: str,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:

    allowed = {

        x.strip().lower()

        for x
        in qtl_types.split(",")

        if x.strip()
    }

    if not allowed:

        raise RuntimeError(
            "No QTL types requested."
        )

    # -------------------------------------------------------------------------
    # Read local eQTL Catalogue metadata
    # -------------------------------------------------------------------------

    catalogue = (
        core.load_eqtl_catalogue(
            root,
            False,
        )
    )

    selector_args = (
        SimpleNamespace(
            catalog_study_regex=study_regex,
            catalog_all_studies=False,
        )
    )

    contexts = (
        core.select_catalogue_contexts(
            catalogue,
            selector_args,
            allowed,
        )
    )

    if contexts.empty:

        raise RuntimeError(
            "No eQTL Catalogue contexts "
            "matched:\n"
            f"study_regex={study_regex!r}\n"
            f"qtl_types={qtl_types!r}"
        )

    dense_root = (

        root
        / "resources"
        / "coloc"
        / "eqtl_catalogue"
        / "dense"

    )

    artifacts = []

    registry = []

    # -------------------------------------------------------------------------
    # Helper
    # -------------------------------------------------------------------------

    def add_artifact(
        dataset_id: str,
        study_label: str,
        tissue_label: str,
        qtl_type: str,
        kind: str,
        url: str,
        local: Path,
    ) -> None:

        url = clean_url(
            url
        )

        if not url:
            return

        artifacts.append(
            {
                "DATASET_ID":
                    dataset_id,

                "STUDY_LABEL":
                    study_label,

                "TISSUE_LABEL":
                    tissue_label,

                "QTL_TYPE":
                    qtl_type,

                "ARTIFACT_TYPE":
                    kind,

                "URL":
                    https_if_ebi_ftp(
                        url
                    ),

                "LOCAL_PATH":
                    str(
                        local
                    ),
            }
        )

    # -------------------------------------------------------------------------
    # Every eQTL Catalogue context
    # -------------------------------------------------------------------------

    for _, row in contexts.iterrows():

        dataset_id = str(
            row.get(
                "dataset_id",
                "",
            )
        ).strip()

        study_label = str(
            row.get(
                "study_label",
                "",
            )
        ).strip()

        tissue_label = str(
            row.get(
                "tissue_label",
                row.get(
                    "tissue",
                    "",
                ),
            )
        ).strip()

        qtl_type = str(
            row.get(
                "QTL_TYPE",
                "",
            )
        ).strip()

        dense_url = clean_url(
            row.get(
                "ftp_path",
                "",
            )
        )

        lbf_url = clean_url(
            row.get(
                "ftp_lbf_path",
                "",
            )
        )

        cs_url = clean_url(
            row.get(
                "ftp_cs_path",
                "",
            )
        )

        if (
            not dataset_id
            or not dense_url
        ):

            continue

        # ---------------------------------------------------------------------
        # Dense summary statistics
        # ---------------------------------------------------------------------

        dense_name = url_basename(
            dense_url
        )

        dense_local = (

            dense_root
            / dataset_id
            / dense_name

        )

        index_url = (
            dense_url
            + ".tbi"
        )

        index_local = Path(
            str(
                dense_local
            )
            + ".tbi"
        )

        # ---------------------------------------------------------------------
        # Provider SuSiE files
        # ---------------------------------------------------------------------

        lbf_local = None

        cs_local = None

        if lbf_url:

            lbf_local = (
                provider.provider_local_path(
                    root,
                    dataset_id,
                    lbf_url,
                )
            )

        if cs_url:

            cs_local = (
                provider.provider_local_path(
                    root,
                    dataset_id,
                    cs_url,
                )
            )

        # ---------------------------------------------------------------------
        # Add each FILE as an independent task
        # ---------------------------------------------------------------------

        add_artifact(
            dataset_id,
            study_label,
            tissue_label,
            qtl_type,
            "DENSE_SUMSTATS",
            dense_url,
            dense_local,
        )

        add_artifact(
            dataset_id,
            study_label,
            tissue_label,
            qtl_type,
            "DENSE_INDEX",
            index_url,
            index_local,
        )

        if (
            lbf_url
            and lbf_local is not None
        ):

            add_artifact(
                dataset_id,
                study_label,
                tissue_label,
                qtl_type,
                "PROVIDER_LBF",
                lbf_url,
                lbf_local,
            )

        if (
            cs_url
            and cs_local is not None
        ):

            add_artifact(
                dataset_id,
                study_label,
                tissue_label,
                qtl_type,
                "PROVIDER_CS",
                cs_url,
                cs_local,
            )

        # ---------------------------------------------------------------------
        # Context registry
        # ---------------------------------------------------------------------

        registry.append(
            {
                "DATASET_ID":
                    dataset_id,

                "STUDY_LABEL":
                    study_label,

                "TISSUE_LABEL":
                    tissue_label,

                "QTL_TYPE":
                    qtl_type,

                "DENSE_URL":
                    https_if_ebi_ftp(
                        dense_url
                    ),

                "DENSE_PATH":
                    str(
                        dense_local
                    ),

                "INDEX_URL":
                    https_if_ebi_ftp(
                        index_url
                    ),

                "INDEX_PATH":
                    str(
                        index_local
                    ),

                "PROVIDER_LBF_URL":
                    (
                        https_if_ebi_ftp(
                            lbf_url
                        )
                        if lbf_url
                        else ""
                    ),

                "PROVIDER_LBF_PATH":
                    str(
                        lbf_local or ""
                    ),

                "PROVIDER_CS_URL":
                    (
                        https_if_ebi_ftp(
                            cs_url
                        )
                        if cs_url
                        else ""
                    ),

                "PROVIDER_CS_PATH":
                    str(
                        cs_local or ""
                    ),
            }
        )

    if not registry:

        raise RuntimeError(
            "No downloadable cell-type "
            "QTL contexts were produced."
        )

    # -------------------------------------------------------------------------
    # Remove duplicate file targets
    # -------------------------------------------------------------------------

    artifacts_df = (
        pd.DataFrame(
            artifacts
        )
        .drop_duplicates(
            subset=[
                "LOCAL_PATH"
            ],
            keep="first",
        )
        .reset_index(
            drop=True
        )
    )

    registry_df = (
        pd.DataFrame(
            registry
        )
        .drop_duplicates(
            subset=[
                "DATASET_ID"
            ],
            keep="first",
        )
        .reset_index(
            drop=True
        )
    )

    # -------------------------------------------------------------------------
    # Determine which files actually need downloading
    # -------------------------------------------------------------------------

    artifacts_df[
        "VALID_NOW"
    ] = [

        (
            "YES"
            if artifact_valid(
                kind,
                Path(
                    path
                ),
            )
            else "NO"
        )

        for kind, path
        in zip(
            artifacts_df[
                "ARTIFACT_TYPE"
            ],
            artifacts_df[
                "LOCAL_PATH"
            ],
        )
    ]

    pending = (
        artifacts_df[
            artifacts_df[
                "VALID_NOW"
            ].eq(
                "NO"
            )
        ]
        .copy()
        .reset_index(
            drop=True
        )
    )

    if not pending.empty:

        pending.insert(
            0,
            "TASK_ID",
            range(
                1,
                len(
                    pending
                )
                + 1,
            ),
        )

    else:

        pending = (
            pd.DataFrame(
                columns=[
                    "TASK_ID",
                    *artifacts_df.columns,
                ]
            )
        )

    return (
        pending,
        registry_df,
    )


# =============================================================================
# Full resource verification
# =============================================================================

def verify_registry(
    registry_file: Path,
    verbose: bool = True,
) -> bool:

    registry = pd.read_csv(
        registry_file,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    rows = []

    tabix = shutil.which(
        "tabix"
    )

    for _, row in registry.iterrows():

        dense = Path(
            row[
                "DENSE_PATH"
            ]
        )

        index = Path(
            row[
                "INDEX_PATH"
            ]
        )

        lbf_url = clean_url(
            row.get(
                "PROVIDER_LBF_URL",
                "",
            )
        )

        cs_url = clean_url(
            row.get(
                "PROVIDER_CS_URL",
                "",
            )
        )

        lbf = (
            Path(
                row[
                    "PROVIDER_LBF_PATH"
                ]
            )
            if lbf_url
            else None
        )

        cs = (
            Path(
                row[
                    "PROVIDER_CS_PATH"
                ]
            )
            if cs_url
            else None
        )

        dense_ok = artifact_valid(
            "DENSE_SUMSTATS",
            dense,
        )

        index_ok = artifact_valid(
            "DENSE_INDEX",
            index,
        )

        lbf_ok = (
            artifact_valid(
                "PROVIDER_LBF",
                lbf,
            )
            if lbf is not None
            else False
        )

        cs_ok = (
            artifact_valid(
                "PROVIDER_CS",
                cs,
            )
            if cs is not None
            else False
        )

        # ---------------------------------------------------------------------
        # Confirm data/index pair is actually readable by tabix
        # ---------------------------------------------------------------------

        index_pair_ok = False

        if (
            dense_ok
            and index_ok
        ):

            if tabix:

                proc = subprocess.run(

                    [
                        tabix,
                        "-l",
                        str(
                            dense
                        ),
                    ],

                    text=True,

                    capture_output=True,
                )

                index_pair_ok = (
                    proc.returncode
                    == 0
                    and bool(
                        proc.stdout.strip()
                    )
                )

            else:

                # Existence validation only.
                index_pair_ok = True

        # ---------------------------------------------------------------------
        # Provider status
        # ---------------------------------------------------------------------

        if not lbf_url:

            provider_status = (
                "NO_LBF_URL_IN_METADATA"
            )

        elif not cs_url:

            provider_status = (
                "NO_CS_URL_IN_METADATA"
            )

        elif (
            lbf_ok
            and cs_ok
        ):

            provider_status = (
                "READY"
            )

        else:

            provider_status = (
                "MISSING_OR_INVALID"
            )

        ready = (

            dense_ok
            and index_ok
            and index_pair_ok
            and lbf_ok
            and cs_ok

        )

        rows.append(
            {
                "DATASET_ID":
                    row[
                        "DATASET_ID"
                    ],

                "STUDY_LABEL":
                    row[
                        "STUDY_LABEL"
                    ],

                "TISSUE_LABEL":
                    row[
                        "TISSUE_LABEL"
                    ],

                "QTL_TYPE":
                    row[
                        "QTL_TYPE"
                    ],

                "DENSE_OK":
                    dense_ok,

                "INDEX_OK":
                    index_ok,

                "TABIX_PAIR_OK":
                    index_pair_ok,

                "PROVIDER_STATUS":
                    provider_status,

                "READY_FOR_STEP11B":
                    ready,
            }
        )

    result = pd.DataFrame(
        rows
    )

    status_file = (
        registry_file
        .with_name(
            "celltype_qtl_resource_status.tsv"
        )
    )

    result.to_csv(
        status_file,
        sep="\t",
        index=False,
    )

    if verbose:

        banner(
            "CELL-TYPE QTL RESOURCE VERIFICATION"
        )

        print(
            result.to_string(
                index=False
            )
        )

        print()

        print(
            f"Contexts          : "
            f"{len(result)}"
        )

        print(
            f"Ready Step11B     : "
            f"{int(result['READY_FOR_STEP11B'].sum())}"
        )

        print(
            f"Not ready         : "
            f"{int((~result['READY_FOR_STEP11B']).sum())}"
        )

        print(
            f"Status table      : "
            f"{status_file}"
        )

    all_ready = (

        bool(
            len(
                result
            )
        )

        and bool(
            result[
                "READY_FOR_STEP11B"
            ].all()
        )

    )

    marker = (
        registry_file
        .with_name(
            ".complete"
        )
    )

    if all_ready:

        marker.write_text(

            json.dumps(
                {
                    "STATUS":
                        "COMPLETE",

                    "N_CONTEXTS":
                        len(
                            result
                        ),

                    "VERIFIED_UNIX":
                        time.time(),
                },

                indent=2,
            )

            + "\n",

            encoding="utf-8",
        )

        print()
        print(
            "CELL-TYPE QTL RESOURCE SET: COMPLETE"
        )

    else:

        if marker.exists():

            marker.unlink()

        print()
        print(
            "CELL-TYPE QTL RESOURCE SET: INCOMPLETE"
        )

    return all_ready


# =============================================================================
# SLURM
# =============================================================================

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
        str(
            script
        )
    )

    out = (
        subprocess.check_output(
            cmd,
            text=True,
        )
        .strip()
    )

    return (
        out.split(
            ";",
            1,
        )[0]
    )


def submit_downloads(
    root: Path,
    pending: pd.DataFrame,
    registry_file: Path,
    partition: str,
    walltime: str,
    mem: str,
    cpus: int,
    array_limit: int,
) -> None:

    # -------------------------------------------------------------------------
    # Nothing missing
    # -------------------------------------------------------------------------

    if pending.empty:

        print(
            "[CACHE] No missing artifacts."
        )

        verify_registry(
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
        / "celltype_qtl"
        / slug(
            work.name
        )

    )

    logs.mkdir(
        parents=True,
        exist_ok=True,
    )

    script_dir = (
        work
        / "slurm"
    )

    script_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    py = Path(
        sys.executable
    ).resolve()

    this_script = Path(
        __file__
    ).resolve()

    job_ids = []

    total = len(
        pending
    )

    array_limit = max(
        1,
        int(
            array_limit
        ),
    )

    # -------------------------------------------------------------------------
    # Split arrays if >1000 tasks
    # -------------------------------------------------------------------------

    for batch_no, start in enumerate(

        range(
            1,
            total + 1,
            array_limit,
        ),

        start=1,
    ):

        end = min(
            total,
            start
            + array_limit
            - 1,
        )

        n_tasks = (
            end
            - start
            + 1
        )

        offset = (
            start
            - 1
        )

        script = (

            script_dir
            / (
                f"download_batch_"
                f"{batch_no:03d}.sh"
            )

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

{shlex.quote(str(py))} {shlex.quote(str(this_script))} \\
    --worker \\
    --manifest {shlex.quote(str(manifest))} \\
    --task-id "$GLOBAL_TASK"
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
            f"[SUBMIT] batch {batch_no}: "
            f"tasks {start}-{end} "
            f"-> job {job_id}"
        )

    # -------------------------------------------------------------------------
    # Verification job
    # -------------------------------------------------------------------------

    verify_script = (

        script_dir
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

{shlex.quote(str(py))} {shlex.quote(str(this_script))} \\
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
        + ":".join(
            job_ids
        )
    )

    verify_job = sbatch(
        verify_script,
        dependency=dependency,
    )

    banner(
        "CELL-TYPE QTL DOWNLOADS SUBMITTED"
    )

    print(
        f"Missing artifacts      : "
        f"{total}"
    )

    print(
        f"Download job(s)        : "
        f"{', '.join(job_ids)}"
    )

    print(
        f"Verification job       : "
        f"{verify_job}"
    )

    print(
        f"Array limit            : "
        f"{array_limit}"
    )

    print(
        "Array throttle         : NONE"
    )

    print(
        "ONE ARRAY TASK         : ONE FILE"
    )

    print(
        f"Manifest               : "
        f"{manifest}"
    )

    print(
        f"Registry               : "
        f"{registry_file}"
    )

    print(
        f"Logs                   : "
        f"{logs}"
    )


# =============================================================================
# CLI
# =============================================================================

def arguments():

    parser = argparse.ArgumentParser(

        description=(
            "Download all cell-type QTL resources "
            "required by GWAS2m Step11B."
        )

    )

    parser.add_argument(
        "--root",
        default=".",
    )

    parser.add_argument(
        "--study-regex",
        default=r"^Jerber_2021$",
    )

    parser.add_argument(
        "--qtl-types",
        default="eQTL",
    )

    parser.add_argument(
        "--partition",
        default="general",
    )

    parser.add_argument(
        "--time",
        default="24:00:00",
    )

    parser.add_argument(
        "--mem",
        default="2G",
    )

    parser.add_argument(
        "--cpus",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--array-limit",
        type=int,
        default=DEFAULT_ARRAY_LIMIT,
    )

    parser.add_argument(
        "--plan-only",
        action="store_true",
    )

    # -------------------------------------------------------------------------
    # Internal SLURM worker modes
    # -------------------------------------------------------------------------

    parser.add_argument(
        "--worker",
        action="store_true",
    )

    parser.add_argument(
        "--manifest",
        default="",
    )

    parser.add_argument(
        "--task-id",
        type=int,
    )

    parser.add_argument(
        "--verify-only",
        action="store_true",
    )

    parser.add_argument(
        "--registry",
        default="",
    )

    return parser.parse_args()


# =============================================================================
# Main
# =============================================================================

def main() -> None:

    args = arguments()

    root = Path(
        args.root
    ).resolve()

    # -------------------------------------------------------------------------
    # Worker
    # -------------------------------------------------------------------------

    if args.worker:

        if (
            not args.manifest
            or args.task_id is None
        ):

            raise SystemExit(
                "--worker requires "
                "--manifest and --task-id"
            )

        worker(
            Path(
                args.manifest
            ),
            int(
                args.task_id
            ),
        )

        return

    # -------------------------------------------------------------------------
    # Verification worker
    # -------------------------------------------------------------------------

    if args.verify_only:

        if not args.registry:

            raise SystemExit(
                "--verify-only requires --registry"
            )

        ok = verify_registry(
            Path(
                args.registry
            )
        )

        raise SystemExit(
            0
            if ok
            else 2
        )

    # -------------------------------------------------------------------------
    # Controller
    # -------------------------------------------------------------------------

    banner(
        "GWAS2m SETUP 11B -- "
        "CELL-TYPE QTL RESOURCE DOWNLOADER"
    )

    print(
        f"Version       : {VERSION}"
    )

    print(
        f"Root          : {root}"
    )

    print(
        f"Study regex   : {args.study_regex}"
    )

    print(
        f"QTL types     : {args.qtl_types}"
    )

    print(
        f"Partition     : {args.partition}"
    )

    print(
        f"Array limit   : {args.array_limit}"
    )

    print(
        "Throttle      : NONE"
    )

    print(
        "Task unit     : ONE FILE"
    )

    # -------------------------------------------------------------------------
    # Build plan
    # -------------------------------------------------------------------------

    pending, registry = build_plans(

        root,

        args.study_regex,

        args.qtl_types,

    )

    study_slug = slug(
        args.study_regex
        .strip(
            "^$"
        )
    )

    if not study_slug:

        study_slug = (
            "celltype_study"
        )

    work = (

        root
        / "resources"
        / "coloc"
        / "eqtl_catalogue"
        / "celltype_sets"
        / study_slug

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

    # -------------------------------------------------------------------------
    # Print plan
    # -------------------------------------------------------------------------

    banner(
        "DOWNLOAD PLAN"
    )

    print(
        f"Catalogue contexts     : "
        f"{len(registry)}"
    )

    print(
        f"Missing file artifacts : "
        f"{len(pending)}"
    )

    if not pending.empty:

        counts = (
            pending[
                "ARTIFACT_TYPE"
            ]
            .value_counts()
        )

        print()

        for kind, count in counts.items():

            print(
                f"{kind:18s}: "
                f"{count}"
            )

        print()

        print(
            "First pending artifacts:"
        )

        print(

            pending[
                [
                    "TASK_ID",
                    "DATASET_ID",
                    "TISSUE_LABEL",
                    "ARTIFACT_TYPE",
                    "LOCAL_PATH",
                ]
            ]
            .head(
                30
            )
            .to_string(
                index=False
            )

        )

    print()

    print(
        f"Registry: "
        f"{registry_file}"
    )

    # -------------------------------------------------------------------------
    # Plan only
    # -------------------------------------------------------------------------

    if args.plan_only:

        print()
        print(
            "[PLAN ONLY] Nothing submitted."
        )

        return

    # -------------------------------------------------------------------------
    # Submit
    # -------------------------------------------------------------------------

    submit_downloads(

        root=root,

        pending=pending,

        registry_file=registry_file,

        partition=args.partition,

        walltime=args.time,

        mem=args.mem,

        cpus=args.cpus,

        array_limit=args.array_limit,

    )


if __name__ == "__main__":

    main()