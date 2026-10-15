#!/usr/bin/env python3
"""
GWAS2m Step 07D Setup
Download compact ENCODE regulatory annotations for GRCh38.

Strategy
--------
1. Query ENCODE REST API dynamically.
2. Keep released GRCh38 BED/bed.gz regulatory interval files.
3. Never download FASTQ/BAM/bigWig for variant annotation.
4. Prefer processed peak files.
5. Download:
      - candidate cis-Regulatory Elements (cCRE)
      - DNase-seq peaks
      - ATAC-seq peaks
      - H3K27ac peaks
      - H3K4me1 peaks
      - H3K4me3 peaks
      - CTCF peaks
6. Check available disk space before downloading.
7. Resume partial downloads.
8. Write complete provenance manifest.

Usage
-----
python Setup_07D_Download_ENCODE.py --plan-only
python Setup_07D_Download_ENCODE.py
python Setup_07D_Download_ENCODE.py --core-only
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlencode

try:
    import requests
except ImportError:
    sys.exit(
        "ERROR: requests is required.\n"
        "Install with:\n"
        "  python -m pip install requests"
    )


ROOT = Path("/data/ascher02/uqmmune1/Splice2/GWAS2m")

RESOURCE_ROOT = ROOT / "resources" / "encode"
MANIFEST_DIR = RESOURCE_ROOT / "manifests"

CCRE_DIR = RESOURCE_ROOT / "ccre"
PEAK_DIR = RESOURCE_ROOT / "peaks"

API = "https://www.encodeproject.org"

ASSEMBLY = "GRCh38"

# Keep substantial free space because this filesystem is shared
# with GWAS/QTL/fine-mapping resources.
DEFAULT_RESERVE_GB = 50.0

HEADERS = {
    "accept": "application/json",
    "User-Agent": "GWAS2m-ENCODE/1.0"
}


def banner(title: str) -> None:
    print()
    print("=" * 120)
    print(title)
    print("=" * 120)


def human_bytes(n: int | float) -> str:
    n = float(n)

    units = ["B", "KiB", "MiB", "GiB", "TiB"]

    for unit in units:
        if abs(n) < 1024.0:
            return f"{n:,.2f} {unit}"
        n /= 1024.0

    return f"{n:,.2f} PiB"


def disk_status(path: Path) -> tuple[int, int, int]:
    path.mkdir(parents=True, exist_ok=True)

    usage = shutil.disk_usage(path)

    return usage.total, usage.used, usage.free


def get_json(url: str, params: dict[str, Any] | None = None) -> dict:
    for attempt in range(1, 6):
        try:
            r = requests.get(
                url,
                params=params,
                headers=HEADERS,
                timeout=120
            )

            r.raise_for_status()

            return r.json()

        except Exception as exc:
            print(
                f"[WARN] API request failed "
                f"(attempt {attempt}/5): {exc}"
            )

            if attempt == 5:
                raise

            time.sleep(attempt * 3)

    raise RuntimeError("unreachable")


def normalize_href(obj: dict) -> str | None:
    """
    ENCODE File objects normally expose href such as:
    /files/ENCFF.../@@download/ENCFF....bed.gz
    """

    href = obj.get("href")

    if href:
        if href.startswith("http"):
            return href

        return API + href

    accession = obj.get("accession") or obj.get("title")

    if not accession:
        return None

    submitted = obj.get("submitted_file_name", "")

    if submitted.endswith(".bed.gz"):
        suffix = ".bed.gz"
    elif submitted.endswith(".bed"):
        suffix = ".bed"
    else:
        suffix = ".bed.gz"

    return (
        f"{API}/files/{accession}/@@download/"
        f"{accession}{suffix}"
    )


def biosample_text(obj: dict) -> str:
    values = []

    for key in (
        "simple_biosample_summary",
        "biosample_summary",
        "biosample_ontology"
    ):
        value = obj.get(key)

        if isinstance(value, str):
            values.append(value)

        elif isinstance(value, dict):
            for k in (
                "term_name",
                "classification",
                "organ_slims",
                "system_slims"
            ):
                v = value.get(k)

                if isinstance(v, str):
                    values.append(v)

                elif isinstance(v, list):
                    values.extend(map(str, v))

    return " | ".join(values)


def query_files(params: dict[str, Any]) -> list[dict]:
    query = {
        "type": "File",
        "status": "released",
        "assembly": ASSEMBLY,
        "file_format": "bed",
        "limit": "all",
        "format": "json",
        "frame": "embedded",
    }

    query.update(params)

    banner("ENCODE API QUERY")

    print("Parameters:")

    for key, value in query.items():
        print(f"  {key:25s}: {value}")

    data = get_json(
        f"{API}/search/",
        params=query
    )

    graph = data.get("@graph", [])

    print(f"\nReturned objects: {len(graph):,}")

    return graph


def is_peak_file(obj: dict) -> bool:
    fmt = str(obj.get("file_format", "")).lower()

    if fmt != "bed":
        return False

    output = str(obj.get("output_type", "")).lower()
    fmt_type = str(obj.get("file_format_type", "")).lower()

    useful = (
        "peak" in output
        or "candidate cis-regulatory" in output
        or "narrowpeak" in fmt_type
        or "broadpeak" in fmt_type
    )

    return useful


def choose_preferred(files: list[dict]) -> list[dict]:
    """
    Remove obvious duplicate representations.

    Preference:
       preferred_default
       conservative IDR thresholded peaks
       optimal IDR thresholded peaks
       replicated peaks
       pseudoreplicated peaks
       peaks
    """

    ranking = {
        "conservative idr thresholded peaks": 0,
        "optimal idr thresholded peaks": 1,
        "idr thresholded peaks": 2,
        "replicated peaks": 3,
        "pseudoreplicated peaks": 4,
        "candidate cis-regulatory elements": 5,
        "peaks": 6,
    }

    grouped: dict[tuple, list[dict]] = {}

    for obj in files:

        if not is_peak_file(obj):
            continue

        dataset = str(obj.get("dataset", ""))
        target = ""

        target_obj = obj.get("target")

        if isinstance(target_obj, dict):
            target = str(
                target_obj.get("label")
                or target_obj.get("name")
                or ""
            )

        assay = str(
            obj.get("assay_title")
            or obj.get("assay_term_name")
            or ""
        )

        biosample = biosample_text(obj)

        key = (
            dataset,
            assay,
            target,
            biosample
        )

        grouped.setdefault(key, []).append(obj)

    selected = []

    for _, candidates in grouped.items():

        def score(obj: dict):
            output = str(
                obj.get("output_type", "")
            ).lower()

            preferred = bool(
                obj.get("preferred_default", False)
            )

            return (
                0 if preferred else 1,
                ranking.get(output, 50),
                int(obj.get("file_size") or 10**18)
            )

        selected.append(
            sorted(candidates, key=score)[0]
        )

    return selected


def classify(obj: dict) -> str:
    output = str(
        obj.get("output_type", "")
    ).lower()

    assay = str(
        obj.get("assay_title")
        or obj.get("assay_term_name")
        or ""
    ).lower()

    target = ""

    target_obj = obj.get("target")

    if isinstance(target_obj, dict):
        target = str(
            target_obj.get("label")
            or target_obj.get("name")
            or ""
        ).lower()

    text = " ".join([
        output,
        assay,
        target,
        biosample_text(obj).lower()
    ])

    if "candidate cis-regulatory" in text:
        return "cCRE"

    if "dnase" in text:
        return "DNase"

    if "atac" in text:
        return "ATAC"

    if "h3k27ac" in text:
        return "H3K27ac"

    if "h3k4me1" in text:
        return "H3K4me1"

    if "h3k4me3" in text:
        return "H3K4me3"

    if "ctcf" in text:
        return "CTCF"

    if "chip-seq" in text:
        return "TF_CHIP"

    return "OTHER"


def relevant_category(obj: dict) -> bool:
    return classify(obj) in {
        "cCRE",
        "DNase",
        "ATAC",
        "H3K27ac",
        "H3K4me1",
        "H3K4me3",
        "CTCF",
    }


def destination(obj: dict) -> Path:
    accession = (
        obj.get("accession")
        or obj.get("title")
        or "UNKNOWN"
    )

    category = classify(obj)

    href = normalize_href(obj) or ""

    if href.endswith(".bed"):
        suffix = ".bed"
    else:
        suffix = ".bed.gz"

    if category == "cCRE":
        directory = CCRE_DIR
    else:
        directory = PEAK_DIR / category

    directory.mkdir(
        parents=True,
        exist_ok=True
    )

    return directory / f"{accession}{suffix}"


def existing_complete(obj: dict, path: Path) -> bool:
    if not path.exists():
        return False

    expected = int(
        obj.get("file_size") or 0
    )

    if expected <= 0:
        return path.stat().st_size > 0

    return path.stat().st_size == expected


def download(
    url: str,
    path: Path,
    expected_size: int
) -> None:

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    part = Path(str(path) + ".part")

    existing = (
        part.stat().st_size
        if part.exists()
        else 0
    )

    headers = dict(HEADERS)

    if existing > 0:
        headers["Range"] = f"bytes={existing}-"

    print(
        f"[DOWNLOAD] {path.name}\n"
        f"           expected : {human_bytes(expected_size)}\n"
        f"           existing : {human_bytes(existing)}"
    )

    with requests.get(
        url,
        headers=headers,
        stream=True,
        timeout=300
    ) as r:

        # Server ignored Range request.
        if existing and r.status_code == 200:
            existing = 0
            mode = "wb"

        else:
            r.raise_for_status()
            mode = "ab" if existing else "wb"

        with open(part, mode) as fh:

            for chunk in r.iter_content(
                chunk_size=8 * 1024 * 1024
            ):
                if chunk:
                    fh.write(chunk)

    final_size = part.stat().st_size

    if expected_size and final_size != expected_size:
        raise RuntimeError(
            f"Size mismatch for {path.name}: "
            f"expected={expected_size:,}, "
            f"observed={final_size:,}"
        )

    part.replace(path)

    print(
        f"[OK] {path} "
        f"({human_bytes(final_size)})"
    )


def write_manifest(files: list[dict]) -> Path:
    MANIFEST_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    path = MANIFEST_DIR / "encode_grch38_manifest.tsv"

    columns = [
        "accession",
        "category",
        "dataset",
        "assay",
        "target",
        "biosample",
        "output_type",
        "file_format",
        "file_format_type",
        "file_size",
        "md5sum",
        "url",
        "local_path",
    ]

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t"
        )

        writer.writeheader()

        for obj in files:
            target_obj = obj.get("target")

            if isinstance(target_obj, dict):
                target = (
                    target_obj.get("label")
                    or target_obj.get("name")
                    or ""
                )
            else:
                target = ""

            writer.writerow({
                "accession":
                    obj.get("accession")
                    or obj.get("title", ""),

                "category":
                    classify(obj),

                "dataset":
                    obj.get("dataset", ""),

                "assay":
                    obj.get("assay_title")
                    or obj.get("assay_term_name")
                    or "",

                "target":
                    target,

                "biosample":
                    biosample_text(obj),

                "output_type":
                    obj.get("output_type", ""),

                "file_format":
                    obj.get("file_format", ""),

                "file_format_type":
                    obj.get("file_format_type", ""),

                "file_size":
                    obj.get("file_size", 0),

                "md5sum":
                    obj.get("md5sum", ""),

                "url":
                    normalize_href(obj) or "",

                "local_path":
                    destination(obj),
            })

    return path


def discover_encode_files() -> list[dict]:

    banner("DISCOVERING ENCODE REGULATORY FILES")

    all_files = []

    # cCRE
    all_files.extend(
        query_files({
            "output_type":
                "candidate Cis-Regulatory Elements"
        })
    )

    # Accessibility.
    for assay in (
        "DNase-seq",
        "ATAC-seq",
    ):
        all_files.extend(
            query_files({
                "assay_title": assay
            })
        )

    # Histone marks / CTCF.
    for target in (
        "H3K27ac",
        "H3K4me1",
        "H3K4me3",
        "CTCF",
    ):
        all_files.extend(
            query_files({
                "target.label": target
            })
        )

    # Deduplicate accession.
    by_accession = {}

    for obj in all_files:
        accession = (
            obj.get("accession")
            or obj.get("title")
        )

        if accession:
            by_accession[accession] = obj

    files = list(by_accession.values())

    print(
        f"\nUnique ENCODE file objects: "
        f"{len(files):,}"
    )

    files = [
        x for x in files
        if relevant_category(x)
        and is_peak_file(x)
    ]

    print(
        f"Useful regulatory BED objects: "
        f"{len(files):,}"
    )

    files = choose_preferred(files)

    print(
        f"After duplicate/representation filtering: "
        f"{len(files):,}"
    )

    return files


def print_category_summary(files: list[dict]) -> None:

    summary: dict[str, list[dict]] = {}

    for obj in files:
        summary.setdefault(
            classify(obj),
            []
        ).append(obj)

    banner("ENCODE DOWNLOAD SUMMARY")

    grand_total = 0

    for category in sorted(summary):
        subset = summary[category]

        size = sum(
            int(x.get("file_size") or 0)
            for x in subset
        )

        grand_total += size

        print(
            f"{category:15s} "
            f"{len(subset):6,d} files   "
            f"{human_bytes(size):>12s}"
        )

    print("-" * 60)

    print(
        f"{'TOTAL':15s} "
        f"{len(files):6,d} files   "
        f"{human_bytes(grand_total):>12s}"
    )


def main() -> None:

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--plan-only",
        action="store_true",
        help="Query ENCODE and report files/space without downloading."
    )

    parser.add_argument(
        "--core-only",
        action="store_true",
        help="Download only cCRE files."
    )

    parser.add_argument(
        "--reserve-gb",
        type=float,
        default=DEFAULT_RESERVE_GB,
        help="Minimum disk space that must remain after download."
    )

    args = parser.parse_args()

    banner("GWAS2m STEP 07D — ENCODE REGULATORY RESOURCE SETUP")

    print(f"Root       : {ROOT}")
    print(f"Resources  : {RESOURCE_ROOT}")
    print(f"Assembly   : {ASSEMBLY}")
    print("BigWig     : NEVER DOWNLOADED")
    print("FASTQ/BAM  : NEVER DOWNLOADED")
    print("Mode       :", "PLAN ONLY" if args.plan_only else "DOWNLOAD")

    files = discover_encode_files()

    if args.core_only:
        files = [
            x for x in files
            if classify(x) == "cCRE"
        ]

    print_category_summary(files)

    manifest = write_manifest(files)

    print(f"\nManifest: {manifest}")

    required = 0
    already_present = 0
    missing_files = []

    for obj in files:
        path = destination(obj)

        size = int(
            obj.get("file_size") or 0
        )

        if existing_complete(obj, path):
            already_present += size

        else:
            required += size
            missing_files.append(obj)

    total, used, free = disk_status(RESOURCE_ROOT)

    reserve = int(
        args.reserve_gb
        * 1024 ** 3
    )

    banner("DISK SPACE CHECK")

    print(
        f"Filesystem capacity : {human_bytes(total)}"
    )
    print(
        f"Filesystem used     : {human_bytes(used)}"
    )
    print(
        f"Filesystem free     : {human_bytes(free)}"
    )
    print()
    print(
        f"Already downloaded  : {human_bytes(already_present)}"
    )
    print(
        f"Still required      : {human_bytes(required)}"
    )
    print(
        f"Safety reserve      : {human_bytes(reserve)}"
    )
    print(
        f"Missing files       : {len(missing_files):,}"
    )

    if required + reserve > free:

        print()
        print("[STOP] INSUFFICIENT DISK SPACE")
        print(
            f"Need free space of at least: "
            f"{human_bytes(required + reserve)}"
        )
        print(
            f"Currently available: "
            f"{human_bytes(free)}"
        )

        sys.exit(2)

    print()
    print("[OK] Sufficient disk space.")

    if args.plan_only:
        print()
        print(
            "[PLAN ONLY] No files downloaded."
        )
        return

    banner("DOWNLOADING ENCODE FILES")

    for i, obj in enumerate(
        missing_files,
        start=1
    ):

        accession = (
            obj.get("accession")
            or obj.get("title")
        )

        category = classify(obj)

        path = destination(obj)

        url = normalize_href(obj)

        size = int(
            obj.get("file_size") or 0
        )

        print()
        print(
            f"[{i:,}/{len(missing_files):,}] "
            f"{category} {accession}"
        )

        if not url:
            print("[SKIP] No download URL.")
            continue

        try:
            download(
                url,
                path,
                size
            )

        except Exception as exc:
            print(
                f"[ERROR] {accession}: {exc}"
            )

    banner("FINAL STATUS")

    complete = 0
    missing = 0

    for obj in files:

        if existing_complete(
            obj,
            destination(obj)
        ):
            complete += 1
        else:
            missing += 1

    print(f"Complete : {complete:,}")
    print(f"Missing  : {missing:,}")

    if missing:
        print(
            "\nRe-run this script to resume "
            "missing downloads."
        )
        sys.exit(1)

    print()
    print(
        "[SUCCESS] ENCODE regulatory resources ready."
    )


if __name__ == "__main__":
    main()