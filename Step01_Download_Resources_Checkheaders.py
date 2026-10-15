#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2m STEP 01 - CHECK DOWNLOADED RESOURCES (readable? corrupt? print headers)
===============================================================================

Checks every file planned by Step00_Download_Resources.py
(setup_manifests/download_resources*_plan.tsv) and prints its header:

    .parquet          footer read with pyarrow   -> rows + column names
    .vcf.gz           gzip readable              -> the #CHROM line
    .gz (gtf/tsv/fa)  gzip readable              -> first line
    .tbi              tabix magic "TBI\\1"
    .tar / .tar.gz    archive readable           -> first members
    .bw               bigWig magic
    .zip              archive readable           -> member names
    .fa / .fai / text readable                   -> first line

Also reports MISSING files and files whose size differs from the size the
server reported (= partial download). Large archives that were deleted after a
verified extraction are reported OK with their extracted contents.

    envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py            # fast (headers only)
    envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py --full     # + read every compressed file to the end
    envs/pipeline/bin/python Step01_Download_Resources_Checkheaders.py --group opentargets

Run it with envs/pipeline/bin/python (has pyarrow). The fast check is fine on a
login node; --full reads every downloaded byte, so submit that one as a job.

Report: setup_manifests/download_check_headers.tsv
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


ROOT = Path.cwd().resolve()
PLAN = ROOT / "setup_manifests" / "download_resources_plan.tsv"
REPORT = ROOT / "setup_manifests" / "download_check_headers.tsv"

parser = argparse.ArgumentParser(description="Check that downloaded resources are readable and print headers.")
parser.add_argument("--full", action="store_true",
                    help="Also read every compressed file / archive to the end (slow, catches any corruption).")
parser.add_argument("--group", default=None,
                    help="Only check this group, e.g. genome, 1000g, gtex, vep, opentargets, eqtl_catalogue.")
parser.add_argument("--threads", type=int, default=8, help="Files checked in parallel (default 8).")
args = parser.parse_args()

PIGZ = shutil.which("pigz") or str(ROOT / "envs" / "pipeline" / "bin" / "pigz")
PIGZ = PIGZ if Path(PIGZ).exists() else None


# =============================================================================
# ONE CHECK PER FILE TYPE  ->  returns the header text, raises if unreadable
# =============================================================================

def first_line_text(path: Path) -> str:
    opener = gzip.open if path.name.endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.readline().rstrip("\n")


def vcf_header(path: Path) -> str:
    meta = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                meta += 1
                continue
            if line.startswith("#CHROM"):
                samples = len(line.rstrip("\n").split("\t")) - 9
                return f"{meta} ## lines | {samples} samples | " + "\t".join(line.split("\t")[:10])
            raise ValueError("no #CHROM header line")
    raise ValueError("file ended inside the header")


def tbi_header(path: Path) -> str:
    with gzip.open(path, "rb") as handle:
        magic = handle.read(4)
    if magic != b"TBI\x01":
        raise ValueError(f"not a tabix index (magic {magic!r})")
    return "tabix index (TBI\\1)"


def parquet_header(path: Path) -> str:
    try:
        import pyarrow.parquet as pq
    except ImportError:  # no pyarrow: at least check both parquet magic markers
        with open(path, "rb") as handle:
            head = handle.read(4)
            handle.seek(-4, 2)
            tail = handle.read(4)
        if head != b"PAR1" or tail != b"PAR1":
            raise ValueError("parquet magic bytes missing")
        return "parquet (magic OK; install pyarrow for columns)"
    meta = pq.ParquetFile(path).metadata
    names = [meta.schema.column(i).name for i in range(meta.num_columns)]
    return f"rows={meta.num_rows:,} cols={meta.num_columns}: " + ", ".join(names[:12]) + (" ..." if len(names) > 12 else "")


def bigwig_header(path: Path) -> str:
    with open(path, "rb") as handle:
        magic = handle.read(4)
    if magic not in (b"\x26\xfc\x8f\x88", b"\x88\x8f\xfc\x26"):
        raise ValueError(f"not a bigWig file (magic {magic!r})")
    return "bigWig"


def zip_header(path: Path) -> str:
    import zipfile

    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
    return f"zip, {len(names)} files: " + ", ".join(names[:5])


def tar_header(path: Path) -> str:
    members = []
    with tarfile.open(path, "r|*") as archive:          # streaming: only reads the start
        for member in archive:
            members.append(member.name)
            if len(members) == 5:
                break
    if not members:
        raise ValueError("empty archive")
    return "members: " + ", ".join(members)


def full_test(path: Path) -> None:
    """Read the whole file (only with --full)."""
    name = path.name
    if name.endswith((".tar", ".tar.gz", ".tgz")):
        result = subprocess.run(["tar", "-tf", str(path)], stdout=subprocess.DEVNULL,
                                stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise ValueError(f"tar -t failed: {result.stderr.strip()[-200:]}")
    elif name.endswith((".gz", ".bgz", ".tbi")):
        if PIGZ:
            result = subprocess.run([PIGZ, "-t", str(path)], stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                raise ValueError(f"pigz -t failed: {result.stderr.strip()[-200:]}")
        else:
            with gzip.open(path, "rb") as handle:
                while handle.read(1 << 24):
                    pass
    elif name.endswith(".parquet"):
        import pyarrow.parquet as pq
        pq.read_table(path)


def header_of(path: Path) -> str:
    name = path.name
    if name.endswith(".parquet"):
        return parquet_header(path)
    if name.endswith(".tbi"):
        return tbi_header(path)
    if name.endswith((".tar", ".tar.gz", ".tgz")):
        return tar_header(path)
    if name.endswith((".vcf.gz", ".vcf.bgz")):
        return vcf_header(path)
    if name.endswith((".bw", ".bigwig", ".bigWig")):
        return bigwig_header(path)
    if name.endswith(".zip"):
        return zip_header(path)
    return first_line_text(path)


def extracted_summary(directory: Path) -> str:
    files = [p for p in directory.rglob("*") if p.is_file() and not p.name.startswith(".")]
    names = ", ".join(str(p.relative_to(directory)) for p in files[:3])
    return f"archive extracted to {directory.relative_to(ROOT)}: {len(files)} files ({names}{' ...' if len(files) > 3 else ''})"


# =============================================================================
# CHECK ONE FILE
# =============================================================================

def check(row: dict) -> dict:
    path = ROOT / row["DEST"]
    expected = int(row["SIZE_BYTES"]) if row.get("SIZE_BYTES", "").isdigit() else None
    result = {"GROUP": row["GROUP"], "FILE": row["DEST"], "SIZE_BYTES": "", "STATUS": "", "HEADER": ""}

    if not path.is_file():
        # Large archives are deleted after a verified extraction: check the extracted copy.
        post_dir, marker = row.get("POST_DIR", "-"), row.get("POST_MARKER", "-")
        if post_dir not in ("", "-") and (ROOT / post_dir / marker).is_file():
            result["STATUS"] = "OK"
            result["HEADER"] = extracted_summary(ROOT / post_dir)
            return result
        result["STATUS"] = "MISSING"
        return result
    size = path.stat().st_size
    result["SIZE_BYTES"] = size
    if size == 0:
        result["STATUS"] = "EMPTY"
        return result
    if Path(str(path) + ".aria2").exists() or (expected and size != expected):
        result["STATUS"] = "PARTIAL"
        result["HEADER"] = f"size {size:,} != expected {expected:,}" if expected and size != expected \
            else "aria2 control file present (download in progress)"
        return result
    try:
        result["HEADER"] = header_of(path)
        if args.full:
            full_test(path)
        result["STATUS"] = "OK"
    except Exception as exc:
        result["STATUS"] = "CORRUPT"
        result["HEADER"] = f"{type(exc).__name__}: {exc}"
    return result


def planned_files() -> list[dict]:
    # One plan per generated script (e.g. a separate meQTL plan writes its own).
    plans = sorted(PLAN.parent.glob("download_resources*_plan.tsv"))
    if not plans:
        sys.exit(f"ERROR: no plan in {PLAN.parent}.\n       Run envs/pipeline/bin/python Step00_Download_Resources.py first.")
    rows: dict[str, dict] = {}
    for plan in plans:
        with open(plan, encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows[row["DEST"]] = row
    rows = list(rows.values())
    extra = []
    for row in rows:
        # The FASTA is decompressed + indexed by the download job: check those too.
        if row["DEST"].endswith(".fa.gz"):
            fasta = row["DEST"][:-3]
            extra += [dict(row, DEST=fasta, SIZE_BYTES=""), dict(row, DEST=fasta + ".fai", SIZE_BYTES="")]
    rows += extra
    if args.group:
        rows = [r for r in rows if r["GROUP"] == args.group or r["GROUP"].startswith(args.group + "/")]
    return rows


# =============================================================================
# MAIN
# =============================================================================

def main() -> int:
    rows = planned_files()
    print("=" * 100)
    print(f"GWAS2m - CHECK DOWNLOADED RESOURCES  ({len(rows)} files, {'FULL read' if args.full else 'headers only'})")
    print("=" * 100)

    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
        results = list(pool.map(check, rows))

    for r in results:
        size = f"{r['SIZE_BYTES'] / 1e6:,.1f} MB" if r["SIZE_BYTES"] else "-"
        print(f"[{r['STATUS']:<7}] {size:>12}  {r['FILE']}")
        if r["HEADER"]:
            print(f"            {r['HEADER'][:160]}")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with open(REPORT, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["STATUS", "GROUP", "SIZE_BYTES", "FILE", "HEADER"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)

    counts: dict[str, int] = {}
    for r in results:
        counts[r["STATUS"]] = counts.get(r["STATUS"], 0) + 1
    print("=" * 100)
    print("SUMMARY : " + "  ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"Report  : {REPORT}")
    bad = len(results) - counts.get("OK", 0)
    if bad:
        print(f"{bad} file(s) not OK. Re-plan and resubmit the missing ones:")
        print("  envs/pipeline/bin/python Step00_Download_Resources.py && sbatch Step00_Download_Resources.sh")
    else:
        print("All files readable.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
