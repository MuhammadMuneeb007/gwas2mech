#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 08 - Local annotation of fine-mapped variants
===============================================================================

    # plan + SLURM array (one task per fine-mapped study; >1000 studies share tasks)
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR
    sbatch Step08_Annotate_parkinson_s_disease_european.sh

    # study 2 (the number the Step02 planner prints) on this machine
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --index 2

    # which local resources are present, progress, per-study statistics
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --status
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --stats

    # one variant (GRCh38, CHR:POS:REF:ALT)
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --test-variant 4:89822321:G:A

    # verify everything annotation needs (files + indexes, VEP cache, tools, Python modules,
    # Step06 input for this panel, API keys); exits 1 if anything is missing
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --check-resources
    # ... and also test the web APIs (login node with internet)
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --check-resources --online screen,regulomedb

    # ---- 1000 Genomes 30x LD panel (Steps 04-06 with --reference-panel 1000G30x) ----
    # reads 06_finemapping/<p>/<a>_1000G30x/, writes 08_annotation/<p>/<a>_1000G30x/
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x
    sbatch Step08_Annotate_parkinson_s_disease_european_1000G30x.sh
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --index 2
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --status
    envs/pipeline/bin/python Step08_Annotate_Finemapped_Variants.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

LOCAL FIRST. Every annotation comes from versioned files downloaded by
Step00_Download_Resources.py (--groups vep genome annotation cadd); no web
service is needed, so the step runs on compute nodes without internet, has no
rate limits, annotates ALL fine-mapped variants (not only credible sets) and
is reproducible:

    Ensembl VEP 116, offline cache  consequences on every transcript, the most
        + GRCh38 FASTA              severe / MANE / canonical one, gene, biotype,
                                    exon/intron, HGVS c./p., amino acids, protein
                                    domains, UniProt, SIFT, PolyPhen, regulatory
                                    features, rsIDs, gnomAD v4 / 1000G frequencies
    CADD v1.7 (SNVs + indels)       raw score, PHRED
    AlphaMissense (hg38)            missense pathogenicity + class
    ClinVar VCF (dated release)     clinical significance, review status, condition
    ENCODE SCREEN cCRE V4           overlapping candidate cis-regulatory element + class
    phyloP 100-way (UCSC)           base-wise conservation
    gnomAD v4.1 constraint          LOEUF, pLI, missense Z of the affected gene

OPTIONAL ONLINE (--online screen,regulomedb), credible-set variants only:
    SCREEN GraphQL  cell-type-specific cCRE activity, cCRE -> gene links
                    (needs SCREEN_API_KEY, see "API keys" below)
    RegulomeDB      regulatory probability / rank (no key)
Network failures are recorded as API_ERROR, never as "no hit".

API keys (only for --online screen; the default local run needs no key and no internet):
    SCREEN_API_KEY is read, in this order, from
      1. the environment:           export SCREEN_API_KEY=...
      2. config/api_keys.env:       SCREEN_API_KEY=...   (git-ignored, one KEY=value per line;
                                    template: config/api_keys.example.env)
    Never paste a key into this script: it is tracked by git.

Inputs (Step06):
    06_finemapping/<p>/<a>/<GCST>/<GCST>_finemapped_variants.tsv.gz
    (LOCUS_ID CHR REFERENCE_POS REFERENCE_REF REFERENCE_ALT PIP IN_95_CREDIBLE_SET ...)

Outputs (08_annotation/<p>/<a>/):
    <GCST>/<GCST>_annotated_variants.tsv.gz       every variant: fine-mapping + annotation
    <GCST>/<GCST>_credible_set_variants.tsv        95% credible-set variants, by locus / PIP
    <GCST>/<GCST>_locus_annotation_summary.tsv     per locus: top variant, gene, coding /
                                                   LoF / missense / ClinVar / cCRE evidence,
                                                   PIP-weighted functional composition
    <GCST>/annotation_provenance.json              resources, versions, counts per source
    annotation_stats.tsv                           (--stats)

Alleles: annotation is for the 1000G REF/ALT pair matched in Step06. Z_ALT /
BETA_ALT give the GWAS effect of the ALT allele (Step06 orients Z to REF).

SLURM settings come from the #SBATCH header of HPC_CPU.sh; --partition --qos
--account --cpus --mem --time --max-parallel override them.
===============================================================================
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

try:
    import requests
except ImportError:
    requests = None


VERSION = "6.0.0-local"
ASSEMBLY = "GRCh38"
SPECIES = "homo_sapiens"
PROJECT_ROOT = Path.cwd().resolve()
MAX_ARRAY_TASKS = 1000

# ---- web APIs (only with --online; the default local annotation calls none) --
SCREEN_GRAPHQL = "https://screen.api.wenglab.org/graphql"       # needs SCREEN_API_KEY
REGULOMEDB_API = "https://api.regulomedb.org/search"            # public, no key
REGULOMEDB_FALLBACK = "https://regulomedb.org/regulome-search/"  # public, no key
# Keys: environment first, then this git-ignored file (KEY=value lines). Never put a key in this script.
API_KEYS_FILE = PROJECT_ROOT / "config" / "api_keys.env"
API_KEY_NAMES = {"screen": "SCREEN_API_KEY"}                    # --online provider -> key variable


def read_api_keys_file(path: Path = API_KEYS_FILE) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.removeprefix("export ").split("=", 1)
        values[key.strip()] = value.strip().strip("'\"")
    return values


def api_key(name: str) -> tuple[str, str]:
    """(value, source) of an API key: environment, then config/api_keys.env; ('', '') if unset."""
    value = os.environ.get(name, "").strip()
    if value:
        return value, "environment"
    value = read_api_keys_file().get(name, "").strip()
    return (value, API_KEYS_FILE.relative_to(PROJECT_ROOT).as_posix()) if value else ("", "")


SCREEN_API_KEY, SCREEN_API_KEY_SOURCE = api_key("SCREEN_API_KEY")

# ---- local resources (Step00_Download_Resources.py) -------------------------
RES = PROJECT_ROOT / "resources"
VEP_RELEASE = "116"
VEP_CACHE_DIR = RES / "vep" / "cache"
FASTA = RES / "genome" / "GRCh38" / "GRCh38.primary_assembly.genome.fa"
CADD_SNV = RES / "annotation" / "cadd" / "v1.7" / "whole_genome_SNVs.tsv.gz"
CADD_INDEL = RES / "annotation" / "cadd" / "v1.7" / "gnomad.genomes.r4.0.indel.tsv.gz"
ALPHAMISSENSE = RES / "annotation" / "alphamissense" / "AlphaMissense_hg38.tsv.gz"
CLINVAR = RES / "annotation" / "clinvar" / "clinvar.vcf.gz"
CCRE_BED = RES / "annotation" / "encode_screen" / "GRCh38-cCREs.V4.bed"
CONSTRAINT = RES / "annotation" / "gnomad" / "gnomad.v4.1.constraint_metrics.tsv"
PHYLOP = RES / "annotation" / "conservation" / "hg38.phyloP100way.bw"

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

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

CONSEQUENCE_RANK = {
    "transcript_ablation": 0,
    "splice_acceptor_variant": 1,
    "splice_donor_variant": 2,
    "stop_gained": 3,
    "frameshift_variant": 4,
    "stop_lost": 5,
    "start_lost": 6,
    "transcript_amplification": 7,
    "inframe_insertion": 8,
    "inframe_deletion": 9,
    "missense_variant": 10,
    "protein_altering_variant": 11,
    "splice_region_variant": 12,
    "incomplete_terminal_codon_variant": 13,
    "start_retained_variant": 14,
    "stop_retained_variant": 15,
    "synonymous_variant": 16,
    "coding_sequence_variant": 17,
    "mature_miRNA_variant": 18,
    "5_prime_UTR_variant": 19,
    "3_prime_UTR_variant": 20,
    "non_coding_transcript_exon_variant": 21,
    "intron_variant": 22,
    "NMD_transcript_variant": 23,
    "non_coding_transcript_variant": 24,
    "upstream_gene_variant": 25,
    "downstream_gene_variant": 26,
    "TFBS_ablation": 27,
    "TFBS_amplification": 28,
    "TF_binding_site_variant": 29,
    "regulatory_region_ablation": 30,
    "regulatory_region_amplification": 31,
    "feature_elongation": 32,
    "regulatory_region_variant": 33,
    "feature_truncation": 34,
    "intergenic_variant": 35,
}


# =============================================================================
# GENERAL UTILITIES
# =============================================================================

def banner(text: str) -> None:
    print()
    print("=" * 112)
    print(text)
    print("=" * 112)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_text(value: Any) -> str:
    x = str(value).lower()
    x = re.sub(r"[^a-z0-9]+", " ", x)
    x = re.sub(r"\s+", " ", x)
    return x.strip()


def slugify(value: Any) -> str:
    return normalize_text(value).replace(" ", "_")


def canonical_ancestry(value: str) -> tuple[str, str]:
    key = normalize_text(value)
    if key not in ANCESTRY_ALIASES:
        raise ValueError(
            f"Unsupported ancestry {value!r}. Use EUR, AFR, EAS, SAS or AMR."
        )
    return ANCESTRY_ALIASES[key]


def safe_float(value: Any) -> float:
    try:
        x = float(value)
        return x if np.isfinite(x) else np.nan
    except Exception:
        return np.nan


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def compact_json(value: Any) -> str:
    if value in (None, "", [], {}):
        return ""
    try:
        return json.dumps(value, separators=(",", ":"), sort_keys=True)
    except Exception:
        return str(value)


def join_unique(values: Iterable[Any], sep: str = ";", limit: int | None = None) -> str:
    seen = set()
    out = []
    for value in values:
        if value is None:
            continue
        if isinstance(value, float) and np.isnan(value):
            continue
        if isinstance(value, (list, tuple, set)):
            items = value
        else:
            items = [value]
        for item in items:
            if item is None:
                continue
            s = str(item).strip()
            if not s or s.lower() in {"nan", "none", "null"}:
                continue
            if s not in seen:
                seen.add(s)
                out.append(s)
                if limit is not None and len(out) >= limit:
                    return sep.join(out)
    return sep.join(out)


def primary_rsid(*values: Any) -> str:
    for value in values:
        m = re.search(r"\brs\d+\b", str(value), flags=re.I)
        if m:
            return m.group(0).lower()
    return ""


def variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    c = str(chrom).replace("chr", "").replace("CHR", "").strip()
    return f"{c}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def functional_category(consequences: Iterable[str], biotype: str = "") -> str:
    terms = {str(x).lower() for x in consequences if str(x).strip()}
    bt = str(biotype).lower()

    if any(
        t in terms
        for t in {
            "splice_acceptor_variant",
            "splice_donor_variant",
            "splice_region_variant",
            "splice_donor_5th_base_variant",
            "splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant",
        }
    ):
        return "SPLICING"

    if any(
        t in terms
        for t in {
            "transcript_ablation",
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "transcript_amplification",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant",
            "coding_sequence_variant",
            "synonymous_variant",
            "start_retained_variant",
            "stop_retained_variant",
            "incomplete_terminal_codon_variant",
        }
    ):
        return "CODING"

    if "5_prime_utr_variant" in terms or "3_prime_utr_variant" in terms:
        return "UTR"

    if any("regulatory" in t or "tf_binding" in t or "tfbs" in t for t in terms):
        return "REGULATORY"

    if (
        any("non_coding_transcript" in t for t in terms)
        or any(
            x in bt
            for x in (
                "lncrna",
                "lincrna",
                "mirna",
                "snrna",
                "snorna",
                "rrna",
                "ncrna",
                "antisense",
                "processed_transcript",
            )
        )
    ):
        return "NONCODING_RNA"

    if "intron_variant" in terms:
        return "INTRONIC"
    if "upstream_gene_variant" in terms:
        return "UPSTREAM"
    if "downstream_gene_variant" in terms:
        return "DOWNSTREAM"
    if "intergenic_variant" in terms:
        return "INTERGENIC"
    return "OTHER"


def request_session() -> requests.Session:
    if requests is None:
        raise RuntimeError(
            "Python requests is required. Install it in your pipeline environment."
        )
    s = requests.Session()
    s.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": f"GWAS2Mech-Step08/{VERSION}",
        }
    )
    return s


def request_with_retry(
    session: requests.Session,
    method: str,
    url: str,
    *,
    tries: int = 6,
    timeout: int = 180,
    **kwargs: Any,
) -> requests.Response:
    last: Exception | None = None

    for attempt in range(tries):
        try:
            r = session.request(method, url, timeout=timeout, **kwargs)

            if r.status_code == 429:
                retry = r.headers.get("Retry-After")
                try:
                    wait = float(retry) if retry else min(120.0, 2.0 ** attempt)
                except Exception:
                    wait = min(120.0, 2.0 ** attempt)
                print(f"[429] {url} — sleeping {wait:.1f}s", flush=True)
                time.sleep(wait)
                continue

            if r.status_code in {500, 502, 503, 504}:
                wait = min(60.0, 2.0 ** attempt)
                print(
                    f"[HTTP {r.status_code}] transient server error — sleeping {wait:.1f}s",
                    flush=True,
                )
                time.sleep(wait)
                continue

            r.raise_for_status()
            return r

        except Exception as exc:
            last = exc
            if attempt == tries - 1:
                break
            wait = min(60.0, 2.0 ** attempt)
            print(
                f"[WARN] request failed attempt {attempt + 1}/{tries}: "
                f"{type(exc).__name__}: {exc}; sleeping {wait:.1f}s",
                flush=True,
            )
            time.sleep(wait)

    raise RuntimeError(f"Request failed after {tries} attempts: {url}\n{last}")


def atomic_write_gzip_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.{uuid.uuid4().hex}.tmp")
    with gzip.open(tmp, "wt", encoding="utf-8") as fh:
        json.dump(obj, fh, separators=(",", ":"))
    os.replace(tmp, path)


def read_gzip_json(path: Path) -> Any:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    df.to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, path)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def shell(cmd: list[str]) -> str:
    cmd = [str(x) for x in cmd]
    print("$ " + " ".join(shlex.quote(x) for x in cmd), flush=True)
    p = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(p.stdout, end="")
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {p.returncode}: {' '.join(cmd)}"
        )
    return p.stdout


# =============================================================================
# INPUT / PLANNING
# =============================================================================

def load_finemapped(path: Path) -> pd.DataFrame:
    x = read_tsv(path)

    required = [
        "LOCUS_ID",
        "CHR",
        "REFERENCE_POS",
        "REFERENCE_REF",
        "REFERENCE_ALT",
        "PIP",
    ]
    missing = [c for c in required if c not in x.columns]
    if missing:
        raise RuntimeError(f"{path} missing required columns: {missing}")

    x["CHR"] = (
        x["CHR"]
        .astype(str)
        .str.replace(r"^chr", "", regex=True, case=False)
        .str.strip()
    )
    x["CHR"] = pd.to_numeric(x["CHR"], errors="coerce")
    x["REFERENCE_POS"] = pd.to_numeric(x["REFERENCE_POS"], errors="coerce")
    x["PIP"] = pd.to_numeric(x["PIP"], errors="coerce")

    x = x.dropna(
        subset=["CHR", "REFERENCE_POS", "REFERENCE_REF", "REFERENCE_ALT", "PIP"]
    ).copy()
    x = x[x["CHR"].between(1, 22, inclusive="both")].copy()

    x["CHR"] = x["CHR"].astype(int)
    x["REFERENCE_POS"] = x["REFERENCE_POS"].astype(int)
    x["REFERENCE_REF"] = x["REFERENCE_REF"].astype(str).str.upper().str.strip()
    x["REFERENCE_ALT"] = x["REFERENCE_ALT"].astype(str).str.upper().str.strip()
    x = x[(x["REFERENCE_REF"] != "") & (x["REFERENCE_ALT"] != "")].copy()

    x["VARIANT_KEY"] = [
        variant_key(c, p, r, a)
        for c, p, r, a in zip(
            x["CHR"],
            x["REFERENCE_POS"],
            x["REFERENCE_REF"],
            x["REFERENCE_ALT"],
        )
    ]

    def get_rsid(row: pd.Series) -> str:
        return primary_rsid(
            row.get("rsID", ""),
            row.get("SNPID", ""),
            row.get("REFERENCE_ID", ""),
        )

    x["PRIMARY_RSID"] = x.apply(get_rsid, axis=1)

    if "IN_95_CREDIBLE_SET" in x.columns:
        s = x["IN_95_CREDIBLE_SET"].fillna("").astype(str).str.lower().str.strip()
        x["IN_95_CREDIBLE_SET"] = s.isin({"1", "true", "yes", "y", "t"})
    elif "CREDIBLE_SET" in x.columns:
        s = x["CREDIBLE_SET"].fillna("").astype(str).str.lower().str.strip()
        x["IN_95_CREDIBLE_SET"] = ~s.isin({"", ".", "na", "nan", "none", "false", "0"})
    else:
        x["IN_95_CREDIBLE_SET"] = False

    x["PIP_GE_0_01"] = x["PIP"] >= 0.01
    x["PIP_GE_0_10"] = x["PIP"] >= 0.10
    x["PIP_GE_0_50"] = x["PIP"] >= 0.50
    x["PIP_GE_0_90"] = x["PIP"] >= 0.90

    return x


# =============================================================================
# VEP JSON -> one row per variant (most severe / MANE / canonical transcript)
# =============================================================================

def transcript_rank(tc: dict[str, Any]) -> tuple:
    mane = 0 if tc.get("mane_select") else 1
    canonical = 0 if tc.get("canonical") in (1, True, "1") else 1
    impact = IMPACT_RANK.get(str(tc.get("impact", "")).upper(), 9)

    consequences = tc.get("consequence_terms", []) or []
    c_rank = min(
        [CONSEQUENCE_RANK.get(str(c), 99) for c in consequences] or [99]
    )

    appris = str(tc.get("appris", "") or "")
    appris_rank = 0 if "principal" in appris.lower() else 1

    return mane, canonical, appris_rank, impact, c_rank


def collect_key_values(obj: Any, wanted_tokens: Iterable[str]) -> dict[str, list[str]]:
    wanted = tuple(str(x).lower() for x in wanted_tokens)
    found: dict[str, list[str]] = {}

    def walk(x: Any) -> None:
        if isinstance(x, dict):
            for k, v in x.items():
                lk = str(k).lower()
                if any(token in lk for token in wanted):
                    if not isinstance(v, (dict, list)):
                        found.setdefault(str(k), [])
                        sv = str(v)
                        if sv not in found[str(k)]:
                            found[str(k)].append(sv)
                    elif isinstance(v, list) and all(
                        not isinstance(y, (dict, list)) for y in v
                    ):
                        found.setdefault(str(k), [])
                        for y in v:
                            sy = str(y)
                            if sy not in found[str(k)]:
                                found[str(k)].append(sy)
                walk(v)
        elif isinstance(x, list):
            for y in x:
                walk(y)

    walk(obj)
    return found


def extract_ensembl_vep_record(item: dict[str, Any], api_id: str, vkey: str) -> dict[str, Any]:
    out: dict[str, Any] = {
        "VARIANT_KEY": vkey,
        "API_ID": api_id,
        "ENSEMBL_VEP_STATUS": "OK",
        "ENSEMBL_VEP_MOST_SEVERE": item.get("most_severe_consequence", ""),
        "ENSEMBL_VEP_ID": item.get("id", ""),
        "ENSEMBL_VEP_ASSEMBLY": item.get("assembly_name", ""),
        "ENSEMBL_VEP_SEQ_REGION": item.get("seq_region_name", ""),
        "ENSEMBL_VEP_START": item.get("start", ""),
        "ENSEMBL_VEP_END": item.get("end", ""),
        "ENSEMBL_VEP_ALLELE_STRING": item.get("allele_string", ""),
    }

    transcripts = [
        x
        for x in (item.get("transcript_consequences", []) or [])
        if isinstance(x, dict)
    ]

    all_consequences = []
    all_genes = []
    all_symbols = []
    all_transcripts = []
    all_biotypes = []

    for tc in transcripts:
        all_consequences.extend(tc.get("consequence_terms", []) or [])
        all_genes.append(tc.get("gene_id", ""))
        all_symbols.append(tc.get("gene_symbol", ""))
        all_transcripts.append(tc.get("transcript_id", ""))
        all_biotypes.append(tc.get("biotype", ""))

    out["ENSEMBL_ALL_CONSEQUENCES"] = join_unique(all_consequences)
    out["ENSEMBL_ALL_GENES"] = join_unique(all_genes)
    out["ENSEMBL_ALL_SYMBOLS"] = join_unique(all_symbols)
    out["ENSEMBL_ALL_TRANSCRIPTS"] = join_unique(all_transcripts)
    out["ENSEMBL_ALL_BIOTYPES"] = join_unique(all_biotypes)
    out["ENSEMBL_N_TRANSCRIPT_CONSEQUENCES"] = len(transcripts)

    best: dict[str, Any] = {}
    if transcripts:
        best = sorted(transcripts, key=transcript_rank)[0]

    best_cons = best.get("consequence_terms", []) or []
    out["ENSEMBL_FUNCTIONAL_CATEGORY"] = functional_category(
        best_cons, best.get("biotype", "")
    )
    out["ENSEMBL_BEST_CONSEQUENCE"] = join_unique(best_cons)
    out["ENSEMBL_BEST_IMPACT"] = best.get("impact", "")
    out["ENSEMBL_GENE_ID"] = best.get("gene_id", "")
    out["ENSEMBL_GENE_SYMBOL"] = best.get("gene_symbol", "")
    out["ENSEMBL_TRANSCRIPT_ID"] = best.get("transcript_id", "")
    out["ENSEMBL_TRANSCRIPT_BIOTYPE"] = best.get("biotype", "")
    out["ENSEMBL_EXON"] = best.get("exon", "")
    out["ENSEMBL_INTRON"] = best.get("intron", "")
    out["ENSEMBL_HGVSC"] = best.get("hgvsc", "")
    out["ENSEMBL_HGVSP"] = best.get("hgvsp", "")
    out["ENSEMBL_PROTEIN_ID"] = best.get("protein_id", "")
    out["ENSEMBL_PROTEIN_START"] = best.get("protein_start", "")
    out["ENSEMBL_PROTEIN_END"] = best.get("protein_end", "")
    out["ENSEMBL_AMINO_ACIDS"] = best.get("amino_acids", "")
    out["ENSEMBL_CODONS"] = best.get("codons", "")
    out["ENSEMBL_CANONICAL"] = best.get("canonical", "")
    out["ENSEMBL_MANE_SELECT"] = best.get("mane_select", "")
    out["ENSEMBL_MANE_PLUS_CLINICAL"] = best.get("mane_plus_clinical", "")
    out["ENSEMBL_APPRIS"] = best.get("appris", "")
    out["ENSEMBL_TSL"] = best.get("tsl", "")
    out["ENSEMBL_SIFT_PREDICTION"] = best.get("sift_prediction", "")
    out["ENSEMBL_SIFT_SCORE"] = best.get("sift_score", "")
    out["ENSEMBL_POLYPHEN_PREDICTION"] = best.get("polyphen_prediction", "")
    out["ENSEMBL_POLYPHEN_SCORE"] = best.get("polyphen_score", "")

    # High-value VEP plugin/prediction fields.
    wanted = (
        "alphamissense",
        "am_pathogenicity",
        "am_class",
        "cadd",
        "revel",
        "protvar",
        "spliceai",
        "blosum",
        "clinpred",
        "eve",
        "loeuf",
        "mavedb",
        "opentarget",
        "swissprot",
        "uniprot",
        "uniparc",
        "domain",
    )
    kv = collect_key_values(transcripts, wanted)
    for key, values in sorted(kv.items()):
        col = "ENSEMBL_PLUGIN_" + re.sub(r"[^A-Za-z0-9]+", "_", key).upper()
        out[col] = join_unique(values)

    # Co-located variants give rsIDs, ClinVar significance and frequency evidence.
    colocated = [
        x for x in (item.get("colocated_variants", []) or []) if isinstance(x, dict)
    ]
    out["ENSEMBL_COLOCATED_IDS"] = join_unique(
        [x.get("id", "") for x in colocated]
    )
    out["ENSEMBL_CLINICAL_SIGNIFICANCE"] = join_unique(
        [
            sig
            for x in colocated
            for sig in (x.get("clin_sig", []) or [])
        ]
    )
    out["ENSEMBL_FREQUENCIES_JSON"] = join_unique(
        [
            compact_json(x.get("frequencies", {}))
            for x in colocated
            if x.get("frequencies")
        ]
    )
    out["ENSEMBL_PHENOTYPE_VARIANTS"] = join_unique(
        [
            x.get("id", "")
            for x in colocated
            if x.get("phenotype_or_disease") in (1, True, "1")
        ]
    )

    regs = [
        x
        for x in (item.get("regulatory_feature_consequences", []) or [])
        if isinstance(x, dict)
    ]
    out["ENSEMBL_REGULATORY_FEATURE_IDS"] = join_unique(
        [x.get("regulatory_feature_id", "") for x in regs]
    )
    out["ENSEMBL_REGULATORY_CONSEQUENCES"] = join_unique(
        [c for x in regs for c in (x.get("consequence_terms", []) or [])]
    )
    out["ENSEMBL_REGULATORY_BIOTYPES"] = join_unique(
        [x.get("biotype", "") for x in regs]
    )

    return out


# =============================================================================
# OPTIONAL ONLINE PROVIDERS: SCREEN (GraphQL) and RegulomeDB
# =============================================================================

def screen_post(
    session: requests.Session,
    query: str,
    api_key: str,
) -> dict[str, Any]:
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    r = request_with_retry(
        session,
        "POST",
        SCREEN_GRAPHQL,
        headers=headers,
        json={"query": query},
    )
    payload = r.json()
    if payload.get("errors"):
        raise RuntimeError(
            "SCREEN GraphQL error: "
            + compact_json(payload["errors"])[:5000]
        )
    return payload.get("data", {}) or {}


def screen_coords_literal(batch: pd.DataFrame) -> str:
    coords = []
    for r in batch.itertuples(index=False):
        pos = int(r.REFERENCE_POS)
        coords.append(
            f'{{chromosome:"chr{int(r.CHR)}",start:{max(0, pos - 1)},end:{pos}}}'
        )
    return ",".join(coords)


def screen_query_overlaps(
    session: requests.Session,
    batch: pd.DataFrame,
    api_key: str,
) -> list[dict[str, Any]]:
    coords = screen_coords_literal(batch)
    query = f"""
    query {{
      getmaxZScoresQuery(
        assembly: "grch38",
        coordinates: [{coords}]
      ) {{
        accession
        chromosome
        start
        stop
        ccre_group
        atac_max_zscore
        ctcf_max_zscore
        dnase_max_zscore
        h3k27ac_max_zscore
        h3k4me3_max_zscore
        nearestgenes {{
          distance
          gene
        }}
      }}
    }}
    """
    data = screen_post(session, query, api_key)
    x = data.get("getmaxZScoresQuery", []) or []
    return x if isinstance(x, list) else []


def screen_query_ccre_detail(
    session: requests.Session,
    accession: str,
    api_key: str,
) -> dict[str, Any]:
    # A single GraphQL request retrieves multiple rich evidence layers for one cCRE.
    query = f"""
    query {{
      signals: getcCREZScoresQuery(
        assembly: "grch38",
        accession: ["{accession}"],
        include_biosample_details: true
      ) {{
        chromosome
        start
        stop
        accession
        atac_max_zscore
        ctcf_max_zscore
        dnase_max_zscore
        h3k27ac_max_zscore
        h3k4me3_max_zscore
        zscores
      }}

      linked: linkedGenesQuery(
        assembly: "grch38",
        accession: ["{accession}"]
      ) {{
        accession
        p_val
        gene
        geneid
        genetype
        method
        effectsize
        assay
        celltype
        experiment_accession
        tissue
        variantid
        source
        slope
        score
        displayname
      }}

      computational: ComputationalGeneLinksQuery(
        accession: ["{accession}"]
      ) {{
        genename
        accession
        geneid
        genetype
        method
        celltype
        score
        methodregion
        fileaccession
      }}

      tfpeaks: getcCRERPeaksQuery(
        accession: ["{accession}"]
      ) {{
        rpeak_id
        rpeak_chromosome
        rpeak_start
        rpeak_stop
        tf
        metadata
        strand
      }}

      silencers: silencersQuery(
        accession: ["{accession}"]
      ) {{
        silencer_studies
      }}

      dynamic: dynamicEnhancersQuery(
        accession: ["{accession}"]
      ) {{
        celltype
        accession
      }}
    }}
    """
    return screen_post(session, query, api_key)


def summarize_screen_detail(detail: dict[str, Any], active_z: float) -> dict[str, Any]:
    linked = detail.get("linked", []) or []
    comp = detail.get("computational", []) or []
    tfpeaks = detail.get("tfpeaks", []) or []
    silencers = detail.get("silencers", []) or []
    dynamic = detail.get("dynamic", []) or []
    signals = detail.get("signals", []) or []

    experimental_genes = []
    linked_displaynames = []
    computational_genes = []
    celltypes = []
    tissues = []
    methods = []
    assays = []
    experiment_accessions = []
    file_accessions = []
    tfs = []
    active_biosamples = []
    active_tissues = []
    active_assays = []
    active_lifestages = []

    for x in linked:
        if not isinstance(x, dict):
            continue
        experimental_genes.append(x.get("gene", ""))
        linked_displaynames.append(x.get("displayname", ""))
        celltypes.append(x.get("celltype", ""))
        tissues.append(x.get("tissue", ""))
        methods.append(x.get("method", ""))
        assays.append(x.get("assay", ""))
        experiment_accessions.append(x.get("experiment_accession", ""))

    for x in comp:
        if not isinstance(x, dict):
            continue
        computational_genes.append(x.get("genename", ""))
        celltypes.append(x.get("celltype", ""))
        methods.append(x.get("method", ""))
        file_accessions.append(x.get("fileaccession", ""))

    for x in tfpeaks:
        if isinstance(x, dict):
            tfs.append(x.get("tf", ""))

    signal_rows = []
    for s in signals:
        if not isinstance(s, dict):
            continue
        z = s.get("zscores", []) or []
        if isinstance(z, str):
            try:
                z = json.loads(z)
            except Exception:
                z = []

        if isinstance(z, list):
            for row in z:
                if not isinstance(row, (list, tuple)) or len(row) < 9:
                    continue
                signal_rows.append(row)

                score = safe_float(row[8])
                if not np.isfinite(score) or score < active_z:
                    continue

                # SCREEN documentation:
                # [experiment_accession, file_accession, assay, biosample_value,
                #  biosample_name, organ/tissue, sampletype, lifestage, zscore]
                experiment_accessions.append(row[0])
                file_accessions.append(row[1])
                active_assays.append(row[2])
                active_biosamples.append(row[4])
                active_tissues.append(row[5])
                active_lifestages.append(row[7])

    return {
        "SCREEN_LINKED_GENES_EXPERIMENTAL": join_unique(experimental_genes),
        "SCREEN_LINKED_DISPLAY_NAMES": join_unique(linked_displaynames, limit=100),
        "SCREEN_LINKED_GENES_COMPUTATIONAL": join_unique(computational_genes),
        "SCREEN_CELL_TYPES": join_unique(celltypes, limit=100),
        "SCREEN_TISSUES": join_unique(tissues, limit=100),
        "SCREEN_LINK_METHODS": join_unique(methods),
        "SCREEN_LINK_ASSAYS": join_unique(assays),
        "SCREEN_EXPERIMENT_ACCESSIONS": join_unique(experiment_accessions, limit=200),
        "SCREEN_FILE_ACCESSIONS": join_unique(file_accessions, limit=200),
        "SCREEN_TFS": join_unique(tfs, limit=100),
        "SCREEN_ACTIVE_BIOSAMPLES": join_unique(active_biosamples, limit=100),
        "SCREEN_ACTIVE_TISSUES": join_unique(active_tissues, limit=100),
        "SCREEN_ACTIVE_ASSAYS": join_unique(active_assays),
        "SCREEN_ACTIVE_LIFESTAGES": join_unique(active_lifestages),
        "SCREEN_N_SIGNAL_ROWS": len(signal_rows),
        "SCREEN_SILENCER_JSON": compact_json(silencers),
        "SCREEN_DYNAMIC_ENHANCERS_JSON": compact_json(dynamic),
    }


def screen_cache_file(cache_root: Path, accession: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", accession)
    return cache_root / f"{safe}.json.gz"


def get_screen_detail_cached(
    accession: str,
    api_key: str,
    cache_root: Path,
    active_z: float,
    force: bool,
) -> tuple[str, dict[str, Any], str]:
    path = screen_cache_file(cache_root, accession)

    if path.exists() and not force:
        try:
            detail = read_gzip_json(path)
            return accession, summarize_screen_detail(detail, active_z), "CACHE"
        except Exception:
            pass

    session = request_session()
    try:
        detail = screen_query_ccre_detail(session, accession, api_key)
        atomic_write_gzip_json(path, detail)
        return accession, summarize_screen_detail(detail, active_z), "API"
    except Exception as exc:
        return accession, {}, f"API_ERROR:{type(exc).__name__}:{exc}"


def worker_screen(batch: pd.DataFrame, args: argparse.Namespace, output_file: Path) -> pd.DataFrame:
    api_key = SCREEN_API_KEY
    if not api_key:
        return pd.DataFrame(
            [
                {
                    "VARIANT_KEY": r.VARIANT_KEY,
                    "SCREEN_STATUS": "SCREEN_API_KEY_NOT_SET",
                }
                for r in batch.itertuples(index=False)
            ]
        )

    session = request_session()

    try:
        ccres = screen_query_overlaps(session, batch, api_key)
    except Exception as exc:
        err = f"{type(exc).__name__}: {exc}"
        return pd.DataFrame(
            [
                {
                    "VARIANT_KEY": r.VARIANT_KEY,
                    "SCREEN_STATUS": "API_ERROR",
                    "SCREEN_ERROR": err,
                }
                for r in batch.itertuples(index=False)
            ]
        )

    # Determine which returned cCREs overlap each one-base variant interval.
    var_hits: dict[str, list[dict[str, Any]]] = {
        str(v): [] for v in batch["VARIANT_KEY"]
    }

    for c in ccres:
        if not isinstance(c, dict):
            continue

        chrom = str(c.get("chromosome", "")).replace("chr", "")
        start = safe_int(c.get("start"), -1)
        stop = safe_int(c.get("stop"), -1)

        if start < 0 or stop < 0:
            continue

        candidates = batch[
            (batch["CHR"].astype(int).astype(str) == chrom)
        ]

        for r in candidates.itertuples(index=False):
            # GWAS position is 1-based. SCREEN intervals are treated as BED-like
            # [start, stop), so the one-base variant interval is [pos-1, pos).
            v0 = int(r.REFERENCE_POS) - 1
            if start <= v0 < stop:
                var_hits[str(r.VARIANT_KEY)].append(c)

    accessions = sorted(
        {
            str(c.get("accession"))
            for hits in var_hits.values()
            for c in hits
            if c.get("accession")
        }
    )

    cache_root = output_file.parent.parent.parent / "cache" / "screen_ccre"
    cache_root.mkdir(parents=True, exist_ok=True)

    detail_map: dict[str, dict[str, Any]] = {}
    detail_status: dict[str, str] = {}

    workers = max(1, min(int(args.screen_detail_workers), 6))

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [
            ex.submit(
                get_screen_detail_cached,
                acc,
                api_key,
                cache_root,
                args.screen_active_z,
                args.force,
            )
            for acc in accessions
        ]

        for fut in as_completed(futs):
            acc, summary, status = fut.result()
            detail_map[acc] = summary
            detail_status[acc] = status

    rows = []

    for r in batch.itertuples(index=False):
        vkey = str(r.VARIANT_KEY)
        hits = var_hits.get(vkey, [])

        if not hits:
            rows.append(
                {
                    "VARIANT_KEY": vkey,
                    "SCREEN_STATUS": "NO_HIT",
                }
            )
            continue

        access = [str(x.get("accession", "")) for x in hits if x.get("accession")]
        details = [detail_map.get(a, {}) for a in access]

        row: dict[str, Any] = {
            "VARIANT_KEY": vkey,
            "SCREEN_STATUS": "OK",
            "SCREEN_CCRE_IDS": join_unique(access),
            "SCREEN_CCRE_GROUPS": join_unique(
                [x.get("ccre_group", "") for x in hits]
            ),
            "SCREEN_ATAC_MAX_Z": np.nanmax(
                [safe_float(x.get("atac_max_zscore")) for x in hits]
            )
            if hits else np.nan,
            "SCREEN_CTCF_MAX_Z": np.nanmax(
                [safe_float(x.get("ctcf_max_zscore")) for x in hits]
            )
            if hits else np.nan,
            "SCREEN_DNASE_MAX_Z": np.nanmax(
                [safe_float(x.get("dnase_max_zscore")) for x in hits]
            )
            if hits else np.nan,
            "SCREEN_H3K27AC_MAX_Z": np.nanmax(
                [safe_float(x.get("h3k27ac_max_zscore")) for x in hits]
            )
            if hits else np.nan,
            "SCREEN_H3K4ME3_MAX_Z": np.nanmax(
                [safe_float(x.get("h3k4me3_max_zscore")) for x in hits]
            )
            if hits else np.nan,
            "SCREEN_NEAREST_GENES": join_unique(
                [
                    g.get("gene", "")
                    for x in hits
                    for g in (x.get("nearestgenes", []) or [])
                    if isinstance(g, dict)
                ]
            ),
        }

        merge_fields = [
            "SCREEN_LINKED_GENES_EXPERIMENTAL",
            "SCREEN_LINKED_GENES_COMPUTATIONAL",
            "SCREEN_CELL_TYPES",
            "SCREEN_TISSUES",
            "SCREEN_LINK_METHODS",
            "SCREEN_LINK_ASSAYS",
            "SCREEN_EXPERIMENT_ACCESSIONS",
            "SCREEN_FILE_ACCESSIONS",
            "SCREEN_TFS",
            "SCREEN_ACTIVE_BIOSAMPLES",
            "SCREEN_ACTIVE_TISSUES",
            "SCREEN_ACTIVE_ASSAYS",
            "SCREEN_ACTIVE_LIFESTAGES",
        ]

        for field in merge_fields:
            row[field] = join_unique([d.get(field, "") for d in details], limit=200)

        row["SCREEN_N_SIGNAL_ROWS"] = sum(
            safe_int(d.get("SCREEN_N_SIGNAL_ROWS"), 0) for d in details
        )
        row["SCREEN_SILENCER_JSON"] = join_unique(
            [d.get("SCREEN_SILENCER_JSON", "") for d in details]
        )
        row["SCREEN_DYNAMIC_ENHANCERS_JSON"] = join_unique(
            [d.get("SCREEN_DYNAMIC_ENHANCERS_JSON", "") for d in details]
        )
        row["SCREEN_DETAIL_STATUS"] = join_unique(
            [detail_status.get(a, "") for a in access]
        )

        rows.append(row)

    return pd.DataFrame(rows)


# =============================================================================
# REGULOMEDB
# =============================================================================

def compact_regulomedb(
    raw: Any,
    expected_ref: str,
    expected_alt: str,
) -> dict[str, Any]:
    """Select the allele-specific RegulomeDB record for this input variant.

    RegulomeDB coordinate queries can return multiple alleles at the same
    position.  Never silently combine their scores.  Prefer an exact GRCh38
    REF/ALT match.  If only the reversed pair exists, expose the evidence but
    mark the input as REF_ALT_SWAPPED so downstream code does not treat it as
    a trusted allele-specific match.
    """
    if not isinstance(raw, dict):
        return {}

    expected_ref = str(expected_ref).upper().strip()
    expected_alt = str(expected_alt).upper().strip()
    variants = raw.get("variants", []) or []
    variants = [v for v in variants if isinstance(v, dict)]

    exact = None
    reversed_pair = None
    for v in variants:
        ref = str(v.get("ref", "")).upper().strip()
        alt = str(v.get("alt", "")).upper().strip()
        if ref == expected_ref and alt == expected_alt:
            exact = v
            break
        if ref == expected_alt and alt == expected_ref:
            reversed_pair = v

    if exact is not None:
        selected = exact
        allele_match = "EXACT"
    elif reversed_pair is not None:
        selected = reversed_pair
        allele_match = "REF_ALT_SWAPPED"
    else:
        selected = None
        allele_match = "NO_EXACT_ALLELE"

    selected = selected or {}
    score = selected.get("regulome_score", {}) or {}
    features = selected.get("features", {}) or {}

    # Some API versions provide query-level supporting evidence too. Preserve
    # it separately but keep the allele-specific selected record authoritative.
    query_features = raw.get("features", []) or []
    supporting = raw.get("total_supporting_datasets")
    if supporting is None:
        if isinstance(query_features, (dict, list)):
            supporting = len(query_features)
        else:
            supporting = 0

    return {
        "allele_match": allele_match,
        "expected_ref": expected_ref,
        "expected_alt": expected_alt,
        "matched_ref": selected.get("ref", ""),
        "matched_alt": selected.get("alt", ""),
        "regulome_score": score,
        "features": features,
        "frequency": selected.get("freq", {}) or {},
        "gencode_category": selected.get("gencode_category", ""),
        "hgvs": selected.get("hgvs", ""),
        "spdi": selected.get("spdi", ""),
        "rsids": selected.get("rsids", []) or [],
        "query_features": query_features,
        "variants": variants,
        "nearby_snps": raw.get("nearby_snps", []) or [],
        "supporting_dataset_count": safe_int(supporting, 0),
        "notifications": raw.get("notifications", {}),
    }


def query_one_regulomedb(
    row: Any,
    force: bool,
    cache_root: Path,
) -> dict[str, Any]:
    vkey = str(row.VARIANT_KEY)
    chrom = int(row.CHR)
    pos = int(row.REFERENCE_POS)
    ref = str(row.REFERENCE_REF).upper().strip()
    alt = str(row.REFERENCE_ALT).upper().strip()
    region = f"chr{chrom}:{max(0, pos - 1)}-{pos}"

    # Cache is allele-specific; a coordinate can contain multiple alternate
    # alleles and must not reuse a compacted record for another allele.
    key = hashlib.sha256(vkey.encode()).hexdigest()[:24]
    cache = cache_root / f"{key}.json.gz"

    if cache.exists() and not force:
        try:
            data = read_gzip_json(cache)
            allele_match = str(data.get("allele_match", ""))
            if allele_match == "EXACT":
                status = "OK"
            elif allele_match == "REF_ALT_SWAPPED":
                status = "ALLELE_MISMATCH"
            else:
                status = "NO_HIT"
            return {
                "VARIANT_KEY": vkey,
                "REGULOMEDB_STATUS": status,
                "REGULOMEDB_SOURCE": "CACHE",
                **regulomedb_columns(data),
            }
        except Exception:
            pass

    session = request_session()
    params = {"regions": region, "genome": "GRCh38"}
    last: Exception | None = None

    for endpoint in (REGULOMEDB_API, REGULOMEDB_FALLBACK):
        try:
            r = request_with_retry(
                session,
                "GET",
                endpoint,
                params=params,
                timeout=120,
            )
            raw = r.json()
            data = compact_regulomedb(raw, ref, alt)
            atomic_write_gzip_json(cache, data)

            allele_match = str(data.get("allele_match", ""))
            if allele_match == "EXACT":
                status = "OK"
            elif allele_match == "REF_ALT_SWAPPED":
                status = "ALLELE_MISMATCH"
            else:
                status = "NO_HIT"

            return {
                "VARIANT_KEY": vkey,
                "REGULOMEDB_STATUS": status,
                "REGULOMEDB_SOURCE": endpoint,
                **regulomedb_columns(data),
            }

        except Exception as exc:
            last = exc

    return {
        "VARIANT_KEY": vkey,
        "REGULOMEDB_STATUS": "API_ERROR",
        "REGULOMEDB_ERROR": f"{type(last).__name__}: {last}",
    }


def regulomedb_columns(data: dict[str, Any]) -> dict[str, Any]:
    score = data.get("regulome_score", {}) or {}
    ranking = ""
    probability = ""
    tissue = ""

    if isinstance(score, dict):
        ranking = score.get("ranking", score.get("rank", ""))
        probability = score.get("probability", "")
        tissue = compact_json(score.get("tissue_specific_scores", {}))
    elif score:
        ranking = str(score)

    features = data.get("features", {}) or {}
    true_features = []
    if isinstance(features, dict):
        true_features = [
            k for k, v in features.items()
            if v is True or str(v).lower() in {"true", "1", "yes"}
        ]

    return {
        "REGULOMEDB_ALLELE_MATCH": data.get("allele_match", ""),
        "REGULOMEDB_EXPECTED_REF": data.get("expected_ref", ""),
        "REGULOMEDB_EXPECTED_ALT": data.get("expected_alt", ""),
        "REGULOMEDB_MATCHED_REF": data.get("matched_ref", ""),
        "REGULOMEDB_MATCHED_ALT": data.get("matched_alt", ""),
        "REGULOMEDB_RANKING": ranking,
        "REGULOMEDB_PROBABILITY": probability,
        "REGULOMEDB_TISSUE_SCORES_JSON": tissue,
        "REGULOMEDB_GENCODE_CATEGORY": data.get("gencode_category", ""),
        "REGULOMEDB_HGVS": data.get("hgvs", ""),
        "REGULOMEDB_SPDI": data.get("spdi", ""),
        "REGULOMEDB_RSIDS": join_unique(data.get("rsids", []) or []),
        "REGULOMEDB_TRUE_FEATURES": join_unique(true_features),
        "REGULOMEDB_FEATURE_FLAGS_JSON": compact_json(features),
        "REGULOMEDB_FREQUENCIES_JSON": compact_json(data.get("frequency", {})),
        "REGULOMEDB_SUPPORTING_DATASETS": data.get("supporting_dataset_count", 0),
        "REGULOMEDB_QUERY_FEATURES_JSON": compact_json(data.get("query_features", [])),
        "REGULOMEDB_ALL_POSITION_ALLELES_JSON": compact_json(data.get("variants", [])),
    }


def worker_regulomedb(
    batch: pd.DataFrame,
    args: argparse.Namespace,
    output_file: Path,
) -> pd.DataFrame:
    cache_root = output_file.parent.parent.parent / "cache" / "regulomedb"
    cache_root.mkdir(parents=True, exist_ok=True)

    workers = max(1, min(int(args.regulomedb_workers), 4))
    rows = []

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [
            ex.submit(
                query_one_regulomedb,
                r,
                args.force,
                cache_root,
            )
            for r in batch.itertuples(index=False)
        ]

        for fut in as_completed(futs):
            rows.append(fut.result())

    return pd.DataFrame(rows)


# =============================================================================
# LOCAL RESOURCES
# =============================================================================

def tool(name: str) -> str | None:
    local = PROJECT_ROOT / "envs" / "pipeline" / "bin" / name
    return str(local) if local.exists() else shutil.which(name)


def resource_table() -> list[dict[str, Any]]:
    """Every local resource Step08 can use, with presence and the download group."""
    rows = [
        ("VEP cache", VEP_CACHE_DIR / SPECIES / f"{VEP_RELEASE}_{ASSEMBLY}", "vep", True),
        ("GRCh38 FASTA (+.fai)", FASTA, "genome", True),
        ("CADD v1.7 SNVs (+.tbi)", CADD_SNV, "cadd", False),
        ("CADD v1.7 indels (+.tbi)", CADD_INDEL, "cadd", False),
        ("AlphaMissense (+.tbi)", ALPHAMISSENSE, "annotation", False),
        ("ClinVar VCF (+.tbi)", CLINVAR, "annotation", False),
        ("ENCODE cCRE V4", CCRE_BED, "annotation", False),
        ("phyloP 100-way bigWig", PHYLOP, "annotation", False),
        ("gnomAD v4.1 constraint", CONSTRAINT, "annotation", False),
    ]
    out = []
    for name, path, group, needed in rows:
        index = {".gz": ".tbi", ".fa": ".fai"}.get(path.suffix)
        problem = resource_problem(path, index)
        out.append({"RESOURCE": name, "PATH": str(path), "GROUP": group, "REQUIRED_FOR_VEP": needed,
                    "PRESENT": not problem, "PROBLEM": problem})
    return out


def resource_problem(path: Path, index: str | None) -> str:
    """'' when the resource is usable, otherwise what is wrong with it."""
    if not path.exists():
        return "missing"
    if path.is_dir():
        if not (path / "info.txt").is_file():
            return "VEP cache incomplete (no info.txt)"
        chroms = [d for d in path.iterdir() if d.is_dir()]
        return "" if len(chroms) >= 22 else f"VEP cache incomplete ({len(chroms)} chromosome folders)"
    if path.stat().st_size == 0:
        return "empty file"
    if index:
        idx = Path(str(path) + index)
        if not idx.exists():
            return f"no {index} index"
    return ""


def size_text(path: Path) -> str:
    if not path.exists():
        return ""
    if path.is_dir():
        return "dir"
    n = float(path.stat().st_size)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return ""


def python_module_ok(name: str) -> bool:
    import importlib.util
    return importlib.util.find_spec(name) is not None


def probe_api(provider: str) -> tuple[bool, str]:
    """One small live request per provider (only with --check-resources --online)."""
    if requests is None:
        return False, "python requests not installed"
    try:
        if provider == "screen":
            if not SCREEN_API_KEY:
                return False, "SCREEN_API_KEY not set"
            r = requests.post(SCREEN_GRAPHQL, timeout=30, json={"query": "{ __typename }"},
                              headers={"Authorization": f"Bearer {SCREEN_API_KEY}",
                                       "Content-Type": "application/json", "Accept": "application/json"})
        elif provider == "regulomedb":
            r = requests.get(REGULOMEDB_API, timeout=30, headers={"Accept": "application/json"},
                             params={"regions": "chr12:40340399-40340400", "genome": "GRCh38"})
        else:
            return False, "unknown provider (use screen, regulomedb)"
    except Exception as exc:
        return False, f"{type(exc).__name__}: {str(exc)[:120]}"
    if r.status_code in (401, 403):
        return False, f"HTTP {r.status_code}: key rejected"
    return r.ok, f"HTTP {r.status_code}"


def print_resources() -> list[str]:
    """Print the local data files and tools; return the problems found."""
    banner("STEP 08 - LOCAL ANNOTATION RESOURCES")
    missing_groups, problems = set(), []
    for r in resource_table():
        state = "OK     " if r["PRESENT"] else "MISSING"
        print(f"  {state}  {r['RESOURCE']:<28} [{r['GROUP']:<10}] {size_text(Path(r['PATH'])):>9}  {r['PATH']}"
              + (f"  <- {r['PROBLEM']}" if r["PROBLEM"] else ""))
        if not r["PRESENT"]:
            missing_groups.add(r["GROUP"])
            problems.append(f"{r['RESOURCE']}: {r['PROBLEM']}")
    vep = tool("vep")
    print(f"  {'OK     ' if vep else 'MISSING'}  {'vep executable':<28} [install   ] {'':>9}  {vep or ''}")
    if not vep:
        problems.append("vep executable: not found (Step00_Installation.py)")
    for module, use in (("pysam", "CADD / AlphaMissense / ClinVar lookups"), ("pyBigWig", "phyloP"),
                        ("requests", "--online only")):
        ok = python_module_ok(module)
        print(f"  {'OK     ' if ok else 'MISSING'}  {'python ' + module:<28} [install   ] {'':>9}  {use}")
        if not ok and module != "requests":
            problems.append(f"python {module}: not installed (Step00_Installation.py)")
    if missing_groups:
        print("\nDownload the missing groups with:")
        print("  envs/pipeline/bin/python Step00_Download_Resources.py --groups "
              + " ".join(sorted(missing_groups)) + " --output Step00_Download_Step08.sh")
        print("  sbatch Step00_Download_Step08.sh")
    return problems


def check_resources_mode(args) -> None:
    """--check-resources: is everything there to annotate this phenotype / ancestry / panel?"""
    P = project_paths(args)
    problems = print_resources()

    banner(f"STEP 08 - INPUT FROM STEP06 - {P['phenotype']} / {P['label']} - LD panel {P['panel']}")
    ready, skipped = finemapped_studies(P)
    print(f"  Fine-mapping folder : {P['finemap']}  ({'exists' if P['finemap'].exists() else 'MISSING'})")
    print(f"  Ready to annotate   : {len(ready)} studies  {' '.join(ready[:10])}{' ...' if len(ready) > 10 else ''}")
    for accession, why in skipped:
        print(f"  Skipped             : {accession}  ({why})")
    if not ready:
        problems.append(f"no fine-mapped studies under {P['finemap']} - run Step06"
                        f"{panel_option(P['panel'])} first")
    print(f"  Results will go to  : {P['out']}")

    banner("STEP 08 - WEB APIs (optional, only used with --online)")
    print("  The default annotation is fully local: no API and no internet needed.")
    print(f"  SCREEN_API_KEY      : {'SET (from ' + SCREEN_API_KEY_SOURCE + ')' if SCREEN_API_KEY else 'not set'}"
          f"  - needed only for --online screen")
    print(f"                        set it with: export SCREEN_API_KEY=...  or  SCREEN_API_KEY=... in "
          f"{API_KEYS_FILE.relative_to(PROJECT_ROOT).as_posix()}")
    print("  RegulomeDB          : public, no key needed")
    print(f"  Endpoints           : SCREEN {SCREEN_GRAPHQL}")
    print(f"                        RegulomeDB {REGULOMEDB_API}")
    providers = [x.strip().lower() for x in (args.online or "").split(",") if x.strip()]
    for provider in providers:
        ok, detail = probe_api(provider)
        print(f"  {'OK     ' if ok else 'FAILED '}  live test {provider:<12} {detail}")
        if not ok:
            problems.append(f"--online {provider}: {detail}")
    if not providers:
        print("  (add --online screen,regulomedb to test the APIs live from this node)")

    banner("STEP 08 - CHECK RESULT")
    if problems:
        print(f"NOT READY - {len(problems)} problem(s):")
        for problem in problems:
            print(f"  - {problem}")
        raise SystemExit(1)
    print("READY - every annotation resource, tool and Step06 input is present.")
    base = (f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {shlex.quote(P['phenotype'])} "
            f"--ancestry {P['code']}{panel_option(P['panel'])}")
    print(f"NEXT:\n  {base}")


def tabix_contig(contigs: Iterable[str], chrom: int) -> str | None:
    contigs = set(contigs)
    for name in (str(chrom), f"chr{chrom}"):
        if name in contigs:
            return name
    return None


# =============================================================================
# 1. ENSEMBL VEP (offline cache)
# =============================================================================

def write_vep_input(universe: pd.DataFrame, path: Path) -> None:
    rows = universe.sort_values(["CHR", "REFERENCE_POS"])
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in rows.itertuples(index=False):
            handle.write(f"{r.CHR}\t{r.REFERENCE_POS}\t{r.VARIANT_KEY}\t{r.REFERENCE_REF}\t{r.REFERENCE_ALT}\t.\t.\t.\n")


def write_synonyms(path: Path) -> Path | None:
    """VEP input uses 1..22 like the cache; the GENCODE FASTA uses chr1..: map them."""
    fai = Path(str(FASTA) + ".fai")
    if not fai.exists():
        return None
    first = fai.read_text(encoding="utf-8", errors="replace").split("\t", 1)[0]
    if not first.startswith("chr"):
        return None
    names = [str(c) for c in range(1, 23)] + ["X", "Y"]
    path.write_text("".join(f"{c}\tchr{c}\n" for c in names) + "MT\tchrM\n", encoding="utf-8")
    return path


def vep_command(vcf: Path, out_json: Path, threads: int, args) -> list[str]:
    vep = args.vep_path or tool("vep")
    if not vep:
        raise RuntimeError("vep not found. Run: python Step00_Installation.py")
    command = [vep, "--input_file", str(vcf), "--format", "vcf", "--output_file", str(out_json), "--json",
               "--force_overwrite", "--no_stats", "--species", SPECIES, "--assembly", ASSEMBLY,
               "--fork", str(max(1, threads)), "--buffer_size", "5000",
               "--hgvs", "--symbol", "--canonical", "--mane", "--biotype", "--numbers", "--domains",
               "--uniprot", "--protein", "--variant_class", "--regulatory", "--sift", "b", "--polyphen", "b",
               "--check_existing", "--af", "--af_gnomade", "--af_gnomadg", "--max_af", "--appris", "--tsl",
               "--ccds", "--gene_phenotype"]
    if args.vep_gtf:
        # Testing only: gene models from a bgzipped + tabix-indexed GTF instead of the
        # cache. Frequencies, rsIDs, SIFT/PolyPhen, UniProt, domains and regulatory
        # features exist only in the cache, so those options are dropped here.
        cache_only = {"--af", "--af_gnomade", "--af_gnomadg", "--max_af", "--check_existing", "--regulatory",
                      "--domains", "--uniprot", "--appris", "--tsl", "--ccds", "--gene_phenotype", "--mane"}
        keep, skip_next = [], False
        for x in command:
            if skip_next:
                skip_next = False
            elif x in ("--sift", "--polyphen"):
                skip_next = True
            elif x not in cache_only:
                keep.append(x)
        keep += ["--gtf", str(args.vep_gtf), "--fasta", str(FASTA)]
        synonyms = write_synonyms(out_json.with_name("chr_synonyms.txt"))
        return keep + (["--synonyms", str(synonyms)] if synonyms else [])
    cache = VEP_CACHE_DIR / SPECIES / f"{VEP_RELEASE}_{ASSEMBLY}"
    if not cache.exists():
        raise RuntimeError(f"VEP cache missing: {cache}\nRun: Step00_Download_Resources.py --groups vep")
    if not FASTA.exists():
        raise RuntimeError(f"GRCh38 FASTA missing: {FASTA}\nRun: Step00_Download_Resources.py --groups genome")
    command += ["--offline", "--cache", "--dir_cache", str(VEP_CACHE_DIR), "--cache_version", VEP_RELEASE,
                "--fasta", str(FASTA)]
    synonyms = write_synonyms(out_json.with_name("chr_synonyms.txt"))
    if synonyms:
        command += ["--synonyms", str(synonyms)]
    return command


def frequency_columns(item: dict[str, Any], alt: str) -> dict[str, Any]:
    """gnomAD v4 / 1000G allele frequencies of the ALT allele (co-located variants)."""
    out: dict[str, Any] = {}
    rsids = []
    for colocated in item.get("colocated_variants", []) or []:
        if not isinstance(colocated, dict):
            continue
        cid = str(colocated.get("id", ""))
        if cid.startswith("rs"):
            rsids.append(cid)
        freqs = (colocated.get("frequencies") or {}).get(alt) or {}
        for key, column in (("gnomadg", "GNOMADG_AF"), ("gnomadg_nfe", "GNOMADG_NFE_AF"),
                            ("gnomade", "GNOMADE_AF"), ("gnomade_nfe", "GNOMADE_NFE_AF"),
                            ("af", "KG_AF"), ("eur", "KG_EUR_AF"), ("max_af", "MAX_AF")):
            if key in freqs and column not in out:
                out[column] = freqs[key]
        if "max_af_pops" in freqs and "MAX_AF_POPS" not in out:
            out["MAX_AF_POPS"] = ",".join(freqs["max_af_pops"]) if isinstance(freqs["max_af_pops"], list) \
                else freqs["max_af_pops"]
    out["VEP_RSIDS"] = join_unique(rsids)
    return out


def annotate_vep(universe: pd.DataFrame, work: Path, threads: int, args) -> tuple[pd.DataFrame, dict]:
    vcf, out_json = work / "vep_input.vcf", work / "vep_output.json"
    write_vep_input(universe, vcf)
    command = vep_command(vcf, out_json, threads, args)
    print("$ " + " ".join(shlex.quote(x) for x in command), flush=True)
    log = work / "vep.log"
    with open(log, "w", encoding="utf-8") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode != 0 or not out_json.exists():
        tail = log.read_text(encoding="utf-8", errors="replace")[-3000:]
        raise RuntimeError(f"VEP failed (exit {result.returncode}). Log: {log}\n{tail}")

    alt_of = dict(zip(universe["VARIANT_KEY"], universe["REFERENCE_ALT"]))
    rows, version = [], ""
    with open(out_json, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            item = json.loads(line)
            key = str(item.get("id", ""))
            record = extract_ensembl_vep_record(item, key, key)
            record.update(frequency_columns(item, alt_of.get(key, "")))
            rows.append(record)
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        if "ensembl-vep" in line.lower() or "VEP version" in line:
            version = line.strip()
            break
    table = pd.DataFrame(rows)
    if not table.empty:
        table = table.drop(columns=[c for c in ("API_ID",) if c in table.columns])
        table = table.rename(columns={c: "VEP_" + c[len("ENSEMBL_VEP_"):] if c.startswith("ENSEMBL_VEP_")
                                      else c.replace("ENSEMBL_", "VEP_", 1) for c in table.columns
                                      if c.startswith("ENSEMBL_")})
    info = {"source": "Ensembl VEP " + (f"--gtf {args.vep_gtf} (test)" if args.vep_gtf else f"{VEP_RELEASE} offline cache"),
            "annotated": int(len(table)), "log_version_line": version}
    return table, info


# =============================================================================
# 2-4. TABIX LOOKUPS: CADD, ALPHAMISSENSE, CLINVAR
# =============================================================================

def tabix_lookup(universe: pd.DataFrame, path: Path, parse, label: str) -> tuple[pd.DataFrame, dict]:
    """parse(fields, ref, alt) -> dict or None for every record at the variant position."""
    import pysam

    if not path.exists() or not Path(str(path) + ".tbi").exists():
        return pd.DataFrame({"VARIANT_KEY": universe["VARIANT_KEY"], f"{label}_STATUS": "NOT_AVAILABLE"}), \
            {"source": str(path), "status": "NOT_AVAILABLE"}
    rows, hits = [], 0
    with pysam.TabixFile(str(path)) as tbx:
        contigs = tbx.contigs
        for r in universe.itertuples(index=False):
            contig = tabix_contig(contigs, r.CHR)
            out = {"VARIANT_KEY": r.VARIANT_KEY, f"{label}_STATUS": "NO_HIT"}
            if contig is not None:
                try:
                    lines = tbx.fetch(contig, int(r.REFERENCE_POS) - 1, int(r.REFERENCE_POS))
                except ValueError:
                    lines = []
                for line in lines:
                    fields = line.split("\t")
                    try:
                        if int(fields[1]) != int(r.REFERENCE_POS):
                            continue
                    except (ValueError, IndexError):
                        continue
                    value = parse(fields, r.REFERENCE_REF, r.REFERENCE_ALT)
                    if value:
                        out.update(value)
                        out[f"{label}_STATUS"] = "OK"
                        hits += 1
                        break
            rows.append(out)
    return pd.DataFrame(rows), {"source": str(path), "status": "OK", "hits": hits}


def parse_cadd(fields, ref, alt):
    # #Chrom Pos Ref Alt RawScore PHRED
    if len(fields) >= 6 and fields[2] == ref and fields[3] == alt:
        return {"CADD_RAW": float(fields[4]), "CADD_PHRED": float(fields[5])}
    return None


def annotate_cadd(universe: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    snv = (universe["REFERENCE_REF"].str.len() == 1) & (universe["REFERENCE_ALT"].str.len() == 1)
    a, info_a = tabix_lookup(universe[snv], CADD_SNV, parse_cadd, "CADD")
    b, info_b = tabix_lookup(universe[~snv], CADD_INDEL, parse_cadd, "CADD")
    status = "OK" if "OK" in (info_a["status"], info_b["status"]) else "NOT_AVAILABLE"
    return pd.concat([a, b], ignore_index=True), {"snv": info_a, "indel": info_b, "version": "CADD v1.7",
                                                  "status": status}


def parse_alphamissense(fields, ref, alt):
    # #CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class
    if len(fields) >= 10 and fields[2] == ref and fields[3] == alt:
        return {"AM_PATHOGENICITY": float(fields[8]), "AM_CLASS": fields[9].strip(),
                "AM_UNIPROT": fields[5], "AM_TRANSCRIPT": fields[6], "AM_PROTEIN_VARIANT": fields[7]}
    return None


def parse_clinvar(fields, ref, alt):
    # VCF: CHROM POS ID REF ALT QUAL FILTER INFO
    if len(fields) < 8 or fields[3] != ref or alt not in fields[4].split(","):
        return None
    info = dict(x.split("=", 1) for x in fields[7].split(";") if "=" in x)
    return {"CLINVAR_VARIATION_ID": fields[2],
            "CLINVAR_SIGNIFICANCE": info.get("CLNSIG", "").replace("_", " "),
            "CLINVAR_REVIEW_STATUS": info.get("CLNREVSTAT", "").replace("_", " "),
            "CLINVAR_CONDITIONS": info.get("CLNDN", "").replace("_", " "),
            "CLINVAR_ALLELE_ID": info.get("ALLELEID", "")}


def clinvar_release() -> str:
    if not CLINVAR.exists():
        return ""
    with gzip.open(CLINVAR, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("##fileDate="):
                return line.strip().split("=", 1)[1]
            if not line.startswith("##"):
                break
    return ""


# =============================================================================
# 5-7. INTERVALS / TRACKS: cCRE, phyloP, gene constraint
# =============================================================================

_CCRE_CACHE: dict | None = None


def load_ccres() -> dict | None:
    global _CCRE_CACHE
    if _CCRE_CACHE is not None or not CCRE_BED.exists():
        return _CCRE_CACHE
    bed = pd.read_csv(CCRE_BED, sep="\t", header=None, usecols=[0, 1, 2, 4, 5],
                      names=["CHR", "START", "END", "ACCESSION", "CLASS"], dtype={"CHR": str})
    bed["CHR"] = bed["CHR"].str.replace("chr", "", regex=False)
    _CCRE_CACHE = {}
    for chrom, group in bed.groupby("CHR"):
        group = group.sort_values("START")
        _CCRE_CACHE[chrom] = (group["START"].to_numpy(), group["END"].to_numpy(),
                              group["ACCESSION"].to_numpy(), group["CLASS"].to_numpy())
    return _CCRE_CACHE


def annotate_ccre(universe: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    ccres = load_ccres()
    if ccres is None:
        return pd.DataFrame({"VARIANT_KEY": universe["VARIANT_KEY"], "CCRE_STATUS": "NOT_AVAILABLE"}), \
            {"source": str(CCRE_BED), "status": "NOT_AVAILABLE"}
    rows = []
    for r in universe.itertuples(index=False):
        out = {"VARIANT_KEY": r.VARIANT_KEY, "CCRE_ACCESSION": "", "CCRE_CLASS": "", "IN_CCRE": False}
        data = ccres.get(str(r.CHR))
        if data is not None:
            starts, ends, acc, cls = data
            i = np.searchsorted(starts, r.REFERENCE_POS - 1, side="right") - 1   # BED 0-based, half-open
            if i >= 0 and r.REFERENCE_POS <= ends[i]:
                out.update({"CCRE_ACCESSION": acc[i], "CCRE_CLASS": cls[i], "IN_CCRE": True})
        rows.append(out)
    table = pd.DataFrame(rows)
    return table, {"source": str(CCRE_BED), "status": "OK", "in_ccre": int(table["IN_CCRE"].sum())}


def annotate_phylop(universe: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    if not PHYLOP.exists():
        return pd.DataFrame({"VARIANT_KEY": universe["VARIANT_KEY"], "PHYLOP_STATUS": "NOT_AVAILABLE"}), \
            {"source": str(PHYLOP), "status": "NOT_AVAILABLE"}
    import pyBigWig
    rows = []
    with pyBigWig.open(str(PHYLOP)) as bw:
        chroms = bw.chroms()
        for r in universe.itertuples(index=False):
            name = tabix_contig(chroms.keys(), r.CHR)
            value = np.nan
            if name is not None:
                start = int(r.REFERENCE_POS) - 1
                end = start + max(1, len(str(r.REFERENCE_REF)))
                vals = [v for v in bw.values(name, start, min(end, chroms[name])) if v == v]
                value = float(np.max(vals)) if vals else np.nan
            rows.append({"VARIANT_KEY": r.VARIANT_KEY, "PHYLOP_100WAY": value})
    return pd.DataFrame(rows), {"source": str(PHYLOP), "status": "OK"}


_CONSTRAINT_CACHE: pd.DataFrame | None = None


def gene_constraint() -> pd.DataFrame | None:
    """LOEUF / pLI / missense Z per Ensembl gene (canonical, else MANE transcript)."""
    global _CONSTRAINT_CACHE
    if _CONSTRAINT_CACHE is not None or not CONSTRAINT.exists():
        return _CONSTRAINT_CACHE
    c = pd.read_csv(CONSTRAINT, sep="\t", dtype=str,
                    usecols=lambda x: x in ("gene", "gene_id", "canonical", "mane_select", "lof.oe_ci.upper",
                                            "lof.pLI", "mis.z_score", "lof.oe"))
    c = c[c["gene_id"].fillna("").str.startswith("ENSG")].copy()
    c["_rank"] = np.where(c["canonical"].str.lower() == "true", 0,
                          np.where(c["mane_select"].str.lower() == "true", 1, 2))
    c = c.sort_values("_rank").drop_duplicates("gene_id")
    _CONSTRAINT_CACHE = pd.DataFrame({
        "GENE_ID_NOVERSION": c["gene_id"].str.split(".").str[0],
        "GNOMAD_LOEUF": pd.to_numeric(c["lof.oe_ci.upper"], errors="coerce"),
        "GNOMAD_PLI": pd.to_numeric(c.get("lof.pLI"), errors="coerce"),
        "GNOMAD_MIS_Z": pd.to_numeric(c.get("mis.z_score"), errors="coerce"),
    })
    return _CONSTRAINT_CACHE


# =============================================================================
# OPTIONAL ONLINE (credible-set variants only)
# =============================================================================

def annotate_online(universe: pd.DataFrame, cs_keys: set, P: dict, accession: str, args) -> tuple[pd.DataFrame, dict]:
    providers = [x.strip() for x in (args.online or "").split(",") if x.strip()]
    if not providers:
        return pd.DataFrame({"VARIANT_KEY": universe["VARIANT_KEY"]}), {}
    batch = universe[universe["VARIANT_KEY"].isin(cs_keys)].copy()
    anchor = P["out"] / "_online" / accession / "shard.tsv"   # cache: <out>/cache/<provider>
    anchor.parent.mkdir(parents=True, exist_ok=True)
    online_args = argparse.Namespace(force=False, regulomedb_workers=2, screen_detail_workers=2,
                                     screen_active_z=args.screen_active_z)
    result, info = pd.DataFrame({"VARIANT_KEY": universe["VARIANT_KEY"]}), {}
    for provider in providers:
        if batch.empty:
            break
        if provider == "screen":
            table = worker_screen(batch, online_args, anchor)
        elif provider == "regulomedb":
            table = worker_regulomedb(batch, online_args, anchor)
        else:
            raise SystemExit(f"Unknown --online provider {provider!r} (screen, regulomedb)")
        info[provider] = {"variants_queried": int(len(batch))}
        if not table.empty and "VARIANT_KEY" in table.columns:
            result = result.merge(table.drop_duplicates("VARIANT_KEY"), on="VARIANT_KEY", how="left")
    return result, info


# =============================================================================
# ONE STUDY
# =============================================================================

# Columns the summaries rely on; created empty when a source is skipped / missing.
ANNOTATION_COLUMNS = {
    "VEP_BEST_CONSEQUENCE": "", "VEP_FUNCTIONAL_CATEGORY": "NOT_ANNOTATED", "VEP_GENE_SYMBOL": "",
    "VEP_GENE_ID": "", "VEP_HGVSC": "", "VEP_HGVSP": "", "CADD_PHRED": np.nan, "CADD_RAW": np.nan,
    "AM_PATHOGENICITY": np.nan, "AM_CLASS": "", "CLINVAR_SIGNIFICANCE": "", "CCRE_ACCESSION": "",
    "CCRE_CLASS": "", "IN_CCRE": False, "PHYLOP_100WAY": np.nan, "GNOMADG_AF": np.nan,
    "GNOMADG_NFE_AF": np.nan, "GNOMAD_LOEUF": np.nan,
}


def ensure_columns(m: pd.DataFrame) -> pd.DataFrame:
    for column, default in ANNOTATION_COLUMNS.items():
        if column not in m.columns:
            m[column] = default
    m["VEP_FUNCTIONAL_CATEGORY"] = m["VEP_FUNCTIONAL_CATEGORY"].fillna("NOT_ANNOTATED")
    return m


LOF_TERMS = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
             "start_lost", "stop_lost", "transcript_ablation"}


def derive_flags(m: pd.DataFrame) -> pd.DataFrame:
    cons = m.get("VEP_BEST_CONSEQUENCE", pd.Series("", index=m.index)).fillna("").astype(str)
    m["IS_LOF"] = cons.apply(lambda s: bool(set(s.split(";")) & LOF_TERMS))
    m["IS_MISSENSE"] = cons.str.contains("missense_variant")
    # missense / LoF count as coding even when a splice-region term ranks first
    m["IS_CODING"] = m["VEP_FUNCTIONAL_CATEGORY"].astype(str).eq("CODING") | m["IS_LOF"] | m["IS_MISSENSE"]
    m["IS_SPLICE"] = m["VEP_FUNCTIONAL_CATEGORY"].astype(str).eq("SPLICING")
    m["CADD_GE_20"] = pd.to_numeric(m.get("CADD_PHRED"), errors="coerce") >= 20
    m["AM_LIKELY_PATHOGENIC"] = m.get("AM_CLASS", pd.Series("", index=m.index)).fillna("").eq("likely_pathogenic")
    sig = m.get("CLINVAR_SIGNIFICANCE", pd.Series("", index=m.index)).fillna("").str.lower()
    m["CLINVAR_PATHOGENIC"] = sig.str.contains("pathogenic") & ~sig.str.contains("conflicting")
    m["PHYLOP_CONSERVED"] = pd.to_numeric(m.get("PHYLOP_100WAY"), errors="coerce") >= 2.0
    evidence = []
    for _, r in m.iterrows():
        e = []
        if r["IS_LOF"]:
            e.append("LoF")
        elif r["IS_MISSENSE"]:
            e.append("missense")
        elif r["IS_SPLICE"]:
            e.append("splice")
        if r["AM_LIKELY_PATHOGENIC"]:
            e.append("AlphaMissense likely pathogenic")
        if r["CADD_GE_20"]:
            e.append("CADD>=20")
        if r["CLINVAR_PATHOGENIC"]:
            e.append("ClinVar P/LP")
        if str(r.get("IN_CCRE", "")) == "True":
            e.append(f"cCRE {r.get('CCRE_CLASS', '')}")
        if r["PHYLOP_CONSERVED"]:
            e.append("conserved")
        evidence.append("; ".join(e))
    m["FUNCTIONAL_EVIDENCE"] = evidence
    return m


def locus_table(m: pd.DataFrame) -> pd.DataFrame:
    rows = []
    category = m.get("VEP_FUNCTIONAL_CATEGORY", pd.Series("NOT_ANNOTATED", index=m.index)).fillna("NOT_ANNOTATED")
    for locus, g in m.groupby("LOCUS_ID", sort=True):
        cs = g[g["IN_95_CREDIBLE_SET"]]
        top = g.sort_values("PIP", ascending=False).iloc[0]
        pip = g["PIP"].clip(lower=0)
        total = float(pip.sum()) or 1.0
        weights = pip.groupby(category.loc[g.index]).sum() / total
        rows.append({
            "LOCUS_ID": locus,
            "CHR": int(top["CHR"]),
            "N_VARIANTS": int(len(g)),
            "N_CS_VARIANTS": int(len(cs)),
            "TOP_VARIANT": top["VARIANT_KEY"],
            "TOP_RSID": top.get("PRIMARY_RSID", ""),
            "TOP_PIP": round(float(top["PIP"]), 4),
            "TOP_GENE": top.get("VEP_GENE_SYMBOL", ""),
            "TOP_CONSEQUENCE": top.get("VEP_BEST_CONSEQUENCE", ""),
            "TOP_EVIDENCE": top.get("FUNCTIONAL_EVIDENCE", ""),
            "CS_GENES": join_unique(cs.get("VEP_GENE_SYMBOL", [])),
            "CS_CODING_VARIANTS": join_unique(cs.loc[cs["IS_CODING"], "VARIANT_KEY"]) if len(cs) else "",
            "CS_LOF_OR_MISSENSE": join_unique(cs.loc[cs["IS_LOF"] | cs["IS_MISSENSE"], "VARIANT_KEY"]) if len(cs) else "",
            "CS_CLINVAR_PATHOGENIC": join_unique(cs.loc[cs["CLINVAR_PATHOGENIC"], "VARIANT_KEY"]) if len(cs) else "",
            "CS_MAX_CADD_PHRED": pd.to_numeric(cs.get("CADD_PHRED"), errors="coerce").max() if len(cs) else np.nan,
            "CS_N_IN_CCRE": int((cs.get("IN_CCRE", pd.Series(dtype=object)).astype(str) == "True").sum()),
            "CS_CCRE_CLASSES": join_unique(cs.get("CCRE_CLASS", [])),
            **{f"PIP_FRACTION_{k}": round(float(v), 4) for k, v in weights.items()},
        })
    return pd.DataFrame(rows)


def annotate_universe(universe: pd.DataFrame, work: Path, threads: int, args, P=None, accession="",
                      cs_keys=frozenset()) -> tuple[pd.DataFrame, dict]:
    provenance: dict[str, Any] = {}
    merged = universe.copy()
    steps = []
    if not args.no_vep:
        steps.append(("vep", lambda: annotate_vep(universe, work, threads, args)))
    steps += [("cadd", lambda: annotate_cadd(universe)),
              ("alphamissense", lambda: tabix_lookup(
                  universe[(universe["REFERENCE_REF"].str.len() == 1) & (universe["REFERENCE_ALT"].str.len() == 1)],
                  ALPHAMISSENSE, parse_alphamissense, "AM")),
              ("clinvar", lambda: tabix_lookup(universe, CLINVAR, parse_clinvar, "CLINVAR")),
              ("ccre", lambda: annotate_ccre(universe)),
              ("phylop", lambda: annotate_phylop(universe))]
    for name, run in steps:
        started = time.time()
        table, info = run()
        info["seconds"] = round(time.time() - started, 1)
        provenance[name] = info
        print(f"  {name:<14} {info.get('status', 'OK'):<14} {info['seconds']:>7}s", flush=True)
        if not table.empty:
            merged = merged.merge(table.drop_duplicates("VARIANT_KEY"), on="VARIANT_KEY", how="left")
    provenance["clinvar"]["release"] = clinvar_release()

    constraint = gene_constraint()
    if constraint is not None and "VEP_GENE_ID" in merged.columns:
        merged["GENE_ID_NOVERSION"] = merged["VEP_GENE_ID"].fillna("").astype(str).str.split(".").str[0]
        merged = merged.merge(constraint, on="GENE_ID_NOVERSION", how="left").drop(columns="GENE_ID_NOVERSION")
        provenance["gnomad_constraint"] = {"source": str(CONSTRAINT), "status": "OK"}
    else:
        provenance["gnomad_constraint"] = {"source": str(CONSTRAINT),
                                           "status": "NOT_AVAILABLE" if constraint is None else "NO_VEP_GENES"}

    if P is not None and args.online:
        table, info = annotate_online(universe, set(cs_keys), P, accession, args)
        provenance["online"] = info
        merged = merged.merge(table.drop_duplicates("VARIANT_KEY"), on="VARIANT_KEY", how="left")
    return merged, provenance


def annotate_study(P: dict, accession: str, args, threads: int) -> None:
    out_dir = P["out"] / accession
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / f"{accession}_annotated_variants.tsv.gz"
    provenance_file = out_dir / "annotation_provenance.json"
    if final.exists() and provenance_file.exists() and not args.force:
        print(f"{accession}: already annotated (use --force to redo)")
        return

    source = P["finemap"] / accession / f"{accession}_finemapped_variants.tsv.gz"
    study = load_finemapped(source)
    study["STUDY_ACCESSION"] = accession
    if args.scope == "cs":
        study = study[study["IN_95_CREDIBLE_SET"]]
    elif args.scope == "pip":
        study = study[study["PIP"] >= args.min_pip]
    if "Z_REF" in study.columns:
        study["Z_ALT"] = -pd.to_numeric(study["Z_REF"], errors="coerce")
    if "BETA_REF" in study.columns:
        study["BETA_ALT"] = -pd.to_numeric(study["BETA_REF"], errors="coerce")

    universe = study[["VARIANT_KEY", "PRIMARY_RSID", "CHR", "REFERENCE_POS", "REFERENCE_REF", "REFERENCE_ALT"]] \
        .drop_duplicates("VARIANT_KEY").reset_index(drop=True)
    banner(f"STEP 08 - {accession}: {len(universe):,} unique variants ({args.scope}), {threads} threads")

    work = out_dir / "work"
    work.mkdir(exist_ok=True)
    cs_keys = set(study.loc[study["IN_95_CREDIBLE_SET"], "VARIANT_KEY"])
    annotation, provenance = annotate_universe(universe, work, threads, args, P, accession, cs_keys)

    keep = [c for c in annotation.columns if c not in ("PRIMARY_RSID", "CHR", "REFERENCE_POS",
                                                        "REFERENCE_REF", "REFERENCE_ALT")]
    master = study.merge(annotation[keep], on="VARIANT_KEY", how="left")
    master = derive_flags(ensure_columns(master))
    front = ["STUDY_ACCESSION", "LOCUS_ID", "VARIANT_KEY", "PRIMARY_RSID", "CHR", "REFERENCE_POS",
             "REFERENCE_REF", "REFERENCE_ALT", "PIP", "IN_95_CREDIBLE_SET", "CREDIBLE_SET", "P", "Z_ALT",
             "BETA_ALT", "VEP_GENE_SYMBOL", "VEP_BEST_CONSEQUENCE", "VEP_FUNCTIONAL_CATEGORY",
             "VEP_HGVSC", "VEP_HGVSP", "CADD_PHRED", "AM_PATHOGENICITY", "AM_CLASS", "CLINVAR_SIGNIFICANCE",
             "CCRE_CLASS", "PHYLOP_100WAY", "GNOMADG_AF", "GNOMADG_NFE_AF", "GNOMAD_LOEUF",
             "FUNCTIONAL_EVIDENCE"]
    front = [c for c in front if c in master.columns]
    master = master[front + [c for c in master.columns if c not in front]]
    master = master.sort_values(["LOCUS_ID", "PIP"], ascending=[True, False], kind="stable")

    master.to_csv(final.with_suffix(".tmp"), sep="\t", index=False, compression="gzip")
    os.replace(final.with_suffix(".tmp"), final)
    master[master["IN_95_CREDIBLE_SET"]].to_csv(out_dir / f"{accession}_credible_set_variants.tsv",
                                               sep="\t", index=False)
    loci = locus_table(master)
    loci.to_csv(out_dir / f"{accession}_locus_annotation_summary.tsv", sep="\t", index=False)

    provenance.update({
        "study": accession, "version": VERSION, "scope": args.scope, "n_variants": int(len(master)),
        "n_unique_variants": int(len(universe)), "n_credible_set_variants": int(master["IN_95_CREDIBLE_SET"].sum()),
        "n_loci": int(master["LOCUS_ID"].nunique()), "input": str(source),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "resources": resource_table(),
    })
    provenance_file.write_text(json.dumps(provenance, indent=2, default=str), encoding="utf-8")
    (out_dir / "ANNOTATION_FAILED.txt").unlink(missing_ok=True)
    if not args.keep_work:
        shutil.rmtree(work, ignore_errors=True)
    print(f"{accession}: {len(master):,} variants, {len(loci)} loci -> {final}")


# =============================================================================
# PATHS / STUDIES / SLURM
# =============================================================================

DEFAULT_PANEL = "1000G"


def results_slug(ancestry_slug: str, panel: str = DEFAULT_PANEL) -> str:
    """Results folder per LD reference panel: <ancestry> (1000G) or <ancestry>_<panel>."""
    return ancestry_slug if panel == DEFAULT_PANEL else f"{ancestry_slug}_{panel}"


def panel_option(panel: str) -> str:
    return "" if panel == DEFAULT_PANEL else f" --reference-panel {panel}"


def project_paths(args) -> dict:
    code, label = canonical_ancestry(args.ancestry)
    p_slug, a_slug = slugify(args.phenotype), slugify(label)
    panel = args.reference_panel
    r_slug = results_slug(a_slug, panel)
    # a_slug: Step02 study numbers (shared by every panel); r_slug: this panel's results.
    return {"phenotype": args.phenotype.strip(), "code": code, "label": label, "p_slug": p_slug,
            "a_slug": a_slug, "r_slug": r_slug, "panel": panel,
            "finemap": PROJECT_ROOT / "06_finemapping" / p_slug / r_slug,
            "out": PROJECT_ROOT / "08_annotation" / p_slug / r_slug}


def study_for_index(phenotype_slug: str, ancestry_slug: str, index: int) -> str:
    """Accession of study N in the Step02 download manifest (the ARRAY_ID the
    Step02 planner prints) - the same N in Steps 02, 05, 06, 07 and 08."""
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


def finemapped_studies(P: dict) -> tuple[list[str], list[tuple[str, str]]]:
    ready, skipped = [], []
    if not P["finemap"].exists():
        return ready, skipped
    for d in sorted(p for p in P["finemap"].iterdir() if p.is_dir()):
        summary = d / "finemapping_summary.json"
        variants = d / f"{d.name}_finemapped_variants.tsv.gz"
        if not summary.exists():
            continue
        info = json.loads(summary.read_text(encoding="utf-8"))
        if info.get("STATUS") == "NO_SIGNIFICANT_LOCI" or int(info.get("N_LOCI_SUCCESS", 0) or 0) == 0:
            skipped.append((d.name, "no fine-mapped loci"))
        elif not variants.exists() or variants.stat().st_size < 100:
            skipped.append((d.name, "fine-mapped variant file missing"))
        else:
            ready.append(d.name)
    return ready, skipped


PER_JOB_KEYS = {"job-name", "output", "error", "array"}


def read_slurm_template(path):
    values = {}
    if Path(path).is_file():
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            match = re.match(r"^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
            if match and (match.group(2) or "").strip():
                values[match.group(1)] = match.group(2).strip()
    return values


def resolve_slurm(args):
    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": "8", "mem": "32G", "time": "12:00:00"}
    values.update({k: v for k, v in read_slurm_template(PROJECT_ROOT / args.template).items()
                   if k not in PER_JOB_KEYS})
    for key, value in (("partition", args.partition), ("qos", args.qos), ("account", args.account),
                       ("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time)):
        if value not in (None, ""):
            values[key] = str(value)
    args.slurm_values = values
    args.cpus = int(values["cpus-per-task"])
    return args


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


def passthrough_options(args) -> list[str]:
    options = ["--scope", args.scope, "--min-pip", str(args.min_pip)]
    if args.online:
        options += ["--online", args.online]
    if args.no_vep:
        options.append("--no-vep")
    if args.force:
        options.append("--force")
    if args.reference_panel != DEFAULT_PANEL:
        options += ["--reference-panel", args.reference_panel]
    return options


# =============================================================================
# MODES
# =============================================================================

def planner_mode(args) -> None:
    P = project_paths(args)
    print_resources()
    ready, skipped = finemapped_studies(P)
    banner(f"STEP 08 PLANNER - {P['phenotype']} / {P['label']} ({P['code']}) - LD panel {P['panel']}")
    for accession, reason in skipped:
        print(f"SKIPPED  {accession}: {reason}")
    if not ready:
        raise SystemExit(f"No fine-mapped studies under {P['finemap']} - run Step06 first"
                         f"{' with' + panel_option(P['panel']) if P['panel'] != DEFAULT_PANEL else ''}.")
    if not args.no_vep and not args.vep_gtf and not all(
            r["PRESENT"] for r in resource_table() if r["REQUIRED_FOR_VEP"]):
        raise SystemExit("VEP cache / FASTA missing (see above), or add --no-vep to annotate without VEP.")

    P["out"].mkdir(parents=True, exist_ok=True)
    n_tasks = min(len(ready), MAX_ARRAY_TASKS)
    manifest = pd.DataFrame({"TASK_ID": [i % n_tasks + 1 for i in range(len(ready))], "STUDY_ACCESSION": ready})
    manifest_file = P["out"] / "annotation_manifest.tsv"
    manifest.to_csv(manifest_file, sep="\t", index=False)

    log_dir = PROJECT_ROOT / "logs" / "step08_annotation" / P["p_slug"] / P["r_slug"]
    log_dir.mkdir(parents=True, exist_ok=True)
    header = sbatch_header_block(args, f"ANNOT_{P['p_slug']}_{P['code'].lower()}",
                                 f"{log_dir}/annot.%A_%a.out", f"{log_dir}/annot.%A_%a.err", n_tasks)
    python = PROJECT_ROOT / "envs" / "pipeline" / "bin" / "python"
    command = [str(python), str(Path(__file__).resolve()), "--phenotype", P["phenotype"], "--ancestry", P["code"],
               "--worker-manifest", str(manifest_file)] + passthrough_options(args)
    bash_file = PROJECT_ROOT / f"Step08_Annotate_{P['p_slug']}_{P['r_slug']}.sh"
    bash_file.write_text(f"""#!/bin/bash --login
{header}

set -euo pipefail
cd {shlex.quote(str(PROJECT_ROOT))}
export PATH={shlex.quote(str(PROJECT_ROOT / 'envs' / 'pipeline' / 'bin'))}:"$PATH"
{'# SCREEN_API_KEY: exported in the environment or read from config/api_keys.env' if 'screen' in (args.online or '') else ''}

{' '.join(shlex.quote(x) for x in command)}
""", encoding="utf-8", newline="\n")
    bash_file.chmod(0o755)

    for accession in ready:
        print(f"READY    {accession}")
    base = (f"envs/pipeline/bin/python {Path(__file__).name} --phenotype {shlex.quote(P['phenotype'])} "
            f"--ancestry {P['code']}{panel_option(P['panel'])}")
    print(f"\nStudies: {len(ready)}  (array 1-{n_tasks})\nManifest: {manifest_file}\n\nRUN:\n  sbatch {bash_file.name}")
    print(f"THEN:\n  {base} --status\n  {base} --stats")


def worker_mode(args) -> None:
    P = project_paths(args)
    task = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    manifest = pd.read_csv(args.worker_manifest, sep="\t", dtype=str)
    selected = manifest[pd.to_numeric(manifest["TASK_ID"]) == task]
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", args.cpus))
    failed = []
    for accession in selected["STUDY_ACCESSION"]:
        try:
            annotate_study(P, accession, args, threads)
        except Exception as error:
            (P["out"] / accession).mkdir(parents=True, exist_ok=True)
            (P["out"] / accession / "ANNOTATION_FAILED.txt").write_text(
                f"{type(error).__name__}: {error}\n\n{traceback.format_exc()}", encoding="utf-8")
            print(f"[FAILED] {accession}: {error}", file=sys.stderr)
            failed.append(accession)
    if failed:
        raise SystemExit(f"Annotation failed for: {', '.join(failed)}")


def index_mode(args) -> None:
    P = project_paths(args)
    accession = study_for_index(P["p_slug"], P["a_slug"], args.index)
    ready, skipped = finemapped_studies(P)
    if accession not in ready:
        reason = dict(skipped).get(accession, "Step06 not run for this study")
        raise SystemExit(f"Study {args.index} ({accession}) cannot be annotated: {reason}. Run:\n"
                         f"  envs/pipeline/bin/python Step06_FineMap_GWAS_Loci.py --phenotype "
                         f"{shlex.quote(P['phenotype'])} --ancestry {P['code']} --index {args.index}"
                         f"{panel_option(P['panel'])}")
    try:
        annotate_study(P, accession, args, int(os.environ.get("SLURM_CPUS_PER_TASK", args.cpus)))
    except Exception as error:
        (P["out"] / accession).mkdir(parents=True, exist_ok=True)
        (P["out"] / accession / "ANNOTATION_FAILED.txt").write_text(
            f"{type(error).__name__}: {error}\n\n{traceback.format_exc()}", encoding="utf-8")
        raise


def status_mode(args) -> None:
    P = project_paths(args)
    ready, skipped = finemapped_studies(P)
    banner(f"STEP 08 STATUS - {P['phenotype']} / {P['label']} ({P['code']}) - LD panel {P['panel']}")
    counts: dict[str, int] = {}
    print(f"{'STUDY':<15}{'STATE':<10}{'VARIANTS':>10}{'LOCI':>6}  DETAIL")
    for accession in ready + [a for a, _ in skipped]:
        d = P["out"] / accession
        provenance = d / "annotation_provenance.json"
        if accession not in ready:
            state, detail, n, loci = "SKIPPED", dict(skipped)[accession], "", ""
        elif provenance.exists():
            info = json.loads(provenance.read_text(encoding="utf-8"))
            state, n, loci = "DONE", f"{info.get('n_variants', 0):,}", info.get("n_loci", "")
            missing = [k for k, v in info.items() if isinstance(v, dict) and v.get("status") == "NOT_AVAILABLE"]
            detail = f"scope {info.get('scope')}" + (f"; not available: {', '.join(missing)}" if missing else "")
        elif (d / "ANNOTATION_FAILED.txt").exists():
            state, n, loci = "FAILED", "", ""
            detail = (d / "ANNOTATION_FAILED.txt").read_text(encoding="utf-8", errors="replace").splitlines()[0][:80]
        else:
            state, detail, n, loci = "PENDING", "", "", ""
        counts[state] = counts.get(state, 0) + 1
        print(f"{accession:<15}{state:<10}{n!s:>10}{loci!s:>6}  {detail}")
    print("\nStudies : " + "  ".join(f"{k}={v}" for k, v in sorted(counts.items())))


def stats_mode(args) -> None:
    P = project_paths(args)
    banner(f"STEP 08 STATISTICS - {P['phenotype']} / {P['label']} ({P['code']}) - LD panel {P['panel']}")
    rows = []
    for d in sorted(p for p in P["out"].iterdir() if p.is_dir()) if P["out"].exists() else []:
        final = d / f"{d.name}_annotated_variants.tsv.gz"
        if not final.exists():
            continue
        m = pd.read_csv(final, sep="\t", low_memory=False)
        cs = m[m["IN_95_CREDIBLE_SET"].astype(str) == "True"]
        loci = pd.read_csv(d / f"{d.name}_locus_annotation_summary.tsv", sep="\t").fillna("")
        b = lambda frame, col: int((frame.get(col, pd.Series(dtype=str)).astype(str) == "True").sum())  # noqa: E731
        category = cs.get("VEP_FUNCTIONAL_CATEGORY", pd.Series(dtype=str)).fillna("NOT_ANNOTATED")
        pip_by_class = (cs["PIP"].groupby(category).sum() / max(cs["PIP"].sum(), 1e-12)).round(3)
        row = {"STUDY_ACCESSION": d.name, "N_VARIANTS": len(m), "N_LOCI": len(loci), "N_CS_VARIANTS": len(cs),
               "CS_CODING": b(cs, "IS_CODING"), "CS_LOF": b(cs, "IS_LOF"), "CS_MISSENSE": b(cs, "IS_MISSENSE"),
               "CS_SPLICE": b(cs, "IS_SPLICE"), "CS_AM_LIKELY_PATHOGENIC": b(cs, "AM_LIKELY_PATHOGENIC"),
               "CS_CADD_GE_20": b(cs, "CADD_GE_20"), "CS_CLINVAR_PATHOGENIC": b(cs, "CLINVAR_PATHOGENIC"),
               "CS_IN_CCRE": b(cs, "IN_CCRE"), "CS_CONSERVED": b(cs, "PHYLOP_CONSERVED"),
               **{f"CS_PIP_FRACTION_{k}": v for k, v in pip_by_class.items()}}
        rows.append(row)
        print(f"\n{d.name}: {len(m):,} variants in {len(loci)} loci; {len(cs):,} in 95% credible sets")
        print(f"   credible-set variants: coding {row['CS_CODING']} (LoF {row['CS_LOF']}, missense "
              f"{row['CS_MISSENSE']}), splice {row['CS_SPLICE']}, AlphaMissense likely pathogenic "
              f"{row['CS_AM_LIKELY_PATHOGENIC']}, CADD>=20 {row['CS_CADD_GE_20']}, ClinVar P/LP "
              f"{row['CS_CLINVAR_PATHOGENIC']}, in cCRE {row['CS_IN_CCRE']}, conserved {row['CS_CONSERVED']}")
        print("   PIP mass in credible sets by class: "
              + ", ".join(f"{k} {v:.0%}" for k, v in pip_by_class.sort_values(ascending=False).items()))
        top = loci.sort_values("TOP_PIP", ascending=False).head(8)
        for _, r in top.iterrows():
            print(f"   {r['LOCUS_ID']}: {r['TOP_VARIANT']} {r.get('TOP_RSID', '') or ''} PIP {r['TOP_PIP']} "
                  f"{r.get('TOP_GENE', '') or ''} {r.get('TOP_CONSEQUENCE', '') or ''} | {r.get('TOP_EVIDENCE', '') or ''}")
    if not rows:
        print("No annotated studies yet (see --status).")
        return
    output = P["out"] / "annotation_stats.tsv"
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)
    print(f"\nTable: {output}")


def test_variant_mode(args) -> None:
    match = re.fullmatch(r"(?:chr)?(\d{1,2})[:_-](\d+)[:_-]([ACGTN]+)[:_-]([ACGTN]+)", args.test_variant.strip(), re.I)
    if not match:
        raise SystemExit("--test-variant must look like 4:89822321:G:A (GRCh38)")
    chrom, pos, ref, alt = int(match.group(1)), int(match.group(2)), match.group(3).upper(), match.group(4).upper()
    universe = pd.DataFrame([{"VARIANT_KEY": variant_key(chrom, pos, ref, alt), "PRIMARY_RSID": "", "CHR": chrom,
                              "REFERENCE_POS": pos, "REFERENCE_REF": ref, "REFERENCE_ALT": alt}])
    work = PROJECT_ROOT / "08_annotation" / "_test_variant" / slugify(universe.at[0, "VARIANT_KEY"])
    work.mkdir(parents=True, exist_ok=True)
    print_resources()
    banner(f"STEP 08 TEST VARIANT {universe.at[0, 'VARIANT_KEY']}")
    table, _ = annotate_universe(universe, work, args.cpus, args)
    for key, value in table.iloc[0].items():
        if str(value) not in ("", "nan", "None"):
            print(f"  {key:<40} {str(value)[:110]}")


def arguments(argv=None):
    parser = argparse.ArgumentParser(description="Local annotation of Step06 fine-mapped variants.")
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--index", type=int, default=None,
                       help="Study N of the Step02 download manifest: annotate it on this machine.")
    modes.add_argument("--status", action="store_true", help="READ-ONLY: done / failed / pending per study.")
    modes.add_argument("--stats", action="store_true", help="Credible-set functional summary per study.")
    modes.add_argument("--check-resources", action="store_true",
                       help="Verify every annotation file (+ index), the VEP cache, tools, Python modules, the "
                            "Step06 input of this panel and the API keys; exit 1 if anything is missing. "
                            "With --online, also test the web APIs live.")
    modes.add_argument("--test-variant", default=None, help="Annotate one GRCh38 variant, e.g. 4:89822321:G:A.")
    modes.add_argument("--worker-manifest", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--reference-panel", default=DEFAULT_PANEL,
                        help="LD reference panel of Steps 04-06: 1000G (default), 1000G30x or a custom name. "
                             "Other panels read 06_finemapping/<p>/<a>_<PANEL>/ and write 08_annotation/<p>/<a>_<PANEL>/.")
    parser.add_argument("--scope", choices=["all", "cs", "pip"], default="all",
                        help="all fine-mapped variants (default), 95%% credible sets only, or PIP >= --min-pip.")
    parser.add_argument("--min-pip", type=float, default=0.01)
    parser.add_argument("--online", default="",
                        help="Optional web lookups for credible-set variants: screen,regulomedb.")
    parser.add_argument("--screen-active-z", type=float, default=1.64)
    parser.add_argument("--no-vep", action="store_true", help="Skip VEP (CADD / AlphaMissense / ClinVar / ... only).")
    parser.add_argument("--vep-gtf", default=None, help=argparse.SUPPRESS)  # testing: GTF gene models, no cache
    parser.add_argument("--vep-path", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--force", action="store_true", help="Re-annotate studies that are already done.")
    parser.add_argument("--keep-work", action="store_true", help="Keep the VEP input / output / log files.")
    slurm = parser.add_argument_group("SLURM (default: the #SBATCH header of HPC_CPU.sh)")
    slurm.add_argument("--template", default="HPC_CPU.sh")
    slurm.add_argument("--partition", default=None)
    slurm.add_argument("--qos", default=None)
    slurm.add_argument("--account", default=None)
    slurm.add_argument("--cpus", type=int, default=None)
    slurm.add_argument("--mem", "--memory", dest="mem", default=None)
    slurm.add_argument("--time", default=None)
    slurm.add_argument("--max-parallel", type=int, default=0)
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.reference_panel):
        parser.error(f"--reference-panel {args.reference_panel!r}: use letters, digits, '_', '-' or '.'")
    return resolve_slurm(args)


def main() -> None:
    args = arguments()
    if args.check_resources:
        check_resources_mode(args)
    elif args.test_variant:
        test_variant_mode(args)
    elif args.status:
        status_mode(args)
    elif args.stats:
        stats_mode(args)
    elif args.worker_manifest:
        worker_mode(args)
    elif args.index is not None:
        index_mode(args)
    else:
        planner_mode(args)


if __name__ == "__main__":
    main()
