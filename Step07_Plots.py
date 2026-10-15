#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
GWAS2Mech STEP 07 - Publication figures for GWAS loci and fine-mapping
===============================================================================

    envs/pipeline/bin/python Step07_Plots.py --phenotype "Parkinson's disease" --ancestry EUR
    envs/pipeline/bin/python Step07_Plots.py --phenotype "Parkinson's disease" --ancestry EUR --formats pdf,svg,png
    envs/pipeline/bin/python Step07_Plots.py --phenotype "Parkinson's disease" --ancestry EUR --status

    # as a SLURM job (settings from HPC_CPU.sh)
    envs/pipeline/bin/python Step07_Plots.py --phenotype "Parkinson's disease" --ancestry EUR --slurm
    sbatch Step07_Plots_parkinson_s_disease_european.sh

Inputs (nothing is recomputed):
    Step02  02_summary_stats/<p>/<a>/qc/<GCST>/<GCST>_GRCh38_QC.tsv.gz, GWAS_study_metadata.tsv
    Step05  05_ld_clumping/<p>/<a>/<GCST>/lead_variants.tsv
    Step06  06_finemapping/<p>/<a>/<GCST>/ (loci, regional GWAS, SuSiE PIPs / credible
            sets, LD matrices, locus QC)
    genes   resources/gencode/*/gencode.v*.annotation.gtf.gz (Step00 download), canonical
            transcripts (Ensembl_canonical) for the gene tracks

Figures (07_plots/<p>/<a>/):
    <GCST>/manhattan             genome-wide association, independent leads labelled
                                 with nearest gene, fine-mapped loci shaded
    <GCST>/ideogram              chromosomes 1-22 with every locus (karyoploteR)
    <GCST>/loci/<L>_region       per locus, stacked panels on one genomic axis:
                                   association (-log10 P, coloured by r2 with the index
                                   variant), SuSiE PIP (95% credible sets), gene models,
                                   LD triangle (r2, 1000G panel)
    <GCST>/all_loci_multipage.pdf
    combined/ideogram_all_gwas   every locus of every GWAS, one band per GWAS
    combined/locus_by_gwas       loci x GWAS matrix: min P and max PIP (ComplexHeatmap)
    combined/shared_loci/<chr>_<pos>_<gene>
                                 loci found in >= 2 GWAS: one association panel per
                                 GWAS on a shared axis, PIP of every GWAS, genes
    combined/locus_summary.tsv   the table behind the matrix

R packages (installed by Step00_Installation.py): ggplot2, patchwork, ggrastr,
ggrepel, scales, svglite, ComplexHeatmap, karyoploteR, data.table. Dense point
layers are rasterised (ggrastr); text, lines and gene models stay vector.
Gene tracks use the pipeline's local GENCODE release (offline, same annotation
as every other step) rather than online Ensembl annotation.
===============================================================================
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path.cwd().resolve()

ANCESTRY_ALIASES = {
    "eur": ("EUR", "European"), "european": ("EUR", "European"),
    "afr": ("AFR", "African"), "african": ("AFR", "African"),
    "eas": ("EAS", "East Asian"), "east asian": ("EAS", "East Asian"),
    "sas": ("SAS", "South Asian"), "south asian": ("SAS", "South Asian"),
    "amr": ("AMR", "Hispanic or Latin American"), "hispanic": ("AMR", "Hispanic or Latin American"),
    "latino": ("AMR", "Hispanic or Latin American"),
    "hispanic or latin american": ("AMR", "Hispanic or Latin American"),
}

R_PACKAGES = ["ggplot2", "patchwork", "ggrastr", "ggrepel", "scales", "svglite",
              "ComplexHeatmap", "karyoploteR", "data.table"]


# =============================================================================
# HELPERS
# =============================================================================

def banner(text):
    print()
    print("=" * 88)
    print(text)
    print("=" * 88)


def normalize_text(value):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", str(value).lower())).strip()


def slugify(value):
    return normalize_text(value).replace(" ", "_")


def resolve_ancestry(value):
    key = normalize_text(value)
    if key not in ANCESTRY_ALIASES:
        raise SystemExit(f"Unsupported ancestry {value!r}: use EUR, AFR, EAS, SAS or AMR")
    return ANCESTRY_ALIASES[key]


def clean_chr(value):
    try:
        chrom = int(float(str(value).strip().lower().replace("chr", "")))
    except ValueError:
        return None
    return chrom if 1 <= chrom <= 22 else None


def read_tsv(path, **kwargs):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def project_paths(args):
    code, label = resolve_ancestry(args.ancestry)
    p_slug, a_slug = slugify(args.phenotype), slugify(label)
    out = PROJECT_ROOT / "07_plots" / p_slug / a_slug
    return {
        "phenotype": args.phenotype.strip(), "code": code, "label": label,
        "p_slug": p_slug, "a_slug": a_slug,
        "qc": PROJECT_ROOT / "02_summary_stats" / p_slug / a_slug / "qc",
        "metadata": PROJECT_ROOT / "01_gwas_catalog" / p_slug / a_slug / "GWAS_study_metadata.tsv",
        "clump": PROJECT_ROOT / "05_ld_clumping" / p_slug / a_slug,
        "finemap": PROJECT_ROOT / "06_finemapping" / p_slug / a_slug,
        "out": out,
        "data": out / "data",
    }


# =============================================================================
# GENE MODELS (local GENCODE; canonical transcripts)
# =============================================================================

def gene_models():
    """gene / transcript / exon rows of protein-coding + lncRNA genes, canonical
    transcripts only (Ensembl_canonical); cached next to the GTF."""
    gtfs = sorted((PROJECT_ROOT / "resources" / "gencode").glob("*/gencode.v*.annotation.gtf.gz"))
    if not gtfs:
        return pd.DataFrame(columns=["FEATURE", "CHR", "START", "END", "STRAND", "GENE", "TYPE"])
    gtf = gtfs[-1]
    cache = gtf.with_name(gtf.name.replace(".gtf.gz", ".step07_models.tsv.gz"))
    if cache.exists():
        return pd.read_csv(cache, sep="\t")

    print(f"Reading gene models from {gtf} (cached once)")
    rows = []
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "transcript", "exon"):
                continue
            chrom = clean_chr(f[0])
            if chrom is None:
                continue
            attributes = f[8]
            gtype = re.search(r'gene_type "([^"]+)"', attributes)
            name = re.search(r'gene_name "([^"]+)"', attributes)
            if not gtype or not name or gtype.group(1) not in ("protein_coding", "lncRNA"):
                continue
            if f[2] != "gene" and 'tag "Ensembl_canonical"' not in attributes:
                continue
            rows.append((f[2], chrom, int(f[3]), int(f[4]), f[6], name.group(1), gtype.group(1)))

    table = pd.DataFrame(rows, columns=["FEATURE", "CHR", "START", "END", "STRAND", "GENE", "TYPE"])
    table.to_csv(cache, sep="\t", index=False, compression="gzip")
    return table


def nearest_gene(genes, chrom, pos):
    """Protein-coding gene containing pos, else the closest one ('' if none)."""
    if genes.empty or pos is None or pd.isna(pos):
        return ""
    g = genes[(genes["FEATURE"] == "gene") & (genes["CHR"] == int(chrom)) & (genes["TYPE"] == "protein_coding")]
    if g.empty:
        return ""
    distance = np.maximum(0, np.maximum(g["START"].to_numpy() - pos, pos - g["END"].to_numpy()))
    i = int(np.argmin(distance))
    return g.iloc[i]["GENE"] if distance[i] <= 1_000_000 else ""


# =============================================================================
# DATA PREPARATION (one compact folder per study for the R script)
# =============================================================================

def study_for_index(phenotype_slug, ancestry_slug, index):
    """Accession of study N in the Step02 download manifest (the ARRAY_ID the
    Step02 planner prints) - the same N in Step02, Step05, Step06 and Step07."""
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


def lead_table(P, accession):
    leads = read_tsv(P["clump"] / accession / "lead_variants.tsv")
    if leads.empty:
        return pd.DataFrame(columns=["CHR", "POS", "ID", "P"])
    chr_col = next(c for c in ("#CHROM", "CHROM", "CHR") if c in leads.columns)
    out = pd.DataFrame({"CHR": leads[chr_col].map(clean_chr), "POS": pd.to_numeric(leads["POS"], errors="coerce"),
                        "ID": leads["ID"].astype(str), "P": pd.to_numeric(leads["P"], errors="coerce")})
    return out.dropna(subset=["CHR", "POS"]).astype({"CHR": int, "POS": int})


def manhattan_input(P, accession, out_file):
    """CHR POS P: every variant with P < 1e-3 plus a fixed 10% of the rest; lambda GC on all."""
    qc_file = P["qc"] / accession / f"{accession}_GRCh38_QC.tsv.gz"
    if not qc_file.exists():
        return None
    from scipy.stats import chi2

    rng = np.random.default_rng(7)
    kept, all_p = [], []
    for chunk in pd.read_csv(qc_file, sep="\t", usecols=lambda c: c in ("CHR", "POS", "P"),
                             chunksize=2_000_000, low_memory=False):
        chunk["CHR"] = chunk["CHR"].map(clean_chr)
        chunk = chunk.dropna(subset=["CHR", "POS", "P"])
        chunk = chunk[chunk["P"].between(0, 1)]
        all_p.append(chunk["P"].to_numpy(dtype=float))
        keep = (chunk["P"] < 1e-3) | (rng.random(len(chunk)) < 0.10)
        kept.append(chunk.loc[keep, ["CHR", "POS", "P"]].astype({"CHR": int, "POS": int}))
    data = pd.concat(kept, ignore_index=True)
    data.to_csv(out_file, sep="\t", index=False)
    p = np.clip(np.concatenate(all_p), 1e-300, 1.0)
    return float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)), int(len(p))


def ld_row(locus_dir, code, ids, variant_id):
    ld_file = locus_dir / f"{code}_LD.unphased.vcor1.bin"
    if not ld_file.exists() or variant_id not in ids:
        return None
    matrix = np.memmap(ld_file, dtype="<f4", mode="r", shape=(len(ids), len(ids)))
    return np.asarray(matrix[ids.index(variant_id)], dtype=float)


def prepare_locus(P, accession, locus, genes, out_dir, ld_max):
    """variants.tsv (+ ld_variants.tsv / ld_pairs.tsv) for one locus; returns summary dict."""
    study_dir = P["finemap"] / accession
    locus_id = str(locus["LOCUS_ID"])
    locus_dir = study_dir / "loci" / locus_id
    regional = read_tsv(study_dir / "regional_gwas" / f"{locus_id}.tsv.gz")
    fm = read_tsv(locus_dir / "susie" / "susie_variants.tsv")
    if regional.empty and fm.empty:
        return None

    keys = ["POS", "EA", "NEA"]
    base = regional[[c for c in ("POS", "EA", "NEA", "P", "SNPID") if c in regional.columns]].copy()
    if not fm.empty:
        extra = fm[[c for c in keys + ["REFERENCE_ID", "PIP", "PIP_L1", "CREDIBLE_SET", "IN_95_CREDIBLE_SET"]
                    if c in fm.columns]]
        variants = base.merge(extra, on=keys, how="outer") if not base.empty else fm.copy()
        if "P" not in variants or variants["P"].isna().any():
            variants = variants.merge(fm[keys + ["P", "SNPID"]], on=keys, how="left", suffixes=("", "_FM"))
            variants["P"] = variants["P"].fillna(variants.get("P_FM"))
            variants["SNPID"] = variants["SNPID"].fillna(variants.get("SNPID_FM"))
    else:
        variants = base.assign(REFERENCE_ID=np.nan, PIP=np.nan, PIP_L1=np.nan, CREDIBLE_SET="")

    variants["P"] = pd.to_numeric(variants["P"], errors="coerce")
    variants = variants.dropna(subset=["POS", "P"]).copy()
    variants["POS"] = variants["POS"].astype(int)
    variants["IN_LD"] = variants["REFERENCE_ID"].notna()
    variants["CREDIBLE_SET"] = variants.get("CREDIBLE_SET", "").fillna("").astype(str)

    # r2 with the index variant (smallest P among variants in the LD matrix)
    vars_file = locus_dir / f"{P['code']}_LD.unphased.vcor1.bin.vars"
    ids = [x.strip() for x in vars_file.read_text().splitlines() if x.strip()] if vars_file.exists() else []
    variants["R2"] = np.nan
    variants["INDEX"] = False
    in_ld = variants[variants["IN_LD"]]
    if not in_ld.empty:
        index_row = in_ld.loc[in_ld["P"].idxmin()]
        variants.loc[index_row.name, "INDEX"] = True
        r = ld_row(locus_dir, P["code"], ids, str(index_row["REFERENCE_ID"]))
        if r is not None:
            variants["R2"] = variants["REFERENCE_ID"].map(dict(zip(ids, r ** 2)))

    variants["LABEL"] = variants["SNPID"].fillna(variants.get("REFERENCE_ID")).astype(str)
    out_dir.mkdir(parents=True, exist_ok=True)
    columns = ["POS", "P", "LABEL", "PIP", "PIP_L1", "CREDIBLE_SET", "IN_LD", "R2", "INDEX"]
    variants[[c for c in columns if c in variants.columns]].to_csv(out_dir / "variants.tsv", sep="\t", index=False)

    # LD triangle: credible-set members + strongest associations (<= ld_max), by position
    if ids:
        cs = variants[variants["CREDIBLE_SET"] != ""]
        top = variants[variants["IN_LD"]].sort_values("P").head(ld_max)
        chosen = pd.concat([cs[cs["IN_LD"]], top]).drop_duplicates("REFERENCE_ID").head(ld_max)
        chosen = chosen[chosen["REFERENCE_ID"].isin(ids)].sort_values("POS").reset_index(drop=True)
        if len(chosen) >= 2:
            index = [ids.index(x) for x in chosen["REFERENCE_ID"]]
            matrix = np.memmap(locus_dir / f"{P['code']}_LD.unphased.vcor1.bin", dtype="<f4", mode="r",
                               shape=(len(ids), len(ids)))
            sub = np.asarray(matrix[np.ix_(index, index)], dtype=float) ** 2
            i, j = np.triu_indices(len(index), k=1)
            pd.DataFrame({"I": i + 1, "J": j + 1, "R2": np.round(sub[i, j], 4)}).to_csv(
                out_dir / "ld_pairs.tsv", sep="\t", index=False)
            pd.DataFrame({"I": np.arange(1, len(chosen) + 1), "POS": chosen["POS"].to_numpy()}).to_csv(
                out_dir / "ld_variants.tsv", sep="\t", index=False)

    best = variants.sort_values(["PIP", "P"], ascending=[False, True]).iloc[0] \
        if variants["PIP"].notna().any() else variants.sort_values("P").iloc[0]
    return {
        "TOP_POS": int(best["POS"]),
        "TOP_LABEL": str(best["LABEL"]),
        "MIN_P": float(variants["P"].min()),
        "MAX_PIP": float(variants["PIP"].max()) if variants["PIP"].notna().any() else np.nan,
        "NEAREST_GENE": nearest_gene(genes, int(locus["CHR"]), int(best["POS"])),
    }


def merge_windows(frame):
    """Overlapping (CHR, START, END) windows -> combined locus ids."""
    frame = frame.sort_values(["CHR", "START"]).reset_index(drop=True)
    ids, current, end, chrom = [], 0, -1, None
    for _, row in frame.iterrows():
        if row["CHR"] != chrom or row["START"] > end:
            current += 1
            chrom, end = row["CHR"], row["END"]
        else:
            end = max(end, row["END"])
        ids.append(f"C{current:03d}")
    frame["COMBINED_ID"] = ids
    return frame


def prepare(P, args):
    banner(f"STEP 07 - PREPARING FIGURE DATA - {P['phenotype']} / {P['label']} ({P['code']})")
    data = P["data"]
    if data.exists():
        shutil.rmtree(data)
    data.mkdir(parents=True)

    meta = read_tsv(P["metadata"], dtype=str).fillna("")
    meta = {r["GCST_ID"]: r for _, r in meta.iterrows()} if not meta.empty else {}
    manifest = read_tsv(P["finemap"] / "finemapping_manifest.tsv", dtype=str)
    finemapped = {r["STUDY_ACCESSION"]: r for _, r in manifest.iterrows()} if not manifest.empty else {}

    studies = sorted(p.name for p in P["clump"].iterdir() if (p / "lead_variants.tsv").exists()) \
        if P["clump"].exists() else []
    if args.index is not None:
        wanted = study_for_index(P["p_slug"], P["a_slug"], args.index)
        if wanted not in studies:
            raise SystemExit(f"Study {args.index} ({wanted}) has no Step05 result yet. Run:\n"
                             f"  envs/pipeline/bin/python Step05_Clump_GWAS_Loci.py --phenotype "
                             f"{shlex.quote(P['phenotype'])} --ancestry {P['code']} --index {args.index}")
        studies = [wanted]
    if not studies:
        raise SystemExit(f"No Step05 results under {P['clump']} - run Step05 (and Step06) first.")

    genes = gene_models()
    if genes.empty:
        print("No GENCODE GTF under resources/gencode/ - figures have no gene track / gene names "
              "(Step00_Download_Resources.py --groups genome).")

    study_rows, locus_rows, windows = [], [], []
    for accession in studies:
        out = data / accession
        out.mkdir(parents=True)
        result = manhattan_input(P, accession, out / "manhattan.tsv")
        if result is None:
            print(f"{accession}: no Step02 QC file - skipped")
            continue
        lambda_gc, n_variants = result

        leads = lead_table(P, accession)
        leads["GENE"] = [nearest_gene(genes, c, p) for c, p in zip(leads["CHR"], leads["POS"])]
        leads["LABEL"] = np.where(leads["GENE"] != "", leads["GENE"] + " " + leads["ID"], leads["ID"])
        leads.to_csv(out / "leads.tsv", sep="\t", index=False)

        study_dir = P["finemap"] / accession
        summary = read_tsv(study_dir / f"{accession}_finemapping_summary.tsv")
        loci = read_tsv(study_dir / "loci_definition.tsv")
        summary_json = study_dir / "finemapping_summary.json"
        done = not loci.empty and summary_json.exists()
        if done and accession not in finemapped:   # fine-mapped with Step06 --index
            finemapped[accession] = json.loads(summary_json.read_text(encoding="utf-8"))

        if done:
            if not summary.empty:
                keep = [c for c in ("LOCUS_ID", "N_CREDIBLE_SETS", "LOCUS_QC", "LOCUS_QC_NOTES", "ESTIMATE_S_RSS")
                        if c in summary.columns]
                loci = loci.merge(summary[keep], on="LOCUS_ID", how="left")
            n_loci = 0
            for _, locus in loci.iterrows():
                if args.max_loci and n_loci >= args.max_loci:
                    break
                info = prepare_locus(P, accession, locus, genes, out / "loci" / str(locus["LOCUS_ID"]), args.ld_max)
                if info is None:
                    continue
                n_loci += 1
                row = {"GCST": accession, "LOCUS_ID": locus["LOCUS_ID"], "CHR": int(locus["CHR"]),
                       "START": int(locus["LOCUS_START"]), "END": int(locus["LOCUS_END"]),
                       "N_CS": locus.get("N_CREDIBLE_SETS", ""), "LOCUS_QC": locus.get("LOCUS_QC", ""),
                       "QC_NOTES": locus.get("LOCUS_QC_NOTES", ""), "S_RSS": locus.get("ESTIMATE_S_RSS", ""),
                       **info}
                locus_rows.append(row)
                windows.append(row)
        else:
            # Not fine-mapped (excluded / pending): loci = Step05 leads +- 1 Mb, no PIP.
            for chrom, group in leads.groupby("CHR"):
                for _, lead in group.iterrows():
                    windows.append({"GCST": accession, "LOCUS_ID": f"lead_{lead['ID']}", "CHR": int(chrom),
                                    "START": max(1, int(lead["POS"]) - 1_000_000),
                                    "END": int(lead["POS"]) + 1_000_000, "TOP_POS": int(lead["POS"]),
                                    "TOP_LABEL": lead["ID"], "MIN_P": float(lead["P"]), "MAX_PIP": np.nan,
                                    "NEAREST_GENE": lead["GENE"], "N_CS": "", "LOCUS_QC": "", "QC_NOTES": "",
                                    "S_RSS": ""})

        m = meta.get(accession, {})
        f = finemapped.get(accession, {})
        study_rows.append({
            "GCST": accession,
            "PHENOTYPE_TYPE": m.get("PHENOTYPE_TYPE", "") if len(m) else "",
            "N_USED": f.get("GWAS_N", "") if len(f) else "",
            "N_SOURCE": f.get("GWAS_N_SOURCE", "") if len(f) else "not fine-mapped",
            "LAMBDA_GC": round(lambda_gc, 3),
            "N_VARIANTS": n_variants,
            "N_LEADS": len(leads),
            "FINEMAPPED": bool(done),
        })
        print(f"{accession}: {n_variants:,} variants, lambda {lambda_gc:.3f}, {len(leads)} leads, "
              f"{'fine-mapped' if done else 'NOT fine-mapped (Manhattan / ideogram only)'}")

    pd.DataFrame(study_rows).to_csv(data / "studies.tsv", sep="\t", index=False)
    all_windows = pd.DataFrame(windows)
    if all_windows.empty:
        all_windows = pd.DataFrame(columns=["GCST", "LOCUS_ID", "CHR", "START", "END", "TOP_POS", "TOP_LABEL",
                                            "MIN_P", "MAX_PIP", "NEAREST_GENE"])
    all_windows = merge_windows(all_windows) if not all_windows.empty else all_windows.assign(COMBINED_ID=[])
    all_windows.to_csv(data / "loci.tsv", sep="\t", index=False)

    # combined loci: name by nearest gene of the strongest signal
    combined = []
    for cid, group in all_windows.groupby("COMBINED_ID"):
        best = group.sort_values("MIN_P").iloc[0]
        mid = (group["START"].min() + group["END"].max()) / 2
        gene = best["NEAREST_GENE"] if isinstance(best["NEAREST_GENE"], str) else ""
        name = f"chr{int(best['CHR'])}:{mid / 1e6:.1f}Mb" + (f" {gene}" if gene else "")
        combined.append({"COMBINED_ID": cid, "CHR": int(best["CHR"]), "START": int(group["START"].min()),
                         "END": int(group["END"].max()), "NAME": name, "N_GWAS": group["GCST"].nunique(),
                         "FILE_STEM": f"chr{int(best['CHR'])}_{int(mid)}" + (f"_{gene}" if gene else "")})
    pd.DataFrame(combined, columns=["COMBINED_ID", "CHR", "START", "END", "NAME", "N_GWAS", "FILE_STEM"]).to_csv(
        data / "combined_loci.tsv", sep="\t", index=False)

    region_genes = genes
    if not genes.empty and not all_windows.empty:
        keep = np.zeros(len(genes), dtype=bool)
        for _, w in all_windows.iterrows():
            keep |= (genes["CHR"].to_numpy() == w["CHR"]) & (genes["END"].to_numpy() >= w["START"]) & \
                    (genes["START"].to_numpy() <= w["END"])
        region_genes = genes[keep]
    region_genes.to_csv(data / "genes.tsv", sep="\t", index=False)

    lead_genes = pd.DataFrame(columns=["GCST", "CHR", "POS", "LABEL"])
    for accession in studies:
        f = data / accession / "leads.tsv"
        if f.exists():
            lg = read_tsv(f)
            if not lg.empty:
                lead_genes = pd.concat([lead_genes, lg.assign(GCST=accession)[["GCST", "CHR", "POS", "LABEL"]]])
    lead_genes.to_csv(data / "all_leads.tsv", sep="\t", index=False)
    return study_rows, all_windows


# =============================================================================
# R FIGURES
# =============================================================================

R_CODE = r'''options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]; out_dir <- args[2]
formats <- strsplit(args[3], ",")[[1]]; dpi <- as.integer(args[4]); title_prefix <- args[5]
ld_window_bp <- as.numeric(args[6])

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggrepel); library(data.table); library(scales)
})
has <- function(p) requireNamespace(p, quietly = TRUE)
HAS_RASTR <- has("ggrastr"); HAS_KARYO <- has("karyoploteR"); HAS_CHM <- has("ComplexHeatmap"); HAS_SVG <- has("svglite")

GRCH38 <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
            145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
            101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468)
OFFSET <- c(0, cumsum(as.numeric(GRCH38)))[1:22]
R2_COLOURS <- c(">=0.8" = "#d7191c", "0.6-0.8" = "#fdae61", "0.4-0.6" = "#4dac26",
                "0.2-0.4" = "#8fbfe0", "<0.2" = "#2c3e91", "not in LD panel" = "#cfcfcf")
CS_COLOURS <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#17becf")
GWAS_COLOURS <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d",
                  "#1f78b4", "#b2df8a", "#fb9a99", "#cab2d6", "#666666")

theme_pub <- theme_classic(base_size = 8) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.5, colour = "grey30"),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold", hjust = 0))
raster <- function(layer) if (HAS_RASTR) ggrastr::rasterise(layer, dpi = 300) else layer
mb <- function(x) sprintf("%.2f", x / 1e6)
readf <- function(path) if (file.exists(path) && file.info(path)$size > 0) fread(path) else data.table()

save_gg <- function(p, base, w, h) {
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  for (f in formats) {
    file <- paste0(base, ".", f)
    if (f == "svg" && !HAS_SVG) next
    if (f == "svg") ggsave(file, p, width = w, height = h, units = "in", device = svglite::svglite)
    else if (f == "png") ggsave(file, p, width = w, height = h, units = "in", dpi = dpi)
    else ggsave(file, p, width = w, height = h, units = "in", device = f)
  }
  invisible(base)
}

save_base <- function(base, w, h, draw) {
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  for (f in formats) {
    file <- paste0(base, ".", f)
    if (f == "pdf") pdf(file, width = w, height = h, useDingbats = FALSE)
    else if (f == "svg") { if (!HAS_SVG) next; svglite::svglite(file, width = w, height = h) }
    else if (f == "png") png(file, width = w, height = h, units = "in", res = dpi)
    else next
    tryCatch(draw(), finally = dev.off())
  }
}

# ---------------------------------------------------------------------------- data
studies <- readf(file.path(data_dir, "studies.tsv"))
loci_all <- readf(file.path(data_dir, "loci.tsv"))
combined <- readf(file.path(data_dir, "combined_loci.tsv"))
genes <- readf(file.path(data_dir, "genes.tsv"))
all_leads <- readf(file.path(data_dir, "all_leads.tsv"))
gwas_col <- setNames(rep(GWAS_COLOURS, length.out = nrow(studies)), studies$GCST)
cat("Studies:", nrow(studies), " loci:", nrow(loci_all), " combined loci:", nrow(combined), "\n")

# ---------------------------------------------------------------------------- Manhattan
plot_manhattan <- function(gcst, info) {
  m <- readf(file.path(data_dir, gcst, "manhattan.tsv"))
  leads <- readf(file.path(data_dir, gcst, "leads.tsv"))
  loci <- loci_all[GCST == gcst]
  m[, X := POS + OFFSET[CHR]][, LOGP := -log10(pmax(P, 1e-300))]
  ymax <- max(c(m$LOGP, 10)) * 1.12
  p <- ggplot(m, aes(X, LOGP))
  if (nrow(loci)) p <- p + geom_rect(data = loci, inherit.aes = FALSE, fill = "#fdebd0",
                                     aes(xmin = START + OFFSET[CHR], xmax = END + OFFSET[CHR], ymin = 0, ymax = Inf))
  p <- p + raster(geom_point(aes(colour = factor(CHR %% 2)), size = 0.35, show.legend = FALSE)) +
    scale_colour_manual(values = c("0" = "#2b6c9e", "1" = "#8fb8d9")) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "#c0392b", linewidth = 0.3) +
    geom_hline(yintercept = -log10(1e-5), linetype = "dotted", colour = "grey55", linewidth = 0.25)
  if (nrow(leads)) {
    leads[, X := POS + OFFSET[CHR]][, LOGP := -log10(pmax(P, 1e-300))]
    p <- p + geom_point(data = leads, shape = 21, fill = "#c0392b", colour = "black", size = 1.4, stroke = 0.2) +
      geom_text_repel(data = head(leads[order(P)], 25), aes(label = LABEL), size = 2, fontface = "italic",
                      min.segment.length = 0, segment.size = 0.2, max.overlaps = 50, box.padding = 0.25)
  }
  p + scale_x_continuous(breaks = OFFSET + GRCH38 / 2, labels = 1:22, expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    labs(x = "Chromosome (GRCh38)", y = expression(-log[10](italic(P))),
         title = sprintf("%s  %s", title_prefix, gcst),
         subtitle = sprintf("%s | N used %s | %s variants | lambda GC %.3f | %d independent leads | shaded: fine-mapping windows",
                            info$PHENOTYPE_TYPE, info$N_USED, format(info$N_VARIANTS, big.mark = ","),
                            info$LAMBDA_GC, info$N_LEADS)) +
    theme_pub + theme(axis.text.x = element_text(size = 6))
}

# ---------------------------------------------------------------------------- ideograms
genome_gr <- function() GenomicRanges::GRanges(paste0("chr", 1:22), IRanges::IRanges(1, GRCH38))

draw_ideogram <- function(tracks, title, labels) {
  groups <- unique(tracks$GROUP)
  if (HAS_KARYO) {
    pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
    pp$data1height <- max(60, 22 * length(groups)); pp$ideogramheight <- 12
    kp <- tryCatch(karyoploteR::plotKaryotype(genome = "hg38", chromosomes = paste0("chr", 1:22),
                                              plot.type = 1, plot.params = pp, main = title, cex = 0.6),
                   error = function(e) karyoploteR::plotKaryotype(genome = genome_gr(), plot.type = 1,
                                              plot.params = pp, chromosomes = "all", main = title, cex = 0.6))
    # locus bands in the lower 45% of each data panel (one band per GWAS),
    # labels in the upper half, so nothing reaches the chromosome above
    n <- length(groups)
    for (i in seq_along(groups)) {
      t <- tracks[GROUP == groups[i]]
      karyoploteR::kpRect(kp, chr = paste0("chr", t$CHR), x0 = t$START, x1 = t$END, y0 = 0, y1 = 1,
                          r0 = 0.45 * (i - 1) / n, r1 = 0.45 * i / n - 0.04 / n, col = t$COLOUR, border = NA)
    }
    if (nrow(labels)) karyoploteR::kpPlotMarkers(kp, chr = paste0("chr", labels$CHR), x = labels$POS,
                                                 labels = labels$LABEL, cex = 0.42, r0 = 0.5, r1 = 0.98,
                                                 marker.parts = c(0.15, 0.75, 0.1), line.color = "grey50",
                                                 ignore.chromosome.ends = TRUE)
    if (n > 1) legend("bottomright", legend = groups, fill = unique(tracks$COLOUR), cex = 0.55, bty = "n")
  } else {
    plot.new(); title(main = paste(title, "(karyoploteR not installed)"))
  }
}

# ---------------------------------------------------------------------------- gene track
gene_track <- function(chr, start, end) {
  g <- genes[CHR == chr & END >= start & START <= end]
  if (!nrow(g)) return(list(plot = ggplot() + annotate("text", x = (start + end) / 2, y = 0, size = 2,
                                                       label = "no genes / no GENCODE annotation") +
                              coord_cartesian(xlim = c(start, end), expand = FALSE) + theme_void(), rows = 1))
  spans <- g[FEATURE == "transcript"]
  spans <- rbind(spans, g[FEATURE == "gene" & !(GENE %in% spans$GENE)])
  if (nrow(spans[TYPE == "protein_coding"]) > 25) spans <- spans[TYPE == "protein_coding"]
  spans <- spans[order(START)]
  pad <- (end - start) * 0.06
  row_end <- c(); spans[, ROW := 0L]
  for (k in seq_len(nrow(spans))) {
    free <- which(row_end < spans$START[k] - pad)
    if (length(free)) { r <- free[1] } else { r <- length(row_end) + 1 }
    row_end[r] <- spans$END[k] + nchar(spans$GENE[k]) * (end - start) / 120
    spans$ROW[k] <- r
  }
  exons <- merge(g[FEATURE == "exon"], spans[, .(GENE, ROW)], by = "GENE")
  spans[, COLOUR := ifelse(TYPE == "protein_coding", "#1f4e79", "#8c8c8c")]
  exons[, COLOUR := ifelse(TYPE == "protein_coding", "#1f4e79", "#8c8c8c")]
  spans[, LABEL := paste0(GENE, ifelse(STRAND == "+", " >", " <"))]
  p <- ggplot() +
    geom_segment(data = spans, aes(x = pmax(START, start), xend = pmin(END, end), y = -ROW, yend = -ROW,
                                   colour = COLOUR), linewidth = 0.35) +
    geom_rect(data = exons, aes(xmin = pmax(START, start), xmax = pmin(END, end), ymin = -ROW - 0.18,
                                ymax = -ROW + 0.18, fill = COLOUR), colour = NA) +
    geom_text(data = spans, aes(x = pmin(pmax((START + END) / 2, start), end), y = -ROW + 0.42, label = LABEL,
                                colour = COLOUR), size = 1.9, fontface = "italic") +
    scale_colour_identity() + scale_fill_identity() +
    coord_cartesian(xlim = c(start, end), ylim = c(-max(spans$ROW) - 0.5, -0.3), expand = FALSE) +
    labs(y = "Genes", x = NULL) + theme_pub +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  list(plot = p, rows = max(spans$ROW))
}

# ---------------------------------------------------------------------------- LD triangle
ld_triangle <- function(dir, start, end) {
  v <- readf(file.path(dir, "ld_variants.tsv")); pairs <- readf(file.path(dir, "ld_pairs.tsv"))
  if (nrow(v) < 2 || !nrow(pairs)) return(NULL)
  pos <- v$POS; n <- length(pos)
  pairs <- pairs[pos[J] - pos[I] <= ld_window_bp]
  if (!nrow(pairs)) return(NULL)
  mids <- (pos[-1] + pos[-n]) / 2
  left <- c(pos[1] - (pos[2] - pos[1]) / 2, mids); right <- c(mids, pos[n] + (pos[n] - pos[n - 1]) / 2)
  corners <- function(u, w) list(x = (u + w) / 2, y = -(w - u) / 2)
  a <- corners(left[pairs$I], left[pairs$J]); b <- corners(right[pairs$I], left[pairs$J])
  c <- corners(right[pairs$I], right[pairs$J]); d <- corners(left[pairs$I], right[pairs$J])
  id <- seq_len(nrow(pairs))
  poly <- data.table(ID = rep(id, 4), R2 = rep(pairs$R2, 4),
                     X = c(a$x, b$x, c$x, d$x), Y = c(a$y, b$y, c$y, d$y))
  setorder(poly, ID)
  ggplot(poly, aes(X, Y, group = ID, fill = R2)) + geom_polygon(colour = NA) +
    scale_fill_gradientn(colours = c("white", "#fcbba1", "#ef3b2c", "#67000d"), limits = c(0, 1),
                         name = expression(r^2)) +
    coord_cartesian(xlim = c(start, end), expand = FALSE) +
    scale_x_continuous(labels = mb) +
    labs(x = NULL, y = sprintf("LD r2
(%d var., <=%g kb)", n, ld_window_bp / 1000)) + theme_pub +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
}

# ---------------------------------------------------------------------------- per-locus figure
r2_bin <- function(d) {
  b <- cut(d$R2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), right = FALSE,
           labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">=0.8"))
  b <- as.character(b); b[is.na(b)] <- "not in LD panel"
  factor(b, levels = names(R2_COLOURS))
}

plot_locus <- function(gcst, locus) {
  dir <- file.path(data_dir, gcst, "loci", locus$LOCUS_ID)
  d <- readf(file.path(dir, "variants.tsv"))
  if (!nrow(d)) return(NULL)
  start <- locus$START; end <- locus$END; chr <- locus$CHR
  d[, LOGP := -log10(pmax(P, 1e-300))][, R2BIN := r2_bin(d)]
  d[, CS := sub(";.*", "", CREDIBLE_SET)][is.na(CS), CS := ""]
  idx <- d[INDEX == TRUE]
  ymax <- max(d$LOGP, -log10(5e-8)) * 1.15

  pa <- ggplot(d[order(R2BIN != "not in LD panel", R2)], aes(POS, LOGP)) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "#c0392b", linewidth = 0.25) +
    geom_point(aes(fill = R2BIN), shape = 21, colour = "grey25", size = 1.25, stroke = 0.12) +
    geom_point(data = idx, shape = 23, fill = "#7b1fa2", colour = "black", size = 2.8, stroke = 0.3) +
    geom_text_repel(data = idx, aes(label = LABEL), size = 2.2, min.segment.length = 0, segment.size = 0.2) +
    scale_fill_manual(values = R2_COLOURS, drop = FALSE, name = expression(r^2 ~ "with index")) +
    coord_cartesian(xlim = c(start, end), ylim = c(0, ymax), expand = FALSE) +
    labs(x = NULL, y = expression(-log[10](italic(P))),
         title = sprintf("%s  %s  %s  chr%d:%s-%s Mb%s", title_prefix, gcst, locus$LOCUS_ID, chr, mb(start), mb(end),
                         ifelse(is.na(locus$NEAREST_GENE) || locus$NEAREST_GENE == "", "", paste0("  (", locus$NEAREST_GENE, ")"))),
         subtitle = paste(strwrap(sprintf("Locus QC: %s%s | LD/z mismatch s = %s | credible sets: %s",
                            locus$LOCUS_QC, ifelse(is.na(locus$QC_NOTES) || locus$QC_NOTES == "", "", paste0(" - ", locus$QC_NOTES)),
                            ifelse(is.na(suppressWarnings(as.numeric(locus$S_RSS))), "NA", sprintf("%.2f", as.numeric(locus$S_RSS))),
                            ifelse(is.na(locus$N_CS), 0, locus$N_CS)), 120), collapse = "
")) +
    theme_pub + theme(axis.text.x = element_blank(), legend.position = "right")

  fm <- d[!is.na(PIP)]
  cs_names <- sort(unique(fm$CS[fm$CS != ""]))
  cs_col <- setNames(rep(CS_COLOURS, length.out = length(cs_names)), cs_names)
  fm[, CSLAB := ifelse(CS == "", "not in 95% CS", sub(".*_", "", CS))]
  lab_col <- c(setNames(cs_col, sub(".*_", "", names(cs_col))), "not in 95% CS" = "grey70")
  best <- fm[CS != "", .SD[which.max(PIP)], by = CS]
  pb <- ggplot(fm[order(CS != "")], aes(POS, PIP)) +
    geom_segment(data = fm[CS != ""], aes(xend = POS, y = 0, yend = PIP, colour = CSLAB), linewidth = 0.25) +
    geom_point(aes(colour = CSLAB), size = ifelse(fm[order(CS != "")]$CS != "", 1.5, 0.8)) +
    geom_text_repel(data = best, aes(label = LABEL, colour = CSLAB), size = 2.1, show.legend = FALSE,
                    min.segment.length = 0, segment.size = 0.2) +
    scale_colour_manual(values = lab_col, name = "95% credible set") +
    coord_cartesian(xlim = c(start, end), ylim = c(0, 1.08), expand = FALSE) +
    labs(x = NULL, y = "PIP (SuSiE)") + theme_pub + theme(axis.text.x = element_blank())

  gt <- gene_track(chr, start, end)
  ld <- ld_triangle(dir, start, end)
  bottom_axis <- theme(axis.text.x = element_text(size = 6))
  if (is.null(ld)) {
    fig <- pa / pb / (gt$plot + scale_x_continuous(labels = mb) + bottom_axis + labs(x = sprintf("chr%d (Mb, GRCh38)", chr))) +
      plot_layout(heights = c(3, 1.6, 0.5 + 0.35 * gt$rows))
  } else {
    fig <- pa / pb / gt$plot / (ld + labs(x = sprintf("chr%d (Mb, GRCh38)", chr)) + bottom_axis) +
      plot_layout(heights = c(3, 1.6, 0.5 + 0.35 * gt$rows, 1.6))
  }
  fig
}

# ---------------------------------------------------------------------------- shared loci across GWAS
plot_shared <- function(cl) {
  members <- loci_all[COMBINED_ID == cl$COMBINED_ID & !grepl("^lead_", LOCUS_ID)]
  if (length(unique(members$GCST)) < 2) return(NULL)
  rows <- list()
  for (k in seq_len(nrow(members))) {
    d <- readf(file.path(data_dir, members$GCST[k], "loci", members$LOCUS_ID[k], "variants.tsv"))
    if (nrow(d)) rows[[k]] <- d[, GCST := members$GCST[k]]
  }
  d <- rbindlist(rows, fill = TRUE)
  d <- d[POS >= cl$START & POS <= cl$END]
  d[, LOGP := -log10(pmax(P, 1e-300))][, R2BIN := r2_bin(d)]
  pa <- ggplot(d[order(R2BIN != "not in LD panel", R2)], aes(POS, LOGP)) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "#c0392b", linewidth = 0.25) +
    geom_point(aes(fill = R2BIN), shape = 21, colour = "grey25", size = 1.1, stroke = 0.12) +
    geom_point(data = d[INDEX == TRUE], shape = 23, fill = "#7b1fa2", size = 2.4, stroke = 0.3) +
    facet_grid(GCST ~ ., scales = "free_y") +
    scale_fill_manual(values = R2_COLOURS, drop = FALSE, name = expression(r^2 ~ "with index")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +   # room for the index diamond
    coord_cartesian(xlim = c(cl$START, cl$END)) +
    labs(x = NULL, y = expression(-log[10](italic(P))),
         title = sprintf("%s  shared locus %s  chr%d:%s-%s Mb  (%d GWAS)", title_prefix, cl$NAME, cl$CHR,
                         mb(cl$START), mb(cl$END), length(unique(d$GCST)))) +
    theme_pub + theme(axis.text.x = element_blank())
  # one open symbol per GWAS so coincident PIPs (e.g. both = 1) stay visible
  fm <- d[!is.na(PIP)][order(CREDIBLE_SET != "")]
  gwas_shape <- setNames(rep(c(21, 24, 22, 23, 25), length.out = length(gwas_col)), names(gwas_col))
  pb <- ggplot(fm, aes(POS, PIP, colour = GCST, shape = GCST)) +
    geom_point(data = fm[CREDIBLE_SET == ""], size = 0.7, alpha = 0.6) +
    geom_point(data = fm[CREDIBLE_SET != ""], aes(fill = GCST), size = 2.4, stroke = 0.5, alpha = 0.55) +
    scale_colour_manual(values = gwas_col, name = "GWAS") +
    scale_fill_manual(values = gwas_col, guide = "none") +
    scale_shape_manual(values = gwas_shape, name = "GWAS") +
    coord_cartesian(xlim = c(cl$START, cl$END), ylim = c(0, 1.08), expand = FALSE) +
    labs(x = NULL, y = "PIP (large = in 95% CS)") + theme_pub + theme(axis.text.x = element_blank())
  gt <- gene_track(cl$CHR, cl$START, cl$END)
  pa / pb / (gt$plot + scale_x_continuous(labels = mb) + theme(axis.text.x = element_text(size = 6)) +
               labs(x = sprintf("chr%d (Mb, GRCh38)", cl$CHR))) +
    plot_layout(heights = c(1.7 * length(unique(d$GCST)), 1.4, 0.5 + 0.35 * gt$rows))
}

# ============================================================================ run
for (s in seq_len(nrow(studies))) {
  info <- studies[s]; gcst <- info$GCST
  base <- file.path(out_dir, gcst)
  cat("\n==", gcst, "==\n")
  save_gg(plot_manhattan(gcst, info), file.path(base, "manhattan"), 11, 4.2)
  cat("  manhattan\n")

  loci <- loci_all[GCST == gcst]
  if (nrow(loci)) {
    tracks <- loci[, .(CHR, START, END, GROUP = gcst, COLOUR = "#c0392b")]
    labels <- loci[order(MIN_P)][1:min(.N, 40), .(CHR, POS = TOP_POS,
                   LABEL = ifelse(is.na(NEAREST_GENE) | NEAREST_GENE == "", TOP_LABEL, NEAREST_GENE))]
    save_base(file.path(base, "ideogram"), 9, 6,
              function() draw_ideogram(tracks, sprintf("%s %s: %d loci", title_prefix, gcst, nrow(loci)), labels))
    cat("  ideogram\n")
  }

  fine <- loci[!grepl("^lead_", LOCUS_ID)]
  pages <- list()
  for (k in seq_len(nrow(fine))) {
    fig <- tryCatch(plot_locus(gcst, fine[k]), error = function(e) { message("  ", fine$LOCUS_ID[k], ": ", conditionMessage(e)); NULL })
    if (is.null(fig)) next
    save_gg(fig, file.path(base, "loci", paste0(fine$LOCUS_ID[k], "_region")), 7.2, 7.8)
    pages[[length(pages) + 1]] <- fig
  }
  if (length(pages)) {
    pdf(file.path(base, "all_loci_multipage.pdf"), width = 7.2, height = 7.8, onefile = TRUE, useDingbats = FALSE)
    for (fig in pages) print(fig)
    dev.off()
    cat("  ", length(pages), "locus figures + all_loci_multipage.pdf\n")
  }
}

# ---------------------------------------------------------------------------- combined
cdir <- file.path(out_dir, "combined")
dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
if (nrow(loci_all)) {
  tracks <- loci_all[, .(CHR, START, END, GROUP = GCST, COLOUR = gwas_col[GCST])]
  labels <- merge(combined, loci_all[, .(MIN_P = min(MIN_P)), by = COMBINED_ID], by = "COMBINED_ID")
  labels <- labels[order(MIN_P)][1:min(.N, 45), .(CHR, POS = (START + END) / 2, LABEL = sub("^chr[0-9]+:[0-9.]+Mb ?", "", NAME))]
  labels[LABEL == "", LABEL := "."]
  save_base(file.path(cdir, "ideogram_all_gwas"), 10, 6.5,
            function() draw_ideogram(tracks, sprintf("%s: loci of %d GWAS", title_prefix, nrow(studies)), labels))
  cat("\ncombined ideogram\n")

  summary <- merge(loci_all, combined[, .(COMBINED_ID, NAME)], by = "COMBINED_ID")
  cells <- summary[, .(MIN_P = min(MIN_P), MAX_PIP = suppressWarnings(max(MAX_PIP, na.rm = TRUE)),
                       N_CS = suppressWarnings(sum(as.numeric(N_CS), na.rm = TRUE)),
                       TOP = TOP_LABEL[which.min(MIN_P)]), by = .(COMBINED_ID, NAME, GCST)]
  cells[!is.finite(MAX_PIP), MAX_PIP := NA]
  fwrite(cells, file.path(cdir, "locus_summary.tsv"), sep = "\t")

  if (HAS_CHM && nrow(combined)) {
    rn <- combined[order(CHR, START)]$NAME; cn <- studies$GCST
    mat <- function(v) { m <- matrix(NA_real_, length(rn), length(cn), dimnames = list(rn, cn))
                         for (i in seq_len(nrow(cells))) m[cells$NAME[i], cells$GCST[i]] <- cells[[v]][i]; m }
    mp <- -log10(pmax(mat("MIN_P"), 1e-300)); mpip <- mat("MAX_PIP"); mcs <- mat("N_CS")
    top <- max(c(mp, 10), na.rm = TRUE)
    h1 <- ComplexHeatmap::Heatmap(mp, name = "-log10 P", na_col = "grey95",
            heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
                                        labels_gp = grid::gpar(fontsize = 6)),
            col = circlize::colorRamp2(c(0, -log10(5e-8), top), c("white", "#fdd49e", "#b30000")),
            cluster_rows = FALSE, cluster_columns = FALSE, column_title = "min P (locus)",
            row_names_gp = grid::gpar(fontsize = 6.5), column_names_gp = grid::gpar(fontsize = 6.5),
            column_title_gp = grid::gpar(fontsize = 7.5, fontface = "bold"), border = TRUE,
            rect_gp = grid::gpar(col = "white", lwd = 0.5))
    h2 <- ComplexHeatmap::Heatmap(mpip, name = "max PIP", na_col = "grey95",
            heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
                                        labels_gp = grid::gpar(fontsize = 6)),
            col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#9ecae1", "#08519c")),
            cluster_rows = FALSE, cluster_columns = FALSE, column_title = "max PIP (n credible sets)",
            row_names_gp = grid::gpar(fontsize = 6.5), column_names_gp = grid::gpar(fontsize = 6.5),
            column_title_gp = grid::gpar(fontsize = 7.5, fontface = "bold"), border = TRUE,
            rect_gp = grid::gpar(col = "white", lwd = 0.5),
            cell_fun = function(j, i, x, y, w, h, fill) if (!is.na(mcs[i, j]) && mcs[i, j] > 0)
              grid::grid.text(mcs[i, j], x, y, gp = grid::gpar(fontsize = 6)))
    height <- max(3, 0.17 * length(rn) + 1.6); width <- max(4.5, 0.55 * length(cn) * 2 + 3)
    save_base(file.path(cdir, "locus_by_gwas"), width, height,
              function() ComplexHeatmap::draw(h1 + h2, merge_legend = TRUE,
                                              column_title = sprintf("%s: %d loci x %d GWAS", title_prefix, length(rn), length(cn)),
                                              column_title_gp = grid::gpar(fontsize = 8, fontface = "bold")))
    cat("locus x GWAS matrix\n")
  }

  shared <- combined[N_GWAS >= 2]
  n_shared <- 0
  for (k in seq_len(nrow(shared))) {
    fig <- tryCatch(plot_shared(shared[k]), error = function(e) { message("  shared ", shared$NAME[k], ": ", conditionMessage(e)); NULL })
    if (is.null(fig)) next
    n_members <- length(unique(loci_all[COMBINED_ID == shared$COMBINED_ID[k]]$GCST))
    save_gg(fig, file.path(cdir, "shared_loci", shared$FILE_STEM[k]), 7.2, 2.6 + 1.5 * n_members)
    n_shared <- n_shared + 1
  }
  cat(n_shared, "shared-locus comparison figures\n")
}
cat("\nFigures in", out_dir, "\n")
'''


def run_figures(P, args):
    rscript = PROJECT_ROOT / "envs" / "pipeline" / "bin" / "Rscript"
    rscript = str(rscript) if rscript.exists() else shutil.which("Rscript")
    if not rscript:
        raise SystemExit("Rscript not found. Run: python Step00_Installation.py")

    check = subprocess.run(
        [rscript, "-e", "pk <- c(" + ",".join(f'"{p}"' for p in R_PACKAGES) + "); "
         "m <- pk[!vapply(pk, requireNamespace, logical(1), quietly=TRUE)]; cat(m, sep=' ')"],
        capture_output=True, text=True)
    missing = check.stdout.split()
    required = {"ggplot2", "patchwork", "ggrepel", "scales", "data.table"}
    if required & set(missing):
        raise SystemExit(f"Missing R packages: {' '.join(sorted(required & set(missing)))}\n"
                         "Run: python Step00_Installation.py")
    if missing:
        print(f"Optional R packages missing (reduced figures): {' '.join(missing)}")

    script = P["out"] / "step07_plots.R"
    script.write_text(R_CODE, encoding="utf-8", newline="\n")
    title = f"{P['phenotype']} ({P['code']})"
    banner("STEP 07 - DRAWING FIGURES (R)")
    result = subprocess.run([rscript, str(script), str(P["data"]), str(P["out"]), args.formats,
                             str(args.dpi), title, str(args.ld_window_kb * 1000)])
    if result.returncode != 0:
        raise SystemExit(f"R figure script failed (exit {result.returncode}): {script}")

    (P["out"] / "step07_run.json").write_text(json.dumps({
        "phenotype": P["phenotype"], "ancestry": P["code"], "formats": args.formats, "dpi": args.dpi,
        "max_loci": args.max_loci, "ld_max": args.ld_max}, indent=2), encoding="utf-8")


# =============================================================================
# --status
# =============================================================================

def status_mode(P):
    banner(f"STEP 07 STATUS - {P['phenotype']} / {P['label']} ({P['code']})")
    studies = read_tsv(P["data"] / "studies.tsv")
    loci = read_tsv(P["data"] / "loci.tsv")
    if studies.empty:
        print("No figures yet. Run without --status.")
        return
    print(f"{'STUDY':<15}{'MANHATTAN':<11}{'IDEOGRAM':<10}{'LOCUS FIGS':>11}  MULTIPAGE")
    for accession in studies["GCST"]:
        base = P["out"] / accession
        n_loci = int(((loci["GCST"] == accession) & ~loci["LOCUS_ID"].astype(str).str.startswith("lead_")).sum()) \
            if not loci.empty else 0
        n_figs = len({p.stem for p in (base / "loci").glob("*_region.*")}) if (base / "loci").exists() else 0
        print(f"{accession:<15}{'yes' if list(base.glob('manhattan.*')) else 'no':<11}"
              f"{'yes' if list(base.glob('ideogram.*')) else 'no':<10}{f'{n_figs}/{n_loci}':>11}  "
              f"{'yes' if (base / 'all_loci_multipage.pdf').exists() else 'no'}")
    combined = P["out"] / "combined"
    print()
    print(f"combined ideogram   : {'yes' if list(combined.glob('ideogram_all_gwas.*')) else 'no'}")
    print(f"locus x GWAS matrix : {'yes' if list(combined.glob('locus_by_gwas.*')) else 'no'}")
    shared = {p.stem for p in (combined / "shared_loci").glob("*.*")} if (combined / "shared_loci").exists() else set()
    print(f"shared-locus figures: {len(shared)}")
    print(f"\nFolder: {P['out']}")


# =============================================================================
# SLURM (one job; settings from HPC_CPU.sh)
# =============================================================================

def read_slurm_template(path):
    values = {}
    if Path(path).is_file():
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            match = re.match(r"^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$", line)
            if match and (match.group(2) or "").strip():
                values[match.group(1)] = match.group(2).strip()
    return values


def write_slurm(P, args):
    values = {"nodes": "1", "ntasks": "1", "cpus-per-task": "4", "mem": "32G", "time": "12:00:00"}
    values.update({k: v for k, v in read_slurm_template(PROJECT_ROOT / args.template).items()
                   if k not in {"job-name", "output", "error", "array"}})
    for key, value in (("partition", args.partition), ("qos", args.qos), ("account", args.account),
                       ("cpus-per-task", args.cpus), ("mem", args.mem), ("time", args.time)):
        if value not in (None, ""):
            values[key] = str(value)
    log_dir = PROJECT_ROOT / "logs" / "step07_plots" / P["p_slug"] / P["a_slug"]
    log_dir.mkdir(parents=True, exist_ok=True)
    values.update({"job-name": f"PLOTS_{P['p_slug']}_{P['code'].lower()}"[:100],
                   "output": f"{log_dir}/plots.%j.out", "error": f"{log_dir}/plots.%j.err"})
    order = ["job-name", "account", "partition", "qos", "nodes", "ntasks", "cpus-per-task", "mem", "time",
             "output", "error"]
    header = "\n".join(f"#SBATCH --{k}={values[k]}" for k in order if values.get(k))
    header += "".join(f"\n#SBATCH --{k}={v}" for k, v in values.items() if k not in order and v)

    command = [str(PROJECT_ROOT / "envs" / "pipeline" / "bin" / "python"), str(Path(__file__).resolve()),
               "--phenotype", P["phenotype"], "--ancestry", P["code"], "--formats", args.formats,
               "--dpi", str(args.dpi), "--max-loci", str(args.max_loci), "--ld-max", str(args.ld_max),
               *(["--index", str(args.index)] if args.index is not None else []),
               "--ld-window-kb", str(args.ld_window_kb)]
    bash_file = PROJECT_ROOT / f"Step07_Plots_{P['p_slug']}_{P['a_slug']}.sh"
    bash_file.write_text(f"""#!/bin/bash --login
{header}

set -euo pipefail
cd {shlex.quote(str(PROJECT_ROOT))}
export PATH={shlex.quote(str(PROJECT_ROOT / 'envs' / 'pipeline' / 'bin'))}:"$PATH"

{' '.join(shlex.quote(x) for x in command)}
""", encoding="utf-8", newline="\n")
    bash_file.chmod(0o755)
    print(f"Generated: {bash_file}\nRUN:\n  sbatch {bash_file.name}")


# =============================================================================
# MAIN
# =============================================================================

def arguments(argv=None):
    parser = argparse.ArgumentParser(description="Publication figures for Step05 loci and Step06 fine-mapping.")
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)
    parser.add_argument("--index", type=int, default=None,
                        help="Only study N of the Step02 download manifest (the number Step02 prints).")
    parser.add_argument("--formats", default="pdf,png",
                        help="Comma-separated: pdf, svg, png (default pdf,png).")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300).")
    parser.add_argument("--max-loci", type=int, default=0,
                        help="Locus figures per study (0 = all, strongest loci first in Step06 order).")
    parser.add_argument("--ld-max", type=int, default=120,
                        help="Variants in the LD triangle (credible sets + strongest; default 120).")
    parser.add_argument("--ld-window-kb", type=int, default=500,
                        help="LD triangle shows pairs at most this far apart (default 500 kb).")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--status", action="store_true", help="READ-ONLY: which figures exist.")
    modes.add_argument("--slurm", action="store_true", help="Write a SLURM script instead of plotting here.")
    slurm = parser.add_argument_group("SLURM (default: the #SBATCH header of HPC_CPU.sh)")
    slurm.add_argument("--template", default="HPC_CPU.sh")
    slurm.add_argument("--partition", default=None)
    slurm.add_argument("--qos", default=None)
    slurm.add_argument("--account", default=None)
    slurm.add_argument("--cpus", type=int, default=None)
    slurm.add_argument("--mem", default=None)
    slurm.add_argument("--time", default=None)
    args = parser.parse_args(argv)
    bad = set(args.formats.split(",")) - {"pdf", "svg", "png"}
    if bad:
        parser.error(f"unknown format(s): {', '.join(sorted(bad))}")
    return args


def main():
    args = arguments()
    P = project_paths(args)
    if args.status:
        status_mode(P)
        return
    if args.slurm:
        write_slurm(P, args)
        return
    prepare(P, args)
    run_figures(P, args)
    status_mode(P)


if __name__ == "__main__":
    main()
