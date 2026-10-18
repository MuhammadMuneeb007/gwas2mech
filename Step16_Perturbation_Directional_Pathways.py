#!/usr/bin/env python3
"""GWAS2Mech Step16 v0.1: public perturbation discovery and signed responses.

First run inventories local *real* scPerturb .h5ad experiments against the
Step15 CORE_STRONG gene list. A conservative manually reviewed manifest is
required for differential expression and pathway-response inference.

Commands from GWAS2Mech root:
  envs/perturb/bin/python Step16_Perturbation_Directional_Pathways.py \
    --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --inventory
  # Inspect/edit Step16_Analysis_Manifest_Template.tsv, especially biological replicates.
  envs/perturb/bin/python Step16_Perturbation_Directional_Pathways.py \
    --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x \
    --run-manifest /path/to/reviewed_manifest.tsv --gmt /path/to/full_reactome.gmt
  envs/perturb/bin/python Step16_Perturbation_Directional_Pathways.py \
    --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

Biological requirements:
  * A contrast is performed ONLY when manually confirmed control/target
    labels and >=2 biological replicates per condition are available.
  * Technical batches, cells and gRNAs are not automatically treated as
    independent biological replicates.
  * Count matrices: pseudobulk summed by biological replicate, library-size
    normalised (log2 CPM+1), exploratory Welch's t test over replicates.
  * Log-expression matrices: mean expression by replicate, exploratory
    Welch's t test; differences are *not* log-fold changes.
  * Gene-set enrichment uses ALL pathways from supplied complete GMT file(s),
    correcting across eligible pathways. Step14 significant pathways are
    annotated only after the full-library GSEA (avoids FDR only on selected hits).
  * NES is a direction of transcriptional response (perturbed vs control), NOT
    proof of biochemical pathway activation, direct regulation, risk-allele
    direction or causality.
  * Noncoding genes can be screened for perturbations even without STRING PPI.
  * Lack of a tested gene or incomplete metadata => NOT TESTED, not no effect.

This stage does not suggest or rank drugs and is not GWAS2Therapeutics.
"""
from __future__ import annotations

import argparse
import gc
import gzip
import json
import re
import warnings
import zipfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

VERSION = "0.1.0-discovery-replicate-gsea"
CONTROL_EXACT = {"control", "ctrl", "non-targeting", "nontargeting", "non targeting", "negative control",
                 "non-targeting control", "nontargeting control", "ntc", "scramble", "scrambled", "vehicle"}
TARGET_CANDIDATES = ["perturbation_name", "perturbation", "gene_target", "target_gene", "perturbation_target",
                     "target", "gene", "guide_target", "gene_name", "knockdown_gene", "perturbation_genes"]
ID_RE = re.compile(r"ENSG\d+", re.I)


def slug(s):
    return re.sub(r"[^a-z0-9]+", "_", str(s).lower()).strip("_")


def clean(s):
    if pd.isna(s):
        return ""
    x = str(s).strip()
    return "" if x.lower() in ("nan", "na", "none", "null", "", "<na>") else x


def clean_gene(s):
    s = clean(s)
    match = ID_RE.search(s)
    return match.group(0).upper() if match else s.upper()


def is_control(s):
    s = re.sub(r"[_\s]+", " ", clean(s).lower()).strip()
    return s in CONTROL_EXACT


def get_paths(args):
    root = Path(args.root).expanduser().resolve()
    ancestry = {"eur":"european", "afr":"african", "eas":"east_asian",
                "sas":"south_asian", "amr":"hispanic_or_latin_american"}.get(args.ancestry.lower())
    if not ancestry:
        raise SystemExit("--ancestry must be EUR, AFR, EAS, SAS or AMR")
    suffix = "" if args.reference_panel == "1000G" else "_" + args.reference_panel
    rel = Path(slug(args.phenotype)) / (ancestry + suffix)
    return root, root / "15_network" / rel, root / "16_perturbation" / rel


def load_step15(step15):
    handoff = step15 / "Step15_Gene_Pathway_Handoff.tsv"
    audit = step15 / "Step15_Gene_Mapping_Audit.tsv"
    if not handoff.exists():
        raise FileNotFoundError(f"Step15 handoff missing: {handoff}")
    h = pd.read_csv(handoff, sep="\t", dtype=str).fillna("")
    if not {"GENE_ID", "GENE_SYMBOL"}.issubset(h):
        raise ValueError(f"Missing required GENE_ID and GENE_SYMBOL in {handoff}")
    if audit.exists():
        a = pd.read_csv(audit, sep="\t", dtype=str).fillna("")
        a = a.loc[a["IS_GWAS2M_PRIMARY"].astype(str).str.upper().eq("YES")].copy()
        a["GENE_ID"] = a["GENE_ID"].map(clean_gene)
        a = a.drop_duplicates("GENE_ID")
        genes = a[["GENE_ID", "GENE_SYMBOL"]].copy()
    else:
        h["GENE_ID"] = h["GENE_ID"].map(clean_gene)
        genes = h[["GENE_ID", "GENE_SYMBOL"]].drop_duplicates("GENE_ID")
    genes = genes.loc[genes.GENE_ID.str.match(r"ENSG\d+$")].copy()
    if genes.empty:
        raise ValueError("No Ensembl gene IDs available in Step15 handoff")
    # Resolve symbols from handoff when audit symbol is absent, while preserving ID.
    symmap = dict(zip(h.GENE_ID.map(clean_gene), h.GENE_SYMBOL.map(clean)))
    genes["GENE_SYMBOL"] = [clean(row.GENE_SYMBOL) or symmap.get(row.GENE_ID, "") for row in genes.itertuples()]
    h["GENE_ID"] = h["GENE_ID"].map(clean_gene)
    for c in ["SOURCE", "TERM_ID", "TERM_NAME", "FDR"]:
        if c not in h:
            h[c] = ""
    pathways = h.loc[h.TERM_ID.ne("") | h.TERM_NAME.ne(""), ["SOURCE", "TERM_ID", "TERM_NAME", "FDR"]].drop_duplicates()
    return genes.reset_index(drop=True), pathways.reset_index(drop=True)


def dataset_paths(root, data_root):
    if data_root:
        base = Path(data_root).expanduser()
        if not base.is_absolute():
            base = root / base
    else:
        base = root / "resources" / "perturbation" / "scperturb"
    return base, sorted(base.rglob("*.h5ad")) if base.exists() else []


def load_adata(path):
    try:
        import anndata as ad
    except ImportError:
        raise SystemExit("Missing anndata: use envs/perturb/bin/python (pip install anndata)")
    return ad.read_h5ad(path, backed="r")


def inventory(root, step15, out, args):
    genes, pathways = load_step15(step15)
    base, paths = dataset_paths(root, args.data_root)
    if args.limit_files:
        paths = paths[:args.limit_files]
    targets = {}
    for row in genes.itertuples(index=False):
        targets[row.GENE_ID.upper()] = (row.GENE_ID, clean(row.GENE_SYMBOL))
        if clean(row.GENE_SYMBOL):
            targets[row.GENE_SYMBOL.upper()] = (row.GENE_ID, clean(row.GENE_SYMBOL))
    recs, coverage, template = [], [], []
    print(f"[Step16] Scanning {len(paths)} local .h5ad datasets below {base}", flush=True)
    for idx, path in enumerate(paths, 1):
        drec = {"DATASET_PATH": str(path), "DATASET": path.name, "STATUS": "", "N_CELLS": 0,
                "N_MEASURED_GENES": 0, "N_MATCHED_PRIMARY_TARGETS": 0, "BEST_PERTURBATION_COLUMN": "",
                "INFERRED_CONTROL_LABEL": "", "MATCHED_TARGETS": "", "POSSIBLE_REPLICATE_COLUMNS": "",
                "CANDIDATE_PERTURBATION_COLUMNS": "", "ERROR": ""}
        adata = None
        try:
            adata = load_adata(path)
            drec["N_CELLS"], drec["N_MEASURED_GENES"] = adata.shape
            cols = list(adata.obs.columns)
            cand = [c for c in cols if c.lower() in TARGET_CANDIDATES or
                    any(w in c.lower() for w in ("perturb", "target", "guide", "knockdown"))]
            drec["CANDIDATE_PERTURBATION_COLUMNS"] = ";".join(cand[:30])
            drec["POSSIBLE_REPLICATE_COLUMNS"] = ";".join([c for c in cols if any(t in c.lower() for t in
                ("replicate", "donor", "patient", "sample", "batch"))][:30])
            best = None
            for col in cand:
                series = adata.obs[col]
                if isinstance(series.dtype, pd.CategoricalDtype) or series.dtype == object or pd.api.types.is_string_dtype(series):
                    vc = series.astype(str).value_counts(dropna=True)
                    hits = [(label, targets[clean_gene(label)]) for label in vc.index
                            if clean_gene(label) in targets]
                    controls = [label for label in vc.index if is_control(label)]
                    if not hits:
                        continue
                    score = (len({t[0] for _, t in hits}), sum(int(vc[label]) for label, _ in hits), bool(controls))
                    if best is None or score > best[0]:
                        best = (score, col, hits, controls, vc)
            if best:
                _, col, hits, controls, vc = best
                control = str(max(controls, key=lambda x: vc[x])) if controls else ""
                drec.update({"BEST_PERTURBATION_COLUMN": col,
                             "INFERRED_CONTROL_LABEL": control,
                             "N_MATCHED_PRIMARY_TARGETS": len({g[0] for _, g in hits}),
                             "MATCHED_TARGETS": ";".join(sorted({g[1] or g[0] for _,g in hits}))})
                for label, (gid,symbol) in hits:
                    coverage.append({"DATASET": path.name, "DATASET_PATH": str(path),
                                     "GENE_ID": gid, "GENE_SYMBOL": symbol, "TARGET_VALUE": str(label),
                                     "TARGET_CELLS": int(vc[label]), "PERTURBATION_COLUMN": col,
                                     "CONTROL_VALUE": control, "CONTROL_CELLS": int(vc.get(control, 0)),
                                     "EXPERIMENTAL_CONTRAST_READY": "NO_METADATA_REVIEW_REQUIRED"})
                    template.append({"INCLUDE": "NO", "DATASET_PATH": str(path), "GENE_ID": gid,
                                     "GENE_SYMBOL": symbol, "PERTURBATION_COLUMN": col,
                                     "CONTROL_VALUE": control, "TARGET_VALUE": str(label),
                                     "REPLICATE_COLUMN": "", "CELL_TYPE_COLUMN": "", "CELL_TYPE_VALUE": "",
                                     "PERTURBATION_MODE": "REVIEW_REQUIRED",
                                     "EXPRESSION_LAYER": "X", "EXPRESSION_KIND": "REVIEW_REQUIRED",
                                     "GENE_SYMBOL_COLUMN": "var_names", "NOTES": "Review counts/log-expression, biological replicates, controls and cell type before INCLUDE=YES"})
            drec["STATUS"] = "READ_OK"
        except Exception as e:
            drec["STATUS"] = "UNREADABLE_OR_PARTIAL"
            drec["ERROR"] = f"{type(e).__name__}: {str(e)[:180]}"
        finally:
            if adata is not None:
                try: adata.file.close()
                except Exception: pass
                del adata
                gc.collect()
        recs.append(drec)
        print(f"[{idx:02d}/{len(paths):02d}] {path.name}: {drec['STATUS']} "
              f"matched={drec['N_MATCHED_PRIMARY_TARGETS']}", flush=True)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(recs, columns=list(recs[0]) if recs else ["DATASET_PATH","DATASET","STATUS"]).to_csv(
        out/"Step16_Dataset_Inventory.tsv",sep="\t",index=False)
    covcols = ["DATASET", "DATASET_PATH", "GENE_ID", "GENE_SYMBOL", "TARGET_VALUE", "TARGET_CELLS",
               "PERTURBATION_COLUMN", "CONTROL_VALUE", "CONTROL_CELLS", "EXPERIMENTAL_CONTRAST_READY"]
    pd.DataFrame(coverage,columns=covcols).to_csv(out/"Step16_Perturbation_Coverage.tsv",sep="\t",index=False)
    tcols = ["INCLUDE", "DATASET_PATH", "GENE_ID", "GENE_SYMBOL", "PERTURBATION_COLUMN", "CONTROL_VALUE",
             "TARGET_VALUE", "REPLICATE_COLUMN", "CELL_TYPE_COLUMN", "CELL_TYPE_VALUE", "PERTURBATION_MODE", "EXPRESSION_LAYER",
             "EXPRESSION_KIND", "GENE_SYMBOL_COLUMN", "NOTES"]
    pd.DataFrame(template,columns=tcols).to_csv(out/"Step16_Analysis_Manifest_Template.tsv",sep="\t",index=False)
    matched = {r["GENE_ID"] for r in coverage}
    pd.DataFrame([{**row._asdict(),"FOUND_IN_LOCAL_SC_PERTURB": "YES" if row.GENE_ID in matched else "NO",
                   "N_DATASETS": len({r["DATASET_PATH"] for r in coverage if r["GENE_ID"]==row.GENE_ID})}
                  for row in genes.itertuples(index=False)]).to_csv(out/"Step16_Gene_Coverage_Summary.tsv",sep="\t",index=False)
    summary = {"version": VERSION,"mode":"inventory","n_files_scanned":len(paths),
               "n_files_readable":sum(x["STATUS"]=="READ_OK" for x in recs),
               "n_genes_primary":len(genes),"n_genes_with_exact_target_label":len(matched),
               "n_matching_dataset_gene_pairs":len(coverage),"n_step14_pathway_terms":len(pathways),
               "note":"Exact label matching is discovery only; not evidence of usable controls or replicates."}
    (out/"Step16_summary.json").write_text(json.dumps(summary,indent=2))
    print(f"[DONE] {summary}; files in {out}")
    return summary


@contextmanager
def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield handle
    elif path.suffix == ".zip":
        import io
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.lower().endswith(".gmt")]
            if len(names) != 1:
                raise ValueError(f"GMT archive must contain exactly one .gmt file: {path} ({names})")
            with z.open(names[0]) as raw:
                with io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as handle:
                    yield handle
    else:
        with path.open("rt", encoding="utf-8", errors="replace") as handle:
            yield handle


def load_gmts(paths):
    libraries = {}
    meta = {}
    for p in paths:
        path = Path(p).expanduser()
        if not path.exists():
            raise FileNotFoundError(f"Full pathway gene-set file not found: {path}")
        library = re.sub(r"\.gmt(?:\.gz|\.zip)?$", "", path.name, flags=re.I)
        with open_text(path) as handle:
            for line in handle:
                parts=line.rstrip("\r\n").split("\t")
                if len(parts)<3: continue
                name,desc=parts[0].strip(),parts[1].strip()
                genes={clean_gene(t) for t in parts[2:] if clean_gene(t)}
                if len(genes)<5: continue
                key=f"{library}::{name}"
                if key in libraries: continue
                libraries[key]=sorted(genes)
                meta[key]={"LIBRARY":library, "NAME":name, "DESCRIPTION":desc}
    if not libraries:
        raise ValueError("No pathway gene sets with >=5 entries found in --gmt files")
    return libraries,meta


def gene_names(adata, col):
    if col == "var_names":
        names = pd.Index(adata.var_names).astype(str).to_numpy()
    elif col in adata.var:
        names = adata.var[col].astype(str).to_numpy()
    else:
        raise ValueError(f"GENE_SYMBOL_COLUMN {col!r} absent. Try var_names or: {list(adata.var.columns)[:30]}")
    return np.array([clean_gene(x) for x in names])


def group_pseudobulk(adata, obs, group_col, replicate_col, expr_layer, kind, min_cells):
    if kind not in {"counts", "lognorm"}:
        raise ValueError("EXPRESSION_KIND must be 'counts' or 'lognorm'")
    X = adata.X if expr_layer.upper() == "X" else adata.layers[expr_layer]
    groups = obs[group_col].astype(str)
    replicas = obs[replicate_col].astype(str)
    labels = sorted({(g,r) for g,r in zip(groups,replicas) if g in {"PERTURB", "CONTROL"} and r not in {"nan", "", "None"}})
    sample_idx = {(g,r):i for i,(g,r) in enumerate(labels)}
    indices = defaultdict(list)
    for i,(g,r) in enumerate(zip(groups,replicas)):
        if (g,r) in sample_idx:
            indices[(g,r)].append(i)
    good = [k for k in labels if len(indices[k])>=min_cells]
    if len([k for k in good if k[0]=="CONTROL"])<2 or len([k for k in good if k[0]=="PERTURB"])<2:
        raise ValueError(f"Fewer than 2 independent replicates per condition with >= {min_cells} cells: "
                         f"{dict(Counter(k[0] for k in good))}")
    idx = {k:i for i,k in enumerate(good)}
    sums = np.zeros((len(good), adata.n_vars), dtype=np.float64)
    n_cells = np.zeros(len(good),dtype=np.int64)
    # Bounded-memory sequential reads: do not load the full experimental matrix.
    chunk_size=256
    group_arr = groups.to_numpy()
    rep_arr = replicas.to_numpy()
    for start in range(0,adata.n_obs,chunk_size):
        stop=min(start+chunk_size,adata.n_obs)
        valid = [j for j in range(start,stop) if (group_arr[j],rep_arr[j]) in idx]
        if not valid: continue
        block = X[start:stop]
        for key in set((group_arr[j],rep_arr[j]) for j in valid):
            row_local = [j-start for j in valid if (group_arr[j],rep_arr[j])==key]
            sub = block[row_local]
            if sparse.issparse(sub):
                v=np.asarray(sub.sum(axis=0)).ravel()
            else:
                v=np.asarray(sub).sum(axis=0)
            sums[idx[key]] += v
            n_cells[idx[key]] += len(row_local)
    if kind=="counts":
        if np.any(sums < 0):
            raise ValueError("Negative values in count matrix. Use EXPRESSION_KIND=lognorm if already normalized.")
        libs=sums.sum(axis=1)
        if np.any(libs<=0): raise ValueError("Zero total counts in a pseudobulk replicate")
        expr=np.log2(1.0+(sums/libs[:,None])*1e6)
    else:
        expr=sums/n_cells[:,None]
    labels_df=pd.DataFrame([{"CONDITION":g,"REPLICATE":r,"N_CELLS":int(n_cells[i])}
                            for i,(g,r) in enumerate(good)])
    return expr, labels_df


def diff_stats(expr, info, names, kind):
    cond=info.CONDITION.to_numpy()
    target=expr[cond=="PERTURB"]
    ctrl=expr[cond=="CONTROL"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        test=stats.ttest_ind(target,ctrl,equal_var=False,axis=0,nan_policy="omit")
    difference=target.mean(axis=0)-ctrl.mean(axis=0)
    t=np.nan_to_num(test.statistic,nan=0.0,posinf=0,neginf=0)
    p=np.nan_to_num(test.pvalue,nan=1.0,posinf=1,neginf=1)
    p=np.clip(p,0,1)
    q=multipletests(p,method="fdr_bh")[1]
    df=pd.DataFrame({"GENE":names,"DIFF_MEAN_PERTURB_MINUS_CONTROL":difference,
                     "WELCH_T":t,"P_VALUE":p,"FDR":q})
    df=df[df.GENE.ne("")].copy()
    df=df.sort_values(["GENE", "P_VALUE"],ascending=[True,True]).drop_duplicates("GENE")
    # Genes with same expression in every replicate receive no direction claim.
    df["RANK_SCORE"] = df.WELCH_T.where(df.WELCH_T.ne(0),np.sign(df.DIFF_MEAN_PERTURB_MINUS_CONTROL)*1e-8)
    return df.sort_values("RANK_SCORE",ascending=False).reset_index(drop=True)


def run_gsea(de, gmt, meta, min_size=10, max_size=500, permutation_num=1000, threads=4):
    try: import gseapy as gp
    except ImportError: raise SystemExit("Missing gseapy: install in envs/perturb (pip install gseapy)")
    available=set(de.GENE)
    usable={key:genes for key,genes in gmt.items()
            if min_size<=len(set(genes)&available)<=max_size}
    if not usable:
        return pd.DataFrame(),0
    ranked=de[["GENE","RANK_SCORE"]].dropna()
    if len(ranked)<20:
        raise ValueError("<20 unique ranked genes; inappropriate for GSEA")
    if float(ranked.RANK_SCORE.max())==float(ranked.RANK_SCORE.min()):
        raise ValueError("All rank scores are equal; no meaningful directional analysis")
    result=gp.prerank(rnk=ranked,gene_sets=usable,permutation_num=permutation_num,
                       min_size=min_size,max_size=max_size,threads=threads,
                       outdir=None,seed=42,verbose=False)
    x=result.res2d.copy()
    if x.empty: return pd.DataFrame(),len(usable)
    columns={"Term":"GMT_KEY","NES":"NES","NOM p-val":"NOMINAL_P","FDR q-val":"GSEA_FDR",
             "Lead_genes":"LEADING_EDGE_GENES"}
    x=x.rename(columns=columns)
    if "GMT_KEY" not in x:
        raise ValueError(f"Unexpected GSEApy results schema: {list(x)}")
    for c in ["NES","NOMINAL_P","GSEA_FDR"]:
        if c not in x: x[c]=np.nan
        x[c]=pd.to_numeric(x[c],errors="coerce")
    x["LIBRARY"]=[meta[k]["LIBRARY"] for k in x.GMT_KEY]
    x["PATHWAY_NAME"]=[meta[k]["NAME"] for k in x.GMT_KEY]
    x["PATHWAY_ID_OR_DESC"]=[meta[k]["DESCRIPTION"] for k in x.GMT_KEY]
    x["TRANSCRIPTIONAL_RESPONSE"] = np.select(
        [(x.GSEA_FDR<=0.05)&(x.NES>0),(x.GSEA_FDR<=0.05)&(x.NES<0)],
        ["UP_AFTER_PERTURBATION","DOWN_AFTER_PERTURBATION"],default="NOT_SIGNIFICANT")
    return x,len(usable)


def annotate_step14(gsea, step14):
    if gsea.empty: return gsea
    x=gsea.copy()
    x["STEP14_SIGNIFICANT_PATHWAY"]="NO"
    x["STEP14_TERM_ID"]=""
    x["STEP14_TERM_NAME"]=""
    def norm(s): return re.sub(r"[^a-z0-9]+"," ",str(s).lower()).strip()
    for idx,row in x.iterrows():
        name=norm(row.PATHWAY_NAME)
        description=norm(row.PATHWAY_ID_OR_DESC)
        possible=step14.loc[(step14.TERM_NAME.map(norm)==name) |
                           (step14.TERM_ID.map(norm)==name) |
                           ((step14.TERM_ID.map(norm)==description)&(description!=""))]
        if not possible.empty:
            x.at[idx,"STEP14_SIGNIFICANT_PATHWAY"]="YES"
            x.at[idx,"STEP14_TERM_ID"]=str(possible.iloc[0].TERM_ID)
            x.at[idx,"STEP14_TERM_NAME"]=str(possible.iloc[0].TERM_NAME)
    return x


def validate_manifest(m):
    necessary={"INCLUDE","DATASET_PATH","GENE_ID","PERTURBATION_COLUMN","CONTROL_VALUE","TARGET_VALUE",
               "REPLICATE_COLUMN","PERTURBATION_MODE","EXPRESSION_LAYER","EXPRESSION_KIND","GENE_SYMBOL_COLUMN"}
    missing=necessary-set(m.columns)
    if missing: raise ValueError(f"Manifest missing columns: {sorted(missing)}")
    return m.loc[m.INCLUDE.astype(str).str.upper().isin(["YES","Y","TRUE","1"])].copy()


def analyze(root, step15, out, args):
    genes, step14 = load_step15(step15)
    if not args.gmt:
        raise SystemExit("For directional pathway inference supply --gmt FULL_LIBRARY.gmt "
                         "(a complete gene-set resource, not the Step15 overlapping genes). "
                         "The Reactome ZIP is currently MISSING; use an existing valid .gmt or finish download.")
    gmt,meta=load_gmts(args.gmt)
    m=pd.read_csv(args.run_manifest,sep="\t",dtype=str).fillna("")
    m=validate_manifest(m)
    if m.empty:
        raise SystemExit("No manifest rows with INCLUDE=YES. Review controls, replicates and expression first.")
    allowed=set(genes.GENE_ID)
    out.mkdir(parents=True,exist_ok=True)
    de_dir=out/"differential_expression"
    de_dir.mkdir(parents=True,exist_ok=True)
    runs=[]; responses=[]
    for i,row in enumerate(m.to_dict("records"),1):
        p=Path(row["DATASET_PATH"])
        if not p.is_absolute(): p=root/p
        gid=clean_gene(row["GENE_ID"])
        key=f"{slug(p.stem)[:65]}_{slug(gid)}_{i:03d}"
        status={"RUN_ID":key,"DATASET_PATH":str(p),"GENE_ID":gid,
                "GENE_SYMBOL":row.get("GENE_SYMBOL",""),"STATUS":"FAILED", "ERROR":"",
                "N_PERTURB_REPLICATES":0,"N_CONTROL_REPLICATES":0,
                "N_GSEA_TESTED_PATHWAYS":0,"N_GSEA_FDR_LE_0_05":0}
        adata=None
        try:
            if gid not in allowed: raise ValueError("Gene is not among Step15 primary 24")
            if not p.is_file(): raise FileNotFoundError(p)
            perturb_col=row["PERTURBATION_COLUMN"]
            rep_col=row["REPLICATE_COLUMN"]
            mode=row["PERTURBATION_MODE"].strip().upper()
            permitted={"CRISPRI", "CRISPRA", "KNOCKOUT", "OVEREXPRESSION", "RNAI", "OTHER_SPECIFIED"}
            if mode not in permitted:
                raise ValueError(f"Explicit PERTURBATION_MODE required: {sorted(permitted)} (received {mode!r})")
            if not rep_col: raise ValueError("REPLICATE_COLUMN required; do not use individual cells as replicates")
            if not row["CONTROL_VALUE"] or not row["TARGET_VALUE"]:
                raise ValueError("CONTROL_VALUE and TARGET_VALUE must be explicitly confirmed")
            if row["TARGET_VALUE"]==row["CONTROL_VALUE"]:
                raise ValueError("Target and control are identical")
            adata=load_adata(p)
            for col in [perturb_col,rep_col]:
                if col not in adata.obs: raise ValueError(f"Missing metadata column: {col}")
            obs=adata.obs.copy()
            select=np.ones(len(obs),dtype=bool)
            if row.get("CELL_TYPE_COLUMN",""):
                cell_col=row["CELL_TYPE_COLUMN"]
                if cell_col not in obs: raise ValueError(f"Missing cell type column: {cell_col}")
                if not row.get("CELL_TYPE_VALUE",""):
                    raise ValueError("CELL_TYPE_VALUE required when CELL_TYPE_COLUMN set")
                select &= obs[cell_col].astype(str).to_numpy()==row["CELL_TYPE_VALUE"]
            labels=obs[perturb_col].astype(str).to_numpy()
            control_mask=select & (labels==row["CONTROL_VALUE"])
            target_mask=select & (labels==row["TARGET_VALUE"])
            if not control_mask.any() or not target_mask.any():
                raise ValueError(f"No cells in contrast target={int(target_mask.sum())}, control={int(control_mask.sum())}")
            obs["__Step16_group"] = np.where(target_mask,"PERTURB",np.where(control_mask,"CONTROL","OTHER"))
            names=gene_names(adata,row["GENE_SYMBOL_COLUMN"])
            expr,info=group_pseudobulk(adata,obs,"__Step16_group",rep_col,
                                      row["EXPRESSION_LAYER"],row["EXPRESSION_KIND"],args.min_cells)
            status["N_PERTURB_REPLICATES"]=int((info.CONDITION=="PERTURB").sum())
            status["N_CONTROL_REPLICATES"]=int((info.CONDITION=="CONTROL").sum())
            info.to_csv(de_dir/f"{key}_replicates.tsv",sep="\t",index=False)
            de=diff_stats(expr,info,names,row["EXPRESSION_KIND"])
            de.to_csv(de_dir/f"{key}_DE.tsv",sep="\t",index=False)
            g,n_eligible=run_gsea(de,gmt,meta,args.gsea_min_size,args.gsea_max_size,args.permutations,args.threads)
            status["N_GSEA_TESTED_PATHWAYS"]=n_eligible
            if not g.empty:
                g=annotate_step14(g,step14)
                g["RUN_ID"]=key
                g["DATASET"]=p.name
                g["GENE_ID"]=gid
                g["GENE_SYMBOL"]=row.get("GENE_SYMBOL","")
                g["CONTRAST"]=f"{row['TARGET_VALUE']} vs {row['CONTROL_VALUE']}"
                g["PERTURBATION_MODE"]=mode
                g["TARGET_ACTIVITY_DIRECTION_INFERRED_FROM_MODE"]=(
                    "EXPECTED_DECREASE_NOT_VALIDATED" if mode in {"CRISPRI", "KNOCKOUT", "RNAI"} else
                    "EXPECTED_INCREASE_NOT_VALIDATED" if mode in {"CRISPRA", "OVEREXPRESSION"} else
                    "NOT_ESTABLISHED")
                g["EXPRESSION_KIND"]=row["EXPRESSION_KIND"]
                g["CELL_TYPE"]=row.get("CELL_TYPE_VALUE","") or "UNSTRATIFIED"
                g["INFERENCE_SCOPE"]="EXPLORATORY_REPLICATE_LEVEL_EXPRESSION_RESPONSE_NOT_CAUSAL_ACTIVATION"
                responses.append(g)
                status["N_GSEA_FDR_LE_0_05"]=int((g.GSEA_FDR<=0.05).sum())
            status["STATUS"]="SUCCESS"
            print(f"[{i}/{len(m)}] {gid}: {status['STATUS']} "
                  f"reps={status['N_PERTURB_REPLICATES']}/{status['N_CONTROL_REPLICATES']} "
                  f"pathways={n_eligible} significant={status['N_GSEA_FDR_LE_0_05']}",flush=True)
        except Exception as exc:
            status["ERROR"]=f"{type(exc).__name__}: {str(exc)[:500]}"
            print(f"[{i}/{len(m)}] {gid}: {status['ERROR']}",flush=True)
        finally:
            if adata is not None:
                try: adata.file.close()
                except Exception: pass
                del adata
            gc.collect()
        runs.append(status)
        pd.DataFrame(runs).to_csv(out/"Step16_Run_Audit.tsv",sep="\t",index=False)
    merged=pd.concat(responses,ignore_index=True) if responses else pd.DataFrame(columns=["GENE_ID","GMT_KEY","NES","GSEA_FDR","STEP14_SIGNIFICANT_PATHWAY","TRANSCRIPTIONAL_RESPONSE"])
    merged.to_csv(out/"Step16_Directional_Pathway_Responses_All.tsv",sep="\t",index=False)
    selected=merged.loc[merged.STEP14_SIGNIFICANT_PATHWAY.eq("YES")].copy()
    selected.to_csv(out/"Step16_Step14_Pathway_Responses.tsv",sep="\t",index=False)
    if not selected.empty:
        # This is a plot-friendly and Step17-friendly summary, not a causal claim.
        selected.loc[:,[c for c in ["RUN_ID","GENE_ID","GENE_SYMBOL","STEP14_TERM_ID","STEP14_TERM_NAME","NES","GSEA_FDR","TRANSCRIPTIONAL_RESPONSE","CELL_TYPE","CONTRAST","PERTURBATION_MODE","TARGET_ACTIVITY_DIRECTION_INFERRED_FROM_MODE"] if c in selected]].to_csv(
            out/"Step16_Gene_Pathway_Direction_Handoff.tsv",sep="\t",index=False)
    summary={"version":VERSION,"mode":"analysis","n_manifest_rows":len(m),
             "n_success":sum(r["STATUS"]=="SUCCESS" for r in runs),
             "n_not_testable_or_failed":sum(r["STATUS"]!="SUCCESS" for r in runs),
             "n_tested_pathway_response_rows":len(merged),"n_step14_pathway_response_rows":len(selected),
             "n_significant_step14_response_rows":int((selected.GSEA_FDR<=0.05).sum()) if not selected.empty else 0,
             "interpretation":"GSEA NES is transcriptional response to experimental perturbation, not risk-allele causal direction or therapeutic effect."}
    (out/"Step16_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(f"[DONE] {summary}; output={out}")
    return summary


def stats_report(out):
    summary=out/"Step16_summary.json"
    if not summary.exists():
        raise SystemExit(f"No Step16 outputs found: {summary}. Run --inventory first.")
    info=json.loads(summary.read_text())
    print("="*100,"\nGWAS2Mech Step16 status\n","="*100,sep="")
    for k,v in info.items(): print(f"{k:<38} {v}")
    for fname,cols in [
        ("Step16_Gene_Coverage_Summary.tsv",["GENE_ID","GENE_SYMBOL","FOUND_IN_LOCAL_SC_PERTURB","N_DATASETS"]),
        ("Step16_Run_Audit.tsv",["GENE_ID","STATUS","ERROR","N_PERTURB_REPLICATES","N_CONTROL_REPLICATES","N_GSEA_TESTED_PATHWAYS"]),
        ("Step16_Step14_Pathway_Responses.tsv",["GENE_ID","STEP14_TERM_NAME","NES","GSEA_FDR","TRANSCRIPTIONAL_RESPONSE","CELL_TYPE"])]:
        p=out/fname
        if not p.exists(): continue
        data=pd.read_csv(p,sep="\t")
        selected=[c for c in cols if c in data]
        print(f"\n{fname}: {len(data)} rows")
        if data.empty or not selected: continue
        if "GSEA_FDR" in data: data=data.sort_values("GSEA_FDR")
        print(data[selected].head(35).to_string(index=False))


def main():
    ap=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description=__doc__)
    ap.add_argument("--root",default=".")
    ap.add_argument("--phenotype",required=True)
    ap.add_argument("--ancestry",required=True)
    ap.add_argument("--reference-panel",default="1000G30x")
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--inventory",action="store_true",help="Index local h5ad perturbations; no pathway statistics")
    mode.add_argument("--run-manifest",help="Use manually reviewed manifest with INCLUDE=YES, replicates and controls")
    mode.add_argument("--stats",action="store_true")
    ap.add_argument("--data-root",help="Override resources/perturbation/scperturb")
    ap.add_argument("--limit-files",type=int,default=0)
    ap.add_argument("--gmt",action="append",default=[],help="Complete full-library GMT, .gmt.gz or .gmt.zip; repeatable")
    ap.add_argument("--min-cells",type=int,default=10,help="Minimum cells per biological replicate")
    ap.add_argument("--gsea-min-size",type=int,default=10)
    ap.add_argument("--gsea-max-size",type=int,default=500)
    ap.add_argument("--permutations",type=int,default=1000)
    ap.add_argument("--threads",type=int,default=4)
    args=ap.parse_args()
    root,step15,out=get_paths(args)
    if args.stats: stats_report(out)
    elif args.inventory: inventory(root,step15,out,args)
    else: analyze(root,step15,out,args)


if __name__=="__main__":
    main()
