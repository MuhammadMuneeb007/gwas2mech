#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GWAS2Mech Step16: public gene-perturbation -> transcriptional pathway responses.

Run from GWAS2Mech root with envs/perturb/bin/python:
    python Step16_Perturbation_Directional_Pathways.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --check-resources
    python Step16_Perturbation_Directional_Pathways.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --slurm
    sbatch Step16_Perturb_parkinson_s_disease_eur_1000G30x.sh
    python Step16_Perturbation_Directional_Pathways.py --phenotype "Parkinson's disease" --ancestry EUR --reference-panel 1000G30x --stats

No manual manifest. Use --prepare-libraries once to cache source-specific GO BP /
WikiPathways / KEGG full GMT libraries before submitting --slurm.
The external library versions are audited, NOT assumed equivalent to Step14
g:Profiler annotation releases. Only stable, same-source term IDs count as matches.
Reads Step15_Gene_Pathway_Handoff.tsv and
Step15_Gene_Mapping_Audit.tsv. Scans locally completed scPerturb AnnData files.
Handles H5AD dense and CSR/CSC sparse matrices with bounded-memory selected-row
processing, and mixed .obs/.var storage incl. pandas categoricals in h5py.

IMPORTANT BIOLOGICAL INTERPRETATION:
    Default analysis is descriptive. Cells, guides and technical batches are NOT
    independent biological replicates. Without verified independent biological
    replicates, a gene-set-permutation FDR is NOT an experiment-level FDR.
    Expression scale is inferred from sampled matrix values and flagged as
    unverified, not presumed to be normalized/raw based on file names.
    No claim about disease risk allele, pathway biochemical activity, causal
    effector genes or therapeutic effect is made.
    LINCS Level-5 is inventoried separately; it is not treated as raw cells.

Large downloads still in progress are skipped if HDF5 cannot open them.
"""
from __future__ import annotations
import argparse
import gc
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import sys
import time
import zipfile
import urllib.request
import urllib.parse
from datetime import datetime, timezone
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

VERSION = "3.0.0-source-matched-gsea"
CONTROL = {"control", "ctrl", "non-targeting", "nontargeting", "non targeting", "ntc",
           "negative control", "non-targeting control", "scramble", "scrambled", "neg ctrl"}
PERTURB_COLS = ["gene", "perturbation", "gene_target", "target_gene", "perturbation_name",
                "target", "perturbation_target", "perturbed_gene", "gene_id"]
CELL_TYPE_COLS = ["celltype", "cell_type", "cell_line", "cell_type_label"]
BIOREP_COLS = ["biological_replicate", "bioreplicate", "donor_id", "donor", "sample_id"]
OUTPUT_COLS = ["DATASET", "DATASET_PATH", "GENE_ID", "GENE_SYMBOL", "TARGET_VALUE",
               "TARGET_CELLS", "CONTROL_VALUE", "CONTROL_CELLS", "PERTURBATION_COLUMN",
               "CELL_TYPE", "EXPERIMENT_MODE"]

def put(df, path, cols=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(df, list): df = pd.DataFrame(df)
    if cols is not None: df = df.reindex(columns=cols)
    df.to_csv(path, sep="\t", index=False)


def label(x):
    if x is None: return ""
    if isinstance(x, bytes): x = x.decode("utf-8", "replace")
    if isinstance(x, (float, np.floating)) and not np.isfinite(x): return ""
    s = str(x).strip()
    return "" if s.lower() in {"nan", "none", "null", "na", "<na>"} else s


def gene_id(x):
    s = label(x)
    m = re.search(r"ENSG\d+", s, re.I)
    return m.group(0).upper() if m else s.upper()


def slug(x): return re.sub(r"[^a-z0-9]+", "_", str(x).lower()).strip("_")


def norm(x): return re.sub(r"[^a-z0-9]+", " ", label(x).lower()).strip()


def is_control(x): return norm(x) in {norm(v) for v in CONTROL}


def paths(args):
    root = Path(args.root).expanduser().resolve()
    ancestry = {"eur":"european", "afr":"african", "eas":"east_asian",
                "sas":"south_asian", "amr":"hispanic_or_latin_american"}.get(args.ancestry.lower())
    if not ancestry: raise SystemExit("Invalid ancestry: expected EUR/AFR/EAS/SAS/AMR")
    rel = Path(slug(args.phenotype)) / (ancestry + ("" if args.reference_panel == "1000G" else "_" + args.reference_panel))
    if args.index is not None:
        manifest = root / "01_gwas_catalog" / slug(args.phenotype) / ancestry / "GWAS_download_manifest.tsv"
        m = pd.read_csv(manifest, sep="\t", dtype=str).fillna("")
        hit = m[pd.to_numeric(m.ARRAY_ID, errors="coerce") == args.index]
        if hit.empty: raise SystemExit(f"Study index {args.index} not found in {manifest}")
        rel /= str(hit.iloc[0].STUDY_ACCESSION)
    return (root, root/"15_network"/rel, root/"16_perturbation"/rel,
            root/"resources"/"perturbation"/"scperturb")


def inputs(step15):
    a = step15 / "Step15_Gene_Mapping_Audit.tsv"
    h = step15 / "Step15_Gene_Pathway_Handoff.tsv"
    if not h.is_file(): raise FileNotFoundError(f"Required Step15 input absent: {h}")
    handoff = pd.read_csv(h, sep="\t", dtype=str).fillna("")
    if not {'GENE_ID', 'GENE_SYMBOL'}.issubset(handoff.columns):
        raise ValueError(f"Step15 handoff requires GENE_ID and GENE_SYMBOL; got {list(handoff.columns)}")
    handoff['GENE_ID'] = handoff.GENE_ID.map(gene_id)
    if a.is_file():
        audit = pd.read_csv(a, sep="\t", dtype=str).fillna("")
        audit = audit[audit.IS_GWAS2M_PRIMARY.str.upper().eq('YES')]
        audit['GENE_ID'] = audit.GENE_ID.map(gene_id)
        genes = audit[['GENE_ID','GENE_SYMBOL']].drop_duplicates('GENE_ID').copy()
    else:
        genes = handoff[['GENE_ID','GENE_SYMBOL']].drop_duplicates('GENE_ID').copy()
    genes = genes[genes.GENE_ID.str.fullmatch(r'ENSG\d+')].reset_index(drop=True)
    symbol_map = handoff[['GENE_ID','GENE_SYMBOL']].drop_duplicates('GENE_ID').set_index('GENE_ID').GENE_SYMBOL.to_dict()
    genes['GENE_SYMBOL'] = [label(s) or label(symbol_map.get(g)) for g,s in genes.itertuples(index=False, name=None)]
    if genes.empty: raise ValueError("No primary Step15 genes found")
    for c in ('TERM_ID','TERM_NAME','SOURCE','FDR'):
        if c not in handoff: handoff[c] = ''
    pathway = handoff[['TERM_ID','TERM_NAME','SOURCE','FDR']].drop_duplicates()
    pathway = pathway[(pathway.TERM_NAME!='')|(pathway.TERM_ID!='')]
    return genes, pathway


# Enrichr supplies full gene-set libraries in GMT format.  They can differ from
# g:Profiler annotations used in Step14; this difference is ALWAYS recorded.
# The library list is intentionally frozen in code so reruns are reproducible.
LIBRARY_CANDIDATES = {
    'GO:BP': ['GO_Biological_Process_2023'],
    'WP': ['WikiPathway_2021_Human', 'WikiPathways_2024_Human'],
    'KEGG': ['KEGG_2021_Human'],
}


def canon_source(x):
    v=label(x).upper().replace(' ','')
    if v in ('GO:BP','GOBP','GO_BIOLOGICAL_PROCESS','GO-BP'): return 'GO:BP'
    if v in ('WIKIPATHWAYS','WIKIPATHWAY','WP'): return 'WP'
    if v in ('KEGG','KEGG_LEGACY','KEGG_MEDICUS'): return 'KEGG'
    if v in ('REAC','REACTOME','REACTOMEPATHWAYS'): return 'REAC'
    return v


def native_term_id(raw, source):
    raw=label(raw)
    if source=='GO:BP':
        m=re.search(r'GO\s*:\s*(\d{7})',raw,re.I)
        return 'GO:'+m.group(1) if m else ''
    if source=='WP':
        m=re.search(r'\bWP\s*[-_: ]?\s*(\d{1,6})\b',raw,re.I)
        return 'WP'+m.group(1) if m else ''
    if source=='KEGG':
        m=re.search(r'(?:\bhsa|\bKEGG\s*:\s*)(\d{5})\b',raw,re.I)
        return 'hsa'+m.group(1) if m else ''
    if source=='REAC':
        m=re.search(r'R-HSA-\d+',raw,re.I)
        return m.group(0).upper() if m else ''
    return ''


def present_sources(step14):
    return sorted({canon_source(x) for x in step14.SOURCE if canon_source(x) in LIBRARY_CANDIDATES})


def library_dir(root,args):
    return Path(args.library_dir).expanduser().resolve() if args.library_dir else root/'resources/pathways/step16_source_matched'


def checksum(p):
    h=hashlib.sha256()
    with p.open('rb') as handle:
        while 1:
            block=handle.read(1024*1024)
            if not block:break
            h.update(block)
    return h.hexdigest()


def count_gmt_rows(path):
    # Lightweight validity check; does NOT load millions of gene-set members.
    if not path.is_file() or path.stat().st_size<100: return 0
    with path.open('rt',encoding='utf8',errors='replace') as handle:
        return sum(1 for line in handle if len(line.rstrip('\r\n').split('\t'))>=4)


def acquire_libraries(root,args,step14,allow_download=False):
    dest=library_dir(root,args); dest.mkdir(parents=True,exist_ok=True)
    result=[]
    for source in present_sources(step14):
        choices=LIBRARY_CANDIDATES[source]
        preferred=choices[0]
        existing=next((dest/(x+'.gmt') for x in choices if count_gmt_rows(dest/(x+'.gmt'))>10),None)
        record={'SOURCE':source,'LIBRARY_NAME':existing.stem if existing else preferred,
                'PATH':str(existing or dest/(preferred+'.gmt')),
                'STATUS':'READY' if existing else 'MISSING',
                'GENESET_ROWS':count_gmt_rows(existing) if existing else 0,
                'SHA256':checksum(existing) if existing else '',
                'DOWNLOAD_URL':'','ERROR':''}
        if existing is None and allow_download:
            err=[]
            for lib in choices:
                url='https://maayanlab.cloud/Enrichr/geneSetLibrary?'+urllib.parse.urlencode(
                    {'mode':'text','libraryName':lib})
                output=dest/(lib+'.gmt');tmp=dest/(lib+'.gmt.partial')
                record['DOWNLOAD_URL']=url
                try:
                    request=urllib.request.Request(url,headers={'User-Agent':'GWAS2Mech-Step16/3.0'})
                    with urllib.request.urlopen(request,timeout=100) as response, tmp.open('wb') as sink:
                        while True:
                            block=response.read(1024*1024)
                            if not block:break
                            sink.write(block)
                    if count_gmt_rows(tmp)<10:
                        raise ValueError('Fetched resource not a valid GMT library (>10 sets required)')
                    tmp.replace(output)
                    record.update(LIBRARY_NAME=lib,PATH=str(output),STATUS='READY',
                                  GENESET_ROWS=count_gmt_rows(output),SHA256=checksum(output),ERROR='')
                    break
                except Exception as e:
                    err.append(f'{lib}: {type(e).__name__}: {e}')
                    tmp.unlink(missing_ok=True)
            if record['STATUS']!='READY':record['ERROR']='; '.join(err)[:900]
        elif existing is None:
            record['DOWNLOAD_URL']='https://maayanlab.cloud/Enrichr/geneSetLibrary?'+urllib.parse.urlencode(
                {'mode':'text','libraryName':preferred})
        result.append(record)
    return result


def load_libraries(records,out):
    # GMT rows: name TAB description TAB members... ; membership MUST use
    # full, independently curated gene sets, never the Step14 intersection genes.
    all_sets={}; meta={}; audit=[]
    for item in records:
        source=item['SOURCE']; file=Path(item['PATH'])
        if item['STATUS']!='READY':
            audit.append({'SOURCE':source,'STATUS':'LIBRARY_MISSING','LIBRARY_NAME':item['LIBRARY_NAME']})
            continue
        n=0; identified=0; duplicate=0
        with file.open('rt',encoding='utf8',errors='replace') as handle:
            for lnum,line in enumerate(handle,1):
                fields=line.rstrip('\n\r').split('\t')
                if len(fields)<4:continue
                display,description=fields[:2]
                ident=native_term_id(display,source) or native_term_id(description,source)
                if not ident:
                    # Name-only libraries can still be screened, but their
                    # results cannot automatically claim a Step14 ID match.
                    ident='NAME:'+hashlib.sha1(display.encode()).hexdigest()[:16]
                else:identified+=1
                key=source+'|'+ident
                members={gene_id(g) for g in fields[2:] if label(g)}
                if len(members)<5:continue
                if key in all_sets:
                    duplicate+=1
                    if all_sets[key]!=sorted(members):
                        audit.append({'SOURCE':source,'STATUS':'CONFLICTING_DUPLICATE_TERM',
                                      'TERM_ID':ident,'TERM_NAME':display,'LIBRARY_NAME':item['LIBRARY_NAME']})
                    continue
                all_sets[key]=sorted(members)
                meta[key]={'PATHWAY_NAME':display,'PATHWAY_DESCRIPTION':description,
                           'PATHWAY_SOURCE':source,'PATHWAY_ID':ident if not ident.startswith('NAME:') else '',
                           'LIBRARY_NAME':item['LIBRARY_NAME'], 'LIBRARY_SHA256':item['SHA256'],
                           'PATHWAY_GENE_SET_SIZE':len(members)}
                n+=1
        audit.append({'SOURCE':source,'STATUS':'LOADED','LIBRARY_NAME':item['LIBRARY_NAME'],
                      'N_LOADED_SETS':n,'N_ROWS_WITH_NATIVE_ID':identified,'N_DUPLICATES':duplicate,
                      'SHA256':item['SHA256']})
    put(audit,out/'Step16_GeneSet_Definition_Audit.tsv')
    return all_sets,meta


def pathway_mapping_audit(meta,step14,out):
    # This reports all of Step14, even terms lacking a gene-set definition.
    keys={(m['PATHWAY_SOURCE'],m['PATHWAY_ID']) for m in meta.values() if m['PATHWAY_ID']}
    bysource=defaultdict(list)
    for m in meta.values(): bysource[m['PATHWAY_SOURCE']].append(m)
    rows=[]
    for r in step14.to_dict('records'):
        source=canon_source(r.get('SOURCE')); ident=native_term_id(r.get('TERM_ID',''),source)
        names=[m for m in bysource[source] if norm(m['PATHWAY_NAME'])==norm(r.get('TERM_NAME',''))]
        exact=(source,ident) in keys if ident else False
        rows.append({'STEP14_SOURCE':r.get('SOURCE',''),'STEP14_TERM_ID':r.get('TERM_ID',''),
          'STEP14_TERM_NAME':r.get('TERM_NAME',''),'STEP14_FDR':r.get('FDR',''),
          'MAPPING_STATUS':'EXACT_SOURCE_NATIVE_ID' if exact else
                           'SAME_SOURCE_NAME_ONLY_NOT_ACCEPTED' if names else 'UNMAPPED',
          'SOURCE':source,'TERM_ID':ident,'N_NAME_ONLY_CANDIDATES':len(names),
          'LIBRARY_NOTE':'Gene-set definitions may differ by enrichment provider/release'})
    put(rows,out/'Step16_Step14_Pathway_Mapping_Audit.tsv')
    return rows


def load_gmt(path):
    """Legacy Reactome loader retained for optional use, not used by default v3."""
    if path is None:return {},{}
    if path.suffix=='.zip':
        with zipfile.ZipFile(path) as z:
            member=next((s for s in z.namelist() if s.lower().endswith('.gmt')),None)
            if not member:raise ValueError(f'No GMT in {path}')
            lines=io.TextIOWrapper(z.open(member),encoding='utf8').read().splitlines()
    elif path.suffix=='.gz':
        with gzip.open(path,'rt') as handle:lines=handle.read().splitlines()
    else:lines=path.read_text().splitlines()
    sets={};meta={}
    for line in lines:
        t=line.split('\t')
        if len(t)<4:continue
        genes=sorted({gene_id(x) for x in t[2:] if label(x)})
        if len(genes)<5:continue
        key='REAC|'+(native_term_id(t[0],'REAC') or t[0])
        sets[key]=genes
        meta[key]={'PATHWAY_NAME':t[0],'PATHWAY_DESCRIPTION':t[1],
                   'PATHWAY_SOURCE':'REAC','PATHWAY_ID':native_term_id(t[0],'REAC'),
                   'LIBRARY_NAME':str(path),'LIBRARY_SHA256':checksum(path),
                   'PATHWAY_GENE_SET_SIZE':len(genes)}
    return sets,meta


def as_str_array(values):
    a=np.asarray(values)
    if a.dtype.kind=='S': return np.char.decode(a,'utf8',errors='replace').astype(str)
    if a.dtype.kind in {'U'}: return a.astype(str)
    if a.dtype.kind=='O': return np.array([label(v) for v in a],dtype=object)
    return np.array([label(v) for v in a],dtype=object)


def column(h5group, name):
    """Read AnnData HDF5 obs/var column; categorical codes -> actual labels."""
    if name not in h5group: raise KeyError(name)
    node = h5group[name]
    if isinstance(node, h5py.Dataset): return as_str_array(node[:])
    if isinstance(node,h5py.Group) and 'codes' in node and 'categories' in node:
        codes=np.asarray(node['codes'][:], dtype=np.int64)
        cats=as_str_array(node['categories'][:])
        result=np.empty(len(codes),dtype=object)
        valid=(codes>=0)&(codes<len(cats))
        result[~valid]=''
        result[valid]=cats[codes[valid]]
        return result
    raise TypeError(f"Unsupported H5AD field {name}: {type(node).__name__}")


def headers(rootgroup):
    return [label(k) for k in rootgroup.keys()]


def matrix_meta(f):
    if 'X' not in f: return 'MISSING_X', (0,0)
    x=f['X']
    shape=tuple(x.shape) if isinstance(x,h5py.Dataset) else tuple(int(v) for v in x.attrs.get('shape',[]))
    enc=label(x.attrs.get('encoding-type', 'array'))
    return enc,shape


def detect_experiment(path, f, aliases, min_cells):
    obs=f['obs']; cols=headers(obs)
    best=None
    for pcol in PERTURB_COLS:
        if pcol not in obs: continue
        try:
            if isinstance(obs[pcol],h5py.Group) and 'categories' in obs[pcol]:
                cats=as_str_array(obs[pcol]['categories'][:]); keys={gene_id(x) for x in cats}
                if not (set(aliases)&keys): continue
            labels=column(obs,pcol)
        except Exception: continue
        counts=Counter(labels)
        controls=[k for k,v in counts.items() if is_control(k) and v>=min_cells]
        matches=[(k,aliases[gene_id(k)],int(v)) for k,v in counts.items()
                 if gene_id(k) in aliases and v>=min_cells]
        if not matches: continue
        score=(bool(controls),len(matches),sum(z[2] for z in matches),-PERTURB_COLS.index(pcol))
        if best is None or score>best[0]: best=(score,pcol,labels,counts,controls,matches)
    return best


def read_index_or_var(f):
    group=f['var']; n=matrix_meta(f)[1][1]
    for col in ('gene_name','gene_symbol','gene_symbols','symbol','_index'):
        if col in group:
            x=column(group,col)
            if len(x)==n: return x,col
    idx=label(group.attrs.get('_index', ''))
    if idx in group:
        x=column(group,idx)
        if len(x)==n: return x,idx
    for col in ('ensembl_id','gene_id'):
        if col in group:
            x=column(group,col)
            if len(x)==n: return x,col
    raise ValueError('No usable gene labels in /var')


def expression_scale(f, n_obs, n_genes):
    """Inspect actual values, never assume all h5ad X is counts."""
    x=f['X']
    try:
        if isinstance(x,h5py.Dataset):
            values=np.asarray(x[:min(32,n_obs), :min(256,n_genes)],dtype=float).ravel()
        elif isinstance(x,h5py.Group) and 'data' in x:
            values=np.asarray(x['data'][:min(100000,len(x['data']))],dtype=float).ravel()
        else: return 'UNSUPPORTED','Unknown H5AD X encoding'
        finite=values[np.isfinite(values)]
        if len(finite)==0 or np.any(finite<0): return 'UNVERIFIED','Unknown/negative matrix scale; use observed-unit differences only'
        positive=finite[finite>0]
        if len(positive)==0: return 'UNVERIFIED','All sampled entries zero; no scale inference'
        fraction=np.mean(np.abs(positive-np.rint(positive))<1e-5)
        if fraction>=0.999:
            return 'COUNT_LIKE_UNVERIFIED',f'{fraction:.3f} fraction sampled positive entries integer-like; verify study provenance'
        return 'NONINTEGER_UNVERIFIED',f'{fraction:.3f} fraction integer-like; no raw-count/log interpretation without study provenance'
    except Exception as e:
        return 'UNVERIFIED',f'Expression scale sampling error: {type(e).__name__}: {e}'


def draw_controls(target_idxs, control_idxs, obs, max_ctrl, seed):
    """Prefer controls from batches containing targets; stratified deterministic subsampling."""
    rng=np.random.default_rng(seed)
    batch='batch' if 'batch' in obs else None
    note='UNPAIRED_CONTROL_SUBSAMPLE'
    if batch:
        b=column(obs,batch)
        t_batches=set(b[target_idxs]); same=control_idxs[np.isin(b[control_idxs],list(t_batches))]
        if len(same)>=10:
            control_idxs=same; note='TECHNICAL_BATCH_MATCHED_NOT_BIOLOGICAL_REPLICATES'
        else: note='NO_SHARED_TECHNICAL_BATCH_CONTROL_SAMPLED'
    if len(control_idxs)>max_ctrl:
        if batch and 'b' in locals() and note.startswith('TECHNICAL'):
            pieces=[]
            uniq,cnt=np.unique(b[target_idxs],return_counts=True)
            for u,n in zip(uniq,cnt):
                ids=control_idxs[b[control_idxs]==u]
                take=min(len(ids), max(1, int(round(max_ctrl * n / cnt.sum()))))
                if take: pieces.extend(rng.choice(ids,size=take,replace=False).tolist())
            control_idxs=np.array(sorted(set(pieces)),dtype=int)
        if len(control_idxs)>max_ctrl:
            control_idxs=np.sort(rng.choice(control_idxs,size=max_ctrl,replace=False))
    return np.sort(control_idxs),note


def extract_selected(f, indices, size_hint=512):
    """Return expression for selected cells only; reads bounded dense blocks or CSC/CSR."""
    x=f['X']; indices=np.array(sorted(set(map(int,indices))),dtype=int)
    if len(indices)==0: raise ValueError('No selected cells')
    if isinstance(x,h5py.Dataset):
        # h5py selects sorted noncontiguous rows without materialising full X.
        blocks=[]
        for batch in np.array_split(indices, max(1,int(np.ceil(len(indices)/size_hint)))):
            if not len(batch): continue
            blocks.append(np.asarray(x[batch.tolist(),:],dtype=np.float32))
        return np.concatenate(blocks,axis=0)
    if isinstance(x,h5py.Group):
        enc=label(x.attrs.get('encoding-type', '')).lower()
        if enc not in {'csr_matrix','csc_matrix'}:
            raise ValueError(f'Unsupported sparse H5AD encoding: {enc}')
        shape=tuple(int(y) for y in x.attrs['shape'])
        data=x['data'][:]; ind=x['indices'][:]; ptr=x['indptr'][:]
        cls=sparse.csr_matrix if enc=='csr_matrix' else sparse.csc_matrix
        full=cls((data,ind,ptr),shape=shape)
        picked=full[indices,:].toarray().astype(np.float32)
        del full,data,ind,ptr
        return picked
    raise ValueError('Unsupported /X type')


def normalize_and_rank(X, which, scale, max_genes=0):
    # Only count-like matrices get library-size normalization (still provenance-unverified).
    X=np.asarray(X,dtype=np.float64)
    if not np.isfinite(X).all(): raise ValueError('NaN/Inf expression values')
    if X.min()<0: raise ValueError('Negative values in expression matrix')
    if scale=='COUNT_LIKE_UNVERIFIED':
        lib=X.sum(axis=1)
        if np.any(lib<=0): raise ValueError('Zero library-size cell in selected expression data')
        X=np.log2(1+1e4*X/lib[:,None]); transform='log2(1+CP10K), COUNT_LIKE_UNVERIFIED'
    else: transform='AS_STORED_EXPRESSION_VALUES_SCALE_UNVERIFIED'
    score=X[which==1].mean(axis=0)-X[which==0].mean(axis=0)
    return score,transform


def do_gsea(ranked, gmt, gmeta, args):
    try: import gseapy as gp
    except ImportError as e:raise RuntimeError('gseapy missing from envs/perturb') from e
    present=set(ranked.GENE)
    usable={k:members for k,members in gmt.items()
            if args.gsea_min_size<=len(set(members)&present)<=args.gsea_max_size}
    if not usable:return pd.DataFrame(),0
    if ranked.RANK_SCORE.nunique()<4:raise ValueError('Too little expression-score variation for GSEA')
    result=gp.prerank(rnk=ranked[['GENE','RANK_SCORE']],gene_sets=usable,
                     permutation_num=args.permutations,min_size=args.gsea_min_size,
                     max_size=args.gsea_max_size,threads=args.threads,
                     seed=args.seed,outdir=None,verbose=False)
    df=result.res2d.rename(columns={'Term':'PATHWAY_KEY','NES':'NES','NOM p-val':'GENE_SET_NOMINAL_P',
                       'FDR q-val':'GENE_SET_PERMUTATION_FDR','Lead_genes':'LEADING_EDGE_GENES'}).copy()
    for c in ('NES','GENE_SET_NOMINAL_P','GENE_SET_PERMUTATION_FDR'):
        if c in df:df[c]=pd.to_numeric(df[c],errors='coerce')
    for c in ('PATHWAY_NAME','PATHWAY_DESCRIPTION','PATHWAY_SOURCE','PATHWAY_ID',
              'LIBRARY_NAME','LIBRARY_SHA256','PATHWAY_GENE_SET_SIZE'):
        df[c]=[gmeta[k].get(c,'') for k in df.PATHWAY_KEY]
    df['RESPONSE_DIRECTION']=np.select([df.NES>0,df.NES<0],
        ['TRANSCRIPTIONAL_SIGNATURE_UP_AFTER_PERTURBATION','TRANSCRIPTIONAL_SIGNATURE_DOWN_AFTER_PERTURBATION'],
        default='NO_DIRECTION')
    df['EXPERIMENT_LEVEL_FDR']=np.nan
    df['EVIDENCE_CLASS']='DESCRIPTIVE_EXPLORATORY_UNREPLICATED_OR_DESIGN_UNVERIFIED'
    df['FDR_SCOPE']='GENE_SET_PERMUTATION_WITHIN_TESTED_LIBRARY_UNREPLICATED'
    return df,len(usable)


def add_pathway_match(df,step14):
    if df.empty:return df
    lookup={}
    same_names=defaultdict(list)
    for r in step14.to_dict('records'):
        src=canon_source(r.get('SOURCE'))
        ident=native_term_id(r.get('TERM_ID',''),src)
        if src and ident:lookup[(src,ident)]=r
        same_names[(src,norm(r.get('TERM_NAME','')))].append(r)
    identifiers=[];names=[];marks=[];types=[];pvalues=[];matched_sources=[]
    for x in df.to_dict('records'):
        src=canon_source(x.get('PATHWAY_SOURCE'))
        ident=native_term_id(x.get('PATHWAY_ID'),src)
        match=lookup.get((src,ident)) if ident else None
        if match is not None:
            markers='YES';kind='EXACT_SOURCE_NATIVE_ID'
        else:
            candidates=same_names[(src,norm(x.get('PATHWAY_NAME')))]
            markers='NO';kind='SAME_SOURCE_NAME_ONLY_NOT_ACCEPTED' if candidates else 'UNMAPPED'
        identifiers.append(label(match.get('TERM_ID')) if match else '')
        names.append(label(match.get('TERM_NAME')) if match else '')
        pvalues.append(label(match.get('FDR')) if match else '')
        matched_sources.append(label(match.get('SOURCE')) if match else '')
        types.append(kind);marks.append(markers)
    df['STEP14_TERM_ID']=identifiers
    df['STEP14_TERM_NAME']=names
    df['STEP14_SOURCE']=matched_sources
    df['STEP14_FDR']=pvalues
    df['STEP14_PATHWAY_MATCH']=marks
    df['PATHWAY_MATCH_METHOD']=types
    return df


def lincs_metadata(root,genes,out):
    """Signature coverage only; Level-5 MODZ needs its own separate analysis."""
    folder=root/'resources/perturbation/lincs'
    files=sorted(folder.rglob('*_sig_info.txt.gz')) if folder.exists() else []
    alias={gene_id(g):g for g in genes.GENE_ID}
    for g,s in genes[['GENE_ID','GENE_SYMBOL']].itertuples(index=False,name=None):
        if label(s): alias[gene_id(s)]=g
    counts=Counter()
    rows=[]
    for f in files:
        try:
            head=pd.read_csv(f,sep='\t',compression='gzip',nrows=1)
            pc=next((x for x in ('pert_iname','pert_name') if x in head.columns),None)
            tc='pert_type' if 'pert_type' in head.columns else None
            if not pc or not tc: raise ValueError('Cannot identify pert_iname/pert_type columns')
            n=0
            for chunk in pd.read_csv(f,sep='\t',compression='gzip',usecols=[pc,tc],dtype=str,chunksize=100000):
                types=chunk[tc].fillna('').str.lower().str.contains(r'trt_sh|trt_oe|trt_xpr|trt_si|trt_crispr|knock|shrna',regex=True)
                ids=chunk.loc[types,pc].fillna('').map(gene_id).map(alias)
                vc=ids.value_counts()
                for k,v in vc.items(): counts[k]+=int(v); n+=int(v)
            rows.append({'FILE':str(f),'STATUS':'READ_OK','GENETIC_SIGNATURES_MATCHED':n})
        except Exception as e: rows.append({'FILE':str(f),'STATUS':'FAILED','ERROR':str(e)})
    put(rows,out/'Step16_LINCS_Metadata_Audit.tsv')
    put([{'GENE_ID':g,'GENE_SYMBOL':s,'N_GENETIC_L1000_SIGNATURES':int(counts[g]),
          'L1000_PATHWAY_DIRECTION_STATUS':'NOT_COMPUTED_BY_SCPERTURB_BACKEND'}
         for g,s in genes[['GENE_ID','GENE_SYMBOL']].itertuples(index=False,name=None)],
         out/'Step16_LINCS_Gene_Coverage.tsv')
    return {'n_lincs_genetic_signatures':sum(counts.values()),
            'n_genes_with_lincs_signatures':sum(v>0 for v in counts.values())}


def scan(root,step15,out,data_dir,args):
    genes, step14=inputs(step15)
    aliases={}
    for gid,s in genes[['GENE_ID','GENE_SYMBOL']].itertuples(index=False,name=None):
        aliases[gene_id(gid)]=(gid,s)
        if label(s): aliases[gene_id(s)]=(gid,s)
    base=Path(args.data_root).expanduser() if args.data_root else data_dir
    if not base.is_absolute():base=root/base
    files=sorted(base.rglob('*.h5ad')) if base.exists() else []
    if args.limit_files: files=files[:args.limit_files]
    audit=[];hits=[];schema=[];plans=[]
    for i,p in enumerate(files,1):
        record={'DATASET':p.name,'DATASET_PATH':str(p),'STATUS':'','N_CELLS':0,'N_GENES':0,
                'N_MATCHED_PRIMARY_TARGETS':0,'ERROR':''}
        opened=False
        try:
            with h5py.File(p,'r') as f:
                opened=True
                enc,shape=matrix_meta(f)
                if len(shape)!=2 or 'obs' not in f or 'var' not in f: raise ValueError(f'Invalid H5AD root, X={shape}')
                record.update(N_CELLS=int(shape[0]),N_GENES=int(shape[1]))
                scale,evidence=expression_scale(f,*shape)
                obscols=headers(f['obs']); varcols=headers(f['var'])
                schema.append({'DATASET':p.name,'X_ENCODING':enc,'N_CELLS':shape[0],'N_GENES':shape[1],
                               'OBS_HEADERS':';'.join(obscols),'VAR_HEADERS':';'.join(varcols),
                               'EXPRESSION_SCALE':scale,'SCALE_EVIDENCE':evidence})
                best=detect_experiment(p,f,aliases,args.min_cells)
                if best:
                    _,pcol,labels,counts,controls,matches=best
                    control=max(controls,key=lambda x:counts[x]) if controls else ''
                    typ=''
                    for cc in CELL_TYPE_COLS:
                        if cc in f['obs']:
                            tt=column(f['obs'],cc)
                            nonempty=Counter(x for x in tt if x)
                            typ=nonempty.most_common(1)[0][0] if nonempty else ''
                            break
                    mode='CRISPRI' if 'crispri' in p.name.lower() else 'CRISPRA' if 'crispra' in p.name.lower() else 'CRISPR_TYPE_UNVERIFIED'
                    record['N_MATCHED_PRIMARY_TARGETS']=len({m[1][0] for m in matches})
                    for value,(gid,sym),n in matches:
                        entry={'DATASET':p.name,'DATASET_PATH':str(p),'GENE_ID':gid,'GENE_SYMBOL':sym,
                               'TARGET_VALUE':label(value),'TARGET_CELLS':int(n),
                               'CONTROL_VALUE':label(control),'CONTROL_CELLS':int(counts.get(control,0)),
                               'PERTURBATION_COLUMN':pcol,'CELL_TYPE':typ,'EXPERIMENT_MODE':mode}
                        hits.append(entry)
                        status='READY_FOR_DESCRIPTIVE' if control else 'SKIP_NO_NEGATIVE_CONTROL'
                        if not control: reason='No unambiguous negative control in same perturbation column'
                        elif n<args.min_cells: status='SKIP_FEW_TARGET_CELLS';reason='Too few perturbed cells'
                        else: reason='No biological-replicate assumption; descriptive contrast only'
                        plans.append({**entry,'AUTO_STATUS':status,'REASON':reason,
                                      'EXPRESSION_SCALE':scale,'SCALE_EVIDENCE':evidence,
                                      'INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED':'NO'})
                record['STATUS']='READ_OK'
        except Exception as e:
            record['STATUS']='METADATA_OR_SCHEMA_FAILED' if opened else 'H5AD_READ_FAILED'
            record['ERROR']=f'{type(e).__name__}: {str(e)[:350]}'
        audit.append(record)
        print(f"[{i:02d}/{len(files):02d}] {p.name}: {record['STATUS']} matched={record['N_MATCHED_PRIMARY_TARGETS']} "+record['ERROR'],flush=True)
    put(audit,out/'Step16_Dataset_Inventory.tsv')
    put(schema,out/'Step16_H5AD_Storage_Audit.tsv')
    put(hits,out/'Step16_Perturbation_Coverage.tsv',OUTPUT_COLS)
    put(plans,out/'Step16_Auto_Analysis_Plan.tsv')
    coverage=[]
    for gid,sym in genes[['GENE_ID','GENE_SYMBOL']].itertuples(index=False,name=None):
        matched=[x for x in hits if x['GENE_ID']==gid]
        coverage.append({'GENE_ID':gid,'GENE_SYMBOL':sym,'FOUND_IN_LOCAL_SC_PERTURB':'YES' if matched else 'NO',
                         'N_DATASETS':len({x['DATASET'] for x in matched}),
                         'N_CONTRASTS':len(matched),
                         'N_TARGET_CELLS':sum(x['TARGET_CELLS'] for x in matched)})
    put(coverage,out/'Step16_Gene_Coverage_Summary.tsv')
    print('\n[PRIMARY GENE COVERAGE]',flush=True)
    for g in coverage:
        print(f"  {g['GENE_SYMBOL'] or g['GENE_ID']:<23} datasets={g['N_DATASETS']:<3} target_cells={g['N_TARGET_CELLS']}",flush=True)
    return genes,step14,plans,{'n_scanned_datasets':len(files),
            'n_readable_datasets':sum(r['STATUS']=='READ_OK' for r in audit),
            'n_genes_scperturb_covered':sum(r['FOUND_IN_LOCAL_SC_PERTURB']=='YES' for r in coverage),
            'n_dataset_gene_matches':len(hits)}


def single_analysis(root,out,plan,gmt,gmeta,step14,args):
    p=Path(plan['DATASET_PATH'])
    gid=plan['GENE_ID']; sym=plan['GENE_SYMBOL']
    record={**{k:plan.get(k,'') for k in OUTPUT_COLS},'STATUS':'FAILED','ERROR':'',
            'N_CELLS_ANALYSED':0,'N_PATHWAYS_TESTED':0,'N_STEP14_MATCHES':0}
    if plan['AUTO_STATUS']!='READY_FOR_DESCRIPTIVE':
        record.update(STATUS='SKIPPED',ERROR=plan['REASON']);return record,None
    try:
        fingerprint=hashlib.sha256(json.dumps({'file':str(p),'mtime':p.stat().st_mtime_ns,'size':p.stat().st_size,
                                                'gid':gid,'target':plan['TARGET_VALUE'],'control':plan['CONTROL_VALUE'],
                                                'version':VERSION,'perms':args.permutations,
                                                'gene_sets_version':args.library_signature,'min_size':args.gsea_min_size,
                                                'max_size':args.gsea_max_size,'seed':args.seed,
                                                'max_target_cells':args.max_target_cells,'max_control_cells':args.max_control_cells},sort_keys=True).encode()).hexdigest()[:12]
        key=slug(p.stem)[:42]+'_'+slug(gid)+'_'+fingerprint
        checkpoint=out/'results_by_contrast'/f'{key}.tsv'
        if checkpoint.is_file() and not args.force:
            df=pd.read_csv(checkpoint,sep='\t')
            record.update(STATUS='SUCCESS_CACHED',N_PATHWAYS_TESTED=len(df),
                          TARGET_EXPRESSION_DIFF=(df['TARGET_EXPRESSION_DIFF'].iloc[0] if 'TARGET_EXPRESSION_DIFF' in df and len(df) else np.nan),
                          N_STEP14_MATCHES=int(df.get('STEP14_PATHWAY_MATCH',pd.Series([])).eq('YES').sum()),
                          CHECKPOINT=str(checkpoint))
            return record,df
        with h5py.File(p,'r') as f:
            obs=f['obs']; labels=column(obs,plan['PERTURBATION_COLUMN'])
            idx_t=np.flatnonzero(labels==plan['TARGET_VALUE'])
            idx_c=np.flatnonzero(labels==plan['CONTROL_VALUE'])
            if len(idx_t)<args.min_cells or len(idx_c)<args.min_cells:
                raise ValueError(f'Insufficient target/control: {len(idx_t)}/{len(idx_c)}')
            if len(idx_t)>args.max_target_cells:
                idx_t=np.sort(np.random.default_rng(42).choice(idx_t,size=args.max_target_cells,replace=False))
                record['TARGET_SAMPLED']=True
            idx_c,note=draw_controls(idx_t,idx_c,obs,args.max_control_cells,42)
            record['CONTROL_SAMPLING']=note
            idx=np.union1d(idx_t,idx_c)
            which=np.isin(idx,idx_t).astype(int)
            X=extract_selected(f,idx)
            record['N_CELLS_ANALYSED']=len(idx)
            feature_names,col_used=read_index_or_var(f)
            if X.shape[1]!=len(feature_names): raise ValueError('Gene labels do not match matrix columns')
            scale,evidence=expression_scale(f,*matrix_meta(f)[1]);record['EXPRESSION_SCALE']=scale
            score,transform=normalize_and_rank(X,which,scale)
            record['EXPRESSION_TRANSFORM']=transform
            # This is an observed, descriptive transcript difference, not
            # independently established guide efficacy or protein activity.
            gene_labels=[gene_id(z) for z in feature_names]
            target_idx=[ii for ii,x in enumerate(gene_labels) if x in {gene_id(sym),gene_id(gid)}]
            target_expr_diff=float(np.nanmean(score[target_idx])) if target_idx else np.nan
            record['TARGET_EXPRESSION_DIFF']=target_expr_diff
            record['TARGET_EXPRESSION_MEASURED']='YES' if target_idx else 'NO'
            de=pd.DataFrame({'GENE':[gene_id(s) for s in feature_names],'GENE_EXPRESSION_DIFF_PERTURB_MINUS_CONTROL':score})
            de=de[de.GENE.ne('')].sort_values('GENE_EXPRESSION_DIFF_PERTURB_MINUS_CONTROL',ascending=False).drop_duplicates('GENE')
            de['RANK_SCORE']=de.GENE_EXPRESSION_DIFF_PERTURB_MINUS_CONTROL
            if len(de)<20: raise ValueError('Too few measured expression features')
            out_de=out/'differential_expression'/f'{key}_descriptive.tsv'
            put(de,out_de)
            record['FEATURE_ID_COLUMN']=col_used
            if not gmt:
                record.update(STATUS='EXPRESSION_ONLY_NO_GMT',ERROR='No compatible source-specific GMT library')
                return record,None
            g,n=do_gsea(de,gmt,gmeta,args)
            record['N_PATHWAYS_TESTED']=int(n)
            if g.empty:
                record.update(STATUS='NO_TESTABLE_PATHWAYS',ERROR='No compatible source-specific sets overlap measured features');return record,None
            g=add_pathway_match(g,step14)
            for k in OUTPUT_COLS: g[k]=plan.get(k,'')
            g['CONTROL_SAMPLING']=note;g['EXPRESSION_SCALE']=scale;g['EXPRESSION_TRANSFORM']=transform
            g['FEATURE_ID_COLUMN']=col_used
            g['N_TARGET_CELLS_ANALYSED']=len(idx_t);g['N_CONTROL_CELLS_ANALYSED']=len(idx_c)
            g['INFERENCE_SCOPE']='DESCRIPTIVE_ONLY_NO_VERIFIED_INDEPENDENT_BIOLOGICAL_REPLICATES'
            g['PERTURBATION_ACTIVITY_DIRECTION']='NOT_VERIFIED_FROM_TARGET_TRANSCRIPT_DIFFERENCE'
            g['TARGET_EXPRESSION_DIFF']=target_expr_diff
            g['TARGET_EXPRESSION_MEASURED']='YES' if target_idx else 'NO'
            g['INDEPENDENT_BIOLOGICAL_REPLICATES_VERIFIED']='NO'
            put(g,checkpoint)
            record.update(STATUS='SUCCESS_DESCRIPTIVE',N_PATHWAYS_TESTED=len(g),
                          N_STEP14_MATCHES=int(g.STEP14_PATHWAY_MATCH.eq('YES').sum()),CHECKPOINT=str(checkpoint))
            return record,g
    except Exception as e:
        record['ERROR']=f'{type(e).__name__}: {str(e)[:400]}'
        return record,None


def output_results(out,genes,step14,plans,report,libraries,args):
    gmt,gmeta=load_libraries(libraries,out)
    mapping=pathway_mapping_audit(gmeta,step14,out)
    if not gmt:
        raise RuntimeError('No full source-specific pathway gene sets loaded. Run --prepare-libraries first.')
    audit=[]; data=[]
    for i,plan in enumerate(plans,1):
        record,df=single_analysis(None,out,plan,gmt,gmeta,step14,args)
        audit.append(record)
        if df is not None:data.append(df)
        put(audit,out/'Step16_Run_Audit.tsv')
        print(f"[CONTRAST {i}/{len(plans)}] {plan['GENE_SYMBOL']} {plan['DATASET']}: "
              f"{record['STATUS']} sets={record.get('N_PATHWAYS_TESTED',0)} {record['ERROR']}",flush=True)
    all_gsea=pd.concat(data,ignore_index=True) if data else pd.DataFrame(columns=[
      'GENE_ID','GENE_SYMBOL','DATASET','PATHWAY_KEY','PATHWAY_SOURCE','PATHWAY_ID','PATHWAY_NAME',
      'NES','GENE_SET_PERMUTATION_FDR','STEP14_PATHWAY_MATCH','PATHWAY_MATCH_METHOD'])
    put(all_gsea,out/'Step16_Directional_Pathway_Responses_All.tsv')
    matched=all_gsea[all_gsea.STEP14_PATHWAY_MATCH.eq('YES')].copy()
    put(matched,out/'Step16_Step14_Pathway_Responses.tsv')
    put(matched,out/'Step16_Gene_Pathway_Direction_Handoff.tsv')
    source_counts=Counter(matched.PATHWAY_SOURCE) if not matched.empty else {}
    summary={'version':VERSION,'n_primary_genes':len(genes),
      'n_step14_pathways':len(step14),'n_step14_pathways_exactly_mappable':sum(
          row['MAPPING_STATUS']=='EXACT_SOURCE_NATIVE_ID' for row in mapping),
      'n_auto_plan_rows':len(plans),
      'n_successful_descriptive_contrasts':sum(x['STATUS'].startswith('SUCCESS') for x in audit),
      'n_failed_or_skipped_contrasts':sum(not x['STATUS'].startswith('SUCCESS') for x in audit),
      'n_exploratory_pathway_rows':len(all_gsea),
      'n_step14_pathway_overlap_rows':len(matched),
      'n_step14_unique_pathway_ids_experimentally_tested':int(matched.STEP14_TERM_ID.nunique()) if len(matched) else 0,
      'n_step14_pathway_overlap_by_source':dict(source_counts),
      'experiment_level_FDR_available':False,
      'pathway_libraries':[{k:v for k,v in row.items() if k in ('SOURCE','LIBRARY_NAME','SHA256','GENESET_ROWS')} for row in libraries],
      'gene_set_provider':'Enrichr source-specific full libraries; may differ from Step14 g:Profiler source releases',
      'gsea_fdr_scope':'Full available GO:BP/WikiPathways/KEGG libraries per contrast; unreplicated gene-set permutations',
      'causal_validation':'NOT_ESTABLISHED',**report}
    tmp=out/'Step16_summary.json.tmp';tmp.write_text(json.dumps(summary,indent=2),encoding='utf8');tmp.replace(out/'Step16_summary.json')
    print('\n[SUMMARY]',json.dumps(summary,indent=2),flush=True)
    if len(matched):
        shown=[x for x in ('GENE_SYMBOL','DATASET','PATHWAY_SOURCE','STEP14_TERM_ID','STEP14_TERM_NAME','NES',
                          'GENE_SET_PERMUTATION_FDR','TARGET_EXPRESSION_DIFF') if x in matched]
        print('\n[EXACT SOURCE-ID MATCHED RESPONSES]\n'+matched[shown].head(30).to_string(index=False),flush=True)
    print('[INTERPRETATION] GSEA NES is a transcriptional enrichment of a descriptive cell contrast. '+
          'A gene-set-permutation FDR is NOT biological replicate-level evidence; risk-allele direction is unresolved.',flush=True)


def stats(out):
    summary=out/'Step16_summary.json'
    if not summary.exists():
        print(f'No completed analysis summary yet: {summary}')
        return
    print(summary.read_text())
    for filename,columns in [
        ('Step16_Gene_Coverage_Summary.tsv',['GENE_SYMBOL','N_DATASETS','N_TARGET_CELLS']),
        ('Step16_Auto_Analysis_Plan.tsv',['GENE_SYMBOL','DATASET','TARGET_CELLS','CONTROL_CELLS','AUTO_STATUS']),
        ('Step16_Run_Audit.tsv',['GENE_SYMBOL','DATASET','STATUS','ERROR','N_PATHWAYS_TESTED','N_STEP14_MATCHES']),
        ('Step16_Step14_Pathway_Responses.tsv',['GENE_SYMBOL','STEP14_TERM_NAME','DATASET','NES','GENE_SET_PERMUTATION_FDR'])]:
        p=out/filename
        if p.is_file():
            try:
                t=pd.read_csv(p,sep='\t')
                print(f'\n{filename} [{len(t)} rows]')
                show=[c for c in columns if c in t.columns]
                if show:print(t[show].head(50).to_string(index=False))
            except Exception as e:print(f'[STATS WARNING] {filename}: {e}')


def status(root,step15,out,data_dir,args=None):
    print(f'GWAS2Mech Step16 {VERSION}\nSCRIPT: {Path(__file__).resolve()}')
    for name,p in [('Step15 handoff',step15/'Step15_Gene_Pathway_Handoff.tsv'),
                   ('Step15 primary audit',step15/'Step15_Gene_Mapping_Audit.tsv'),
                   ('Step16 results',out/'Step16_summary.json')]:
        print(f'{name}: {"OK" if p.is_file() else "MISSING"} {p}')
    if args is not None:
        try:
            _,pathways=inputs(step15)
            for r in acquire_libraries(root,args,pathways,False):
                print(f"PATHWAY_LIBRARY {r['SOURCE']} {r['STATUS']}: {r['PATH']} ({r['GENESET_ROWS']} sets)")
                if r['STATUS']!='READY':print(f"  PREPARE: --prepare-libraries | SOURCE: {r['DOWNLOAD_URL']}")
        except Exception as e:print(f'[PATHWAY RESOURCES] {type(e).__name__}: {e}')
    files=list(data_dir.rglob('*.h5ad')) if data_dir.exists() else []
    print(f'SCPERTURB: {len(files)} visible files in {data_dir}')
    print(f'OUTPUT DIRECTORY: {out}')


def create_slurm(root,args):
    template=root/'HPC_CPU.sh'
    config={'partition':'general','qos':'normal','nodes':'1','ntasks':'1','cpus-per-task':str(args.threads),
            'mem':'64G','time':'24:00:00'}
    if template.exists():
        for line in template.read_text(errors='replace').splitlines():
            m=re.match(r'^\s*#SBATCH\s+--([A-Za-z-]+)(?:[=\s]+(.*))?$',line)
            if m and (m.group(2) or '').strip() and m.group(1) in set(config)|{'account'}:
                config[m.group(1)]=m.group(2).strip()
    for k in ('account','partition','qos','mem','time'):
        if getattr(args,k,None): config[k]=str(getattr(args,k))
    config['cpus-per-task']=str(args.cpus or args.threads)
    name=f'Step16_Perturb_{slug(args.phenotype)}_{args.ancestry.lower()}_{args.reference_panel}'
    config.update({'job-name':name,'output':name+'_%j.out','error':name+'_%j.err'})
    head=['job-name','account','partition','qos','nodes','ntasks','cpus-per-task','mem','time','output','error']
    script=root/(name+'.sh')
    argv=[str(root/'envs/perturb/bin/python'),'-u',str(root/'Step16_Perturbation_Directional_Pathways.py'),
          '--phenotype',args.phenotype,'--ancestry',args.ancestry,'--reference-panel',args.reference_panel,
          '--threads',str(args.threads),'--min-cells',str(args.min_cells),
          '--max-control-cells',str(args.max_control_cells),'--max-target-cells',str(args.max_target_cells),
          '--permutations',str(args.permutations),'--seed',str(args.seed),
          '--gsea-min-size',str(args.gsea_min_size),'--gsea-max-size',str(args.gsea_max_size)]
    if args.library_dir:argv+=['--library-dir',args.library_dir]
    if args.index is not None:argv+=['--index',str(args.index)]
    if args.force:argv.append('--force')
    script.write_text('#!/bin/bash --login\n'+'\n'.join('#SBATCH --'+k+'='+config[k] for k in head if k in config)+
       '\nset -euo pipefail\ncd '+shlex.quote(str(root))+'\n'+
       'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\n'+
       ' '.join(shlex.quote(x) for x in argv)+'\n',encoding='utf8')
    script.chmod(0o755)
    print(f'JOB SCRIPT: {script}\nSUBMIT: sbatch {script.name}')


def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--version',action='version',version='GWAS2Mech Step16 '+VERSION)
    ap.add_argument('--root',default='.')
    ap.add_argument('--phenotype',required=True)
    ap.add_argument('--ancestry',required=True)
    ap.add_argument('--reference-panel',default='1000G')
    ap.add_argument('--index',type=int)
    actions=ap.add_mutually_exclusive_group()
    for a in ('status','stats','check-resources','slurm','inventory-only','prepare-libraries'):
        actions.add_argument('--'+a,action='store_true')
    ap.add_argument('--data-root',default='',help='Optional override; normally auto-discovered')
    ap.add_argument('--library-dir',default='',help='Path to source-specific GMT resources (optional)')
    ap.add_argument('--seed',type=int,default=42)
    ap.add_argument('--limit-files',type=int,default=0)
    ap.add_argument('--min-cells',type=int,default=20)
    ap.add_argument('--max-target-cells',type=int,default=2000)
    ap.add_argument('--max-control-cells',type=int,default=2500)
    ap.add_argument('--gsea-min-size',type=int,default=10)
    ap.add_argument('--gsea-max-size',type=int,default=500)
    ap.add_argument('--permutations',type=int,default=500)
    ap.add_argument('--threads',type=int,default=4)
    ap.add_argument('--force',action='store_true')
    for opt in ('account','partition','qos','mem','time'):ap.add_argument('--'+opt)
    ap.add_argument('--cpus',type=int)
    args=ap.parse_args()
    root,step15,out,data_dir=paths(args)
    print('[Step16] VERSION='+VERSION+' SCRIPT='+str(Path(__file__).resolve()),flush=True)
    if args.status:
        status(root,step15,out,data_dir,args);return
    if args.stats:
        stats(out);return
    try:genes,step14=inputs(step15)
    except Exception as e:raise SystemExit(f'Missing/invalid Step15: {e}')
    if args.prepare_libraries:
        rows=acquire_libraries(root,args,step14,True)
        put(rows,out/'Step16_Pathway_Library_Inventory.tsv')
        for x in rows:print(f"[LIBRARY] {x['SOURCE']} {x['STATUS']} {x['LIBRARY_NAME']} "
                            f"sets={x['GENESET_ROWS']} path={x['PATH']} {x['ERROR']}",flush=True)
        raise SystemExit(0 if all(r['STATUS']=='READY' for r in rows) else 1)
    libraries=acquire_libraries(root,args,step14,False)
    put(libraries,out/'Step16_Pathway_Library_Inventory.tsv')
    missing=[r for r in libraries if r['STATUS']!='READY']
    if args.check_resources:
        status(root,step15,out,data_dir,args)
        checks={'Step15 handoff':(step15/'Step15_Gene_Pathway_Handoff.tsv').is_file(),
                'scPerturb files':data_dir.exists() and bool(list(data_dir.rglob('*.h5ad'))),
                'ALL REQUIRED SOURCE-SPECIFIC GMT':not missing}
        for m in ('gseapy','h5py','scipy','numpy','pandas'):
            try:__import__(m);checks['module '+m]=True
            except ImportError:checks['module '+m]=False
        for k,v in checks.items():print(k,'OK' if v else 'MISSING')
        raise SystemExit(0 if all(checks.values()) else 1)
    if missing and not args.inventory_only:
        print('[ERROR] Missing full gene-set libraries:',flush=True)
        for row in missing:print('  '+row['SOURCE']+': '+row['PATH']+' URL '+row['DOWNLOAD_URL'],flush=True)
        raise SystemExit('Run this command on a network-enabled login node FIRST: '\
           'envs/perturb/bin/python Step16_Perturbation_Directional_Pathways.py '\
           '--phenotype "Parkinson\'s disease" --ancestry EUR --reference-panel 1000G30x --prepare-libraries')
    args.library_signature=hashlib.sha256(json.dumps([{'source':x['SOURCE'],'sha':x['SHA256']} for x in libraries],sort_keys=True).encode()).hexdigest()[:20]
    if args.slurm:
        create_slurm(root,args);return
    genes,step14,plans,report=scan(root,step15,out,data_dir,args)
    report.update(lincs_metadata(root,genes,out))
    if args.inventory_only:
        (out/'Step16_summary.json').write_text(json.dumps({**report,'version':VERSION,'mode':'inventory_only'},indent=2));return
    output_results(out,genes,step14,plans,report,libraries,args)

if __name__=='__main__':main()
