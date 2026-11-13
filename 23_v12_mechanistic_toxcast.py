#!/usr/bin/env python3
"""V12: EPA invitrodb v4.3 assay annotations + optional full mc5/6 results.

Without the optional mc5/6 CSV: annotation inventory ONLY; absolutely no hit claims.
Exact CAS only for NHANES biomarkers; no parent-metabolite proxy.
Requires pandas, numpy, openpyxl, requests, matplotlib (existing v11 environment).
"""
from __future__ import annotations
import argparse, hashlib, io, json, re, sys, zipfile
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import requests

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'/'toxicology'
OUT=ROOT/'results'/'v12_mechanistic_exposomics'
CROSS=ROOT/'results'/'v11_toxicology_bridge'/'v11_verified_chemical_crosswalk.csv'
# Official pages and SHA-256 values from the EPA Clowder file metadata.
SOURCES={
 'annotations':{'id':'68af6bd3e4b02565fc7c3aa8','sha':'eecbdac47fe4b4c3db765d412d15526331c269d36c7d79dde938f59ea60461f9','file':'assay_annotations_invitrodb_v4_3_AUG2025.xlsx'},
 'targets':{'id':'68af6bd3e4b02565fc7c3aa0','sha':'bf45c47dc383c16aea9c4cf92993f7034a3141a252aa02d9f7aa7216ac84c394','file':'assay_target_mappings_invitrodb_v4_3_AUG2025.xlsx'},
}
# Exploratory literal gene/endpoint text matching: pathways NOT automatically confirmed.
TERMS={
 'Inflammatory signaling':r'\b(?:nf[ -]?kappa[ -]?b|nf[ -]?kb|nfkb|rela|relb|tnf|tnfa|il[ -]?1b|il1b|il[ -]?6|il6|tlr[234579]|stat3|jak[123]|cox[ -]?2|ptgs2|n[fk]b)\b',
 'Oxidative stress':r'\b(?:nrf2|nfe2l2|keap1|nqo1|h[m]?ox1|h[o]?1|sod[123]|gpx[1-9]|ros|oxidative|antioxidant|glutathione)\b',
 'Immune response':r'\b(?:immune|immuno|interferon|ifn[ag]|irf[1-9]|cytokine|chemokine|cxcl\d+|ccl\d+|t[ -]?cell|b[ -]?cell|il[ -]?(?:2|4|10|12|17))\b',
}
FALLBACK=[('PFHxS','355-46-4'),('PFNA','375-95-1'),('PFOA','335-67-1'),('PFOS','1763-23-1'),('MBP','131-70-4'),('MEP','2306-33-4'),('MiBP','30833-53-5')]

def norm(s):return re.sub(r'[^a-z0-9]+','_',str(s).strip().lower()).strip('_')
def cas(s):
    if pd.isna(s):return ''
    m=re.fullmatch(r'(\d{2,7})-(\d{2})-(\d)',str(s).strip())
    return '-'.join(m.groups()) if m else ''
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def fetch(name):
    spec=SOURCES[name];DATA.mkdir(parents=True,exist_ok=True);p=DATA/spec['file']
    if p.exists():
        if digest(p)!=spec['sha']:raise ValueError(f'{p} checksum differs from published EPA SHA256; do not use it.')
        print('[cached/verified]',p.name);return p
    urls=[f"https://clowder.edap-cluster.com/api/files/{spec['id']}/blob",f"https://clowder.edap-cluster.com/files/{spec['id']}/blob"]
    last=None
    for url in urls:
        temp=p.with_suffix('.partial')
        try:
            with requests.get(url,timeout=(15,90),stream=True,headers={'User-Agent':'Academic-Exposome-EPA-v12/1.0'}) as r:
                r.raise_for_status()
                with temp.open('wb') as f:
                    for b in r.iter_content(131072):
                        if b:f.write(b)
                        if f.tell()>8_000_000:raise ValueError('Unexpected >8MB annotation file')
            if digest(temp)!=spec['sha']:raise ValueError('SHA mismatch (possibly a website login page)')
            temp.replace(p);print('[download/verified]',p.name);return p
        except Exception as e:last=e;temp.unlink(missing_ok=True)
    raise RuntimeError(f'Unable to download official EPA {name}: {last}\nManual file page: https://clowder.edap-cluster.com/files/{spec["id"]}')

def find_sheet(path):
    """Find a real annotation table; don't silently guess missing keys."""
    xf=pd.ExcelFile(path);best=None
    for sh in xf.sheet_names:
        for header in range(0,7):
            try:sample=pd.read_excel(xf,sheet_name=sh,header=header,nrows=12)
            except Exception:continue
            cc=[norm(x) for x in sample.columns]
            score=4*int('aeid' in cc)+3*int('aenm' in cc)+sum(any(k in c for k in ['gene','target','assay','endpoint','description']) for c in cc)
            if not any(c in cc for c in ('aeid','aenm','assay_component_endpoint_name')):continue
            if best is None or score>best[0]:best=(score,sh,header)
    if not best:raise ValueError(f'Could not find aeid/aenm endpoint key in {path.name}; inspect workbook')
    df=pd.read_excel(xf,sheet_name=best[1],header=best[2]);df.columns=[norm(c) for c in df.columns]
    df=df.loc[:,~df.columns.duplicated()]
    if df.empty:raise ValueError(f'Empty annotation table {path}')
    print('[metadata]',path.name,'sheet',best[1],'rows',len(df),'columns',list(df.columns)[:15]);return df

def endpoint_keys(df):
    e=next((x for x in ('aeid','assay_component_endpoint_id') if x in df),None)
    n=next((x for x in ('aenm','assay_component_endpoint_name','assay_endpoint_name') if x in df),None)
    return e,n

def curated_inventory(anno,targets):
    e,n=endpoint_keys(anno)
    if e is None and n is None:raise ValueError('EPA annotations have neither aeid nor aenm')
    a=anno.copy();a['_key']=a[e].map(str).str.strip() if e else a[n].map(str).str.strip()
    a=a[a['_key'].ne('')&a['_key'].ne('nan')].drop_duplicates('_key')
    te,tn=endpoint_keys(targets)
    if te and e:tk=te
    elif tn and n:tk=tn
    else:tk=None
    if tk:
        t=targets.copy();t['_key']=t[tk].astype(str).str.strip()
        # Retain multiple target annotations instead of one arbitrary record.
        tcols=[c for c in t if c!='_key' and any(z in c for z in ('gene','target','symbol','uniprot','entrez'))]
        if tcols:
            grp=t.groupby('_key',dropna=False)[tcols].agg(lambda z:'; '.join(sorted(set(str(i).strip() for i in z.dropna() if str(i).strip())))[:1000]).reset_index()
            a=a.merge(grp,on='_key',how='left',suffixes=('','_targets'),validate='one_to_one')
    texts=[]
    for _,r in a.iterrows():
        parts=[]
        for c,v in r.items():
            if c=='_key':continue
            if any(x in c for x in ('name','description','gene','symbol','target','biological','pathway','endpoint','function','platform')) and pd.notna(v):parts.append(str(v))
        texts.append(' '.join(parts))
    a['annotation_text']='' # keep outputs compact; don't republish all free-text descriptions
    for label,pat in TERMS.items():a[norm(label)+'_candidate']=[bool(re.search(pat,t,re.I)) for t in texts]
    # Candidate is just a lexical filter and may include false matches, not validation.
    a['candidate_count']=a[[norm(x)+'_candidate' for x in TERMS]].sum(axis=1)
    return a,e,n

def make_cross():
    if CROSS.exists():
        d=pd.read_csv(CROSS,dtype={'casrn':'string'})
        if {'short','casrn'}.issubset(d.columns):
            d=d[['short','casrn']].copy();d['casrn']=d.casrn.map(cas)
            if len(d)==7 and d.casrn.ne('').all() and d.casrn.nunique()==7:return d
            raise ValueError('Existing v11 crosswalk does not contain exactly 7 unique valid CASRNs')
    return pd.DataFrame(FALLBACK,columns=['short','casrn'])

def process_full(csv_path,inv,e,n,chem):
    """Stream full EPA multi-concentration response summary without loading millions of rows."""
    matched=[];total=0;rows_kept=0;cas_set=set(chem.casrn)
    for chunk in pd.read_csv(csv_path,chunksize=65000,low_memory=False):
        chunk.columns=[norm(c) for c in chunk.columns]
        ccol=next((x for x in ('casn','casrn') if x in chunk),None)
        hcol=next((x for x in ('hitcall','hitc') if x in chunk),None)
        acol=next((x for x in ('aeid','aenm') if x in chunk),None)
        if not ccol or not hcol or not acol:raise ValueError(f'Full mc5/6 CSV must contain CAS, hitcall and aeid/aenm. Found: {list(chunk.columns)}')
        total+=len(chunk);mask=chunk[ccol].map(cas).isin(cas_set)
        s=chunk.loc[mask,[ccol,hcol,acol]].copy()
        if s.empty:continue
        s.columns=['casrn','hitcall','endpoint'];s.casrn=s.casrn.map(cas)
        s['hitcall']=pd.to_numeric(s.hitcall,errors='coerce')
        valid=s.hitcall.dropna()
        if ((valid<0)|(valid>1)).any():raise ValueError('Hitcall must be in [0,1]; refusing invalid release/schema')
        s['endpoint']=s.endpoint.astype(str).str.strip();matched.append(s);rows_kept+=len(s)
    if not matched:return pd.DataFrame(),pd.DataFrame(),{'full_rows_scanned':total,'matched_rows':0}
    d=pd.concat(matched,ignore_index=True)
    keycol=e if e and 'aeid' in chunk.columns else n if n and 'aenm' in chunk.columns else None
    # Prefer aeid if present; otherwise endpoint names. Infer from candidate keys (numeric vs text).
    if e is not None and d.endpoint.str.fullmatch(r'\d+(?:\.0)?').mean()>.8:kt=e
    elif n is not None:kt=n
    elif e is not None:kt=e
    else:raise ValueError('No compatible EPA annotation key')
    if kt=='aeid':d['endpoint']=d.endpoint.str.replace(r'\.0$','',regex=True)
    anno=inv.copy();anno['_key']=anno[kt].astype(str).str.strip().str.replace(r'\.0$','',regex=True) if kt in anno else anno['_key']
    keep=['_key']+[norm(t)+'_candidate' for t in TERMS]
    merged=d.merge(anno[keep].drop_duplicates('_key'),left_on='endpoint',right_on='_key',how='left',validate='many_to_one')
    merged['hit_ge_0_9']=np.where(merged.hitcall.notna(),merged.hitcall.ge(.9).astype(float),np.nan)
    # One row per chemical/assay; handle duplicates: duplicates are not new independent endpoints.
    key=['casrn','endpoint'];dup=int(merged.duplicated(key).sum())
    endpoint=merged.groupby(key,as_index=False,dropna=False).agg(hitcall=('hitcall','max'),
        hit_ge_0_9=('hit_ge_0_9','max'),**{norm(t)+'_candidate':(norm(t)+'_candidate','max') for t in TERMS})
    endpoint['annotation_found']=endpoint[[norm(t)+'_candidate' for t in TERMS]].notna().any(axis=1)
    coverage=[]
    for c in chem.itertuples(index=False):
        z=endpoint[endpoint.casrn==c.casrn];tot=len(z);valid=z.hitcall.notna()
        r={'short':c.short,'casrn':c.casrn,'n_assay_endpoints':tot,'n_assays_with_hitcall':int(valid.sum()),'n_active_hitcall_ge_0_9':int((z.hit_ge_0_9==1).sum()),'n_assay_annotation_mapped':int(z.annotation_found.sum())}
        for t in TERMS:
            col=norm(t)+'_candidate';b=z[col].fillna(False).astype(bool)
            r[norm(t)+'_candidate_tested']=int((b&valid).sum());r[norm(t)+'_candidate_active']=int((b&(z.hit_ge_0_9==1)).sum())
        coverage.append(r)
    return endpoint,pd.DataFrame(coverage),{'full_rows_scanned':total,'matched_rows':rows_kept,'duplicate_endpoint_rows':dup,'key_used':kt,'hit_threshold':0.9}

def figures(inv,coverage):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    terms=[norm(t) for t in TERMS]
    if coverage is None or coverage.empty:
        fig,ax=plt.subplots(figsize=(8,4.4));names=list(TERMS)
        vals=[int(inv[t+'_candidate'].sum()) for t in terms]
        ax.barh(names,vals);ax.set_xlabel('Number of assay endpoints tagged by annotation keywords')
        ax.set_title('EPA assay annotation candidates (NOT measured chemical activity)')
        fig.tight_layout();fig.savefig(OUT/'v12_annotation_inventory.png',dpi=170);plt.close(fig)
    else:
        cols=[t+'_candidate_active' for t in terms]
        m=coverage.set_index('short')[cols].astype(float)
        fig,ax=plt.subplots(figsize=(9,4.5));im=ax.imshow(m.to_numpy(),aspect='auto',cmap='Blues')
        ax.set_xticks(range(len(terms)),list(TERMS),rotation=20,ha='right');ax.set_yticks(range(len(m)),m.index)
        for (i,j),v in np.ndenumerate(m.to_numpy()):ax.text(j,i,str(int(v)),ha='center',va='center',color='white' if v>max(1,m.to_numpy().max()/2) else 'black')
        ax.set_title('Assay hit counts (hitcall ≥0.9), candidate pathway annotations')
        fig.colorbar(im,ax=ax,label='Active tested endpoints');fig.tight_layout();fig.savefig(OUT/'v12_candidate_pathway_hit_heatmap.png',dpi=170);plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mc56',type=Path,help='Optional downloaded official mc5-6 CSV; NEVER fetched automatically (EPA ZIP is ~7.5 GB)')
    p.add_argument('--annotations',type=Path,help='Local annotations xlsx (for offline tests or manual download)')
    p.add_argument('--targets',type=Path,help='Local target mappings xlsx')
    a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    paths={k:(getattr(a,k) if getattr(a,k) else fetch(k)) for k in SOURCES}
    audit={'timestamp_utc':datetime.now(timezone.utc).isoformat(),'release':'EPA invitrodb v4.3 (August 2025)',
           'sources':{k:{'path':str(v),'sha256':digest(v),'official_sha256_valid':digest(v)==SOURCES[k]['sha']} for k,v in paths.items()},
           'classification':'Exploratory keyword-based endpoint candidate tags, NOT confirmed pathway effects',
           'parents_not_substituted_for_metabolites':True,'missing_as_inactive':False,'human_risk_score':False}
    if not all(s['official_sha256_valid'] for s in audit['sources'].values()):
        print('[WARNING] Some local annotation files have unverified checksums; provenance recorded')
    annotations=find_sheet(paths['annotations']);targets=find_sheet(paths['targets'])
    inv,e,n=curated_inventory(annotations,targets)
    cols=[c for c in inv if c!='annotation_text']
    inv[cols].to_csv(OUT/'v12_annotated_assay_inventory.csv',index=False)
    summary=pd.DataFrame([{'domain':domain,'n_candidate_endpoints':int(inv[norm(domain)+'_candidate'].sum())} for domain in TERMS])
    summary.to_csv(OUT/'v12_candidate_domain_inventory.csv',index=False)
    cross=make_cross();cross.to_csv(OUT/'v12_exact_cas_crosswalk.csv',index=False)
    coverage=None
    if a.mc56:
        if not a.mc56.exists():raise FileNotFoundError(a.mc56)
        ass,coverage,details=process_full(a.mc56,inv,e,n,cross)
        if not ass.empty:ass.to_csv(OUT/'v12_chemical_assay_hits.csv',index=False)
        coverage.to_csv(OUT/'v12_chemical_domain_coverage_and_hits.csv',index=False)
        audit['assay_status']='FULL_MC56_ANALYZED';audit['full_file_sha256']=digest(a.mc56);audit.update(details)
        print('[V12 FULL]',details,'\n',coverage.to_string(index=False))
    else:
        audit['assay_status']='ANNOTATION_ONLY_NO_CHEMICAL_ACTIVITY_CALCULATED'
        (OUT/'V12_FULL_DATA_NEEDED.txt').write_text(
            'Assay metadata successfully audited. No per-chemical activity is inferred without EPA full mc5-6 results.\n'
            'Download official INVITRODB_SUMMARY.zip (~7.5 GB) from:\n'
            'https://clowder.edap-cluster.com/files/68af6b70e4b02565fc7c3a98\n'
            'Extract mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv only, then run:\n'
            '.venv/bin/python 23_v12_mechanistic_toxcast.py --mc56 /path/to/the/file.csv\n',encoding='utf-8')
        print('[V12] Annotation inventory completed; NO per-chemical mechanistic hit claims without full EPA mc5-6 CSV.')
    figures(inv,coverage)
    (OUT/'v12_provenance_audit.json').write_text(json.dumps(audit,indent=2,default=str),encoding='utf-8')
    print('[V12 saved]',OUT)
if __name__=='__main__':main()
