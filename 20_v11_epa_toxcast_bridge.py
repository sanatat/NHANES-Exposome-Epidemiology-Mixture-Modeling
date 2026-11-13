#!/usr/bin/env python3
"""V11: cross-source NHANES biomarker -> EPA invitrodb v4.3 cytotoxicity bridge.

Only *exact* CAS registry number matches are accepted. Parent plasticizers are
reported in a separate context table and never substituted for urine metabolites.
This is a mechanistic data-availability audit, not a causal risk score.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data' / 'toxicology'
OUT = ROOT / 'results' / 'v11_toxicology_bridge'
POOLED = ROOT / 'data' / 'processed' / 'analysis_joint_pooled_2015_2018.csv'
EPA_FILE = DATA / 'cytotox_invitrodb_v4_3_AUG2025.xlsx'
EPA_URLS = [
    'https://clowder.edap-cluster.com/api/files/68af6bd3e4b02565fc7c3aa4/blob',
    'https://clowder.edap-cluster.com/files/68af6bd3e4b02565fc7c3aa4/blob',
]
EPA_PAGE = 'https://clowder.edap-cluster.com/files/68af6bd3e4b02565fc7c3aa4'
EPA_SHA256 = '3bc444f658c80e4ac9a1189f083c0120365b504e10d0aa4e0f07feda294c08dd'
EPA_RELEASE = 'US EPA ToxCast invitrodb 4.3 (August 2025)'

# Chemical identifiers are CASRNs of *measured* NHANES chemical biomarkers,
# not of parent plasticizers. No unverified DTXSID assignment is fabricated.
CHEMICALS = [
 dict(nhanes_column='LBXPFHS', name='Perfluorohexanesulfonic acid (PFHxS)', short='PFHxS', family='PFAS', specimen='serum', casrn='355-46-4', type='measured substance', source='https://www.ncbi.nlm.nih.gov/books/NBK592148/'),
 dict(nhanes_column='LBXPFNA', name='Perfluorononanoic acid (PFNA)', short='PFNA', family='PFAS', specimen='serum', casrn='375-95-1', type='measured substance', source='https://www.ncbi.nlm.nih.gov/books/NBK592148/'),
 dict(nhanes_column='LBXNFOA', name='n-Perfluorooctanoic acid (linear PFOA)', short='PFOA', family='PFAS', specimen='serum', casrn='335-67-1', type='measured substance', source='https://www.ncbi.nlm.nih.gov/books/NBK592148/'),
 dict(nhanes_column='LBXNFOS', name='n-Perfluorooctanesulfonic acid (linear PFOS)', short='PFOS', family='PFAS', specimen='serum', casrn='1763-23-1', type='measured substance', source='https://www.ncbi.nlm.nih.gov/books/NBK592148/'),
 dict(nhanes_column='URXMBP', name='Mono-n-butyl phthalate (MBP)', short='MBP', family='phthalate metabolite', specimen='urine', casrn='131-70-4', type='measured urinary metabolite', source='https://webbook.nist.gov/cgi/cbook.cgi?ID=131-70-4'),
 dict(nhanes_column='URXMEP', name='Monoethyl phthalate (MEP)', short='MEP', family='phthalate metabolite', specimen='urine', casrn='2306-33-4', type='measured urinary metabolite', source='https://webbook.nist.gov/cgi/cbook.cgi?ID=2306-33-4'),
 dict(nhanes_column='URXMIB', name='Monoisobutyl phthalate (MiBP)', short='MiBP', family='phthalate metabolite', specimen='urine', casrn='30833-53-5', type='measured urinary metabolite', source='https://pubchem.ncbi.nlm.nih.gov/compound/92272'),
]
PARENTS = [
 dict(metabolite='MBP', name='Dibutyl phthalate (DBP)', casrn='84-74-2'),
 dict(metabolite='MEP', name='Diethyl phthalate (DEP)', casrn='84-66-2'),
 dict(metabolite='MiBP', name='Diisobutyl phthalate (DiBP)', casrn='84-69-5'),
]
FLAG_MAP = dict(zip([c['nhanes_column'] for c in CHEMICALS], [
    'LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL','URDMBPLC','URDMEPLC','URDMIBLC']))


def cas_clean(s: object) -> str:
    if pd.isna(s): return ''
    s = str(s).strip().replace('−','-').replace('–','-')
    m = re.fullmatch(r'(\d{2,7})-(\d\d)-(\d)',s)
    return '-'.join(m.groups()) if m else ''


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def verified_xlsx(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 30000 and sha256(path)==EPA_SHA256


def fetch_epa() -> Path:
    DATA.mkdir(parents=True,exist_ok=True)
    if verified_xlsx(EPA_FILE):
        print(f'[cache/sha256 OK] {EPA_FILE} ({EPA_FILE.stat().st_size:,} bytes)')
        return EPA_FILE
    if EPA_FILE.exists():
        raise RuntimeError(f'Existing EPA file checksum mismatch: {EPA_FILE}. Do not overwrite silently.')
    tmp = EPA_FILE.with_suffix('.download')
    error=None
    for attempt,url in enumerate(EPA_URLS, start=1):
        try:
            with requests.get(url,stream=True,timeout=(20,75),headers={'User-Agent':'Exposome-reproducible-academic-research/1.0'}) as r:
                r.raise_for_status()
                with tmp.open('wb') as f:
                    for block in r.iter_content(131072):
                        if block:
                            f.write(block)
                            if f.tell()>20_000_000: raise ValueError('Download larger than expected small XLSX')
            if sha256(tmp)!=EPA_SHA256:
                raise ValueError('EPA checksum mismatch, perhaps HTML/access page instead of the XLSX')
            tmp.replace(EPA_FILE)
            print('[download/sha256 OK] Official EPA invitrodb v4.3 cytotox table')
            return EPA_FILE
        except Exception as exc:
            error=exc
            tmp.unlink(missing_ok=True)
            if attempt<len(EPA_URLS):time.sleep(attempt)
    raise RuntimeError(f'Automatic official file download failed: {error}\n'
                       f'Open {EPA_PAGE}, choose Download, save as {EPA_FILE}.\n'
                       f'Then run the script again. The exact official SHA256 is {EPA_SHA256}')


def norm(s:str)->str:
    return re.sub(r'[^a-z0-9]+','_',str(s).lower().strip()).strip('_')


def load_cytotox(path:Path)->tuple[pd.DataFrame,dict]:
    # Format can vary across public releases; inspect headers across sheets/rows.
    best=None; score=-1; descriptor={}
    xl=pd.ExcelFile(path)
    for sheet in xl.sheet_names:
        for hdr in (0,1,2,3,4):
            try: head=pd.read_excel(xl,sheet_name=sheet,header=hdr,nrows=6)
            except Exception:continue
            cols=[norm(x) for x in head.columns]
            sc=sum(any(x in col for col in cols) for x in ('casn','casrn','dsstox_substance_id','cytotox','ntested','nhit'))
            if sc>score:
                best=(sheet,hdr);score=sc;descriptor={'sheet':sheet,'header_row_zero_indexed':hdr,'preview_columns':cols}
    if best is None or score<2:
        raise ValueError(f'Unexpected EPA Excel schema, best header score={score}; check file manually')
    d=pd.read_excel(xl,sheet_name=best[0],header=best[1])
    d.columns=[norm(x) for x in d.columns]
    d=d.loc[:,~d.columns.duplicated()].copy()
    key=next((c for c in ('casn','casrn','cas','cas_number','cas_registry_number') if c in d),None)
    if key is None:
        raise ValueError(f'EPA file lacks exact CAS column, columns={list(d.columns)}; refuse name-only matching')
    d['exact_cas']=d[key].map(cas_clean)
    d=d[d.exact_cas!=''].copy()
    descriptor.update({'cas_column':key,'rows_exact_cas':len(d),'all_columns':list(d.columns)})
    return d,descriptor


def unique_record(d:pd.DataFrame,cas:str)->tuple[str,pd.Series|None,int]:
    a=d[d.exact_cas==cas]
    if len(a)==0:return 'NO_EXACT_CAS_MATCH',None,0
    # Multiple rows may correspond to substances/forms/test conditions: do not silently select.
    if len(a)>1:return 'AMBIGUOUS_DUPLICATE_CAS',None,len(a)
    return 'EXACT_CAS_MATCH',a.iloc[0],1


def num(obj:object)->float|None:
    try:
        z=float(obj)
        return z if np.isfinite(z) else None
    except (ValueError,TypeError):return None


def extract_epa(rec:pd.Series|None)->dict:
    fields=['casn','chnm','dsstox_substance_id','cytotox_median_um','cytotox_lower_bound_um','ntested','nhit','cytotox_median_log','cytotox_lower_bound_log','burst_pct','cytotox_mad']
    if rec is None:return {f'epa_{x}':None for x in fields}|{'cytotox_1000_uM_possible_default_ceiling':None}
    vals={f'epa_{x}':rec.get(x,None) for x in fields}
    for k in vals:
        if k not in ('epa_casn','epa_chnm','epa_dsstox_substance_id'):
            vals[k]=num(vals[k])
        elif vals[k] is not None and pd.isna(vals[k]):vals[k]=None
    med=vals['epa_cytotox_median_um']
    vals['cytotox_1000_uM_possible_default_ceiling']=bool(med is not None and abs(med-1000)<1e-8)
    return vals


def weighted_geom(values:pd.Series,weights:pd.Series)->float|None:
    x=pd.to_numeric(values,errors='coerce').to_numpy(dtype=float)
    w=pd.to_numeric(weights,errors='coerce').to_numpy(dtype=float)
    good=np.isfinite(x)&np.isfinite(w)&(x>0)&(w>0)
    return float(np.exp(np.average(np.log(x[good]),weights=w[good]))) if good.any() else None


def nhanes_profiles(path:Path)->list[dict]:
    if not path.exists():
        print(f'[WARN] NHANES pooled file missing: {path}. EPA exact-match audit will continue without exposure statistics.')
        return [{"nhanes_column":c['nhanes_column'],'nhanes_n':None,'nhanes_weighted_gmean_ng_mL':None,'nhanes_weighted_pct_below_lod':None} for c in CHEMICALS]
    d=pd.read_csv(path,low_memory=False)
    if 'crp_harmonized_mg_L' not in d or 'survey_weight' not in d:
        raise ValueError('Pooled file required v6 columns missing')
    d=d[d.crp_harmonized_mg_L.gt(0)&d.crp_harmonized_mg_L.le(10)].copy()
    result=[]
    for c in CHEMICALS:
        col=c['nhanes_column'];fl=FLAG_MAP[col]
        if col not in d or fl not in d: raise ValueError(f'Expected {col}/{fl} missing in pooled data')
        x=pd.to_numeric(d[col],errors='coerce');w=pd.to_numeric(d['survey_weight'],errors='coerce')
        flag=pd.to_numeric(d[fl],errors='coerce')
        valid=x.notna()&x.gt(0)&w.gt(0)&w.notna()&flag.isin([0,1])
        v=d.loc[valid]
        lod=float(np.average(flag.loc[valid].to_numpy(float)*100, weights=w.loc[valid].to_numpy(float))) if valid.any() else None
        result.append(dict(nhanes_column=col,nhanes_n=int(valid.sum()),
             nhanes_weighted_gmean_ng_mL=weighted_geom(v[col],v['survey_weight']),
             nhanes_weighted_pct_below_lod=lod))
    return result


def do_cytotox(cytotox:pd.DataFrame, profiles:list[dict])->tuple[pd.DataFrame,pd.DataFrame]:
    profile={r['nhanes_column']:r for r in profiles}
    rows=[]
    for c in CHEMICALS:
        status,rec,n=unique_record(cytotox,c['casrn'])
        rows.append({**c,**profile[c['nhanes_column']], 'match_status':status,'n_epa_rows_matching_cas':n,
                     **extract_epa(rec), 'epa_release':EPA_RELEASE,'epa_file_sha256':EPA_SHA256})
    parents=[]
    for p in PARENTS:
        status,rec,n=unique_record(cytotox,p['casrn'])
        parents.append({**p,'match_status':status,'n_epa_rows_matching_cas':n,
                         **extract_epa(rec),'interpretation':'PARENT PLASTICIZER; NOT MEASURED NHANES URINARY METABOLITE; NEVER TRANSFER THIS BIOACTIVITY TO METABOLITE'})
    return pd.DataFrame(rows),pd.DataFrame(parents)


def optional_mc56(path:Path, rowdf:pd.DataFrame)->pd.DataFrame:
    """Stream only exact CAS records from full EPA mc5/6 CSV, if user supplies it.
    This file is inside a very large official INVITRODB_SUMMARY.zip; never download it automatically.
    """
    requested=set(rowdf.casrn)
    blocks=[];nrows=0; match_n=0
    for chunk in pd.read_csv(path, chunksize=75000,low_memory=False):
        chunk.columns=[norm(c) for c in chunk.columns]
        cascol=next((s for s in ['casn','casrn'] if s in chunk),None)
        hitcol=next((s for s in ['hitc','hitcall'] if s in chunk),None)
        if not cascol or not hitcol:raise ValueError(f'mc5/6 missing CAS or hitcall: {list(chunk.columns)}')
        nrows+=len(chunk)
        chunk['casrn_exact']=chunk[cascol].map(cas_clean)
        keep=chunk[chunk.casrn_exact.isin(requested)].copy()
        if not keep.empty:
            keep['hitcall_value']=pd.to_numeric(keep[hitcol],errors='coerce')
            keep['active_hit_ge_0_9']=keep['hitcall_value'].ge(.9)
            keep['source_release']=EPA_RELEASE
            blocks.append(keep)
            match_n+=len(keep)
    if not blocks:
        print(f'[MC56] Scanned {nrows:,} rows. No exact CAS matches: result remains NA, not zero activity.')
        return pd.DataFrame(columns=['casrn_exact','status','n_sample_endpoint_rows','n_unique_endpoints','n_active_sample_endpoint_rows'])
    full=pd.concat(blocks,ignore_index=True)
    full.to_csv(OUT/'v11_optional_exact_CAS_assay_records.csv',index=False)
    aeid='aeid' if 'aeid' in full else ('aenm' if 'aenm' in full else None)
    records=[]
    for cas,group in full.groupby('casrn_exact'):
        m=group.hitcall_value.notna()
        records.append(dict(casrn_exact=cas,status='EXACT_CAS_ASSAY_ROWS',n_sample_endpoint_rows=len(group),
            n_tested_rows=int(m.sum()), n_unique_endpoints=int(group[aeid].nunique()) if aeid else None,
            n_active_sample_endpoint_rows=int(group.loc[m,'active_hit_ge_0_9'].sum()),
            active_fraction_of_tested_rows=float(group.loc[m,'active_hit_ge_0_9'].mean()) if m.any() else None))
    print(f'[MC56] Scanned {nrows:,} rows, selected {match_n:,} exact-CAS assay rows')
    return pd.DataFrame(records)


def plot_evidence(summary:pd.DataFrame)->None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    df=summary.copy();df['rate']=pd.to_numeric(df.get('epa_nhit'),errors='coerce')/pd.to_numeric(df.get('epa_ntested'),errors='coerce')*100
    fig,axs=plt.subplots(1,2,figsize=(11.6,4.7),layout='constrained')
    x=np.arange(len(df));name=df['short'].to_list()
    vals=df['nhanes_weighted_pct_below_lod'].to_numpy(float)
    axs[0].bar(x,vals);axs[0].set_xticks(x,name);axs[0].set_ylabel('Below LOD (%)');axs[0].set_title('NHANES: survey-weighted detection-limit flags')
    rate=df['rate'].to_numpy(float)
    valid_rate=np.isfinite(rate)
    axs[1].bar(x[valid_rate],rate[valid_rate]);axs[1].set_xticks(x,name);axs[1].set_ylabel('Active burst assays (%)');axs[1].set_title('EPA cytotoxicity burst assay coverage (NOT a risk score)')
    for i,v in enumerate(rate):
        if not np.isfinite(v):axs[1].text(i,0.6,'NA',ha='center',fontsize=10)
    axs[1].set_ylim(bottom=0)
    fig.savefig(OUT/'v11_two_evidence_sources.png',dpi=170)
    plt.close(fig)


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cytotox',type=Path,help='Existing official EPA cytotox XLSX, optionally bypass network and checksum (e.g., for an external export)')
    ap.add_argument('--mc56',type=Path,help='Optional local EPA mc5-6 full CSV (large; not downloaded by default)')
    ap.add_argument('--pooled',type=Path,default=POOLED)
    args=ap.parse_args(argv)
    OUT.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(CHEMICALS).to_csv(OUT/'v11_verified_chemical_crosswalk.csv',index=False)
    if args.cytotox is None:
        path=fetch_epa()
        validated=True
    else:
        path=args.cytotox.resolve()
        validated=sha256(path)==EPA_SHA256
        if not validated:print('[WARNING] Non-official-checksum user file: provenance is NOT verified',file=sys.stderr)
    cdf,format_info=load_cytotox(path)
    profiles=nhanes_profiles(args.pooled)
    summary,parents=do_cytotox(cdf,profiles)
    summary.to_csv(OUT/'v11_exact_metabolites_and_PFAS.csv',index=False)
    parents.to_csv(OUT/'v11_parent_plasticizers_CONTEXT_ONLY.csv',index=False)
    optional_results=None
    if args.mc56:
        if not args.mc56.is_file():raise FileNotFoundError(args.mc56)
        optional_results=optional_mc56(args.mc56,summary)
        optional_results.to_csv(OUT/'v11_optional_assay_coverage.csv',index=False)
    if all(isinstance(r['nhanes_weighted_pct_below_lod'],(float,int)) for r in profiles):
        plot_evidence(summary)
    audit={'retrieved_utc':datetime.now(timezone.utc).isoformat(), 'epa_release':EPA_RELEASE,
           'official_cytotox_file':EPA_PAGE, 'file_sha256':sha256(path), 'verified_against_official_sha256':validated,
           'epa_schema':format_info, 'n_nhanes_biomarkers':7,'n_exact_matched':int(sum(summary.match_status.eq('EXACT_CAS_MATCH'))),
           'n_unmatched_or_ambiguous':int(sum(~summary.match_status.eq('EXACT_CAS_MATCH'))),
           'parent_chemical_proxy_prohibited':True,'impute_missing_assays_as_inactive':False,
           'standalone_human_risk_score_calculated':False, 'full_mc56_provided':args.mc56 is not None,
           'notes':['Exact CAS does not always establish comparable protonation/salt state; inspect chemical/sample QC before mechanism claims.',
                    'Cytotoxicity median at 1000 uM may be a default upper ceiling, not measured potency; flag is heuristic.',
                    'No biological testing in human participants was inferred from in vitro screening.',
                    'Urinary metabolites were not replaced with parent phthalate assay data.',
                    'Missing toxicity profile means not matched/not tested, not harmless.',
                    'Bioactivity data are chemical-level; never treat as individual-level patient measurements.',
                    'Unweighted compound assay proportions do not establish clinical risk.']}
    (OUT/'v11_provenance_and_qc.json').write_text(json.dumps(audit,indent=2,default=str),encoding='utf-8')
    print('\n[V11] EPA exact-identity match and NHANES linkage')
    print(summary[['short','casrn','match_status','nhanes_n','nhanes_weighted_pct_below_lod','epa_ntested','epa_nhit','epa_cytotox_median_um','cytotox_1000_uM_possible_default_ceiling']].to_string(index=False))
    print('\n[V11] Exact matches:',audit['n_exact_matched'],'/ 7; No chemical toxicity risk score was computed.')
    print('[saved]',OUT)
    if args.mc56: print('[assay summary]',len(optional_results),'exact-CAS chemicals')

if __name__=='__main__':main()
