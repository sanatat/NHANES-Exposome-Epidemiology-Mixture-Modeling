#!/usr/bin/env python3
"""V14: reproducible expert-priority panel (EPA ToxCast annotation evidence).

Inputs: V13 shortlist, optionally real EPA v4.3 mc5-6 CSV with casn/casrn,aeid,hitcall.
No per-chemical hits are inferred from assay names or target genes.
"""
import argparse, json, re
from pathlib import Path
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parent
SRC=ROOT/'results/v13_mechanism_evidence_audit/v13_target_gene_supported_shortlist.csv'
OUT=ROOT/'results/v14_mechanistic_panel'
CHEMS=[('PFHxS','355-46-4'),('PFNA','375-95-1'),('PFOA','335-67-1'),('PFOS','1763-23-1'),('MBP','131-70-4'),('MEP','2306-33-4'),('MiBP','30833-53-5')]
# Priority IDs reflect biological interpretability, not effect size or detected hits.
CORE={
210:('Inflammation','IL6 cytokine release readout'),
296:('Inflammation','TNF-alpha cytokine release readout; LPS context'),
284:('Inflammation','CXCL8/IL8 release readout; LPS context'),
286:('Inflammation','CCL2/MCP1 release readout; LPS context'),
1346:('Inflammation','NF-kappaB reporter (not direct protein binding)'),
94:('Inflammation','NF-kappaB cis-regulatory reporter'),
97:('Oxidative stress','ARE/NRF2-response reporter'),
1110:('Oxidative stress','ARE-response reporter; not direct NFE2L2 binding'),
1066:('Oxidative stress','HepaRG NFE2L2-associated readout'),
1070:('Oxidative stress','HepaRG NQO1-associated readout'),
108:('Inflammation','STAT3 cis-regulatory reporter'),
2612:('Inflammation','Human PTGS2 enzyme assay; inhibition/activity is not inflammation induction'),
}

def cas(x):
    if pd.isna(x): return ''
    s=str(x).strip();return s if re.fullmatch(r'\d{2,7}-\d{2}-\d',s) else ''

def normalized_id(v):
    if pd.isna(v):return None
    s=str(v).strip()
    try:return int(float(s)) if float(s).is_integer() else None
    except (ValueError,TypeError):return None

def prepare_panel(inp):
    df=pd.read_csv(inp,low_memory=False)
    need={'aeid','assay_component_endpoint_name','mapped_gene_symbols','domains_matched'}
    if not need.issubset(df.columns):raise ValueError(f'Expected columns: {sorted(need)}; found: {list(df.columns)}')
    df['aeid']=df.aeid.map(normalized_id)
    if df.aeid.isna().any() or df.aeid.duplicated().any():raise ValueError('AEID invalid or not unique in V13 shortlist')
    absent=sorted(set(CORE)-set(df.aeid))
    if absent:raise ValueError(f'Important CORE AEIDs absent in V13 shortlist: {absent}; refusing silent substitution')
    df['tier']=df.aeid.map(lambda x:'CORE_12' if x in CORE else 'EXTENDED_REVIEW')
    df['priority_topic']=df.aeid.map(lambda x: CORE[int(x)][0] if x in CORE else 'Secondary / context')
    df['interpretation_note']=df.aeid.map(lambda x:CORE[int(x)][1] if x in CORE else 'Review primary assay documentation before biological claims')
    return df.sort_values(['tier','aeid']).reset_index(drop=True)

def assay_hits(file_path,panel):
    """Consume only explicit CAS+AEID records. Do not treat missing as inactive."""
    cas_keep=set(c for _,c in CHEMS);aeids=set(map(int,panel.aeid))
    found=[];scanned=0;match=0
    for chunk in pd.read_csv(file_path,chunksize=100000,low_memory=False):
        chunk.columns=[str(c).strip().lower() for c in chunk.columns]
        ccol=next((c for c in ('casrn','casn') if c in chunk),None)
        ecol='aeid' if 'aeid' in chunk else None
        hcol=next((c for c in ('hitcall','hitc') if c in chunk),None)
        if ccol is None or ecol is None or hcol is None:
            raise ValueError(f'Official results need CAS (casn/casrn), AEID, hitcall; found {list(chunk.columns)[:35]}. If it has SPID instead of CAS, supply verified sample-to-CAS crosswalk; no guessing.')
        scanned+=len(chunk)
        chunk['_casrn']=chunk[ccol].map(cas);chunk['_aeid']=chunk[ecol].map(normalized_id)
        m=chunk['_casrn'].isin(cas_keep)&chunk['_aeid'].isin(aeids)
        if not m.any():continue
        cols=['_casrn','_aeid',hcol]
        for col in ('ac50','acc','bmd','spid','sample_id'):
            if col in chunk:cols.append(col)
        z=chunk.loc[m,cols].copy().rename(columns={'_casrn':'casrn','_aeid':'aeid',hcol:'hitcall'})
        z['hitcall']=pd.to_numeric(z.hitcall,errors='coerce')
        if (z.hitcall.dropna().lt(0)|z.hitcall.dropna().gt(1)).any():raise ValueError('Hitcall outside [0,1]; check release and schema.')
        found.append(z);match+=len(z)
    if not found: return pd.DataFrame(),{'rows_scanned':scanned,'matched_records':0}
    x=pd.concat(found,ignore_index=True)
    # Do not convert duplicates or missing hitcalls into negative activity; summarize replicate coverage.
    def summarize(g):
        valid=g.hitcall.dropna()
        ans={'n_records':len(g),'n_valid_hitcall':len(valid),
             'n_positive_hitcall_ge_0_9':int(valid.ge(.9).sum()),
             'max_hitcall':float(valid.max()) if len(valid) else np.nan}
        for col in ('ac50','acc','bmd'):
            if col in g:
                active=(g.hitcall>=.9)
                values=pd.to_numeric(g.loc[active,col],errors='coerce')
                ans['min_'+col+'_raw_active']=float(values.min()) if values.notna().any() else np.nan
        return pd.Series(ans)
    agg=x.groupby(['casrn','aeid'],sort=True).apply(summarize,include_groups=False).reset_index()
    return agg,{'rows_scanned':scanned,'matched_records':match,'matched_chemical_endpoint_pairs':len(agg)}

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('--shortlist',type=Path,default=SRC)
    p.add_argument('--mc56',type=Path,help='Optional official EPA v4.3 mc5-6 results CSV (never downloaded automatically)')
    p.add_argument('--out',type=Path,default=OUT)
    a=p.parse_args(argv)
    if not a.shortlist.exists():raise FileNotFoundError(f'First run v13: {a.shortlist}')
    panel=prepare_panel(a.shortlist);a.out.mkdir(parents=True,exist_ok=True)
    panel.to_csv(a.out/'v14_curated_44_panel.csv',index=False)
    panel[panel.tier=='CORE_12'].to_csv(a.out/'v14_core_12_panel.csv',index=False)
    cross=pd.DataFrame(CHEMS,columns=['chemical','casrn'])
    pairs=cross.assign(_j=1).merge(panel[['aeid','assay_component_endpoint_name','mapped_gene_symbols','domains_matched','tier','interpretation_note']].assign(_j=1),on='_j').drop(columns='_j')
    assert len(pairs)==len(CHEMS)*len(panel)
    audit={'n_chemicals':len(CHEMS),'n_endpoints':len(panel),'n_core':int((panel.tier=='CORE_12').sum()),
           'n_possible_pairs':len(pairs),'interpretation':'Candidate assay metadata; not per-chemical bioactivity',
           'metabolites_not_replaced_by_parents':True,'no_human_hazard_score':True}
    if a.mc56:
        if not a.mc56.is_file():raise FileNotFoundError(a.mc56)
        obs,qc=assay_hits(a.mc56,panel)
        audit.update(qc)
        if len(obs):
            pairs=pairs.merge(obs,how='left',on=['casrn','aeid'],validate='one_to_one')
            pairs['evidence_status']=np.where(pairs.n_records.notna(),
                np.where(pairs.n_valid_hitcall.fillna(0).gt(0),'OBSERVED_HITCALL','OBSERVED_NO_VALID_HITCALL'),
                'NOT_IN_SUPPLIED_RESULTS')
        else:
            pairs['evidence_status']='NOT_IN_SUPPLIED_RESULTS'
        print('[V14] Observed exact-CAS results:',qc)
    else:
        pairs['evidence_status']='NOT_RETRIEVED'
        print('[V14] Annotation-only matrix; no per-chemical activity inferred')
    pairs.to_csv(a.out/'v14_chemical_x_endpoint_evidence.csv',index=False)
    with (a.out/'v14_audit.json').open('w') as f:json.dump(audit,f,indent=2)
    print('[V14] 44 annotated endpoints, 12 core; 7 chemicals x',len(panel),'=',len(pairs),'rows')
    print('[V14] core panel:')
    print(panel.loc[panel.tier=='CORE_12',['aeid','assay_component_endpoint_name','mapped_gene_symbols']].to_string(index=False))
    print('[V14] saved:',a.out)
    return audit,pairs

if __name__=='__main__':main()
