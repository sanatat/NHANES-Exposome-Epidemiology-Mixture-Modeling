#!/usr/bin/env python3
"""V13: audit EPA v4.3 assay candidates; annotation-only, no chemical hits.
Inputs: outputs from 23_v12_mechanistic_toxcast.py, no new download.
"""
import argparse, json, re
from pathlib import Path
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parent
DOMAINS={
    'Inflammatory signaling': ('inflammatory_signaling_candidate',
        {'NFKB1','NFKB2','RELA','RELB','REL','TNF','IL1B','IL6','TLR2','TLR3','TLR4','TLR7','TLR9','STAT3','JAK1','JAK2','JAK3','PTGS2'}),
    'Oxidative stress': ('oxidative_stress_candidate',
        {'NFE2L2','KEAP1','NQO1','HMOX1','SOD1','SOD2','SOD3','GPX1','GPX2','GPX3','GCLC','GCLM'}),
    'Immune response': ('immune_response_candidate',
        {'IFNG','IFNB1','IFNA1','IRF1','IRF2','IRF3','IRF4','IRF5','IRF6','IRF7','IRF8','IRF9','IL2','IL4','IL10','IL12A','IL12B','IL17A','CXCL8','CCL2','CCL5','STAT1','STAT2','TNF','IL6'}),
}

def to_bool(x):
    if pd.isna(x): return False
    if isinstance(x,(bool,np.bool_)):return bool(x)
    return str(x).strip().lower() in ('true','1','yes','t')

def get_target_symbols(row, cols):
    syms=set()
    for col in cols:
        val=row[col]
        if pd.isna(val): continue
        # The EPA target mapping can hold several explicit gene symbols joined by ';'.
        # DO NOT parse broad descriptive text or assay names into putative gene symbols.
        syms |= {s.strip().upper() for s in re.split(r'[;,|]',str(val)) if re.fullmatch(r'[A-Za-z][A-Za-z0-9-]*',s.strip())}
    return syms

def analyze(path,outdir):
    if not path.is_file():raise FileNotFoundError(path)
    df=pd.read_csv(path,low_memory=False)
    flags=[v[0] for v in DOMAINS.values()]
    missing=[x for x in flags if x not in df.columns]
    if missing:raise ValueError(f'Missing v12 classification columns: {missing}')
    idc=next((c for c in ('aeid','aenm','assay_component_endpoint_name') if c in df.columns),None)
    if not idc:raise ValueError('No AEID or endpoint name in inventory')
    df=df[df[idc].notna()].copy()
    df[idc]=df[idc].astype(str).str.strip().str.replace(r'\\.0$','',regex=True)
    df=df[df[idc].ne('')].drop_duplicates(idc,keep='first')
    for flag in flags:df[flag]=df[flag].map(to_bool)
    candidates=df[df[flags].any(axis=1)].copy()
    # Explicit gene symbols, not inferred from keyword occurrences in names/descriptions.
    gene_cols=[c for c in candidates.columns if c in ('official_symbol','gene_symbol') or c.endswith('_official_symbol') or c.endswith('_gene_symbol')]
    if gene_cols:
        sym=candidates.apply(lambda r:get_target_symbols(r,gene_cols),axis=1)
    else:
        sym=pd.Series([set() for _ in range(len(candidates))],index=candidates.index)
    candidates['mapped_gene_symbols']=sym.map(lambda x:';'.join(sorted(x)))
    candidates['domains_matched']=candidates.apply(lambda r:'; '.join(k for k,(f,_) in DOMAINS.items() if r[f]),axis=1)
    candidates['n_domains']=candidates[flags].sum(axis=1)
    for domain,(f,genes) in DOMAINS.items():
        new=domain.lower().replace(' ','_')+'_target_gene_support'
        candidates[new]=[bool(r[f] and (s&genes)) for (_,r),s in zip(candidates.iterrows(),sym)]
    support_cols=[x for x in candidates.columns if x.endswith('_target_gene_support')]
    candidates['target_gene_supported']=candidates[support_cols].any(axis=1)
    candidates['evidence_level']=np.where(candidates['target_gene_supported'],'explicit_target_gene_on_curated_shortlist','keyword_annotation_only')
    namecols=[x for x in ['aeid','aenm','assay_component_endpoint_name','assay_name','assay_desc','official_symbol','official_full_name'] if x in candidates.columns]
    keep=list(dict.fromkeys(namecols+['mapped_gene_symbols','domains_matched','n_domains','target_gene_supported','evidence_level']+flags+support_cols))
    candidates=candidates[keep].sort_values(['target_gene_supported','n_domains'],ascending=[False,False])
    outdir.mkdir(parents=True,exist_ok=True)
    candidates.to_csv(outdir/'v13_all_tagged_endpoints.csv',index=False)
    candidates[candidates.target_gene_supported].to_csv(outdir/'v13_target_gene_supported_shortlist.csv',index=False)
    candidates[~candidates.target_gene_supported].to_csv(outdir/'v13_keyword_only_review.csv',index=False)
    data=[]
    for domain,(f,_) in DOMAINS.items():
        sub=candidates[candidates[f]]
        col=domain.lower().replace(' ','_')+'_target_gene_support'
        data.append({'domain':domain,'n_tagged':len(sub),'n_target_gene_supported':int(sub[col].sum()),'n_annotation_only':int((~sub[col]).sum())})
    domain_stats=pd.DataFrame(data)
    domain_stats.to_csv(outdir/'v13_domain_evidence_summary.csv',index=False)
    subset_counts=[]
    for combo in range(1,8):
        subset=[name for j,name in enumerate(DOMAINS) if (combo>>j)&1]
        mask=np.ones(len(candidates),dtype=bool)
        for j,(_, (f,_)) in enumerate(DOMAINS.items()):mask &= (candidates[f].to_numpy()==bool((combo>>j)&1))
        subset_counts.append({'domain_combination':' + '.join(subset),'n_unique_endpoints':int(mask.sum())})
    overlap=pd.DataFrame(subset_counts)
    overlap.to_csv(outdir/'v13_disjoint_domain_overlap.csv',index=False)
    n=len(candidates)
    report={'n_rows_v12':len(df),'n_unique_tagged_endpoints':n,'n_unique_target_gene_supported':int(candidates.target_gene_supported.sum()),'n_unique_annotation_only':int((~candidates.target_gene_supported).sum()),'n_sum_of_domain_tags':int(candidates.n_domains.sum()),'n_overlaps_extra_tags':int(candidates.n_domains.sum()-n),'gene_symbol_columns':gene_cols,'status':'ANNOTATION_CANDIDATES_ONLY_NO_CHEMICAL_HITS','human_risk_scores_computed':False,'caution':'Curated gene lists are a preliminary triage, not biological pathway confirmation; coverage is not activity.'}
    (outdir/'v13_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,4.4))
    y=np.arange(len(domain_stats))
    ax.barh(y,domain_stats.n_target_gene_supported,label='Explicit target gene on shortlist')
    ax.barh(y,domain_stats.n_annotation_only,left=domain_stats.n_target_gene_supported,label='Keyword annotation only')
    ax.set_yticks(y,domain_stats.domain)
    ax.invert_yaxis();ax.set_xlabel('Tagged assay endpoints (categories may overlap)')
    ax.legend(loc='lower right',fontsize=8);ax.set_title('EPA v4.3 candidate annotations — not chemical bioactivity')
    fig.tight_layout();fig.savefig(outdir/'v13_annotation_evidence.png',dpi=160);plt.close(fig)
    print('[V13] unique tagged endpoints:',n,'; sum of 3 category counts:',int(candidates.n_domains.sum()))
    print('[V13] target-gene shortlist:',report['n_unique_target_gene_supported'],'; annotation-only:',report['n_unique_annotation_only'])
    print('[V13] exact domain stats:\n'+domain_stats.to_string(index=False))
    print('[V13] domain overlap:\n'+overlap.to_string(index=False))
    print('[V13] NO TOXCAST HITS INFERRED; saved',outdir)
    return report, domain_stats, overlap

def main():
    q=argparse.ArgumentParser()
    q.add_argument('--inventory',type=Path,default=ROOT/'results/v12_mechanistic_exposomics/v12_annotated_assay_inventory.csv')
    q.add_argument('--out',type=Path,default=ROOT/'results/v13_mechanism_evidence_audit')
    a=q.parse_args();analyze(a.inventory,a.out)
if __name__=='__main__':main()
