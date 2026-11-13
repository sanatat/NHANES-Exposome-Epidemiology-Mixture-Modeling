#!/usr/bin/env python3
"""v9: Survey-descriptive two-cycle co-exposure atlas (NOT causal inference).

Measures within-cycle weighted Spearman correlations, chemical distributions,
and edges reproducible in both cycles. Never selects features using outcomes.
Requires the v6 pooled CSV; never overwrites earlier results.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'data/processed/analysis_joint_pooled_2015_2018.csv'
OUT = ROOT / 'results/v9_professor_inspired'
EXPOSURES = ['LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB']
FLAGS = ['LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL','URDMBPLC','URDMEPLC','URDMIBLC']
LABELS = {'LBXPFHS':'PFHxS','LBXPFNA':'PFNA','LBXNFOA':'PFOA','LBXNFOS':'PFOS',
          'URXMBP':'MBP','URXMEP':'MEP','URXMIB':'MiBP'}

def weighted_corr(x,y,w):
    x=np.asarray(x,dtype=float); y=np.asarray(y,dtype=float); w=np.asarray(w,dtype=float)
    if len(x)!=len(y) or len(x)!=len(w) or len(x)<3 or not (np.isfinite(x).all() and np.isfinite(y).all() and np.isfinite(w).all()) or np.any(w<=0):
        raise ValueError('Invalid weighted correlation inputs')
    w=w/np.sum(w)
    cx=np.sum(x*w); cy=np.sum(y*w)
    den=np.sqrt(np.sum(w*(x-cx)**2)*np.sum(w*(y-cy)**2))
    return float(np.sum(w*(x-cx)*(y-cy))/den) if den>0 else float('nan')

def weighted_midrank(x,w):
    # Weighted fractional midranks in [0,1]. Ties get the same midpoint.
    x=np.asarray(x,dtype=float); w=np.asarray(w,dtype=float)
    order=np.argsort(x,kind='mergesort'); xs=x[order]; ws=w[order]
    ends=np.r_[np.flatnonzero(xs[1:]!=xs[:-1])+1, len(xs)]
    ranks=np.empty(len(xs),dtype=float); cum=0.0
    for end,start in zip(ends,np.r_[0,ends[:-1]]):
        mass=float(np.sum(ws[start:end])); ranks[start:end]=cum+mass/2; cum+=mass
    inverse=np.empty(len(xs)); inverse[order]=ranks/cum
    return inverse

def weighted_spearman(x,y,w):
    # Weighted Pearson on weighted midranks; descriptive, no design-based CI.
    return weighted_corr(weighted_midrank(x,w),weighted_midrank(y,w),w)

def weighted_quantile(x,w,quantile):
    x=np.asarray(x,float); w=np.asarray(w,float)
    order=np.argsort(x,kind='mergesort')
    c=np.cumsum(w[order])/np.sum(w)
    return float(x[order][np.searchsorted(c,quantile,side='left')])

def build_atlas(d):
    required=['cycle_start','SEQN','survey_weight','crp_harmonized_mg_L','RIDAGEYR']+EXPOSURES+FLAGS
    missing=set(required)-set(d.columns)
    if missing: raise ValueError(f'Missing columns: {sorted(missing)}')
    d=d.loc[(d.RIDAGEYR>=20) & d.crp_harmonized_mg_L.gt(0) & d.crp_harmonized_mg_L.le(10)].copy()
    if not set(d.cycle_start.unique())=={2015,2017}: raise ValueError('Both NHANES cycles required')
    if d.SEQN.duplicated().any(): raise ValueError('Duplicate participants')
    if d[EXPOSURES].isna().any().any() or (d[EXPOSURES]<=0).any().any(): raise ValueError('Bad exposure')
    if not d[FLAGS].isin([0,1]).all().all(): raise ValueError('Bad LOD flag')
    OUT.mkdir(parents=True,exist_ok=True)
    profiles=[]; matrices={}; correlations=[]
    for cycle in [2015,2017]:
        sub=d.loc[d.cycle_start==cycle]
        w=sub.survey_weight.to_numpy(float)
        matrix=np.empty((len(EXPOSURES),len(EXPOSURES)))
        for i,e in enumerate(EXPOSURES):
            profiles.append({'cycle':f'{cycle}-{cycle+1}','chemical':LABELS[e], 'analyte':e,
                'matrix':'serum' if i<4 else 'urine',
                'n':len(sub), 'weighted_median':weighted_quantile(sub[e],w,.5),
                'weighted_p25':weighted_quantile(sub[e],w,.25),
                'weighted_p75':weighted_quantile(sub[e],w,.75),
                'pct_below_LOD_weighted':100*np.sum(w*(sub[FLAGS[i]].to_numpy()==1))/np.sum(w)})
            for j,f in enumerate(EXPOSURES):
                matrix[i,j]=weighted_spearman(sub[e],sub[f],w)
                if j>i:
                    correlations.append({'cycle':f'{cycle}-{cycle+1}',
                        'chemical_a':LABELS[e], 'chemical_b':LABELS[f], 'rho_weighted_spearman':matrix[i,j]})
        matrices[cycle]=matrix
        plt.figure(figsize=(8,6))
        plt.imshow(matrix, vmin=-1,vmax=1,cmap='coolwarm')
        plt.colorbar(label='Weighted Spearman rho (descriptive)')
        plt.xticks(range(7),[LABELS[e] for e in EXPOSURES],rotation=35,ha='right')
        plt.yticks(range(7),[LABELS[e] for e in EXPOSURES])
        plt.title(f'Co-exposure atlas: NHANES {cycle}-{cycle+1} (n={len(sub):,})')
        plt.tight_layout()
        plt.savefig(OUT/f'v9_weighted_correlation_{cycle}.png',dpi=180)
        plt.close()
    # Edges require same direction and |rho|>=0.25 in BOTH cycles (fixed rule).
    edges=[]
    for i in range(len(EXPOSURES)):
        for j in range(i+1,len(EXPOSURES)):
            a,b=matrices[2015][i,j],matrices[2017][i,j]
            edges.append({'chemical_a':LABELS[EXPOSURES[i]],'chemical_b':LABELS[EXPOSURES[j]],
                          'rho_2015':a,'rho_2017':b,
                          'stable_edge_descriptive':bool(np.isfinite(a) and np.isfinite(b) and a*b>0 and min(abs(a),abs(b))>=.25)})
    pd.DataFrame(profiles).to_csv(OUT/'v9_exposure_profiles_by_cycle.csv',index=False)
    pd.DataFrame(correlations).to_csv(OUT/'v9_coexposure_correlations.csv',index=False)
    pd.DataFrame(edges).to_csv(OUT/'v9_coexposure_stable_edges.csv',index=False)
    print('[V9 ATLAS] primary n:',len(d),'cycles:',d.cycle_start.value_counts().to_dict())
    print('[V9 ATLAS] stable edges:',sum(r['stable_edge_descriptive'] for r in edges),'out of',len(edges))
    print('[V9 ATLAS] files in:',OUT)

if __name__=='__main__':
    if not SOURCE.exists(): raise SystemExit('Missing v6 pooled CSV; run 08_build_pooled_2015_2018.py')
    build_atlas(pd.read_csv(SOURCE))
