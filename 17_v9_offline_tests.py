#!/usr/bin/env python3
"""Independent synthetic tests; not a test of scientific findings."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent

def load(name):
    p=ROOT/name
    spec=importlib.util.spec_from_file_location(name.replace('.','_'),p)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

a=load('15_v9_exposure_atlas.py'); p=load('16_v9_crosscycle_prediction.py')
assert abs(a.weighted_corr([1,2,3],[2,4,6],[1,1,1])-1)<1e-12
assert abs(a.weighted_spearman([1,2,3],[3,2,1],[1,1,1])+1)<1e-12
assert p.weighted_mse([1,2],[1,3],[1,3])==.75
rng=np.random.default_rng(42); rows=[]
for cyc in [2015,2017]:
    for i in range(110):
        exposures=np.exp(rng.normal(1,.35,7))
        row={'SEQN':cyc*1000+i,'cycle_start':cyc,'survey_weight':float(rng.uniform(1,4)),
            'pooled_strata':int(cyc*100+i//8),'pooled_psu':int(1+i%2),
            'RIDAGEYR':float(rng.uniform(20,70)),'RIAGENDR':1+(i%2),
            'RIDRETH3':1+i%4,'BMXBMI':float(rng.uniform(18,40)),
            'INDFMPIR':float(rng.uniform(0,5)),'URXUCR':float(rng.uniform(40,160)),
            'crp_harmonized_mg_L':float(np.exp(rng.normal(.2,.5)))}
        for c,x in zip(a.EXPOSURES,exposures): row[c]=float(x)
        for f in a.FLAGS: row[f]=0
        rows.append(row)
d=pd.DataFrame(rows)
t=p.prep_data(d)
assert len(t)==220 and len(t[t.cycle_start==2015])==110
features=p.CONTROLS+p.CATEGORICAL
model=p.build_pipeline('Ridge (covariates)',features,{'alpha':1.})
model.fit(t.loc[:109,features],t.loc[:109,'log_hscrp'],model__sample_weight=t.loc[:109,'survey_weight'])
assert len(model.predict(t.loc[110:,features]))==110
assert len(set(t.loc[t.cycle_start==2015,'SEQN']) & set(t.loc[t.cycle_start==2017,'SEQN']))==0
print('[PASS] Weighted correlation, model pipelines, metadata and nonoverlapping cycles (synthetic only)')
