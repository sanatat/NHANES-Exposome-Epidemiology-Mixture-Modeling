#!/usr/bin/env python3
"""Offline correctness checks. Synthetic data only; not a finding about NHANES."""
import importlib.util
from pathlib import Path
import tempfile
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent

def load(path):
    spec=importlib.util.spec_from_file_location('module_'+path.replace('.','_'),ROOT/path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

v10=load('18_v10_paired_transfer.py')
v9=v10.v9_loader()
assert {'Ridge','ElasticNet','RandomForest'}==set(v10.estimator_name(k) and k for k in ['Ridge','ElasticNet','RandomForest'])
rng=np.random.default_rng(37)
rows=[]
for cycle in (2015,2017):
    for i in range(120):
        chems=np.exp(rng.normal(1,.3,7))
        row={
            'SEQN':cycle*1000+i,'cycle_start':cycle,
            'survey_weight':float(rng.uniform(.5,5)),
            'pooled_strata':int(cycle*100+i//10),
            'pooled_psu':int(1+i%2),
            'RIDAGEYR':float(rng.uniform(20,70)),
            'BMXBMI':float(rng.uniform(18,41)),
            'INDFMPIR':float(rng.uniform(.1,5)),
            'RIAGENDR':1+i%2,'RIDRETH3':1+i%4,
            'URXUCR':float(rng.uniform(50,170)),
        }
        for c,x in zip(v9.EXPOSURES,chems): row[c]=float(x)
        y=.02*(row['RIDAGEYR']-40)+.08*(np.log2(chems[4])-1.5)+rng.normal(0,.4)
        row['crp_harmonized_mg_L']=float(np.exp(y))
        rows.append(row)

with tempfile.TemporaryDirectory() as t:
    table=v10.paired_comparison(pd.DataFrame(rows),v9,output_dir=Path(t),families=('Ridge',),quiet=True)
    assert len(table)==2
    assert table.n_test.eq(120).all()
    assert np.allclose(table['covariates_only_weighted_mse']-table['plus_chemicals_weighted_mse'],
                       table['delta_mse_positive_is_better'])
    preds=pd.read_csv(Path(t)/'v10_matched_algorithm_predictions.csv')
    assert len(preds)==480
    assert preds.SEQN.nunique()==240
    for (_,_,s), group in preds.groupby(['SEQN','test_cycle','model_family']):
        assert set(group.feature_set)=={'covariates','covariates_plus_7_chemicals'}
    assert (Path(t)/'v10_interpretation.txt').exists()
print('[PASS] Matched algorithm, held-out cycles, paired prediction rows, delta MSE (synthetic)')
