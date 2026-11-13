#!/usr/bin/env python3
"""Offline integration smoke test of the NHANES 2015+2017 builder. No internet."""
import importlib.util
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('poolbuilder', ROOT / '08_build_pooled_2015_2018.py')
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fake_data(cycle):
    rng = np.random.default_rng(220 + cycle)
    ids = np.arange(50000 if cycle == 2015 else 60000, 50240 if cycle == 2015 else 60240)
    n = len(ids)
    df = {'DEMO':pd.DataFrame({'SEQN':ids,'RIDAGEYR':rng.integers(20,71,n),
          'RIAGENDR':rng.integers(1,3,n),'RIDRETH3':rng.integers(1,7,n),
          'INDFMPIR':rng.uniform(.2,5,n),'SDMVPSU':np.tile([1,2],n//2),
          'SDMVSTRA':np.repeat(np.arange(1,13),20)}),
          'BMX':pd.DataFrame({'SEQN':ids,'BMXBMI':rng.uniform(18,40,n)}),
          'HSCRP':pd.DataFrame({'SEQN':ids,'LBXHSCRP':rng.uniform(.1,7,n),
              'LBDHRPLC':np.zeros(n)}),
          'ALB_CR':pd.DataFrame({'SEQN':ids,'URXUCR':rng.uniform(25,200,n)})}
    weights = rng.uniform(100,5000,n)
    if cycle == 2015:
        df['HSCRP'].loc[0,'LBXHSCRP'] = 30 # bridge out of range
    pf = pd.DataFrame({'SEQN':ids,'WTSB2YR':weights})
    ph = pd.DataFrame({'SEQN':ids,'WTSB2YR':weights})
    for chem, flag in zip(m.CHEM,m.FLAGS):
        dest = pf if chem.startswith('LBX') else ph
        dest[chem] = rng.uniform(.1,12,n)
        dest[flag] = np.zeros(n)
        # Deliberately inject SAS XPORT 0 artifact; decoder should restore 0.
        dest.loc[0,flag] = m.base.SAS_XPORT_ZERO_ARTIFACT
        dest.loc[1,flag] = 1
    df['PFAS'] = pf; df['PHTHTE'] = ph
    return df

frames = {2015:fake_data(2015),2017:fake_data(2017)}
original_read_sas = pd.read_sas

def read_fake(path, **kwargs):
    stem, suffix = path.stem.rsplit('_',1)
    cycle = 2015 if suffix=='I' else 2017
    return frames[cycle][stem].copy()

with patch.object(m,'download_xpt',side_effect=lambda stem,year:Path(f'/offline/{stem}_{"I" if year==2015 else "J"}.xpt')):
    with patch.object(pd,'read_sas',side_effect=read_fake):
        a,qa = m.read_cycle(2015)
        b,qb = m.read_cycle(2017)
assert len(a)==240 and len(b)==240
assert qa['joint_weights_verified'] and qb['joint_weights_verified']
assert a.loc[a.SEQN==50000,'crp_harmonized_mg_L'].isna().all()
assert abs(a.loc[a.SEQN==50001,'crp_harmonized_mg_L'].iat[0] -
           (0.8695 * frames[2015]['HSCRP'].loc[1,'LBXHSCRP']+0.2954)) < 1e-10
assert b.loc[b.SEQN==60001,'crp_harmonized_mg_L'].iat[0] == frames[2017]['HSCRP'].loc[1,'LBXHSCRP']
assert (a['survey_weight']==a['survey_weight_2yr']/2).all()
assert set(a.pooled_strata).isdisjoint(set(b.pooled_strata))
for flag in m.FLAGS:
    assert a[flag].isin([0,1]).all() and b[flag].isin([0,1]).all()
print('[PASS] 2-cycle merge, XPT zeros, B weights, calibration boundary, strata separation')
