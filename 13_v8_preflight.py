#!/usr/bin/env python3
"""Read-only data QA for v8; does not download/overwrite any data."""
from pathlib import Path
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
p = root / 'data' / 'processed' / 'analysis_joint_pooled_2015_2018.csv'
if not p.is_file():
    raise SystemExit('[ERROR] Missing pooled CSV. Run .venv/bin/python 08_build_pooled_2015_2018.py first.')
d = pd.read_csv(p)
chem = ['LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB']
flags = ['LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL','URDMBPLC','URDMEPLC','URDMIBLC']
need = ['SEQN','cycle_start','survey_weight','pooled_strata','pooled_psu',
    'crp_harmonized_mg_L','RIDAGEYR','RIAGENDR','RIDRETH3','INDFMPIR','BMXBMI','URXUCR',*chem,*flags]
missing = sorted(set(need).difference(d.columns))
if missing: raise SystemExit(f'[ERROR] Missing columns: {missing}')
if d.SEQN.duplicated().any(): raise SystemExit('[ERROR] Duplicate SEQN')
if set(d.cycle_start.dropna().unique()) != {2015,2017}: raise SystemExit('[ERROR] Invalid cycle composition')
if not (np.isfinite(d.survey_weight).all() and d.survey_weight.gt(0).all()):
    raise SystemExit('[ERROR] Invalid survey weights')
for f in flags:
    values = set(d[f].dropna().unique())
    if not values.issubset({0,1}): raise SystemExit(f'[ERROR] Invalid LOD values in {f}: {values}')

# Mirror primary inclusion rules of v6/v8. Doesn't estimate p-values.
cc = d[need].notna().all(axis=1)
m = d.loc[cc & d.RIDAGEYR.ge(20) & d.URXUCR.gt(0) & d.URXMBP.gt(0)
     & d.crp_harmonized_mg_L.gt(0) & d.crp_harmonized_mg_L.le(20)].copy()
pri = m.loc[m.crp_harmonized_mg_L.le(10)].copy()
if len(pri)<1500: raise SystemExit('[ERROR] Fewer than 1500 primary participants')
if not (np.isfinite(pri.URXMBP).all() and np.isfinite(pri.URXUCR).all()):
    raise SystemExit('[ERROR] Invalid raw MBP or creatinine')
ratio = 100000 * pri.URXMBP / pri.URXUCR
if not (ratio.gt(0).all() and np.isfinite(ratio).all()): raise SystemExit('[ERROR] Invalid creatinine ratio')
if not np.allclose(np.log2(ratio), np.log2(pri.URXMBP) - np.log2(pri.URXUCR) + np.log2(100000)):
    raise SystemExit('[ERROR] Creatinine unit conversion failed')
print('[PASS] Pooled input: unique people, cycles, seven exposures, valid LOD flags')
print(f'[PASS] Primary <=10 mg/L: {len(pri):,} people; sensitivity <=20 mg/L: {len(m):,} people')
print(f'[PASS] MBP ng/g creatinine = 100000 * (MBP ng/mL)/(creatinine mg/dL)')
print('[PASS] No files modified; ready for Rscript 12_creatinine_robustness_splines.R')
