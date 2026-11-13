#!/usr/bin/env python3
"""Independent pandas LOD-flag audit for the same NHANES joint cohort.
Run: python3 06_lod_audit_python.py
Optional test: python3 06_lod_audit_python.py --input /path/data.csv --output /path/qc.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
CHEM = ['LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB']
FLAGS = ['LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL','URDMBPLC','URDMEPLC','URDMIBLC']

def audit(df):
    required = ['survey_weight',*FLAGS]
    missing = sorted(set(required).difference(df.columns))
    if missing: raise ValueError(f'Missing columns: {missing}')
    w = pd.to_numeric(df['survey_weight'], errors='coerce')
    valid_w = np.isfinite(w) & (w > 0)
    if not valid_w.all(): raise ValueError('Invalid survey_weight in input')
    out=[]
    for name, flag in zip(CHEM, FLAGS):
        x = pd.to_numeric(df[flag], errors='coerce')
        ok = x.isin([0,1])
        above = ok & x.eq(0)
        below = ok & x.eq(1)
        count0, count1, countok = int(above.sum()), int(below.sum()), int(ok.sum())
        weighted = (100 * w[below].sum() / w[ok].sum()) if countok else float('nan')
        out.append(dict(chemical=name, flag=flag, n_total=len(df), n_detected_0=count0,
                        n_below_LOD_1=count1, n_flag_missing=int(x.isna().sum()),
                        n_flag_invalid=int((~ok & x.notna()).sum()), n_valid_flags=countok,
                        pct_below_LOD_unweighted=100*count1/countok if countok else float('nan'),
                        pct_below_LOD_weighted=weighted))
    return pd.DataFrame(out)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',default=str(ROOT/'data/processed/analysis_joint.csv'))
    p.add_argument('--output',default=str(ROOT/'results/v4_lod_audit_python.csv'))
    a=p.parse_args()
    df=pd.read_csv(a.input)
    q=audit(df)
    o=Path(a.output)
    o.parent.mkdir(parents=True,exist_ok=True)
    q.to_csv(o,index=False)
    print(q.to_string(index=False))
    if (q['n_valid_flags'] != len(df)).any():
        print('[ERROR] LOD flags are invalid; data may contain SAS XPORT encoded zero artifacts.')
        for flag in FLAGS:
            x = pd.to_numeric(df[flag], errors='coerce')
            bad = x.loc[~x.isin([0, 1]) & x.notna()]
            if len(bad): print(f'  {flag} unexpected examples: {bad.unique()[:3]}')
        raise SystemExit('STOP: Fix zero artifacts and rebuild dataset with 01_download_build.py.')
    print('[PASS] All flags are valid CDC 0/1 codes; denominator is correct.')
    print('[saved]',o)
if __name__=='__main__': main()
