#!/usr/bin/env python3
"""Build independent pooled NHANES 2015-16 + 2017-18 exposure analysis.

Uses original CDC XPTs. Checks exposure-specific subsample-B weights, preserves
original v5 data/results, identifies cycles and harmonizes hs-CRP as CDC advises.
No optimization against results or p-values. Explicit outcome scope: inflammation.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import sqlite3
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'data' / 'raw'
PROC = ROOT / 'data' / 'processed'
RESULTS = ROOT / 'results'
NAMES = ['DEMO','BMX','PFAS','PHTHTE','HSCRP','ALB_CR']
CHEM = ['LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB']
FLAGS = ['LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL','URDMBPLC','URDMEPLC','URDMIBLC']
OUT = PROC / 'analysis_joint_pooled_2015_2018.csv'
DB = ROOT / 'data' / 'nhanes_exposome_pooled_2015_2018.sqlite'

spec = importlib.util.spec_from_file_location('base_nhanes', ROOT / '01_download_build.py')
assert spec is not None and spec.loader is not None
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def download_xpt(stem: str, cycle_start: int) -> Path:
    RAW.mkdir(exist_ok=True, parents=True)
    suffix = 'I' if cycle_start == 2015 else 'J'
    fname = f'{stem}_{suffix}.xpt'
    path = RAW / fname
    url = f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle_start}/DataFiles/{fname}'
    if base.verify_binary_xpt(path):
        print(f'[cache] {fname} ({path.stat().st_size:,} bytes)', flush=True)
        return path
    last = None
    for trial in range(1, 5):
        partial = path.with_suffix('.xpt.part')
        try:
            print(f'[download] {url} (attempt {trial})', flush=True)
            with requests.get(url, timeout=(20, 90), stream=True,
                    headers={'User-Agent':'Mozilla/5.0 NHANES-independent-research/1.0'}) as r:
                r.raise_for_status()
                with partial.open('wb') as f:
                    for block in r.iter_content(262144):
                        if block:
                            f.write(block)
            if not base.verify_binary_xpt(partial):
                raise ValueError('Not a valid SAS XPT; check server response')
            partial.replace(path)
            return path
        except (requests.RequestException, ValueError, OSError) as exc:
            last = exc
            partial.unlink(missing_ok=True)
            if trial < 4:
                time.sleep(2*trial)
    raise RuntimeError(f'Cannot download {url}: {last}')


def read_cycle(cycle_start: int) -> tuple[pd.DataFrame, dict]:
    suffix = 'I' if cycle_start == 2015 else 'J'
    inputs: dict[str,pd.DataFrame] = {}
    for stem in NAMES:
        path = download_xpt(stem, cycle_start)
        d = pd.read_sas(path, format='xport', encoding='utf-8')
        d.columns = [str(col).upper() for col in d.columns]
        d = base.correct_sas_xport_zeros(d, f'{stem}_{suffix}')
        # Reuse the validated 0/1 flag assertions from v5.
        base.assert_cdc_lod_flags(d, stem + '_J')
        if 'SEQN' not in d or d.SEQN.isna().any() or d.SEQN.duplicated().any():
            raise ValueError(f'{stem}_{suffix}: SEQN missing or duplicated')
        d['SEQN'] = pd.to_numeric(d['SEQN'], errors='raise').astype('int64')
        inputs[stem + '_J'] = d  # normalize suffix in dict keys only

    _, _, joint, qc = base.assemble(inputs)
    if not qc['joint_weights_verified'] or joint.empty:
        raise ValueError(f'{cycle_start}: PFAS vs phthalate subsample-B weights did not match; no pooling')
    if len(joint) < 100:
        raise ValueError(f'{cycle_start}: too few complete valid joint cases')
    if any(c not in joint for c in CHEM + FLAGS):
        raise ValueError('Missing exposures or LOD flags after merge')
    for flag in FLAGS:
        if not joint[flag].isin([0., 1.]).all():
            raise ValueError(f'{cycle_start}: invalid {flag}')
    joint = joint.copy()
    joint['cycle_start'] = cycle_start
    joint['crp_original_mg_L'] = joint['LBXHSCRP']
    # CDC 2017-18 codebook: map DxC 660i (2015-16) -> Cobas 6000 (2017-18):
    # Cobas = 0.8695*DxC + 0.2954, valid only DxC <= 23 mg/L.
    # Values outside validated calibration range are NOT extrapolated.
    if cycle_start == 2015:
        inside = joint['LBXHSCRP'] <= 23
        joint['crp_harmonized_mg_L'] = np.where(
            inside, 0.8695 * joint['LBXHSCRP'] + 0.2954, np.nan)
    else:
        joint['crp_harmonized_mg_L'] = joint['LBXHSCRP']
    # Two 2-year cycles => four-year weight. v5 already verified the SAME WTSB2YR
    # value for each participant in PFAS & phthalate files.
    joint['survey_weight_2yr'] = joint['survey_weight']
    joint['survey_weight'] = joint['survey_weight_2yr'] / 2.0
    # Nested public-use strata/PSUs across the two 2-year NHANES cycles.
    joint['pooled_strata'] = (joint['SDMVSTRA'].astype(int) +
                               (10000 if cycle_start == 2017 else 0))
    joint['pooled_psu'] = joint['SDMVPSU'].astype(int)
    print(f'[cycle {cycle_start}] n_complete={len(joint)}, '
          f'calibration_outside_range={joint.crp_harmonized_mg_L.isna().sum()}, '
          f'weight_dif_max={qc["joint_weight_max_relative_difference"]:.3g}',
          flush=True)
    return joint, qc


def main() -> None:
    PROC.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    a, qc_i = read_cycle(2015)
    b, qc_j = read_cycle(2017)
    pooled = pd.concat([a, b], ignore_index=True)
    if pooled.SEQN.duplicated().any():
        raise ValueError('Pooled SEQN unexpectedly duplicated')
    if pooled.loc[pooled.survey_weight > 0, 'survey_weight'].isna().any():
        raise ValueError('Invalid survey weights')
    pooled.to_csv(OUT, index=False)
    with sqlite3.connect(DB) as conn:
        pooled.to_sql('analysis_joint_pooled', conn, if_exists='replace', index=False)
        conn.execute('CREATE INDEX IF NOT EXISTS idx_cycle_seqn ON analysis_joint_pooled(cycle_start, SEQN)')
    qc = {
        'cycles':['2015-2016','2017-2018'],
        'n_pooled_complete_pre_hscrp_restriction':int(len(pooled)),
        'n_by_cycle':{str(cycle):int(sum(pooled.cycle_start == cycle)) for cycle in (2015,2017)},
        'n_chronic_hscrp_le10':int(((pooled.crp_harmonized_mg_L > 0) &
            (pooled.crp_harmonized_mg_L <= 10)).sum()),
        'n_hscrp_harmonized_positive_le20':int(((pooled.crp_harmonized_mg_L > 0) &
            (pooled.crp_harmonized_mg_L <= 20)).sum()),
        'n_2015_above_bridge_limit':int(sum((pooled.cycle_start==2015) &
            pooled.crp_harmonized_mg_L.isna())),
        'survey_weight':'WTSB2YR / 2 with person-level PFAS vs phthalate weight match within each cycle',
        'survey_design':'pooled_strata = SDMVSTRA + 10000*indicator(cycle_start==2017); pooled_psu = SDMVPSU',
        'crp_bridge': '2015-2016 <=23 mg/L: 0.8695 * hsCRP + 0.2954; 2017-2018 unchanged',
        'not_for_causal_claims':True,
        'not_for_p_value_optimization':True,
        'source':'CDC NHANES public use 2015-2016 + 2017-2018',
        'cycle_qc':{'2015':qc_i,'2017':qc_j},
    }
    (PROC / 'v6_pooled_qc.json').write_text(json.dumps(qc, indent=2), encoding='utf-8')
    print(f'\n[done] {OUT}\n[done] {DB}')
    print('[QC] pooled subjects:',len(pooled))
    print('[QC] primary hs-CRP <= 10 mg/L:',qc['n_chronic_hscrp_le10'])
    print('[next] Rscript 09_pooled_survey_mixture.R')

if __name__ == '__main__':
    main()
