#!/usr/bin/env python3
"""NHANES 2017-2018 exposure database: download CDC XPT, validate, and make SQLite/CSV.
No synthetic results, no cancer outcome, no assumptions about group sizes.
"""
from __future__ import annotations

import hashlib
import json
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'data' / 'raw'
PROC = ROOT / 'data' / 'processed'
DB = ROOT / 'data' / 'nhanes_exposome_2017_2018.sqlite'
BASE = 'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/'
FILES = ['DEMO_J', 'BMX_J', 'PFAS_J', 'PHTHTE_J', 'HSCRP_J', 'ALB_CR_J']
PFAS = ['LBXPFHS', 'LBXPFNA', 'LBXNFOA', 'LBXNFOS']
PHTH = ['URXMBP', 'URXMEP', 'URXMIB']
# CDC flags: 0 detected, 1 below LOD; values below LOD are provided as LOD/sqrt(2).
PFAS_FLAGS = ['LBDPFHSL', 'LBDPFNAL', 'LBDNFOAL', 'LBDNFOSL']
PHTH_FLAGS = ['URDMBPLC', 'URDMEPLC', 'URDMIBLC']
DEMOGRAPHICS = ['SEQN', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH3', 'INDFMPIR',
                'SDMVPSU', 'SDMVSTRA']
COVARIATES = ['RIDAGEYR', 'RIAGENDR', 'RIDRETH3', 'INDFMPIR', 'BMXBMI']
OUTCOME = 'LBXHSCRP'
# pandas SAS XPORT parser can decode the IBM floating-point all-zero bytes as
# this tiny (but nonzero) number. It is NOT a laboratory measurement.
# See pandas-dev/pandas#30051; affects NHANES comment-code 0 = detected.
SAS_XPORT_ZERO_ARTIFACT = 5.397605346934028e-79


def correct_sas_xport_zeros(df: pd.DataFrame, source: str = '') -> pd.DataFrame:
    """Replace only the exact known XPORT zero artifact, never genuine small values."""
    numeric = df.select_dtypes(include=[np.number]).columns
    if len(numeric) == 0:
        return df
    matches = df.loc[:, numeric].eq(SAS_XPORT_ZERO_ARTIFACT)
    n = int(matches.to_numpy().sum())
    if n:
        df = df.copy()
        df.loc[:, numeric] = df.loc[:, numeric].mask(matches, 0.0)
        print(f'[xport-zero-correction] {source}: restored {n:,} encoded zero values')
    return df


def assert_cdc_lod_flags(df: pd.DataFrame, source: str) -> None:
    """Fail early rather than silently reporting 100% non-detects."""
    columns = (PFAS_FLAGS if source == 'PFAS_J' else
               PHTH_FLAGS if source == 'PHTHTE_J' else [])
    for col in columns:
        if col not in df:
            raise ValueError(f'{source}: flag {col} is missing')
        x = pd.to_numeric(df[col], errors='coerce')
        valid = x.dropna()
        bad = valid.loc[~valid.isin([0.0, 1.0])]
        if len(bad):
            raise ValueError(f'{source}.{col} has unexpected flag values: '
                             f'{bad.unique()[:5]!r}; verify the XPT reader')
        if len(valid) and not (valid == 0).any():
            raise ValueError(f'{source}.{col}: zero detected values: '
                             'check for uncorrected pandas XPT zero artifact')


def verify_binary_xpt(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 1000:
        return False
    with path.open('rb') as fh:
        head = fh.read(256)
    return b'HEADER RECORD' in head and b'LIBRARY' in head


def download(name: str) -> dict:
    RAW.mkdir(parents=True, exist_ok=True)
    path = RAW / f'{name}.xpt'
    url = BASE + path.name
    if verify_binary_xpt(path):
        print(f'[cache] {path.name} ({path.stat().st_size:,} bytes)')
    else:
        if path.exists():
            path.unlink()
        print(f'[download] {url}', flush=True)
        last_error = None
        for attempt in range(1, 5):
            temp = path.with_suffix('.xpt.part')
            try:
                with requests.get(url, headers={
                    'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) NHANES-educational-research/1.0',
                    'Accept': '*/*'}, timeout=(20, 90), stream=True) as response:
                    response.raise_for_status()
                    with temp.open('wb') as fh:
                        for chunk in response.iter_content(chunk_size=1024 * 256):
                            if chunk:
                                fh.write(chunk)
                if not verify_binary_xpt(temp):
                    raise ValueError('Response was not a valid SAS XPT (may be an HTML block page).')
                temp.replace(path)
                last_error = None
                break
            except (requests.RequestException, ValueError, OSError) as exc:
                last_error = exc
                temp.unlink(missing_ok=True)
                print(f'  attempt {attempt}/4 failed: {exc}', file=sys.stderr)
                if attempt < 4:
                    time.sleep(attempt * 2)
        if last_error is not None:
            raise RuntimeError(
                f'Could not download {name}. Download manually from {url} '
                f'and save it to {path}, then re-run. Original error: {last_error}')
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return {'file': path.name, 'url': url, 'bytes': path.stat().st_size,
            'sha256': h.hexdigest()}


def load_table(name: str) -> pd.DataFrame:
    path = RAW / f'{name}.xpt'
    df = pd.read_sas(path, format='xport', encoding='utf-8')
    df.columns = [str(c).upper() for c in df.columns]
    df = correct_sas_xport_zeros(df, source=name)
    assert_cdc_lod_flags(df, name)
    if 'SEQN' not in df.columns:
        raise ValueError(f'{name}: no SEQN identifier; parsed columns={list(df.columns)}')
    if df.SEQN.isna().any() or df.SEQN.duplicated().any():
        raise ValueError(f'{name}: duplicate or null SEQN; unsafe to join')
    df['SEQN'] = pd.to_numeric(df['SEQN'], errors='raise').astype('int64')
    return df


def keep(df: pd.DataFrame, cols: list[str], source: str) -> pd.DataFrame:
    missing = sorted(set(cols) - set(df.columns))
    if missing:
        raise ValueError(f'{source}: missing expected variables {missing}; check CDC codebook')
    out = df[cols].copy()
    for col in out.columns.difference(['SEQN']):
        out[col] = pd.to_numeric(out[col], errors='coerce')
    return out


def assemble(tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    demo = keep(tables['DEMO_J'], DEMOGRAPHICS, 'DEMO_J')
    bmx = keep(tables['BMX_J'], ['SEQN', 'BMXBMI'], 'BMX_J')
    crp = keep(tables['HSCRP_J'], ['SEQN', 'LBXHSCRP', 'LBDHRPLC'], 'HSCRP_J')
    cr = keep(tables['ALB_CR_J'], ['SEQN', 'URXUCR'], 'ALB_CR_J')
    pf = keep(tables['PFAS_J'], ['SEQN', 'WTSB2YR', *PFAS, *PFAS_FLAGS], 'PFAS_J')
    ph = keep(tables['PHTHTE_J'], ['SEQN', 'WTSB2YR', *PHTH, *PHTH_FLAGS], 'PHTHTE_J')
    pf = pf.rename(columns={'WTSB2YR': 'WTSB2YR_PFAS'})
    ph = ph.rename(columns={'WTSB2YR': 'WTSB2YR_PHTH'})

    merged = demo
    for sub in [bmx, crp, cr, pf, ph]:
        merged = merged.merge(sub, on='SEQN', how='left', validate='one_to_one')

    adults = merged.loc[(merged['RIDAGEYR'] >= 20) &
                        (merged['LBXHSCRP'] > 0) &
                        (merged['WTSB2YR_PFAS'] > 0)].copy()
    shared_columns = ['SEQN', *DEMOGRAPHICS[1:], 'BMXBMI', 'LBXHSCRP',
                      'URXUCR', 'WTSB2YR_PFAS', 'WTSB2YR_PHTH',
                      *PFAS, *PHTH, *PFAS_FLAGS, *PHTH_FLAGS]
    pf_needed = ['SEQN', *COVARIATES, 'SDMVSTRA', 'SDMVPSU', 'LBXHSCRP',
                 'WTSB2YR_PFAS', *PFAS]
    pf_eligible = adults[pf_needed].dropna().copy()
    pf_eligible = pf_eligible.rename(columns={'WTSB2YR_PFAS': 'survey_weight'})
    joint = adults.loc[(adults['WTSB2YR_PHTH'] > 0) &
                       (adults['URXUCR'] > 0), shared_columns].dropna().copy()

    # PFAS + phthalates are both environmental subsample B; verify not just the name
    # but person-specific weights before claiming one design-weighted joint model.
    w_pf = joint['WTSB2YR_PFAS'].to_numpy(dtype=float)
    w_ph = joint['WTSB2YR_PHTH'].to_numpy(dtype=float)
    weight_rel_diff = np.abs(w_pf - w_ph) / np.maximum(np.abs(w_pf), 1e-12)
    max_relative_difference = float(weight_rel_diff.max()) if len(weight_rel_diff) else None
    joint_weights_ok = bool(len(joint) >= 100 and max_relative_difference is not None
                            and max_relative_difference < 1e-5)
    if joint_weights_ok:
        joint['survey_weight'] = joint['WTSB2YR_PFAS']
        joint = joint.drop(columns=['WTSB2YR_PFAS', 'WTSB2YR_PHTH'])
    else:
        print('[warning] Joint PFAS+phthalates analysis is disabled: too few complete '
              'cases or person-level weights differ. PFAS-only analysis is valid fallback.')

    checks = {
        'n_demo': len(demo),
        'n_pfas_raw': len(pf),
        'n_phthalates_raw': len(ph),
        'n_adults_pfas_crp_pre_complete_case': len(adults),
        'n_pfas_analysis_complete': len(pf_eligible),
        'n_joint_analysis_complete': len(joint),
        'joint_weight_max_relative_difference': max_relative_difference,
        'joint_weights_verified': joint_weights_ok,
        'joint_analysis_exported': joint_weights_ok,
        'pfas_analytical_exposures': PFAS,
        'joint_analytical_exposures': PFAS + PHTH,
        'crp_units': 'mg/L',
        'PFAS_units': 'ng/mL (serum)',
        'PHTH_units': 'ng/mL (urine)',
        'survey_design': 'WTSB2YR (subsample B), SDMVSTRA, SDMVPSU',
        'below_detection': 'CDC already replaced with LOD/sqrt(2); flags retained in SQLite',
        'cross_sectional_data': True,
    }
    return merged, pf_eligible, joint if joint_weights_ok else pd.DataFrame(), checks


def main():
    PROC.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name in FILES:
        manifest.append(download(name))
    print('\n[read] loading SAS XPT files')
    tables = {name: load_table(name) for name in FILES}
    merged, pfas, joint, checks = assemble(tables)
    if len(pfas) < 100:
        raise RuntimeError('Too few complete PFAS cases (<100); inspect input files and joins.')
    # Overwrite our own derived SQLite database each run. Raw XPT cache is never overwritten.
    DB.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(DB) as conn:
        for name, dataframe in tables.items():
            dataframe.to_sql('raw_' + name.lower(), conn, if_exists='replace', index=False)
        merged.to_sql('merged_all', conn, if_exists='replace', index=False)
        pfas.to_sql('analysis_pfas', conn, if_exists='replace', index=False)
        if len(joint):
            joint.to_sql('analysis_joint', conn, if_exists='replace', index=False)
        else:
            conn.execute('DROP TABLE IF EXISTS analysis_joint')
        conn.execute('CREATE INDEX IF NOT EXISTS idx_merged_seqn ON merged_all(SEQN)')
        conn.execute('CREATE INDEX IF NOT EXISTS idx_pfas_seqn ON analysis_pfas(SEQN)')
        if len(joint):
            conn.execute('CREATE INDEX IF NOT EXISTS idx_joint_seqn ON analysis_joint(SEQN)')
    merged.to_csv(PROC / 'merged_all.csv', index=False)
    pfas.to_csv(PROC / 'analysis_pfas.csv', index=False)
    joint_path = PROC / 'analysis_joint.csv'
    if len(joint):
        joint.to_csv(joint_path, index=False)
    else:
        joint_path.unlink(missing_ok=True)
    report = {'created_utc': datetime.now(timezone.utc).isoformat(),
              'source': 'CDC NHANES 2017-2018',
              'source_files': manifest, **checks}
    (PROC / 'qc_report.json').write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(f'\n[done] SQLite: {DB}')
    print(f'[done] PFAS model: {len(pfas)} complete adults')
    if len(joint):
        print(f'[done] Joint model: {len(joint)} complete adults (subsample B weights verified)')
    else:
        print('[done] Joint model not exported; check qc_report.json')
    print(f'[done] QC: {PROC / "qc_report.json"}')
    print('[next] python 02_exposome_eda.py && Rscript 03_survey_mixture.R')


if __name__ == '__main__':
    main()
