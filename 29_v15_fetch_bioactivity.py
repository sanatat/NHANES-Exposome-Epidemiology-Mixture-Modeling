#!/usr/bin/env python3
"""V15 experimental: retrieve exact chemical x AEID bioactivity from a THIRD-PARTY
mirror of US EPA invitrodb v4.3, using scigantic-comptox (DuckDB/S3).

Data are NOT retrieved from an official EPA endpoint and should be verified
against official v4.3 exports before scientific publication.

Status is NOT_RETRIEVED until live results are available. Missing rows are
NOT_RECORDED_IN_MIRROR, never interpreted as inactive/negative.

US EPA's tcpl documents signed hitc values in [-1,1], where sign can encode
response orientation. Store raw signed hitc, never silently take absolute
value or call a negative sign an inactive assay.
"""
import argparse
import json
import re
import sys
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PANEL = ROOT / 'results/v14_mechanistic_panel/v14_chemical_x_endpoint_evidence.csv'
OUT = ROOT / 'results/v15_measured_bioactivity'
REQUIRED = {'chemical','casrn','aeid','tier'}
CAS_PATTERN = re.compile(r'^\d{2,7}-\d{2}-\d$')


def validate_pairs(frame):
    absent = REQUIRED - set(frame.columns)
    if absent:
        raise ValueError('V14 input missing: ' + repr(sorted(absent)))
    pairs = frame[['chemical','casrn','aeid','tier']].copy()
    pairs['casrn'] = pairs['casrn'].astype(str).str.strip()
    if not pairs['casrn'].map(lambda c: bool(CAS_PATTERN.fullmatch(c))).all():
        raise ValueError('CASRN invalid; refusing fuzzy chemical identity matching')
    a = pd.to_numeric(pairs['aeid'], errors='coerce')
    if a.isna().any() or ((a%1)!=0).any():
        raise ValueError('AEID invalid')
    pairs['aeid'] = a.astype(int)
    if pairs.duplicated(['casrn','aeid']).any():
        raise ValueError('Duplicate CASRN + AEID in source')
    return pairs


def classify_hit(value, orientation_supported=True):
    if pd.isna(value):
        return 'MISSING_HITC'
    h = float(value)
    if not np.isfinite(h) or h < -1 or h > 1:
        raise ValueError(f'Invalid continuous ToxCast hitc: {h}')
    if h >= 0.9:
        return 'HIGH_HITCALL_POSITIVE_DIRECTION'
    if h <= -0.9:
        return 'HIGH_HITCALL_OPPOSITE_DIRECTION' if orientation_supported else 'HIGH_MAGNITUDE_REVIEW'
    return 'BELOW_0_9_MAGNITUDE'


def normalize_real_observations(df, pairs):
    x = df.copy()
    x.columns = [str(c).strip().lower() for c in x.columns]
    ccol = 'casn' if 'casn' in x else 'casrn' if 'casrn' in x else None
    hcol = 'hitc' if 'hitc' in x else 'hitcall' if 'hitcall' in x else None
    if not ccol or 'aeid' not in x or not hcol:
        raise ValueError('Need exact CASRN, AEID, hitc/hitcall. Got '+ repr(x.columns.tolist()))
    x['casrn'] = x[ccol].astype('string').str.strip()
    x['aeid'] = pd.to_numeric(x.aeid, errors='coerce')
    x['hitc_numeric'] = pd.to_numeric(x[hcol], errors='coerce')
    allowed = pairs[['casrn','aeid']]
    x = x.merge(allowed, on=['casrn','aeid'], how='inner', validate='many_to_one')
    if x.empty:
        return x
    # Only signed hitc has documented bidirectional interpretation.
    x['assay_result'] = x.hitc_numeric.map(lambda z: classify_hit(z, orientation_supported=(hcol=='hitc')))
    if 'conc_unit' in x.columns:
        unit = x.conc_unit.astype('string').str.strip().str.lower()
    else:
        unit = pd.Series(pd.NA, index=x.index, dtype='string')
    if 'ac50' in x.columns:
        x['ac50_original'] = pd.to_numeric(x.ac50,errors='coerce')
    else:
        x['ac50_original'] = np.nan
    x['ac50_uM_if_verified'] = x['ac50_original'].where(unit.eq('um'))
    x['ac50_unit_verified'] = unit.eq('um').fillna(False)
    return x


def build_matrix(pairs, raw):
    out = pairs.copy()
    if raw.empty:
        out['retrieval_status'] = 'NOT_RECORDED_IN_SOURCE'
        out['n_records'] = 0
        for c in ['n_valid_hitc','n_high_positive','n_high_opposite','n_below_threshold']:
            out[c] = 0
        out['min_ac50_uM_high_hit'] = np.nan
        return out
    records=[]
    for (casrn, aeid), group in raw.groupby(['casrn','aeid']):
        status=group.assay_result
        pos=status.eq('HIGH_HITCALL_POSITIVE_DIRECTION')
        neg=status.eq('HIGH_HITCALL_OPPOSITE_DIRECTION')
        review=status.eq('HIGH_MAGNITUDE_REVIEW')
        hi=pos | neg | review
        usable=group.loc[hi,'ac50_uM_if_verified'].dropna()
        records.append({
            'casrn':casrn,'aeid':int(aeid),
            'retrieval_status':'OBSERVED_IN_SOURCE',
            'n_records':len(group),
            'n_valid_hitc':int(group.hitc_numeric.notna().sum()),
            'n_high_positive':int(pos.sum()),
            'n_high_opposite':int(neg.sum()),
            'n_high_magnitude_orientation_review':int(review.sum()),
            'n_below_threshold':int(status.eq('BELOW_0_9_MAGNITUDE').sum()),
            'n_missing_hitc':int(status.eq('MISSING_HITC').sum()),
            'min_ac50_uM_high_hit':float(usable.min()) if len(usable) else np.nan,
        })
    summary=pd.DataFrame(records)
    out=out.merge(summary, on=['casrn','aeid'], how='left', validate='one_to_one')
    out['retrieval_status']=out.retrieval_status.fillna('NOT_RECORDED_IN_SOURCE')
    for c in ['n_records','n_valid_hitc','n_high_positive','n_high_opposite','n_high_magnitude_orientation_review','n_below_threshold','n_missing_hitc']:
        out[c]=out[c].fillna(0).astype(int)
    return out


def fetch_from_mirror(pairs):
    try:
        import scigantic_comptox as comptox
    except ImportError as e:
        raise SystemExit('Install optional dependency: .venv/bin/python -m pip install scigantic-comptox') from e
    # Inspect declared view columns rather than guessing the mirror schema.
    schema = comptox.query('SELECT * FROM bioactivity LIMIT 0')
    cols = {str(c).lower() for c in schema.columns}
    if 'casn' not in cols and 'casrn' not in cols:
        raise RuntimeError('Mirror bioactivity view has no CASRN column; exact mapping not available. Columns: '+repr(sorted(cols)))
    if 'aeid' not in cols or not ({'hitc','hitcall'} & cols):
        raise RuntimeError('Mirror view missing AEID or hitcall. Columns: '+repr(sorted(cols)))
    ccol = 'casn' if 'casn' in cols else 'casrn'
    caslist=','.join("'"+c+"'" for c in sorted(pairs.casrn.unique()))
    aeids=','.join(map(str,sorted(pairs.aeid.unique())))
    sql=f'SELECT * FROM bioactivity WHERE {ccol} IN ({caslist}) AND aeid IN ({aeids})'
    return comptox.query(sql), {'source':'THIRD_PARTY_MIRROR_NOT_EPA_ENDPOINT','mirror_package':'scigantic-comptox','query_chemical_count':pairs.casrn.nunique(),'query_aeid_count':pairs.aeid.nunique()}


def run(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-csv', type=Path, help='Optional existing exact-CAS EPA-style CSV, instead of third-party mirror')
    ap.add_argument('--pairs', type=Path, default=PANEL)
    ap.add_argument('--out', type=Path, default=OUT)
    args=ap.parse_args(argv)
    if not args.pairs.is_file():
        raise SystemExit('First run V14; missing '+str(args.pairs))
    pairs = validate_pairs(pd.read_csv(args.pairs,low_memory=False))
    # No output materialized until retrieval is complete; failed live queries never create fake negatives.
    if args.input_csv:
        observations=pd.read_csv(args.input_csv,low_memory=False)
        source={'source':'USER_PROVIDED_LOCAL_CSV','file':str(args.input_csv)}
    else:
        observations,source=fetch_from_mirror(pairs)
    raw=normalize_real_observations(observations,pairs)
    if raw.empty:
        print('[WARNING] No matching exact-CAS/AEID observations in supplied source. This does NOT mean tested negative.')
    result=build_matrix(pairs,raw)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'v15_chemical_x_aeid_matrix.csv').write_text(result.to_csv(index=False))
    (args.out/'v15_matched_raw_assays.csv').write_text(raw.to_csv(index=False))
    core=result[result.tier=='CORE_12'].copy()
    core.to_csv(args.out/'v15_core12_results.csv',index=False)
    audit={**source,'generated_utc':datetime.now(timezone.utc).isoformat(),
           'input_chemicals':int(pairs.casrn.nunique()),'input_endpoints':int(pairs.aeid.nunique()),
           'n_pairs_expected':len(pairs),'n_pairs_observed':int(result.retrieval_status.eq('OBSERVED_IN_SOURCE').sum()),
           'n_pair_absent_from_source':int(result.retrieval_status.eq('NOT_RECORDED_IN_SOURCE').sum()),
           'n_records_exact_chem_endpoint':len(raw),
           'note':'No observed assay record is NOT evidence of a negative result. In vitro data are not human risk scores. Signed hitc is preserved.'}
    (args.out/'v15_audit.json').write_text(json.dumps(audit,indent=2))
    print('\n[V15] Evidence matrix:',len(result),'chemical x endpoint possibilities')
    print('[V15] Observed pairs:',audit['n_pairs_observed'],'; missing from queried source:',audit['n_pair_absent_from_source'])
    print('[V15] Core 12 results:')
    print(core[['chemical','aeid','retrieval_status','n_records','n_high_positive','n_high_opposite']].to_string(index=False))
    print('[V15] saved',args.out)
    return audit

if __name__ == '__main__':
    run()
