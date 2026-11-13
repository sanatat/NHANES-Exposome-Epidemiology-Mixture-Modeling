#!/usr/bin/env python3
"""V20: decode EPA invitrodb numeric mc6 flag *method IDs*.

Inputs: V18 local positive record export (mirrored ToxCast data).
Outputs: record-level, per-ID, core-12 and NRF2 subsets, with provenance/QC.
This is an annotation audit; NEVER equate no/unknown flags with validated assays.

EPA tcpl source: https://usepa.r-universe.dev/tcpl/doc/Introduction_Appendices.html
Sections: Level 6 flags, fit categories, LOEC method.
"""
from __future__ import annotations
import argparse
import json
import re
from collections import Counter
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / 'results/v18_curve_quality_audit/v18_positive_record_quality_fields.csv'
DEFAULT_OUT = ROOT / 'results/v20_numeric_flags'

# IDs are EPA Level 6 method identifiers, NOT bit fields or quality scores.
# See the EPA tcpl vignette table in "Flags (Level 6)".
FLAG_DESCRIPTIONS = {
    5: 'Model directionality questionable',
    6: 'Active only at the highest concentration (single-point)',
    7: 'Active at just one non-highest concentration (single-point)',
    8: 'Inactive despite multiple concentrations above baseline',
    9: 'BMD greater than AC50; high baseline variability',
    10: 'Noisy concentration-response data',
    11: 'Borderline response near cutoff',
    13: 'Average replicates per concentration below 2',
    14: 'Four or fewer concentrations tested',
    15: 'Gain-loss model AC50 range concern',
    17: 'Low efficacy (less than 50%)',
    18: 'AC50 below lowest concentration tested',
    19: 'Cell viability assay with gain-loss winning model',
    20: 'No median response beyond 3x baseline MAD',
}
EMPTY = {'', 'none', 'null', 'nan', '<na>', 'na'}


def parse_flag_ids(raw):
    if raw is None or pd.isna(raw) or str(raw).strip().lower() in EMPTY:
        return [], 'FLAG_INFO_MISSING', ''
    s = str(raw).strip()
    if s in {'0', '0.0', '[]'}:
        return [], 'NO_FLAGS_REPORTED_NOT_VALIDATED', ''
    # Accept explicit comma / semicolon / pipe-separated integer IDs only.
    if not re.fullmatch(r'\s*\d+(?:\.0+)?(?:\s*[,;|]\s*\d+(?:\.0+)?)*\s*', s):
        return [], 'NON_NUMERIC_FORMAT_UNRESOLVED', s
    ids = sorted(set(int(float(tok)) for tok in re.split(r'[,;|]', s)))
    if 0 in ids:
        return [], 'UNEXPECTED_ZERO_WITH_OTHER_IDS', s
    unknown = [x for x in ids if x not in FLAG_DESCRIPTIONS]
    return ids, ('UNKNOWN_FLAG_ID_PRESENT' if unknown else 'NUMERIC_IDS_DECODED'), ','.join(map(str, unknown))


def classify_fit(row):
    cat = pd.to_numeric(row.get('fitc'), errors='coerce')
    modl = str(row.get('modl', '')).strip().lower()
    mtype = pd.to_numeric(row.get('model_type'), errors='coerce')
    if cat == 100 or modl == 'loec' or mtype == 1:
        return 'LOEC_METHOD__NOT_CONVENTIONAL_CURVE_AC50'
    if cat == 41:
        return 'FITC41__AC50_AND_AC95_WITHIN_SCREENED_RANGE__QC_PENDING'
    if pd.isna(cat):
        return 'FIT_CATEGORY_MISSING'
    return 'OTHER_FIT_CATEGORY__REVIEW'


def analyze(df):
    req = {'chemical', 'casrn', 'aeid', 'tier', 'mc6_flags', 'fitc', 'modl', 'model_type'}
    if not req.issubset(df.columns):
        raise ValueError('Missing required V18 fields: ' + repr(sorted(req - set(df.columns))))
    x = df.copy().reset_index(drop=True)
    results = [parse_flag_ids(v) for v in x['mc6_flags']]
    x['mc6_flag_ids_normalized'] = [','.join(map(str, ids)) for ids, _, _ in results]
    x['mc6_flag_decode_status'] = [status for _, status, _ in results]
    x['mc6_unknown_flag_ids_or_text'] = [unknown for _, _, unknown in results]
    x['mc6_flag_definitions'] = [
        '; '.join(f'{code}: {FLAG_DESCRIPTIONS[code]}' for code in ids if code in FLAG_DESCRIPTIONS)
        for ids, _, _ in results
    ]
    x['n_mc6_flag_ids'] = [len(ids) for ids, _, _ in results]
    x['fit_method_interpretation'] = x.apply(classify_fit, axis=1)
    x['potency_interpretation'] = x['fit_method_interpretation'].map(
        lambda s: 'LOEC_CONTEXT__DO_NOT_LABEL_AS_CURVE_AC50'
        if s.startswith('LOEC_METHOD') else
        'AC50_POSSIBLY_INFORMATIVE__CHECK_FLAGS_AND_SOURCE'
        if s.startswith('FITC41') else 'POTENCY_METHOD_REVIEW_REQUIRED'
    )
    x['epa_primary_source_crosscheck'] = 'NOT_PERFORMED__SECONDARY_MIRROR'
    x['tested_concentration_range_check'] = 'NOT_PERFORMED'

    counts = Counter(flag_id for ids, _, _ in results for flag_id in ids)
    freqs = pd.DataFrame([
        {'flag_id': flag_id, 'epa_description': FLAG_DESCRIPTIONS.get(flag_id, 'UNKNOWN_ID_REVIEW_OFFICIAL_METHOD_TABLE'),
         'n_positive_records_with_flag': count}
        for flag_id, count in sorted(counts.items(), key=lambda v: (-v[1], v[0]))
    ], columns=['flag_id', 'epa_description', 'n_positive_records_with_flag'])

    keys = ['chemical', 'casrn', 'aeid', 'tier']
    per_pair = x.groupby(keys, dropna=False).agg(
        n_positive_records=('aeid', 'size'),
        n_decoded=('mc6_flag_decode_status', lambda s: int(s.eq('NUMERIC_IDS_DECODED').sum())),
        n_missing_flag_info=('mc6_flag_decode_status', lambda s: int(s.eq('FLAG_INFO_MISSING').sum())),
        n_unknown_flag_values=('mc6_flag_decode_status', lambda s: int((~s.isin(['NUMERIC_IDS_DECODED', 'FLAG_INFO_MISSING', 'NO_FLAGS_REPORTED_NOT_VALIDATED'])).sum())),
        n_loec_fits=('fit_method_interpretation', lambda s: int(s.str.startswith('LOEC_METHOD').sum())),
        n_fitc41=('fit_method_interpretation', lambda s: int(s.str.startswith('FITC41').sum())),
    ).reset_index().sort_values(['tier', 'chemical', 'aeid'])

    nrf = x[x.chemical.isin(['PFNA', 'PFOS']) & pd.to_numeric(x.aeid, errors='coerce').eq(97)].copy()
    overview = {
        'n_positive_records': int(len(x)),
        'n_chemical_endpoint_pairs': int(len(per_pair)),
        'n_numeric_decoded': int(x.mc6_flag_decode_status.eq('NUMERIC_IDS_DECODED').sum()),
        'n_missing_flag_info': int(x.mc6_flag_decode_status.eq('FLAG_INFO_MISSING').sum()),
        'n_other_format_or_unknown': int((~x.mc6_flag_decode_status.isin(['NUMERIC_IDS_DECODED', 'FLAG_INFO_MISSING', 'NO_FLAGS_REPORTED_NOT_VALIDATED'])).sum()),
        'n_loec_method': int(x.fit_method_interpretation.str.startswith('LOEC_METHOD').sum()),
        'n_fitc41': int(x.fit_method_interpretation.str.startswith('FITC41').sum()),
        'n_nrf2_pfna_pfos_records': int(len(nrf)),
        'source': 'V18 locally curated output of a secondary mirror of EPA invitrodb',
        'official_epa_curve_validation': 'NOT_PERFORMED',
        'caveats': [
            'Flag presence requires inspection, not automatic rejection',
            'Flag absence/missing is not curve validation',
            'LOEC method (fitc 100) must not be reported as a conventional curve-fit AC50',
            'A single sample ID is not independent biological replication',
            'Concentration-response range and original EPA results remain unverified',
        ],
        'reference': 'https://usepa.r-universe.dev/tcpl/doc/Introduction_Appendices.html',
    }
    return x, freqs, per_pair, nrf, overview


def main():
    ap = argparse.ArgumentParser(description='Decode numeric EPA tcpl mc6 caution flag IDs')
    ap.add_argument('--input', type=Path, default=DEFAULT_IN)
    ap.add_argument('--out', type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    if not args.input.is_file():
        raise SystemExit('Missing V18 input: '+str(args.input))
    x, freqs, pairs, nrf, overview = analyze(pd.read_csv(args.input, dtype={'mc6_flags': 'string'}, low_memory=False))
    args.out.mkdir(parents=True, exist_ok=True)
    x.to_csv(args.out / 'v20_positive_records_decoded.csv', index=False)
    freqs.to_csv(args.out / 'v20_flag_frequency.csv', index=False)
    pairs.to_csv(args.out / 'v20_positive_pairs_qc.csv', index=False)
    nrf.to_csv(args.out / 'v20_PFNA_PFOS_NRF2_focus.csv', index=False)
    (args.out / 'v20_audit.json').write_text(json.dumps(overview, indent=2), encoding='utf-8')
    print('[V20] positive records:', overview['n_positive_records'], 'pairs:', overview['n_chemical_endpoint_pairs'])
    print('[V20] decoded:', overview['n_numeric_decoded'], 'missing:', overview['n_missing_flag_info'],
          'unknown/unsupported:', overview['n_other_format_or_unknown'])
    print('[V20] LOEC method:', overview['n_loec_method'], '; fitc 41:', overview['n_fitc41'])
    print('[V20] NUMERIC FLAG FREQUENCY:')
    print(freqs.to_string(index=False) if len(freqs) else '(no flag IDs)')
    print('[V20] PFNA / PFOS NRF2 records:')
    show = ['chemical', 'spid', 'hitc_numeric', 'fitc', 'mc6_flags', 'mc6_flag_definitions']
    print(nrf[[c for c in show if c in nrf]].to_string(index=False))
    print('[V20] All results are screening candidates; official EPA confirmation is pending.')
    print('[V20] saved:', args.out)

if __name__ == '__main__':
    main()
