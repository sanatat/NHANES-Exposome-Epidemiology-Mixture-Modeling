#!/usr/bin/env python3
"""V18: conservative quality metadata audit for positive ToxCast screen hits.

Reads local V15 mirror records and V17 positive review. Does not claim that
missing fit flags mean a reliable curve, and does not rank human risks.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT/'results/v15_measured_bioactivity/v15_matched_raw_assays.csv'
REVIEW = ROOT/'results/v17_positive_hit_review/v17_all44_positive_endpoint_pairs.csv'
OUT = ROOT/'results/v18_curve_quality_audit'
QUALITY_FIELDS = ['fitc', 'flags', 'flag', 'mc6_flags', 'mc6_flag', 'flag_name',
                  'flag_description', 'flag_descriptions', 'n_flags', 'flag_count',
                  'chid_rep', 'model_type', 'modl', 'coff', 'top',
                  'logc_min', 'logc_max', 'logcmin', 'logcmax',
                  'min_conc_um', 'max_conc_um']
FLAG_FIELDS = ['flags', 'flag', 'mc6_flags', 'mc6_flag', 'flag_name',
               'flag_description', 'flag_descriptions', 'n_flags', 'flag_count']
CAUTION_TERMS = ['singlept', 'noise', 'border', 'low.nconc', 'low.nrep',
                 'ac50.lowconc', 'efficacy.50', 'directionality',
                 'viability.gnls', 'bmd.high', 'gnls.lowconc']


def flag_category(value):
    """Returns 'UNKNOWN', 'FLAG_INDICATED', or 'NO_FLAG_REPORTED'.
    Caution: 'NO_FLAG_REPORTED' is *never* synonymous with fit validation.
    """
    if value is None or (not isinstance(value, (str, list, dict, tuple)) and pd.isna(value)):
        return 'UNKNOWN'
    if isinstance(value, (int, float, np.integer, np.floating)):
        return 'FLAG_INDICATED' if float(value)>0 else 'NO_FLAG_REPORTED' if float(value)==0 else 'UNKNOWN'
    s = str(value).strip().lower()
    if s in {'', 'nan', 'none', 'null', 'na', '<na>'}:
        return 'UNKNOWN'
    if s in {'0', '0.0', '[]', '{}', 'false', 'no flags', 'none flagged'}:
        return 'NO_FLAG_REPORTED'
    if s in {'1', 'true'}:
        return 'FLAG_INDICATED'
    if re.fullmatch(r'\d+(\.0+)?',s):
        return 'FLAG_INDICATED' if float(s)>0 else 'NO_FLAG_REPORTED'
    return 'FLAG_INDICATED'


def audit(raw, review):
    r=raw.copy()
    r.columns=[str(c).strip().lower() for c in r.columns]
    req={'casrn','aeid','hitc_numeric'}
    if not req.issubset(r):
        raise ValueError(f'V15 raw missing: {sorted(req-set(r.columns))}')
    p=review.copy()
    if not {'casrn','aeid','chemical','tier'}.issubset(p):
        raise ValueError('V17 review has no exact CAS/AEID/chemical/tier mapping')
    r['casrn']=r.casrn.astype('string').str.strip()
    p['casrn']=p.casrn.astype('string').str.strip()
    r['aeid']=pd.to_numeric(r.aeid,errors='raise').astype(int)
    p['aeid']=pd.to_numeric(p.aeid,errors='raise').astype(int)
    r['hitc_numeric']=pd.to_numeric(r.hitc_numeric,errors='coerce')
    if p.duplicated(['casrn','aeid']).any():
        raise ValueError('V17 duplicate chemical × endpoint pairs')
    if (r['hitc_numeric'].dropna()>1).any() or (r['hitc_numeric'].dropna()<-1).any():
        raise ValueError('Invalid hitc value outside [-1,1]')
    pos=r[r.hitc_numeric>=.9].copy()
    pos=pos.merge(p[['casrn','aeid','chemical','tier']],on=['casrn','aeid'],how='inner',validate='many_to_one',suffixes=('','_v17'))
    if 'chemical_v17' in pos: pos['chemical']=pos.pop('chemical_v17')
    if 'tier_v17' in pos: pos['tier']=pos.pop('tier_v17')
    if len(pos)!=int(pd.to_numeric(p.n_high_positive,errors='raise').sum()):
        raise ValueError('V17 positive counts do not agree with V15 source')
    flag_present=[c for c in FLAG_FIELDS if c in pos.columns]
    quality_present=[c for c in QUALITY_FIELDS if c in pos.columns]
    raw_flag_status=[]
    raw_warning_terms=[]
    for rec in pos.to_dict('records'):
        flags=[flag_category(rec.get(c)) for c in flag_present]
        if 'FLAG_INDICATED' in flags:
            status='REPORTED_FLAG__REVIEW_CURVE'
        elif 'NO_FLAG_REPORTED' in flags:
            status='NO_FLAG_REPORTED__NOT_VALIDATED'
        else:
            status='FLAG_INFO_NOT_AVAILABLE'
        raw_flag_status.append(status)
        sf=' '.join(str(rec.get(c,'')) for c in flag_present).lower()
        raw_warning_terms.append(';'.join(t for t in CAUTION_TERMS if t in sf))
    pos['flag_review_status']=raw_flag_status
    pos['caution_terms_as_text']=raw_warning_terms
    pos['official_epa_validation']='NOT_PERFORMED__THIRD_PARTY_SOURCE'
    pos['concentration_range_check']='NOT_PERFORMED'
    if 'ac50_um_if_verified' in pos.columns:
        pos['ac50_uM_if_verified']=pd.to_numeric(pos['ac50_um_if_verified'],errors='coerce')
    elif 'ac50_uM_if_verified' in pos.columns:
        pos['ac50_uM_if_verified']=pd.to_numeric(pos['ac50_uM_if_verified'],errors='coerce')
    else:
        pos['ac50_uM_if_verified']=np.nan
    if 'chid_rep' in pos:
        cr=pd.to_numeric(pos.chid_rep,errors='coerce')
        pos['representative_sample_status']=np.where(cr.eq(1),'REPORTED_REPRESENTATIVE',np.where(cr.eq(0),'REPORTED_NON_REPRESENTATIVE','UNKNOWN'))
    else:
        pos['representative_sample_status']='UNKNOWN'
    # Only check range when the source columns explicitly say uM; no inferred logc units.
    if {'min_conc_um','max_conc_um'}.issubset(pos.columns):
        lo=pd.to_numeric(pos.min_conc_um,errors='coerce')
        hi=pd.to_numeric(pos.max_conc_um,errors='coerce')
        ac=pos.ac50_uM_if_verified
        ok=lo.notna()&hi.notna()&ac.notna()&(lo>0)&(hi>=lo)
        pos.loc[ok,'concentration_range_check']=np.where((ac[ok]>=lo[ok]) & (ac[ok]<=hi[ok]),'WITHIN_REPORTED_uM_RANGE','OUTSIDE_REPORTED_uM_RANGE__REVIEW')
    cols=['chemical','casrn','aeid','tier']
    all_rows=[]
    for key,g in pos.groupby(cols,dropna=False):
        n=len(g)
        if 'spid' in g.columns: nunique=int(g.spid.dropna().astype(str).nunique())
        else: nunique=None
        all_rows.append({**dict(zip(cols,key)),
           'n_positive_records':n,
           'n_unique_spid_if_available':nunique,
           'n_reported_flags':int(g.flag_review_status.eq('REPORTED_FLAG__REVIEW_CURVE').sum()),
           'n_no_flag_reported_not_validated':int(g.flag_review_status.eq('NO_FLAG_REPORTED__NOT_VALIDATED').sum()),
           'n_flag_info_missing':int(g.flag_review_status.eq('FLAG_INFO_NOT_AVAILABLE').sum()),
           'n_representative_reported':int(g.representative_sample_status.eq('REPORTED_REPRESENTATIVE').sum()),
           'n_representative_unknown':int(g.representative_sample_status.eq('UNKNOWN').sum()),
           'n_with_ac50_uM':int(g.ac50_uM_if_verified.notna().sum()),
           'minimum_screening_ac50_uM':float(g.ac50_uM_if_verified.min()) if g.ac50_uM_if_verified.notna().any() else np.nan,
           'n_outside_explicit_uM_range':int(g.concentration_range_check.eq('OUTSIDE_REPORTED_uM_RANGE__REVIEW').sum()),
           'n_range_not_checked':int(g.concentration_range_check.eq('NOT_PERFORMED').sum()),
           'quality_conclusion':'SCREENING_ONLY__PRIMARY_EPA_QC_PENDING'})
    pair=pd.DataFrame(all_rows)
    if pair.empty:
        pair=pd.DataFrame(columns=cols+['n_positive_records','quality_conclusion'])
    pair=pair.sort_values(['tier','chemical','aeid']).reset_index(drop=True)
    audit={
     'source':'Local V15 third-party ToxCast mirror matched raw records; exact CAS x AEID',
     'n_positive_records':len(pos),'n_candidate_pairs':len(pair),
     'quality_fields_found':quality_present,'flag_fields_found':flag_present,
     'quality_fields_not_found':[x for x in QUALITY_FIELDS if x not in quality_present],
     'n_positive_records_with_flag_indicated':int(pos.flag_review_status.eq('REPORTED_FLAG__REVIEW_CURVE').sum()),
     'n_positive_records_with_flag_info_missing':int(pos.flag_review_status.eq('FLAG_INFO_NOT_AVAILABLE').sum()),
     'official_epa_crosscheck':'NOT_PERFORMED',
     'caveats':['Flag absence is not proof of validity', 'Unverified concentration unit or range is not checked',
                'Sample identifiers are not independent replication', 'Hitc>=0.9 is a screening call not human health risk',
                'Full EPA mc6 QC validation still required']}
    return pos,pair,audit


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--raw',type=Path,default=RAW)
    ap.add_argument('--review',type=Path,default=REVIEW)
    ap.add_argument('--out',type=Path,default=OUT)
    a=ap.parse_args()
    for x in (a.raw,a.review):
        if not x.is_file(): raise SystemExit(f'Missing input file: {x}')
    pos,pairs,audit=audit_fn(pd.read_csv(a.raw,low_memory=False),pd.read_csv(a.review,low_memory=False))
    a.out.mkdir(parents=True,exist_ok=True)
    pos.to_csv(a.out/'v18_positive_record_quality_fields.csv',index=False)
    pairs.to_csv(a.out/'v18_all44_hit_review.csv',index=False)
    pairs[pairs.tier.eq('CORE_12')].to_csv(a.out/'v18_core12_hit_review.csv',index=False)
    (a.out/'v18_qc_schema_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print('[V18] positive records:',audit['n_positive_records'],'candidate pairs:',audit['n_candidate_pairs'])
    print('[V18] QC fields available:',', '.join(audit['quality_fields_found']) or '(none)')
    print('[V18] flags indicated:',audit['n_positive_records_with_flag_indicated'],
          '| missing flag info:',audit['n_positive_records_with_flag_info_missing'])
    print('[V18] CORE_12 quality review:')
    cols=['chemical','aeid','n_positive_records','n_reported_flags','n_flag_info_missing','n_with_ac50_uM','n_range_not_checked']
    print(pairs[pairs.tier.eq('CORE_12')][cols].to_string(index=False))
    print('[V18] All are screening candidates; no verified toxicological risk ranking.')
    print('[V18] Saved:',a.out)

# Helper alias permits straightforward unit testing.
audit_fn=audit
if __name__=='__main__':main()
