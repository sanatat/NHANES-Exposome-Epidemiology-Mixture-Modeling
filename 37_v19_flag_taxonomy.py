#!/usr/bin/env python3
"""V19: transparently audit mc6 caution flag *types* in V18 positive-screening rows.

Input:  results/v18_curve_quality_audit/v18_positive_record_quality_fields.csv
Output: results/v19_flag_taxonomy/*.csv, *.json

This is exploratory. It never labels a record 'verified', 'reliable', or 'safe'.
A flag indicates a review concern, NOT automatic invalidity. The mirror
may not include all flags and official EPA review is still necessary.
"""
from __future__ import annotations
import argparse
import json
import re
from collections import Counter
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_IN = HERE / 'results/v18_curve_quality_audit/v18_positive_record_quality_fields.csv'
DEFAULT_OUT = HERE / 'results/v19_flag_taxonomy'

# Flag vocabulary from EPA tcpl. Text may contain multiple flags or lists.
KNOWN = [
    'modl.directionality.fail', 'low.nrep', 'low.nconc', 'bmd.high',
    'singlept.hit.high', 'singlept.hit.mid', 'multipoint.neg',
    'gnls.lowconc', 'noise', 'border', 'efficacy.50',
    'ac50.lowconc', 'viability.gnls', 'overfit.hit', 'no.med.gt.3bmad',
]
SERIOUS_REVIEW = {
    'modl.directionality.fail', 'singlept.hit.high', 'singlept.hit.mid',
    'noise', 'border', 'ac50.lowconc', 'gnls.lowconc', 'viability.gnls',
}
EMPTY_FLAG_VALUES = {'', 'none', 'null', 'nan', '<na>', '0', '0.0', '[]', '{}', 'false', 'no flags', 'none flagged'}
FLAG_SOURCES = ['mc6_flags', 'flags', 'flag', 'mc6_flag', 'flag_name',
                'flag_description', 'flag_descriptions', 'n_flags', 'flag_count']


def has_flag_value(val):
    if val is None or pd.isna(val) if not isinstance(val, (list, dict)) else False:
        return False
    return str(val).strip().lower() not in EMPTY_FLAG_VALUES


def scan_row(row):
    fields = [c for c in FLAG_SOURCES if c in row.index]
    text = ' ; '.join(str(row[c]) for c in fields if has_flag_value(row[c]))
    tags = sorted({flag for flag in KNOWN if re.search(r'(?<![a-z0-9_.])' + re.escape(flag) + r'(?![a-z0-9_.])', text, re.I)})
    # In the V18 output this status is a pre-existing indication, not a validated assessment.
    status = str(row.get('flag_review_status', '')).strip()
    if status == 'REPORTED_FLAG__REVIEW_CURVE':
        overall = 'KNOWN_CAUTION_TEXT' if tags else 'FLAG_INDICATED_TYPE_UNRESOLVED'
    elif status == 'NO_FLAG_REPORTED__NOT_VALIDATED':
        overall = 'NO_FLAG_REPORTED_NOT_VALIDATED'
    elif status == 'FLAG_INFO_NOT_AVAILABLE':
        overall = 'FLAG_INFO_MISSING'
    else:
        overall = 'UNKNOWN_STATUS_REVIEW'
    # A positive record should never be stamped verified without original EPA review.
    return pd.Series({
        'flag_types_detected': '; '.join(tags),
        'n_known_flag_types': len(tags),
        'review_concern_flag_types': '; '.join(t for t in tags if t in SERIOUS_REVIEW),
        'flag_taxonomy_status': overall,
        'flag_text_raw': text[:2000],
        'official_confirmation': 'NOT_PERFORMED',
    })


def analyze(pos):
    needed={'chemical','casrn','aeid','tier','flag_review_status'}
    if not needed.issubset(pos):
        raise ValueError('Missing V18 required columns: '+repr(sorted(needed-set(pos.columns))))
    x = pd.concat([pos.reset_index(drop=True), pos.apply(scan_row,axis=1).reset_index(drop=True)],axis=1)
    n = len(x)
    assert n == len(pos)
    cnt=Counter()
    for tokens in x.flag_types_detected:
        cnt.update([v.strip() for v in tokens.split(';') if v.strip()])
    stats = pd.DataFrame([{'flag':name,'n_positive_records':cnt.get(name,0)} for name in KNOWN])
    stats=stats.sort_values(['n_positive_records','flag'], ascending=[False,True]).reset_index(drop=True)
    grp=x.groupby(['chemical','aeid','tier'],dropna=False).agg(
        n_positive_records=('aeid','size'),
        n_with_known_flag_type=('flag_taxonomy_status',lambda a: int((a=='KNOWN_CAUTION_TEXT').sum())),
        n_flag_type_unresolved=('flag_taxonomy_status',lambda a: int((a=='FLAG_INDICATED_TYPE_UNRESOLVED').sum())),
        n_missing_flag_info=('flag_taxonomy_status',lambda a: int((a=='FLAG_INFO_MISSING').sum())),
        n_no_flag_reported=('flag_taxonomy_status',lambda a: int((a=='NO_FLAG_REPORTED_NOT_VALIDATED').sum())),
    ).reset_index().sort_values(['tier','chemical','aeid'])
    audit = {
        'positive_records':n,
        'pairs':int(len(grp)),
        'status_counts':{str(k):int(v) for k,v in x.flag_taxonomy_status.value_counts().items()},
        'records_with_any_known_review_concern':int(x.review_concern_flag_types.ne('').sum()),
        'record_count_warning':'Counts reflect recovered screening rows, not independent biological replicates',
        'source':'Existing local V18 results from a non-EPA mirror of invitrodb',
        'official_epa_crosscheck':'NOT_PERFORMED',
        'concentration_range_validation':'NOT_PERFORMED',
        'interpretation':'Flags indicate caution, not proof a hit is invalid; missing flags do not mean a reliable curve.',
    }
    return x, stats, grp, audit


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,default=DEFAULT_IN)
    ap.add_argument('--output',type=Path,default=DEFAULT_OUT)
    args=ap.parse_args()
    if not args.input.is_file():
        raise SystemExit('Missing V18 CSV: '+str(args.input))
    pos = pd.read_csv(args.input, low_memory=False)
    x,stats,grp,audit=analyze(pos)
    args.output.mkdir(parents=True,exist_ok=True)
    x.to_csv(args.output/'v19_positive_records_flag_taxonomy.csv',index=False)
    stats.to_csv(args.output/'v19_flag_frequency.csv',index=False)
    grp.to_csv(args.output/'v19_positive_pairs_review.csv',index=False)
    (args.output/'v19_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print('[V19] positive records',audit['positive_records'],'; chemical-endpoint pairs',audit['pairs'])
    print('[V19] Flag taxonomy status:')
    for name,num in audit['status_counts'].items(): print(' ',name,':',num)
    print('[V19] Recognized EPA caution-flag terms:')
    print(stats[stats.n_positive_records>0].to_string(index=False) if stats.n_positive_records.sum() else '(none recognized — inspect raw flag_text_raw)')
    print('[V19] Core-12 affected pairs:')
    print(grp[grp.tier.eq('CORE_12')].to_string(index=False))
    print('[V19] No official EPA verification / dose-range check performed.')
    print('[V19] Saved:',args.output)

if __name__=='__main__': main()
