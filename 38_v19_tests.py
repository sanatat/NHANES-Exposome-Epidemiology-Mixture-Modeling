#!/usr/bin/env python3
"""Synthetic tests only; never a claim about real ToxCast screening evidence."""
import importlib.util
from pathlib import Path
import pandas as pd
spec=importlib.util.spec_from_file_location('v19',Path(__file__).with_name('37_v19_flag_taxonomy.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
x=pd.DataFrame([
    dict(chemical='PFNA',casrn='375-95-1',aeid=97,tier='CORE_12',flag_review_status='REPORTED_FLAG__REVIEW_CURVE',mc6_flags='["singlept.hit.high", "noise"]'),
    dict(chemical='PFNA',casrn='375-95-1',aeid=97,tier='CORE_12',flag_review_status='REPORTED_FLAG__REVIEW_CURVE',mc6_flags='unknownFlag999'),
    dict(chemical='PFOS',casrn='1763-23-1',aeid=94,tier='CORE_12',flag_review_status='FLAG_INFO_NOT_AVAILABLE',mc6_flags=None),
    dict(chemical='PFOS',casrn='1763-23-1',aeid=97,tier='CORE_12',flag_review_status='NO_FLAG_REPORTED__NOT_VALIDATED',mc6_flags='[]'),
])
a,b,c,d=m.analyze(x)
assert len(a)==4 and len(c)==3
assert b.set_index('flag').loc['noise','n_positive_records']==1
assert b.set_index('flag').loc['singlept.hit.high','n_positive_records']==1
assert d['status_counts']['FLAG_INDICATED_TYPE_UNRESOLVED']==1
assert d['status_counts']['FLAG_INFO_MISSING']==1
assert d['status_counts']['NO_FLAG_REPORTED_NOT_VALIDATED']==1
assert (a.official_confirmation=='NOT_PERFORMED').all()
print('[PASS] V19 known/unknown flag taxonomy, missing vs unflagged, no validation claims')
