#!/usr/bin/env python3
import importlib.util
from pathlib import Path
import pandas as pd
mod_path=Path(__file__).with_name('35_v18_quality_evidence.py')
spec=importlib.util.spec_from_file_location('v18',mod_path)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
raw=pd.DataFrame([
 {'casrn':'355-46-4','aeid':97,'hitc_numeric':0.95,'spid':'x','ac50_uM_if_verified':11,'mc6_flags':'singlept.hit.high','fitc':36},
 {'casrn':'355-46-4','aeid':97,'hitc_numeric':0.96,'spid':'y','ac50_uM_if_verified':10,'mc6_flags':'[]','fitc':30},
 {'casrn':'355-46-4','aeid':97,'hitc_numeric':0.2,'spid':'z','ac50_uM_if_verified':20,'mc6_flags':'[]','fitc':30},
 {'casrn':'375-95-1','aeid':94,'hitc_numeric':.92,'spid':'w','ac50_uM_if_verified':2,'mc6_flags':None,'fitc':20},
])
review=pd.DataFrame([
 {'casrn':'355-46-4','aeid':97,'chemical':'PFHxS','tier':'CORE_12','n_high_positive':2},
 {'casrn':'375-95-1','aeid':94,'chemical':'PFNA','tier':'CORE_12','n_high_positive':1},
])
pos,pairs,audit=m.audit(raw,review)
assert len(pos)==3 and len(pairs)==2
assert audit['n_positive_records_with_flag_indicated']==1
assert audit['n_positive_records_with_flag_info_missing']==1
assert (pos.official_epa_validation=='NOT_PERFORMED__THIRD_PARTY_SOURCE').all()
assert (pos.concentration_range_check=='NOT_PERFORMED').all()
assert m.flag_category('[]')=='NO_FLAG_REPORTED'
assert m.flag_category(None)=='UNKNOWN'
assert m.flag_category('singlept.hit.mid')=='FLAG_INDICATED'
assert pairs.n_with_ac50_uM.sum()==3
print('[PASS] V18: strict screening QC, flag presence != validation, missing != clean, positive-only counts.')
