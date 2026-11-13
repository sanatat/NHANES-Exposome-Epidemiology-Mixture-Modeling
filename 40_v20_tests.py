#!/usr/bin/env python3
"""Self-contained synthetic tests; no download or existing NHANES data required."""
import importlib.util
from pathlib import Path
import pandas as pd

p = Path(__file__).with_name('39_v20_decode_numeric_flags.py')
spec = importlib.util.spec_from_file_location('decode_v20', p)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
assert m.parse_flag_ids('13,14,17')[0] == [13,14,17]
assert m.parse_flag_ids('6,13')[1] == 'NUMERIC_IDS_DECODED'
assert m.parse_flag_ids(None)[1] == 'FLAG_INFO_MISSING'
assert m.parse_flag_ids('bogus')[1] == 'NON_NUMERIC_FORMAT_UNRESOLVED'
assert m.parse_flag_ids('999')[1] == 'UNKNOWN_FLAG_ID_PRESENT'
base = {'casrn':'375-95-1','chemical':'PFNA','tier':'CORE_12','aeid':97}
rows = [
 {**base, 'spid':'a','mc6_flags':'6,13', 'fitc':41,'model_type':2,'modl':'poly2'},
 {**base, 'spid':'b','mc6_flags':pd.NA, 'fitc':41,'model_type':2,'modl':'hill'},
 {**base, 'spid':'c','mc6_flags':'13,14,17', 'fitc':100,'model_type':1,'modl':'loec'},
]
x, freq, pairs, nrf, overview = m.analyze(pd.DataFrame(rows))
assert overview['n_positive_records'] == 3
assert overview['n_numeric_decoded'] == 2
assert overview['n_missing_flag_info'] == 1
assert overview['n_loec_method'] == 1
assert overview['n_fitc41'] == 2
assert len(nrf) == 3
assert 13 in set(freq.flag_id)
assert x.loc[2,'potency_interpretation'] == 'LOEC_CONTEXT__DO_NOT_LABEL_AS_CURVE_AC50'
assert len(pairs) == 1
print('[PASS] V20 numeric flag decoding, unknown/missing, fitc 41 vs LOEC 100, and paired summaries')
