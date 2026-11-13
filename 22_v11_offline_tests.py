#!/usr/bin/env python3
"""Offline synthetic tests: never present synthetic assay metrics as real EPA data."""
import importlib.util, tempfile
from pathlib import Path
import numpy as np
import pandas as pd

p=Path(__file__).resolve().parent/'20_v11_epa_toxcast_bridge.py'
s=importlib.util.spec_from_file_location('v11',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
assert len(m.CHEMICALS)==7
assert {x['nhanes_column'] for x in m.CHEMICALS}==set(m.FLAG_MAP)
assert {x['casrn'] for x in m.PARENTS}.isdisjoint({x['casrn'] for x in m.CHEMICALS})
assert m.cas_clean('131-70-4')=='131-70-4'
assert m.cas_clean('unknown')==''

with tempfile.TemporaryDirectory() as tmp:
    path=Path(tmp)/'synthetic_epa.xlsx'
    sheet=pd.DataFrame({'casn':['131-70-4','2306-33-4','84-74-2','355-46-4'],
        'chnm':['MBP synthetic','MEP synthetic','DBP synthetic','PFHxS synthetic'],
        'cytotox_median_um':[1000.,5.,3.,1000.], 'ntested':[50,60,70,20],'nhit':[0,12,42,0]})
    sheet.to_excel(path,index=False)
    got,meta=m.load_cytotox(path)
    pr=[{'nhanes_column':c['nhanes_column'],'nhanes_n':2568,'nhanes_weighted_gmean_ng_mL':1.0,'nhanes_weighted_pct_below_lod':1.0} for c in m.CHEMICALS]
    sub,parents=m.do_cytotox(got,pr)
    mbp=sub[sub.short=='MBP'].iloc[0]
    assert mbp.match_status=='EXACT_CAS_MATCH'
    assert mbp.epa_cytotox_median_um==1000.0
    assert mbp.cytotox_1000_uM_possible_default_ceiling
    assert sub[sub.short=='MiBP'].iloc[0].match_status=='NO_EXACT_CAS_MATCH'
    assert pd.isna(sub[sub.short=='MiBP'].iloc[0].epa_nhit)
    assert parents[parents.metabolite=='MBP'].iloc[0].epa_ntested==70
    assert sub[sub.short=='MBP'].iloc[0].epa_ntested==50
    pooled=pd.DataFrame({'crp_harmonized_mg_L':[0.5,11.,2.0],
                         'survey_weight':[1.,1.,3.]})
    for c in m.CHEMICALS:
        pooled[c['nhanes_column']]=[1.,2.,3.]
        pooled[m.FLAG_MAP[c['nhanes_column']]]=[0.,1.,1.]
    px=Path(tmp)/'pooled.csv';pooled.to_csv(px,index=False)
    pp=m.nhanes_profiles(px)
    assert pp[0]['nhanes_n']==2
    assert abs(pp[0]['nhanes_weighted_pct_below_lod']-75)<1e-9
print('[PASS] V11: CAS-only identity, parent/metabolite separation, no-match as missing, default toxicity ceiling, pooled NHANES join')
