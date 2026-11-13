#!/usr/bin/env python3
from importlib.machinery import SourceFileLoader
from pathlib import Path
import pandas as pd
import numpy as np

p=Path(__file__).with_name('29_v15_fetch_bioactivity.py')
m=SourceFileLoader('v15',str(p)).load_module()
pairs=pd.DataFrame([{'chemical':'PFHxS','casrn':'355-46-4','aeid':94,'tier':'CORE_12'},
                    {'chemical':'PFNA','casrn':'375-95-1','aeid':94,'tier':'CORE_12'},
                    {'chemical':'MBP','casrn':'131-70-4','aeid':296,'tier':'CORE_12'}])
pairs=m.validate_pairs(pairs)
raw=pd.DataFrame({'casn':['355-46-4','355-46-4','375-95-1','335-67-1'],
                  'aeid':[94,94,94,94],
                  'hitc':[0.99,-0.95,0.0,0.99],
                  'ac50':[12,10,np.nan,1.0],
                  'conc_unit':['uM','uM','uM','uM']})
norm=m.normalize_real_observations(raw,pairs)
assert len(norm)==3
assert norm.assay_result.tolist()==['HIGH_HITCALL_POSITIVE_DIRECTION','HIGH_HITCALL_OPPOSITE_DIRECTION','BELOW_0_9_MAGNITUDE']
mat=m.build_matrix(pairs,norm)
assert len(mat)==3
assert mat.loc[mat.casrn.eq('355-46-4'),'n_records'].item()==2
assert mat.loc[mat.casrn.eq('355-46-4'),'n_high_opposite'].item()==1
assert mat.loc[mat.casrn.eq('131-70-4'),'retrieval_status'].item()=='NOT_RECORDED_IN_SOURCE'
assert mat.loc[mat.casrn.eq('375-95-1'),'n_below_threshold'].item()==1
print('[PASS] V15 exact CAS+AEID, signed hitcall orientation, uM potency, missing != negative, duplicate sample records.')
