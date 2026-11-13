#!/usr/bin/env python3
"""Synthetic-only unit tests for v12; never claim these are EPA results."""
import importlib.util
import tempfile
from pathlib import Path
import pandas as pd
from zipfile import ZipFile

path=Path(__file__).with_name('23_v12_mechanistic_toxcast.py')
spec=importlib.util.spec_from_file_location('v12',path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert m.cas('131-70-4')=='131-70-4' and m.cas('bad')==''
assert len(m.make_cross())==7
with tempfile.TemporaryDirectory() as tmp:
    tmp=Path(tmp)
    ann=pd.DataFrame({'aeid':[101,102,103,104], 'aenm':['TOX21_NFkB','ROS assay','IL6 assay','hormone binding'], 'description':['NF-kB signaling','oxidative ROS stress','cytokine IL6','other']})
    tar=pd.DataFrame({'aeid':[101,102,103,104],'gene_symbol':['NFKB1','NFE2L2','IL6','ESR1']})
    ax=tmp/'ann.xlsx';tx=tmp/'targets.xlsx'
    ann.to_excel(ax,index=False);tar.to_excel(tx,index=False)
    ad=m.find_sheet(ax);td=m.find_sheet(tx);inv,e,n=m.curated_inventory(ad,td)
    assert int(inv['inflammatory_signaling_candidate'].sum())>=2
    assert int(inv['oxidative_stress_candidate'].sum())>=1
    assert int(inv['immune_response_candidate'].sum())>=1
    f=tmp/'mc.csv'
    pd.DataFrame({'casn':['131-70-4','131-70-4','131-70-4','375-95-1','84-74-2'],'aeid':[101,102,103,101,101],'hitcall':[.95,.1,.97,.99,1.0]}).to_csv(f,index=False)
    epi,cov,metadata=m.process_full(f,inv,e,n,m.make_cross())
    mbp=cov[cov.short=='MBP'].iloc[0]
    assert mbp.n_assay_endpoints==3 and mbp.n_active_hitcall_ge_0_9==2
    assert metadata['matched_rows']==4 # parent DBP row must NOT be included
    assert cov[cov.short=='MiBP'].iloc[0].n_assay_endpoints==0
    assert m.TERMS and m.FALLBACK
print('[PASS] V12 synthetic: assay annotation tagging, exact CAS, no parent substitution, hitcall threshold, missing coverage preserved')
