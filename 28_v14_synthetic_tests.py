#!/usr/bin/env python3
"""Test v14 ranking, exact identity, missing vs negatives on fully synthetic input."""
import importlib.util, tempfile
from pathlib import Path
import pandas as pd
spec=importlib.util.spec_from_file_location('v14',Path(__file__).with_name('27_v14_curate_and_bridge.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
with tempfile.TemporaryDirectory() as d:
 d=Path(d)
 ids=sorted(m.CORE)+[9999]
 inp=d/'v13.csv'
 pd.DataFrame([{'aeid':i,'assay_component_endpoint_name':f'A{i}',
 'mapped_gene_symbols':'NFKB1','domains_matched':'Inflammatory signaling'} for i in ids]).to_csv(inp,index=False)
 pan=m.prepare_panel(inp)
 assert len(pan)==len(ids) and sum(pan.tier=='CORE_12')==12
 hits=d/'mc.csv'
 pd.DataFrame([
 {'casn':'355-46-4','aeid':210,'hitcall':0.97,'ac50':2.1},
 {'casn':'355-46-4','aeid':284,'hitcall':0.08,'ac50':30},
 {'casn':'375-95-1','aeid':210,'hitcall':None,'ac50':None},
 {'casn':'84-74-2','aeid':210,'hitcall':0.99,'ac50':1},
 ]).to_csv(hits,index=False)
 agg,qc=m.assay_hits(hits,pan)
 assert qc['matched_records']==3 and len(agg)==3
 assert agg[(agg.casrn=='355-46-4')&(agg.aeid==210)].n_positive_hitcall_ge_0_9.iloc[0]==1
 a,mat=m.main(['--shortlist',str(inp),'--mc56',str(hits),'--out',str(d/'out')])
 assert len(mat)==7*len(ids)
 assert mat[(mat.casrn=='355-46-4')&(mat.aeid==210)].evidence_status.iloc[0]=='OBSERVED_HITCALL'
 assert mat[(mat.casrn=='355-46-4')&(mat.aeid==9999)].evidence_status.iloc[0]=='NOT_IN_SUPPLIED_RESULTS'
 a2,mat2=m.main(['--shortlist',str(inp),'--out',str(d/'out2')])
 assert set(mat2.evidence_status)=={'NOT_RETRIEVED'}
print('[PASS] V14 core panel, strict CAS/AEID, positive vs missing, no parent substitution')
