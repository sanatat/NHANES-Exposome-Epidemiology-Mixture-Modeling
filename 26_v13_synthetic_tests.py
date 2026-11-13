#!/usr/bin/env python3
import importlib.util, tempfile
from pathlib import Path
import pandas as pd
path=Path(__file__).with_name('25_v13_audit_assay_candidates.py')
spec=importlib.util.spec_from_file_location('v13',path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
with tempfile.TemporaryDirectory() as t:
 t=Path(t);p=t/'v12.csv';out=t/'results'
 pd.DataFrame([
  {'aeid':101,'official_symbol':'NFKB1','inflammatory_signaling_candidate':True,'oxidative_stress_candidate':False,'immune_response_candidate':True},
  {'aeid':102,'official_symbol':'NFE2L2','inflammatory_signaling_candidate':False,'oxidative_stress_candidate':True,'immune_response_candidate':False},
  {'aeid':103,'official_symbol':'','inflammatory_signaling_candidate':False,'oxidative_stress_candidate':False,'immune_response_candidate':True},
  {'aeid':104,'official_symbol':'ESR1','inflammatory_signaling_candidate':False,'oxidative_stress_candidate':False,'immune_response_candidate':False},
 ]).to_csv(p,index=False)
 rep,ds,overlap=m.analyze(p,out)
 assert rep['n_unique_tagged_endpoints']==3
 assert rep['n_sum_of_domain_tags']==4
 assert rep['n_overlaps_extra_tags']==1
 assert rep['n_unique_target_gene_supported']==2
 assert rep['n_unique_annotation_only']==1
 assert (out/'v13_annotation_evidence.png').exists()
print('[PASS] V13 disjoint overlap, target-supported vs text-only, no hit claims')
