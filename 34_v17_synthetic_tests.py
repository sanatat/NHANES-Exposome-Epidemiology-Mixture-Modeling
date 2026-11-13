#!/usr/bin/env python3
"""Test V17 against a small self-contained synthetic fixture."""
import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import pandas as pd

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('review_v17', HERE / '33_v17_prioritize_screening_hits.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
rows = [
    dict(chemical='PFHxS', casrn='355-46-4', aeid=97, tier='CORE_12', assay_component_endpoint_name='ATG_NRF2_ARE_CIS',mapped_gene_symbols='NFE2L2',
         n_source_records=2,n_valid_hitc=2,n_high_positive=1,n_high_opposite=0,n_below_positive_threshold=1,n_missing_hitc=0,
         n_unique_sample_ids=2,n_repeated_sample_id_records=0,min_ac50_uM_from_verified_units=97.5,
         n_positive_records_with_verified_ac50=1,evidence_status='MIXED_HIGH_AND_BELOW'),
    dict(chemical='MBP',casrn='131-70-4',aeid=97,tier='CORE_12',assay_component_endpoint_name='ATG_NRF2_ARE_CIS',mapped_gene_symbols='NFE2L2',
         n_source_records=1,n_valid_hitc=1,n_high_positive=0,n_high_opposite=0,n_below_positive_threshold=1,n_missing_hitc=0,
         n_unique_sample_ids=1,n_repeated_sample_id_records=0,min_ac50_uM_from_verified_units=None,
         n_positive_records_with_verified_ac50=0,evidence_status='NO_HIGH_POSITIVE_HITC')]
raw = pd.DataFrame([
    {'casrn':'355-46-4','aeid':97,'hitc_numeric':.94,'ac50_uM_if_verified':97.5,'spid':'s1'},
    {'casrn':'355-46-4','aeid':97,'hitc_numeric':.2,'ac50_uM_if_verified':None,'spid':'s2'},
    {'casrn':'131-70-4','aeid':97,'hitc_numeric':.1,'ac50_uM_if_verified':None,'spid':'s3'},
])
selected = m.select_hits(pd.DataFrame(rows))
assert len(selected)==1
assert selected.iloc[0].positive_fraction_of_valid_records==.5
assert 'mixed_positive' in selected.iloc[0].review_notes
assert 'UNKNOWN__MC6' in selected.iloc[0].flag_review_status
assert len(m.details(raw,selected))==1
with TemporaryDirectory() as folder:
    p = Path(folder)
    pd.DataFrame(rows).to_csv(p/'v16.csv',index=False)
    raw.to_csv(p/'raw.csv',index=False)
    report = m.run(['--v16',str(p/'v16.csv'),'--raw',str(p/'raw.csv'),'--out',str(p/'out')])
    assert report['positive_pairs_core12']==1
    assert report['positive_source_records_all44']==1
    assert (p/'out/v17_core12_positive_review.csv').is_file()
print('[PASS] V17: exact pairs, missing != negative, mixed calls, AC50 units, no unverified confidence claims')
