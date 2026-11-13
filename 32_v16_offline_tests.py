#!/usr/bin/env python3
"""Offline tests for V16. No EPA network access or real results are used."""
import importlib.util
import tempfile
from pathlib import Path
import pandas as pd
import numpy as np

script=Path(__file__).with_name("31_v16_mechanistic_evidence_audit.py")
spec=importlib.util.spec_from_file_location('v16',script)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
with tempfile.TemporaryDirectory() as td:
    td=Path(td)
    matrix=pd.DataFrame([
        ("PFOS","1763-23-1",97,"CORE_12","OBSERVED_IN_SOURCE"),
        ("PFOS","1763-23-1",94,"CORE_12","OBSERVED_IN_SOURCE"),
        ("MBP","131-70-4",97,"CORE_12","OBSERVED_IN_SOURCE"),
        ("MBP","131-70-4",94,"CORE_12","NOT_RECORDED_IN_SOURCE")],
        columns=["chemical","casrn","aeid","tier","retrieval_status"])
    panel=matrix.assign(assay_component_endpoint_name="Test assay",mapped_gene_symbols="NFE2L2")
    raw=pd.DataFrame([
        ("1763-23-1",97,.99,"P1",3.,"uM"),
        ("1763-23-1",97,.2,"P2",4.,"uM"),
        ("1763-23-1",94,-.94,"P3",np.nan,"unknown"),
        ("131-70-4",97,.3,"M1",np.nan,"uM")],
        columns=["casrn","aeid","hitc_numeric","spid","ac50_uM_if_verified","conc_unit"])
    p=m.calculate(matrix,raw,panel)
    by={(r.chemical,r.aeid):r for r in p.itertuples()}
    assert by[('PFOS',97)].evidence_status=='MIXED_HIGH_AND_BELOW'
    assert by[('PFOS',94)].evidence_status=='OPPOSITE_DIRECTION_REVIEW'
    assert by[('MBP',97)].evidence_status=='NO_HIGH_POSITIVE_HITC'
    assert by[('MBP',94)].evidence_status=='NOT_RECORDED_IN_SOURCE'
    assert by[('PFOS',97)].n_unique_sample_ids==2
    assert m.summarize(p).set_index('chemical').loc['PFOS','n_core12_with_positive_hitc']==1
    p.to_csv(td/'pairs.csv',index=False)
    assert (td/'pairs.csv').exists()
    print('[PASS] V16 audited synthetic sample-multiplicity, positive/mixed/opposite, missing vs inactive, and endpoint counts')
