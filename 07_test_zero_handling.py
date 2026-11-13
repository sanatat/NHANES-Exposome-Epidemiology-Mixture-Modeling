#!/usr/bin/env python3
"""Self-contained regression tests: pandas' zero artifact and NHANES LOD flags."""
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd
here = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('nhanes_build', here / '01_download_build.py')
build = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build)
flag = build.PFAS_FLAGS[0]
artifact = build.SAS_XPORT_ZERO_ARTIFACT
raw = pd.DataFrame({flag: [artifact, 1.0, artifact, 0.0, np.nan],
                    'actual_small_measurement': [1e-10, 2e-10, 3e-10, 4e-10, np.nan]})
fixed = build.correct_sas_xport_zeros(raw, 'unit-test')
assert fixed[flag].tolist()[:4] == [0., 1., 0., 0.]
assert fixed['actual_small_measurement'].iloc[0] == 1e-10
assert pd.isna(fixed[flag].iloc[4])
assert raw[flag].iloc[0] == artifact, 'Input should not be mutated'
# Make all PFAS flags consistent and validate the codebook categories.
for f in build.PFAS_FLAGS: fixed[f] = fixed[flag]
build.assert_cdc_lod_flags(fixed, 'PFAS_J')
assert int(fixed[flag].eq(0).sum()) == 3
assert int(fixed[flag].eq(1).sum()) == 1
print('[PASS] NHANES SAS XPORT zero repair: detected=3, below_LOD=1, missing=1')
