# V15: Actual chemical-endpoint observations (experimental third-party mirror)

Your V14 dataset has 7 chemicals x 44 selected endpoints = 308 possible pairs but **no measured hit calls**. V15 queries an independently maintained *third-party S3 mirror* of US EPA invitrodb v4.3 using the `scigantic-comptox` Python package. **Not an official EPA API**; validate important findings against EPA's release before publication. Do not infer human toxicity from screening hits.

## Install / run (existing Linux Mint project)

```bash
cd ~/nhanes_exposome_project
.venv/bin/python -m pip install 'scigantic-comptox'
.venv/bin/python 30_v15_synthetic_tests.py
.venv/bin/python 29_v15_fetch_bioactivity.py
```

This requires internet and access to public S3 from your network. No EPA API key or 7.5GB download is needed if the mirror is reachable. If inaccessible, script exits without inventing results.

**Alternative official local table** (must contain exact `casn` or `casrn`, `aeid`, `hitc` or `hitcall`; no fallback from `SPID`):

```bash
.venv/bin/python 29_v15_fetch_bioactivity.py --input-csv /path/to/epa_filtered.csv
```

## Outputs

- `results/v15_measured_bioactivity/v15_chemical_x_aeid_matrix.csv` (all 308 potential pairs)
- `v15_core12_results.csv` (the core 12 per chemical)
- `v15_matched_raw_assays.csv` (raw rows, no implicit aggregation of chemical samples)
- `v15_audit.json` (provenance and coverage)

## Interpretation

`NOT_RECORDED_IN_SOURCE` is not evidence of inactivity. EPA `hitc` is **signed** in some assay contexts. High positive means strong hitcall in intended direction; high negative means an opposite-direction response and requires endpoint-specific review. Do not treat all negative values as negative assay results. Potency is retained only when concentration units explicitly equal uM. A high hitcall may be affected by cytotoxicity or assay interference; check flags and original protocol.

Source references:
- https://www.epa.gov/comptox-tools/exploring-toxcast-data
- https://pypi.org/project/scigantic-comptox/
- https://usepa.r-universe.dev/tcpl/doc/Introduction_Appendices.html
