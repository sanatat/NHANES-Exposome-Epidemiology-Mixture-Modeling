# V14 — Pre-specified Core Mechanistic Panel and Bioactivity Bridge

This patch builds a **12-endpoint mechanistic priority panel** from the 44 gene-annotated ToxCast v4.3 endpoint candidates in V13. It also builds a 7 chemicals × 44 endpoints matrix. **Without actual EPA results, all chemical × assay rows are explicitly `NOT_RETRIEVED`. No assay activity is inferred.**

## Run (does not download anything)

```bash
cd ~/nhanes_exposome_project
.venv/bin/python 28_v14_synthetic_tests.py
.venv/bin/python 27_v14_curate_and_bridge.py
```

## Optional: Use actual ToxCast v4.3 MC5-6 results

Once you independently obtain the official `mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv` (filename per release documentation), run:

```bash
.venv/bin/python 27_v14_curate_and_bridge.py --mc56 /path/to/mc5-6.csv
```

CSV must provide **CAS (casn or casrn), aeid and hitcall**. If only `spid` is present, a verified EPA mapping of samples to CAS is essential; script stops instead of guessing. The CSV is streamed in chunks and only 7 CAS × 44 AEID records are retained.

`hitcall >= 0.9` counts as an active assay readout, **not causal biological injury**. Results flagged for fitting and cytotoxicity interference require follow-up; raw numeric potency values keep original file units and are NOT assumed to be µM. Multiple samples per chemical are explicitly counted. `NOT_RETRIEVED`, `NOT_IN_SUPPLIED_RESULTS`, and `OBSERVED_NO_VALID_HITCALL` are NOT negative.

Outputs in `results/v14_mechanistic_panel/` include `v14_core_12_panel.csv`, `v14_curated_44_panel.csv`, `v14_chemical_x_endpoint_evidence.csv` and `v14_audit.json`.

Source: US EPA ToxCast invitrodb v4.3 (2025), DOI 10.23645/epacomptox.6062623.v14. Assay metadata from V12, exploratory target shortlist from V13.
