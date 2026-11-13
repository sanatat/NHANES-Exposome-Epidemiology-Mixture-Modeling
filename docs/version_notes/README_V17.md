# V17: Focused ToxCast hit review

Offline review using outputs from V15 and V16. This patch **does not download** EPA data, assign risk scores, claim pathway activation, or calculate p-values.

## Run

```bash
cd ~/nhanes_exposome_project
.venv/bin/python 34_v17_synthetic_tests.py
.venv/bin/python 33_v17_prioritize_screening_hits.py
```

Outputs in `results/v17_positive_hit_review/`:

- `v17_core12_positive_review.csv` — the 17 core12 chemical–endpoint pairs expected from V16 user logs.
- `v17_all44_positive_endpoint_pairs.csv` — all high-hit pairs among the 44 endpoints.
- `v17_positive_source_records.csv` — per-source record details for those pairs.
- `v17_evidence_gaps.json` — the required follow-up verification and caveats.

Numbers of positives in different `spid` samples are not equivalent to independent replicate validation. Curve-fit caution flags (`mc6`), max measured concentration and cytotoxicity interference are not available in the present V16 summarized table. These are **unverified** until linked against a trusted official EPA release. The third-party mirror is a secondary source. Readme status and labels intentionally prevent overclaiming.

## Verification references

- EPA official ToxCast data: https://www.epa.gov/comptox-tools/exploring-toxcast-data
- tcpl flags and assays: https://pmc.ncbi.nlm.nih.gov/articles/PMC10552852/
