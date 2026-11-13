# V16 | ToxCast sample-level QC and biological evidence audit

This patch uses existing v14/v15 outputs and **does not download anything or change earlier results**. It distinguishes *assay coverage* from *activity* and retains per-sample hitcalls without pretending multiple records are independent confirmations.

## Run on Linux Mint

```bash
cd ~/nhanes_exposome_project
.venv/bin/python 32_v16_offline_tests.py
.venv/bin/python 31_v16_mechanistic_evidence_audit.py
```

## Outputs

- `results/v16_bioactivity_quality_audit/v16_per_chemical_endpoint_evidence.csv` — each of 308 candidate CAS–AEID pairs, coverage, high-positive/opposite hits, sample-ID counts when present, and potentcy only when units are verified.
- `v16_chemical_coverage_and_hits.csv` — chemical-level coverage and high-hit counts (descriptive **not** relative human hazard).
- `v16_core12_status.csv` — 84 pairs covering the 12 preselected endpoints and seven chemicals.
- `v16_core12_hitcall_heatmap.png` — labels are *high-positive records / valid hitcalls*, with missing assays clearly distinguished.
- `v16_audit.json` — provenance and method caveats.

## Interpretation safeguards

1. Hitc ≥0.9 is a fitted in-vitro high hitcall, **not** proof of toxicity or causal effect in humans.
2. An unrecorded chemical–assay pair is not an inactive outcome.
3. Multiple records may be distinct preparations/samples for the same chemical and are not independent biological replications; FDA/EPA QC, `mc5_chid` representative-sample designation, and `mc6` flags need primary-source audit.
4. Nrf2/ARE is an indirect antioxidant-response reporter and cannot establish direct molecular binding to NFE2L2.
5. The mirror is third-party. Validate claims against official `invitrodb v4.3` release before publication.
6. AC50 is only carried forward when the original source explicitly verifies `µM` units.

Sources: https://www.epa.gov/comptox-tools/exploring-toxcast-data and https://usepa.r-universe.dev/tcpl/doc/Introduction_Appendices.html
