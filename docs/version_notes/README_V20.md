# V20 — Correct decode of numeric EPA Level 6 quality flag IDs

This patch uses **local v18 positive ToxCast screening records**. It does not redownload NHANES or EPA. The ToxCast source remains a third-party mirror; official EPA curve and concentration-range validation is not performed.

## Run

```bash
cd ~/nhanes_exposome_project
.venv/bin/python 40_v20_tests.py
.venv/bin/python 39_v20_decode_numeric_flags.py
```

## Outputs

Saved to `results/v20_numeric_flags/`:
- `v20_flag_frequency.csv`: counts by numeric flag ID. Multiple IDs per record can be counted, so counts do not sum to the number of records.
- `v20_positive_records_decoded.csv`: all 86 positive screening records with decoded flags, fit interpretation, raw fields.
- `v20_positive_pairs_qc.csv`: summary for 52 positive chemical × endpoint combinations.
- `v20_PFNA_PFOS_NRF2_focus.csv`: PFNA/PFOS endpoint 97 for manual review.
- `v20_audit.json`: totals and limitations.

## Interpretive caveats

- `fitc=41` denotes modeled responses with AC50 and AC95 within screened concentration range, pending curve/flag QC.
- `fitc=100`, `model_type=1`, or `modl=loec` indicates LOEC-based activity. Do **not** describe such an entry as a conventional modeled AC50.
- Flag IDs are method identifiers, **not numeric risk scores**. `13`: fewer than two replicates per concentration on average; `14`: four or fewer tested concentrations; `17`: low efficacy. `6`/`7`: single-point activity.
- Data `hitc >= 0.9` is only a screening positive. No human risk ranking or in-vivo mechanistic confirmation follows.

Official method descriptions: https://usepa.r-universe.dev/tcpl/doc/Introduction_Appendices.html (section *Flags (Level 6)*).
