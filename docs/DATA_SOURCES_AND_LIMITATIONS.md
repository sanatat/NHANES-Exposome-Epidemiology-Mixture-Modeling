# Data, provenance, and scientific guardrails

## Public sources

- **CDC NHANES (2015–2016 and 2017–2018)**: https://wwwn.cdc.gov/nchs/nhanes/ ; source SAS XPT data are downloaded at runtime by `08_build_pooled_2015_2018.py` and (for 2017–2018 only) `01_download_build.py`.
- **EPA CompTox / ToxCast invitrodb v4.3**: https://www.epa.gov/comptox-tools/exploring-toxcast-data ; official cytotoxicity and assay annotation sources are consumed by scripts v11/v12.
- **Optional v15 screening mirror**: `scigantic-comptox` accesses a **third-party** copy of invitrodb. This is not EPA's official validation service.

## Reproducible handling

- The pooled builder preserves the two original NHANES cycles and uses cycle-aware survey strata and appropriate multi-cycle weights. Double-check CDC guidance for each subsample before altering exposure families.
- A SAS XPT encoded-zero correction is applied to laboratory flag fields, with a synthetic test and row-level LOD audit. Do not treat every tiny positive concentration as a zero; inspect `01_download_build.py` and `08_build_pooled_2015_2018.py` for scope.
- Within-cycle effect comparison is exploratory; one period was used during hypothesis development.
- hs-CRP 2015–16 harmonization uses CDC documentation. Calibration and exclusion of out-of-range results are logged by the builder.
- Urinary creatinine enters as a covariate or urinary metabolite normalization in v8 sensitivity specifications; these are not interchangeable without model assumptions.
- Survey PSU/strata/weights matter. The two-year full joint 18-coefficient model exceeded residual design degrees of freedom; prefer pooled v6 models for main descriptive inference.
- Explicitly record data-dependent decisions: hs-CRP ≤10 mg/L primary cut-off, exclusion criteria, and secondary analyses.

## Interpretation warnings

1. **No evidence of statistically significant mixture effects** in the pooled primary models; the nominal MBP p-value did **not** survive FDR across seven chemical tests.
2. **No causality** from cross-sectional NHANES; ToxCast assays are in vitro and do not establish human biological mechanism.
3. `NOT_RECORDED_IN_SOURCE` does not equal an inactive chemical.
4. Sample records (`spid`) are not independent confirmatory experiments; not all high `hitc` values represent valid dose–response curves.
5. Distinguish `fitc=41` modeled responses from `fitc=100` LOEC-based activity. MC6 flag IDs require manual review; no confidence ranking is certified here.
6. The selected LC–MS/MS biomarker concentrations **do not** demonstrate practical processing of untargeted raw LC–HRMS or GC–HRMS signals.

## Data redistribution

Raw CDC, EPA and mirror files are excluded from the Git repository by `.gitignore` to avoid large repository size and ambiguous redistribution rights. Data provenance is stored in generated per-stage audit JSON files locally; users must review those and the original provider terms before reusing derived outputs.
