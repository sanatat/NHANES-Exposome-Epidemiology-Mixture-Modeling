# NHANES Exposome project — v6 four-year survey analysis

This is an **add-on** for the existing v5 project, not a replacement for the original 2017–2018 data. It does not overwrite the original database or result files.

## Run under Linux Mint

From the user's home directory after unzipping the v6 patch:

```bash
cd ~/nhanes_exposome_project
.venv/bin/python 10_test_pooled_builder_offline.py
.venv/bin/python 08_build_pooled_2015_2018.py
Rscript 09_pooled_survey_mixture.R
```

The builder downloads **six additional 2015–2016 CDC SAS XPT files** (`DEMO_I`, `BMX_I`, `PFAS_I`, `PHTHTE_I`, `HSCRP_I`, `ALB_CR_I`) and reads the **six previously downloaded 2017–2018 files**. Files already cached in `data/raw` are reused. It makes:

- `data/processed/analysis_joint_pooled_2015_2018.csv`
- `data/processed/v6_pooled_qc.json`
- `data/nhanes_exposome_pooled_2015_2018.sqlite` (table `analysis_joint_pooled`)

The survey model writes under `results/pooled_2015_2018/`:

- `PRIMARY_pooled_mixtures.csv`
- `PRIMARY_pooled_single_chemical_epidemiology.csv`
- `sensitivity_pooled_mixtures.csv`
- `cohort_and_design_qc.csv`
- `lod_by_cycle.csv`
- `locked_weighted_quartile_cutpoints.csv`
- `METHODS_AND_LIMITATIONS.txt`

## Methodological upgrades

1. Pool two genuine two-year NHANES cycles (2015–2016 and 2017–2018). Derive a four-year environmental subsample weight as `WTSB2YR / 2` **only after testing person-specific PFAS and phthalate weights within each cycle**. Unique survey strata are created per cycle and the model controls for cycle.
2. Apply CDC hs-CRP assay bridge when 2015–2016 original hs-CRP is **at most 23 mg/L**: `2017_scale = 0.8695 * 2015_scale + 0.2954`. Never extrapolate to higher values.
3. Define outcome as log hs-CRP. The primary exploratory low-grade inflammation analysis uses harmonized hs-CRP **<=10 mg/L**. A sensitivity run uses harmonized hs-CRP <=20 mg/L. The <=10 threshold has already been explored in the earlier pilot, so **do not call this preregistered or confirmatory**.
4. Keep the original race/ethnicity, income/poverty, age, sex, BMI, and log urine creatinine covariates, plus cycle. Including BMI may constitute overadjustment under some causal assumptions; it must be assessed against a causal DAG before causal interpretation.
5. Use the `survey` design (`SDMVPSU`, nested `SDMVSTRA`, four-year B subsample weights). Report design/residual degrees of freedom, survey-corrected standard errors and 95% CIs.
6. Construct **population-weighted quartile cutpoints once from the primary pooled analytic cohort**, then lock them for sensitivity analyses; each exposure is quantized into 4 levels. Estimate a **linear-additive joint quantile-shift association**, not nonlinear BKMR.
7. Report all 3 predeclared mixture families (7 exposures / four PFAS / three phthalate metabolites) and all 7 predeclared single-pollutant models. Correct p-values by BH-FDR **within each family of multiple comparisons**, and report nonsignificant results transparently.
8. Keep PFAS concentration in serum and phthalate metabolite concentration in urine; use log urine creatinine as a model covariate (not double-standardizing the urinary biomarkers).

## Limits

These are observational **cross-sectional** data, not a causal experiment or longitudinal early-onset cancer dataset. hs-CRP is an unspecific inflammatory marker. There may be unmeasured confounders, measurement error, selected-subsample effects, and residual nonresponse bias. Extra participants/design degrees of freedom **do not guarantee smaller p-values**, a significant effect, or a clinically meaningful relationship.

`2015–2016` + `2017–2018` is a four-year combination. Do not substitute the `2017–March 2020` **pre-pandemic** series without rewriting the weighting scheme.

The Python integration test uses **synthetic fixture data only** for code logic; no clinical or epidemiological results are derived from that test. Actual CDC download and R fit must be run on the user's computer.

## Official sources

- CDC NHANES 2015–2016 PFAS: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PFAS_I.htm
- CDC 2015–2016 phthalates: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PHTHTE_I.htm
- CDC 2017–2018 hs-CRP assay bridge: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/HSCRP_J.htm
- CDC survey weighting: https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx
- CDC complex sampling design: https://wwwn.cdc.gov/nchs/nhanes/tutorials/SampleDesign.aspx
