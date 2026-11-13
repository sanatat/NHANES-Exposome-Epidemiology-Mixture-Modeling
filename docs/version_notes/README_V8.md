# NHANES v8: Creatinine adjustment, robustness, and restricted cubic splines

**Status:** Exploratory post-hoc sensitivity analysis. This patch does **not** modify the v6 primary models or v7 cycle consistency analysis, and it does **not** search for significance. Original observational data concern hs-CRP, not cancer.

## Inputs

The project directory must already contain `data/processed/analysis_joint_pooled_2015_2018.csv` from v6 (NHANES 2015–16 / 2017–18). No new downloads are required. The input and pre-existing `results/` files are never overwritten by this patch; it writes **new files beginning with `v8_`**.

## Run on Linux Mint

```bash
cd ~/nhanes_exposome_project
.venv/bin/python 14_v8_offline_math_tests.py
.venv/bin/python 13_v8_preflight.py
Rscript 12_creatinine_robustness_splines.R
```

Alternatively: `bash RUN_V8.sh` runs preflight plus R. `survey` is required in R: `sudo apt install r-cran-survey`.

## Scientific design

- **Cohorts:** Use the same complete-case selection and covariates as v6: harmonized hs-CRP >0 and <=10 mg/L for the primary *exploratory* analysis; <=20 mg/L sensitivity. Age >=20.
- **Creatinine:** `URXMBP` urine MBP in `ng/mL`; `URXUCR` creatinine in `mg/dL`. Ratio **ng/g creatinine** = `100000 * URXMBP / URXUCR`. Compare (A) raw MBP plus log creatinine covariate (v6 reference), (B) creatinine-standardized MBP, **without also adjusting for creatinine**, and (C) raw MBP without dilution correction. These represent alternative approaches and partly different estimands, not interchangeable effect estimates.
- **Covariate robustness:** Refit with and without BMI or PIR and with demographic-only covariates. This *does not* determine which adjustment is causally right: that requires an explicit DAG and subject-matter assumptions.
- **Exposure tail robustness:** Winsorize MBP at locked survey-weighted 1st and 99th percentiles, selected without seeing outcomes.
- **Nonlinear exposure–response:** Four-knot restricted cubic spline on log2(MBP), exposure-weighted knot locations at 5/35/65/95 percentiles of the original primary cohort. A **2-degree-of-freedom design-based Wald test** assesses departures from linearity. Plot shows percent changes relative to the weighted median; 95% intervals account for survey design.
- **Design:** Preserve NHANES 2015–2018 two-cycle combined weights (subsample B weight divided by two), pooled strata, and PSU variables. Tests use residual design df from the fitted model.
- **Interpretation:** All findings in v8 are **post-hoc/exploratory** because MBP was highlighted after seeing v6 results. No independent replication or causal claims; do not pick the lowest p-value out of many specifications.

## Files generated

| File | Contents |
|---|---|
| `results/pooled_2015_2018/v8_MBP_robustness_all_specifications.csv` | All 14 sensitivity coefficients, CIs, df, p-values |
| `results/pooled_2015_2018/v8_MBP_RCS_nonlinearity_test.csv` | Nonlinearity Wald p-value and df |
| `results/pooled_2015_2018/v8_MBP_RCS_curve_exploratory.csv` | Dose–response curve with 95% CIs |
| `results/pooled_2015_2018/v8_MBP_exposure_response_spline.png` | Dose–response figure |
| `results/pooled_2015_2018/v8_locked_exposure_cutoffs.csv` | Exposure-only winsor thresholds and spline knots |
| `results/pooled_2015_2018/v8_methods_and_limits.txt` | Interpretive audit trail |

## Method references

- [CDC NHANES multiple-cycle weighting](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx)
- [CDC creatinine units (URXUCR)](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/ALB_CR_I.htm)
- [Survey-weighted regression, R `survey`](https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/svyglm.html)
- [Survey `regTermTest` of nonlinear terms](https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/regTermTest.html)
- [Creatinine correction and covariate dependence of urine dilution](https://pmc.ncbi.nlm.nih.gov/articles/PMC8380693/)

### Reproducibility disclosure

The Python preflight and numerical basis tests can be run without R; the creator's execution environment did not include `Rscript` or CDC samples. The statistical **results must be verified on the user's Linux machine**. This patch does not assert any new statistical significance before that run.
