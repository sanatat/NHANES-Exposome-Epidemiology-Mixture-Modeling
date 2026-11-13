# V10 — Fair matched-model cross-cycle comparison

**Why?** V9 compared a Ridge *covariate-only* baseline with ElasticNet/RandomForest *plus chemicals*, confounding feature value with algorithm choice. V10 fixes that by evaluating **all three model families with exactly the same feature set change**, holding the 2015–16 or 2017–18 cycle out for testing.

## Running (from the existing project folder)

```bash
cd ~/nhanes_exposome_project
.venv/bin/python 19_v10_synthetic_tests.py
.venv/bin/python 18_v10_paired_transfer.py
```

No download or R packages are needed beyond previous v9 dependencies. The script imports `16_v9_crosscycle_prediction.py` and reads the already-built `data/processed/analysis_joint_pooled_2015_2018.csv`; the script **does not modify** them or the v6–v9 result files.

## Outputs

`results/v10_paired_transfer/`

- `v10_matched_algorithm_incremental_value.csv`: 3 estimator families x 2 train/test cycle directions; baseline and +seven chemical weighted MSE, MAE and R², their paired differences and train-only chosen hyperparameters
- `v10_matched_algorithm_predictions.csv`: test-cycle per-person predictions for matched pairs
- `v10_descriptive_summary.csv`: two-direction *descriptive* aggregation, not a pooled inference
- `v10_interpretation.txt`: reporting safeguards

**Positive `delta_mse_positive_is_better`** means improved prediction when chemicals are added *for that very same algorithm*. A negative value means prediction worsened.

## Restrictions and scientific interpretation

- This is a **post-hoc comparison** after examining v9 results, not independent replication or pre-registered hypothesis testing.
- All data preprocessing and hyperparameter selection happen within the training cycle. Public-use strata+PSU labels are used to keep PSU groups together during inner validation.
- Cycle-transfer performance does **not** by itself imply causal associations or clinical utility. Predictive sample weights are not a substitute for survey-design-based standard errors.
- Cross-cycle outcome calibration and the hs-CRP <=10 mg/L cut-off are inherited from v6/v9 and are not re-optimized here.
- Do not choose an algorithm by test-cycle results and then claim an unbiased estimate for the chosen algorithm.
- Inferences about population-level significance require a separate survey-aware analysis.
- This project uses targeted NHANES chemical biomarker data, not raw LC-HRMS or GC-HRMS spectra.
