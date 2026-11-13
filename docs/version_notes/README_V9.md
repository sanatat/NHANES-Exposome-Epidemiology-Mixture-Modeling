# V9 — Utrecht-paper-inspired Exposome Atlas and Cross-cycle Prediction

This **exploratory** patch builds on the real v6–v8 NHANES 2015–2018 pooled dataset. It **does not overwrite** existing outputs or re-download data, and is not cancer research, causal inference, or raw LC-HRMS analysis.

## Research inspiration (not replication)

1. Ohanyan et al. (2022), with Roel Vermeulen, compared LASSO, Random Forest, and ANN for multi-exposure urban exposome prediction. DOI: https://doi.org/10.1016/j.envint.2022.107592
2. Kolijn et al. (2025) emphasized cross-cohort validation of high-dimensional biomarkers. DOI: https://doi.org/10.1038/s41467-025-64534-4
3. Young et al. (2025), including Vermeulen and Kolijn, identified gaps in broad exposure profiling. DOI: https://doi.org/10.1093/jncics/pkae122
4. Dai et al. (2026), including Vermeulen, combined human exposure and toxicity/bioactivity to prioritize chemicals. DOI: https://doi.org/10.1016/j.envint.2026.110110

## Install (once)

```bash
cd ~/nhanes_exposome_project
.venv/bin/python -m pip install 'scikit-learn>=1.5,<2' 'scipy>=1.11,<2'
.venv/bin/python 17_v9_offline_tests.py
.venv/bin/python 15_v9_exposure_atlas.py
.venv/bin/python 16_v9_crosscycle_prediction.py
```

## Output

`results/v9_professor_inspired/`:
- `v9_exposure_profiles_by_cycle.csv`: weighted concentration summaries and below-LOD percentages
- `v9_coexposure_correlations.csv`: weighted Spearman correlations by cycle
- `v9_coexposure_stable_edges.csv`: co-exposure correlation present at |rho|>=0.25 with same sign in both cycles, chosen without reference to outcome
- `v9_weighted_correlation_2015.png` and `_2017.png`: two atlas heatmaps
- `v9_crosscycle_prediction_comparison.csv`: three models, both cycle directions, group-aware train-only tuning and weighted held-out errors
- `v9_crosscycle_predictions.csv`: per-person held-out log-hsCRP predictions
- `v9_prediction_disclaimer.txt`: limitations

## Inference limitations

This is predictive transfer, not cross-cohort replication. Each cycle is a separate NHANES public-use period, but both belong to the same survey system. Do not use prediction scores or differences in weighted MSE as survey-inference p-values. Grouped inner CV uses public design PSU groups, but cannot emulate full replicate-weight uncertainty. The matrix links chemistry from different biospecimens (serum and urine) and does not establish common exposure timing. Excluding CRP>10 mg/L was a post-hoc decision from the previous project. A model that predicts worse after adding chemicals is informative; do not hide results. Predictive benefit does not prove biological causation.

## Next evidence-linked addition (not implemented)

The Vermeulen group has developed a human exposome / toxicity prioritization framework. A defensible extension requires chemically correct EPA CompTox IDs for seven substances **and separate documentation of phthalate biomarkers as metabolites rather than parent plasticizers**. Do not create a toxicity score until identifiers, assay relevance, detection and exposure measurement are validated. See EPA dashboard: https://comptox.epa.gov/dashboard/.

To demonstrate **untargeted LC-HRMS data analysis**, the project needs a different dataset with instrument features, MS/MS spectra, or spectra-to-annotation results. The NHANES targeted assays here do not fulfill this requirement.
