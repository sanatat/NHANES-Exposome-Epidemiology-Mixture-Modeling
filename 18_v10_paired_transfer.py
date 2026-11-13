#!/usr/bin/env python3
"""v10: Fair, paired out-of-cycle prediction comparison for NHANES exposome.

For EACH estimator, compare (a) demographics/creatinine and (b) SAME
estimator/hyperparameter search plus chemicals. Hyperparameter tuning uses
only training-period PSU-grouped folds, and test-period records are untouched
until evaluation. Weighting is prediction-oriented, not design-based inference.
"""
from __future__ import annotations

from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'data/processed/analysis_joint_pooled_2015_2018.csv'
OUT = ROOT / 'results/v10_paired_transfer'


def v9_loader():
    script = ROOT / '16_v9_crosscycle_prediction.py'
    if not script.exists():
        raise FileNotFoundError(f'Missing required v9 script: {script}')
    spec = importlib.util.spec_from_file_location('v9_predictor_reuse', script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def estimator_name(family):
    return {
        'Ridge': 'Ridge (covariates)',
        'ElasticNet': 'ElasticNet (+chemicals)',
        'RandomForest': 'RandomForest (+chemicals)',
    }[family]


def paired_comparison(d, v9, output_dir=OUT, *, families=('Ridge','ElasticNet','RandomForest'), quiet=False):
    """Compute both feature sets using the same model family and held-out people."""
    d = v9.prep_data(d)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    covariates = v9.CONTROLS + v9.CATEGORICAL
    chemicals = ['log2_'+c.lower() for c in v9.EXPOSURES]
    records = []
    predictions = []
    for train_cycle, test_cycle in [(2015,2017), (2017,2015)]:
        train = d[d.cycle_start == train_cycle].reset_index(drop=True)
        test = d[d.cycle_start == test_cycle].reset_index(drop=True)
        assert not set(train.SEQN).intersection(test.SEQN)
        yt = test.log_hscrp.to_numpy(float)
        wt = test.survey_weight.to_numpy(float)
        mean_train = np.average(train.log_hscrp, weights=train.survey_weight)
        constant_mse = v9.weighted_mse(yt, np.full(len(test), mean_train), wt)
        for family in families:
            name = estimator_name(family)
            by_feature_set = {}
            for label, features in [('covariates',covariates), ('covariates_plus_7_chemicals',covariates+chemicals)]:
                cv_err, parameters = v9.select_by_group_cv(train,features,name)
                fitted = v9.build_pipeline(name,features,parameters)
                wtrain = train.survey_weight.to_numpy(float)
                fitted.fit(train[features], train.log_hscrp.to_numpy(float),
                           model__sample_weight=wtrain / wtrain.mean())
                predicted = fitted.predict(test[features])
                by_feature_set[label] = {
                    'mse':v9.weighted_mse(yt,predicted,wt),
                    'mae':v9.weighted_mae(yt,predicted,wt),
                    'cv_mse':cv_err,
                    'params':parameters,
                    'predicted':predicted,
                }
                for seq,truth,pred,weight in zip(test.SEQN.to_numpy(int),yt,predicted,wt):
                    predictions.append({
                        'SEQN':int(seq), 'train_cycle':train_cycle,'test_cycle':test_cycle,
                        'model_family':family,'feature_set':label,'log_hscrp_true':float(truth),
                        'log_hscrp_predicted':float(pred),'test_survey_weight':float(weight)})
            base=by_feature_set['covariates']
            full=by_feature_set['covariates_plus_7_chemicals']
            # Positive improvement = adding chemicals reduces weighted squared error.
            delta_mse = base['mse'] - full['mse']
            result={
                'train_cycle':f'{train_cycle}-{train_cycle+1}',
                'test_cycle':f'{test_cycle}-{test_cycle+1}',
                'model_family':family,'n_train':len(train),'n_test':len(test),
                'covariates_only_weighted_mse':base['mse'],
                'plus_chemicals_weighted_mse':full['mse'],
                'delta_mse_positive_is_better':delta_mse,
                'pct_mse_improvement_positive_is_better':100*delta_mse/base['mse'],
                'delta_mae_positive_is_better':base['mae']-full['mae'],
                'covariates_only_weighted_r2_vs_train_mean':1-base['mse']/constant_mse,
                'plus_chemicals_weighted_r2_vs_train_mean':1-full['mse']/constant_mse,
                'delta_weighted_r2_positive_is_better':delta_mse/constant_mse,
                'train_psu_group_cv_mse_covariates':base['cv_mse'],
                'train_psu_group_cv_mse_plus_chemicals':full['cv_mse'],
                'selected_covariates_hyperparameters':json.dumps(base['params'],sort_keys=True),
                'selected_plus_chemicals_hyperparameters':json.dumps(full['params'],sort_keys=True),
            }
            records.append(result)
            if not quiet:
                flag='IMPROVED' if delta_mse>0 else 'WORSE' if delta_mse<0 else 'SAME'
                print(f'[V10] {train_cycle}->{test_cycle}, {family}: '
                      f'base MSE {base["mse"]:.4f}, chemicals MSE {full["mse"]:.4f}, '
                      f'delta {delta_mse:+.4f} ({flag})',flush=True)
    df=pd.DataFrame(records)
    df.to_csv(output_dir/'v10_matched_algorithm_incremental_value.csv',index=False)
    pd.DataFrame(predictions).to_csv(output_dir/'v10_matched_algorithm_predictions.csv',index=False)
    summary=(df.groupby('model_family',sort=False)
             .agg(n_cycle_transfers=('delta_mse_positive_is_better','size'),
                  transfers_improved=('delta_mse_positive_is_better',lambda x:int((x>0).sum())),
                  average_delta_mse=('delta_mse_positive_is_better','mean'),
                  average_delta_r2=('delta_weighted_r2_positive_is_better','mean'))
             .reset_index())
    # Average of just two cycles is DESCRIPTIVE, not a pooled study estimate.
    summary.to_csv(output_dir/'v10_descriptive_summary.csv',index=False)
    (output_dir/'v10_interpretation.txt').write_text(
        'Matched-family, paired held-out comparison: each algorithm is evaluated with and without 7 chemical biomarkers.\n'
        'Positive delta MSE/R2 means that, within the SAME algorithm and test cycle, adding chemicals improved prediction.\n'
        'Tuning and preprocessing are learned from each training cycle only; PSU groups kept together in inner folds.\n'
        'Weighted predictive metrics are DESCRIPTIVE: do NOT interpret as formal survey-design confidence intervals or causal effects.\n'
        'This analysis is post hoc after inspecting v9 outputs, not an independently pre-registered validation.\n'
        'NHANES 2015-16 and 2017-18 cycles are not independently recruited external cohorts.\n'
        'The complete-case adult <=10 mg/L cohort, clinical covariates and 7 analytes are inherited unchanged from v9.\n'
        'Slight variation in cross-cycle metrics is not proof of added biological signal.\n',
        encoding='utf-8')
    if not quiet: print('[V10] Results:',output_dir)
    return df


if __name__ == '__main__':
    if not SOURCE.exists():
        raise SystemExit('Missing pooled CSV; run 08_build_pooled_2015_2018.py first')
    p=v9_loader()
    paired_comparison(pd.read_csv(SOURCE), p)
