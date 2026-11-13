#!/usr/bin/env python3
"""v9: Cross-cycle predictive transfer, inspired by Utrecht exposome ML research.

Strict 2015->2017 and 2017->2015 held-out transfer. Only inner training
cycle uses GroupKFold (public survey PSU) for hyperparameter choice.
WEIGHTED predictive error is not a survey-robust hypothesis test or effect.
All models/hyperparameters are fixed before seeing results.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge, ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'data/processed/analysis_joint_pooled_2015_2018.csv'
OUT=ROOT/'results/v9_professor_inspired'
EXPOSURES=['LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB']
CONTROLS=['RIDAGEYR','BMXBMI','INDFMPIR','log2_urine_creatinine']
CATEGORICAL=['RIAGENDR','RIDRETH3']

# No cycle variable is used as predictor; each train fold has a single cycle.
# CRP is log-transformed, only measured exposures enter +chemistry models.

def prep_data(d):
    fields=['SEQN','cycle_start','survey_weight','pooled_strata','pooled_psu',
            'RIDAGEYR','BMXBMI','INDFMPIR','RIAGENDR','RIDRETH3','URXUCR',
            'crp_harmonized_mg_L']+EXPOSURES
    missing=set(fields)-set(d)
    if missing: raise ValueError(f'Missing columns {sorted(missing)}')
    d=d.loc[d.RIDAGEYR.ge(20) & d.crp_harmonized_mg_L.gt(0) &
            d.crp_harmonized_mg_L.le(10)].copy()
    if d.SEQN.duplicated().any(): raise ValueError('Duplicate participants')
    if set(d.cycle_start.unique())!={2015,2017}: raise ValueError('Missing cycle')
    if d[fields].isna().any().any(): raise ValueError('Missing data in primary cohort')
    for col in EXPOSURES+['URXUCR','survey_weight']:
        if not np.isfinite(d[col].to_numpy()).all() or d[col].le(0).any():
            raise ValueError(f'Invalid {col}')
    for chem in EXPOSURES:
        d['log2_'+chem.lower()]=np.log2(d[chem].astype(float))
    d['log2_urine_creatinine']=np.log2(d.URXUCR.astype(float))
    d['log_hscrp']=np.log(d.crp_harmonized_mg_L)
    d['psu_group']=d.pooled_strata.astype(str)+'/'+d.pooled_psu.astype(str)
    if min(d.cycle_start.value_counts())<100: raise ValueError('Too few per cycle')
    return d

def weighted_mse(y,yhat,w): return float(np.average((np.asarray(y)-np.asarray(yhat))**2,weights=w))
def weighted_mae(y,yhat,w): return float(np.average(abs(np.asarray(y)-np.asarray(yhat)),weights=w))

def build_pipeline(name,features,pars):
    nums=[f for f in features if f not in CATEGORICAL]
    cats=[f for f in features if f in CATEGORICAL]
    proc=ColumnTransformer([
        ('num',Pipeline([('impute',SimpleImputer(strategy='median')),('scale',StandardScaler())]),nums),
        ('cat',Pipeline([('impute',SimpleImputer(strategy='most_frequent')),('onehot',OneHotEncoder(handle_unknown='ignore',sparse_output=False))]),cats)
    ], sparse_threshold=0)
    if name=='Ridge (covariates)': model=Ridge(alpha=pars['alpha'])
    elif name=='ElasticNet (+chemicals)': model=ElasticNet(alpha=pars['alpha'],l1_ratio=pars['l1_ratio'],max_iter=20000,random_state=42)
    elif name=='RandomForest (+chemicals)': model=RandomForestRegressor(
        n_estimators=80,min_samples_leaf=pars['min_samples_leaf'],max_depth=pars['max_depth'],
        random_state=42,n_jobs=1)
    else: raise ValueError('Unknown model')
    return Pipeline([('prep',proc),('model',model)])

def hyperparameters(name):
    if name=='Ridge (covariates)': return [{'alpha':x} for x in (0.1,10.,100.)]
    if name=='ElasticNet (+chemicals)': return [dict(alpha=a,l1_ratio=r) for a in (0.01,0.1,0.5) for r in (0.2,0.8)]
    if name=='RandomForest (+chemicals)': return [dict(min_samples_leaf=l,max_depth=md) for l,md in [(15,5),(40,6),(40,None)]]
    raise ValueError(name)

def select_by_group_cv(data,features,name):
    x=data[features]; y=data.log_hscrp.to_numpy(float)
    w=data.survey_weight.to_numpy(float); g=data.psu_group.to_numpy()
    if len(set(g))<6: raise ValueError('Not enough distinct public PSUs for 3-fold CV')
    # Normalize weights within folds to keep regularization comparable across folds.
    splitter=GroupKFold(n_splits=3)
    candidates=[]
    for pars in hyperparameters(name):
        errs=[]
        for tr,va in splitter.split(x,y,groups=g):
            model=build_pipeline(name,features,pars)
            wt=w[tr]/np.mean(w[tr]); model.fit(x.iloc[tr],y[tr],model__sample_weight=wt)
            pred=model.predict(x.iloc[va]); errs.append(weighted_mse(y[va],pred,w[va]))
        candidates.append((float(np.mean(errs)),pars))
    return sorted(candidates,key=lambda pair:pair[0])[0]

def run(d):
    d=prep_data(d)
    OUT.mkdir(parents=True,exist_ok=True)
    chem_features=['log2_'+c.lower() for c in EXPOSURES]
    specs=[('Ridge (covariates)',CONTROLS+CATEGORICAL),
           ('ElasticNet (+chemicals)',CONTROLS+CATEGORICAL+chem_features),
           ('RandomForest (+chemicals)',CONTROLS+CATEGORICAL+chem_features)]
    records=[]; preds=[]
    for src,tgt in [(2015,2017),(2017,2015)]:
        train=d.loc[d.cycle_start==src].reset_index(drop=True)
        test=d.loc[d.cycle_start==tgt].reset_index(drop=True)
        mean_train=np.average(train.log_hscrp,weights=train.survey_weight)
        mean_mse=weighted_mse(test.log_hscrp,np.repeat(mean_train,len(test)),test.survey_weight)
        print(f'\n[V9 TRANSFER] Train {src}-{src+1}: {len(train)}; test {tgt}-{tgt+1}: {len(test)}',flush=True)
        baseline_mse=None
        for name,features in specs:
            inner_mse,pars=select_by_group_cv(train,features,name)
            model=build_pipeline(name,features,pars)
            wtrain=train.survey_weight.to_numpy(float)
            model.fit(train[features],train.log_hscrp.to_numpy(),model__sample_weight=wtrain/np.mean(wtrain))
            yhat=model.predict(test[features]); yy=test.log_hscrp.to_numpy(); w=test.survey_weight.to_numpy()
            mse=weighted_mse(yy,yhat,w)
            mae=weighted_mae(yy,yhat,w)
            if name=='Ridge (covariates)': baseline_mse=mse
            record={'train_cycle':f'{src}-{src+1}','test_cycle':f'{tgt}-{tgt+1}',
                'model':name,'n_train':len(train),'n_test':len(test),
                'inner_psu_group_cv_mse':inner_mse,'selected_hyperparams':json.dumps(pars,sort_keys=True),
                'test_weighted_mse':mse,'test_weighted_mae':mae,
                'test_weighted_r2_vs_training_weighted_mean':1-mse/mean_mse,
                'pct_mse_improvement_over_demographics_only':100*(1-mse/baseline_mse),
                'type':'out-of-cycle prediction, not causal or survey design CI'}
            records.append(record)
            print(f' {name}: weighted test MSE={mse:.4f}, R2={record["test_weighted_r2_vs_training_weighted_mean"]:.3f}',flush=True)
            for seq,p,t in zip(test.SEQN.astype(int),yhat,yy):
                preds.append({'SEQN':seq,'train_cycle':src,'test_cycle':tgt,'model':name,
                              'predicted_log_hscrp':p,'observed_log_hscrp':t})
    pd.DataFrame(records).to_csv(OUT/'v9_crosscycle_prediction_comparison.csv',index=False)
    pd.DataFrame(preds).to_csv(OUT/'v9_crosscycle_predictions.csv',index=False)
    (OUT/'v9_prediction_disclaimer.txt').write_text(
        'Pure prediction transfer test; cross-cycle analyses share one national survey, not independent cohorts.\n'
        'Within-cycle CV groups by public PSU; test cycle never influences scaling, imputations or hyperparameter selection.\n'
        'Weights used to fit and calculate test errors. Conventional errors are descriptive, NOT design-based CIs.\n'
        'The 2015-16 hs-CRP calibration bridge and outcome threshold were decided in earlier post-hoc work.\n'
        'Negative R2 or worse errors for chemical models are valid findings, do not suppress them.\n'
        'No p-values are computed; this cannot establish harmful or beneficial chemical effects.\n',encoding='utf-8')
    print('[V9 TRANSFER] files in:',OUT)
    return pd.DataFrame(records)

if __name__=='__main__':
    if not SOURCE.exists(): raise SystemExit('Missing pooled CSV; run 08_build_pooled_2015_2018.py')
    run(pd.read_csv(SOURCE))
