#!/usr/bin/env Rscript
# 2015-2018 pooled NHANES: survey-weighted Exposome/Epidemiology/Mixture.
# Pre-declared scientific estimands. NEVER choose a model by smallest p-value.
# NOT causal inference. Does not process untargeted HRMS raw files.
args <- grep('^--file=', commandArgs(trailingOnly=FALSE), value=TRUE)
if (length(args) != 1) stop('Use: Rscript 09_pooled_survey_mixture.R')
root <- dirname(normalizePath(sub('^--file=', '', args[1])))
if (!requireNamespace('survey', quietly=TRUE)) {
  stop('R survey package required: sudo apt install r-cran-survey')
}
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust')
infile <- file.path(root,'data','processed','analysis_joint_pooled_2015_2018.csv')
outdir <- file.path(root,'results','pooled_2015_2018')
dir.create(outdir, recursive=TRUE,showWarnings=FALSE)
if (!file.exists(infile)) stop('Pooled data missing; run: .venv/bin/python 08_build_pooled_2015_2018.py')

exposures <- c('LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB')
flags <- c('LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL','URDMBPLC','URDMEPLC','URDMIBLC')
needed <- c('SEQN','cycle_start','survey_weight','pooled_strata','pooled_psu',
            'crp_harmonized_mg_L','RIDAGEYR','RIAGENDR','RIDRETH3',
            'INDFMPIR','BMXBMI','URXUCR', exposures,flags)
d0 <- read.csv(infile, check.names=FALSE)
if (!all(needed %in% names(d0))) stop('Missing variables: ',paste(setdiff(needed,names(d0)),collapse=', '))
if (anyDuplicated(d0$SEQN)) stop('Pooled SEQN duplicated')
if (!setequal(unique(d0$cycle_start),c(2015,2017))) stop('Unexpected/missing NHANES cycles')
if (any(!is.finite(d0$survey_weight) | d0$survey_weight<=0)) stop('Weights are missing or invalid')
if (any(!is.finite(d0$URXUCR) | d0$URXUCR<=0)) stop('Urine creatinine missing/invalid')
for (chem in exposures) {
  if (any(!is.finite(d0[[chem]]) | d0[[chem]]<=0)) stop('Invalid exposure: ',chem)
}
for (flag in flags) if (!all(d0[[flag]] %in% c(0,1))) stop('Invalid LOD flag: ',flag)

# Audit each cycle separately and on the pooled file. Counts/percentages must
# be independent of the outcome restriction below.
lod <- list()
for (cyc in c(2015,2017)) for (i in seq_along(exposures)) {
  sub <- d0[d0$cycle_start==cyc,,drop=FALSE]
  f <- as.numeric(sub[[flags[i]]]); w <- sub$survey_weight
  lod[[length(lod)+1]] <- data.frame(
    cycle=paste0(cyc,'-',cyc+1),chemical=exposures[i],n=nrow(sub),
    below_lod=sum(f==1),pct_below_lod_unweighted=100*mean(f==1),
    pct_below_lod_weighted=100*sum(w[f==1])/sum(w))
}
write.csv(do.call(rbind,lod),file.path(outdir,'lod_by_cycle.csv'),row.names=FALSE)

# Primary target: low-grade/chronic inflammation (hs-CRP <=10 mg/L).
# This threshold was explored in v3-v5 and is NOT preregistered/confirmatory.
# Scientific rationale is chronic/low-grade inflammation; report both cohorts.
# Alternative hs-CRP <=20 mg/L is sensitivity only; CDC bridge of 2015 values
# was applied upstream only within its validated instrument range.
base <- d0[complete.cases(d0[,needed,drop=FALSE]) &
           d0$RIDAGEYR>=20 & d0$crp_harmonized_mg_L>0 &
           d0$crp_harmonized_mg_L<=20,,drop=FALSE]
main <- base[base$crp_harmonized_mg_L <= 10,,drop=FALSE]
if (nrow(main)<500) stop('Primary cohort unexpectedly small; check 2015-2018 input')
if (!all(c(2015,2017) %in% main$cycle_start)) stop('Primary cohort needs two survey cycles')

# One and only one set of population-weighted cutpoints, based on the PRIMARY
# analytic cohort. Held fixed in all sensitivity analyses.
weighted_quantile <- function(x,w,probs=c(.25,.50,.75)) {
  if (length(x)!=length(w) || anyNA(x) || anyNA(w) || any(w<=0)) stop('bad weighted quantile input')
  idx <- order(x,seq_along(x)); xs <- x[idx]; ws <- w[idx]
  cs <- cumsum(ws)/sum(ws)
  vapply(probs,function(p) xs[which(cs>=p)[1]],numeric(1))
}
cutoff_rows <- list()
for (i in seq_along(exposures)) {
  nm <- exposures[i]
  q <- weighted_quantile(main[[nm]],main$survey_weight)
  if (length(unique(q))<3) stop('Less than 3 distinct weighted-quartile cut points: ',nm)
  label <- paste0('q_',tolower(nm))
  # findInterval yields 0..3, with concentration ties kept together.
  d0[[label]] <- findInterval(d0[[nm]],q,rightmost.closed=TRUE)
  cutoff_rows[[i]] <- data.frame(chemical=nm,q25=q[1],q50=q[2],q75=q[3],
                                method='weighted quartiles of primary cohort (locked)')
}
write.csv(do.call(rbind,cutoff_rows),file.path(outdir,'locked_weighted_quartile_cutpoints.csv'),row.names=FALSE)

base <- d0[complete.cases(d0[,needed,drop=FALSE]) &
           d0$RIDAGEYR>=20 & d0$crp_harmonized_mg_L>0 &
           d0$crp_harmonized_mg_L<=20,,drop=FALSE]
main <- base[base$crp_harmonized_mg_L<=10,,drop=FALSE]

# All potential confounders kept in main model (not removed to gain p-value).
# Consider causality/overadjustment of BMI separately; income and race are
# population-level sociodemographic covariates, not biological categories.
prepare <- function(d) {
  d$log_crp <- log(d$crp_harmonized_mg_L)
  d$age <- d$RIDAGEYR
  d$bmi <- d$BMXBMI
  d$pir <- d$INDFMPIR
  d$sex <- factor(d$RIAGENDR)
  d$race_ethnicity <- factor(d$RIDRETH3)
  d$cycle <- factor(d$cycle_start)
  d$log_urine_creatinine <- log(d$URXUCR)
  for (e in exposures) d[[paste0('log2_',tolower(e))]] <- log2(d[[e]])
  d
}
main <- prepare(main)
base <- prepare(base)
conf <- c('age','sex','race_ethnicity','pir','bmi','log_urine_creatinine','cycle')
make_design <- function(d) {
  svydesign(ids=~pooled_psu,strata=~pooled_strata,weights=~survey_weight,
            data=d,nest=TRUE)
}

model_test <- function(fit,design,label,model_type,n) {
  coeff <- coef(fit)
  if (anyNA(coeff)) stop('Aliased coefficient: ',label)
  tdf <- degf(design)
  rdf <- tdf - (fit$rank-1L)
  if (rdf <= 0) stop('Insufficient design residual df for ',label,': df=',rdf)
  # 'svyglm' design-based Wald covariance uses independent PSU-strata units.
  list(design_df=tdf,residual_df=rdf)
}

run_mix <- function(d,selected,group,scenario) {
  des <- make_design(d)
  q <- paste0('q_',tolower(selected))
  fit <- svyglm(reformulate(c(q,conf),response='log_crp'),
                design=des, family=gaussian())
  df <- model_test(fit,des,paste(scenario,group), 'mixture',nrow(d))
  beta <- coef(fit)[q]
  v <- vcov(fit)[q,q,drop=FALSE]
  psi <- sum(beta)
  se <- sqrt(sum(v))
  crit <- qt(.975,df$residual_df)
  p <- if (se > 0) 2*pt(-abs(psi/se),df=df$residual_df) else NA_real_
  data.frame(scenario=scenario,exposure_group=group,n=nrow(d),n_exposures=length(q),
             design_df=df$design_df,residual_df=df$residual_df,
             beta_joint_log_crp=psi,se_joint=se,
             percent_change=100*expm1(psi),
             lower95_pct=100*expm1(psi-crit*se),
             upper95_pct=100*expm1(psi+crit*se),p_value=p,
             adjustment=paste(conf,collapse='+'))
}

mixgroups <- list('Seven exposures'=exposures,
                  'PFAS family'=exposures[1:4],
                  'Phthalate family'=exposures[5:7])
# All 3 groups are predeclared and reported, regardless of p-values.
run_groups <- function(d,scenario) {
  do.call(rbind,lapply(names(mixgroups),function(name)
    run_mix(d,mixgroups[[name]],name,scenario)))
}
primary <- run_groups(main,'Primary: harmonized hs-CRP <=10 mg/L')
primary$p_value_BH_3_tests <- p.adjust(primary$p_value,method='BH')
sensitivity <- run_groups(base,'Sensitivity: harmonized hs-CRP <=20 mg/L')
write.csv(primary,file.path(outdir,'PRIMARY_pooled_mixtures.csv'),row.names=FALSE)
write.csv(sensitivity,file.path(outdir,'sensitivity_pooled_mixtures.csv'),row.names=FALSE)

# Independent single-chemical epidemiological models, log2 exposures; all
# predeclared analytes reported with BH-FDR correction. No p-value filtering.
des <- make_design(main)
single <- lapply(exposures,function(chem) {
  name <- paste0('log2_',tolower(chem))
  fit <- svyglm(reformulate(c(name,conf),response='log_crp'),
                design=des, family=gaussian())
  df <- model_test(fit,des,chem,'single-exposure',nrow(main))
  beta <- coef(fit)[name]; se <- sqrt(vcov(fit)[name,name]); crit <- qt(.975,df$residual_df)
  data.frame(chemical=chem,n=nrow(main),design_df=df$design_df,
             residual_df=df$residual_df,beta=as.numeric(beta),se=as.numeric(se),
             percent_change_per_doubling=100*expm1(beta),
             lower95_pct=100*expm1(beta-crit*se),
             upper95_pct=100*expm1(beta+crit*se),
             p_value=2*pt(-abs(beta/se),df=df$residual_df))
})
single <- do.call(rbind,single)
single$p_value_BH_FDR <- p.adjust(single$p_value,method='BH')
write.csv(single,file.path(outdir,'PRIMARY_pooled_single_chemical_epidemiology.csv'),row.names=FALSE)

# Basic model diagnostics; inspect the small effective design df and
# effective sample size (Kish is descriptive only, not survey SE replacement).
meta <- data.frame(scenario=c('primary <=10','sensitivity <=20'),
                   n=c(nrow(main),nrow(base)),
                   cycle_2015_n=c(sum(main$cycle_start==2015),sum(base$cycle_start==2015)),
                   cycle_2017_n=c(sum(main$cycle_start==2017),sum(base$cycle_start==2017)),
                   design_df=c(degf(make_design(main)),degf(make_design(base))),
                   kish_weight_ess=c(sum(main$survey_weight)^2/sum(main$survey_weight^2),
                     sum(base$survey_weight)^2/sum(base$survey_weight^2)))
write.csv(meta,file.path(outdir,'cohort_and_design_qc.csv'),row.names=FALSE)
writeLines(c(
  'NHANES 2015-2018 pooled (2015-16 and 2017-18) survey-weighted analysis.',
  'Outcome: log-transformed harmonized hs-CRP (mg/L), not cancer.',
  'Primary: hs-CRP <=10 mg/L; sensitivity <=20 mg/L.',
  'CDC 2015->2017 bridge: 0.8695 * hs-CRP + 0.2954 when original <=23 mg/L.',
  '2015 samples outside bridge range are excluded; 2017 values unchanged.',
  'All 3 predeclared mixture families reported; all 7 individual exposures reported.',
  'Multiple comparison correction BH across 3 mixtures and 7 single exposures separately.',
  'Quartiles: weighted pooled PRIMARY sample; cuts fixed in sensitivity cohort.',
  'Pool weight: 2-year environmental subsample B WTSB2YR divided by two.',
  'Survey strata: nested cycle-specific SDMVSTRA; PSU: SDMVPSU.',
  paste('Covariates:',paste(conf,collapse=', ')),
  'Urinary creatinine is covariate-adjusted; NOT both ratio-standardized and adjusted.',
  'Mixture model is linear additive quantile-index; interaction/nonlinear effects not evaluated.',
  'No causal claims; cross-sectional exposure and outcome can be confounded/reverse-associated.',
  'BMI may be an intermediate variable under some DAGs; alternative adjustment requires pre-specification.',
  'This expansion is for design precision, NOT to select significant results.',
  'Do not pool 2017-2020 pre-pandemic files with these 2-year weights.',
  '2015-2016 and 2017-2018 use different hs-CRP instrument methods; see CDC HSCRP_J analytic notes.'
),file.path(outdir,'METHODS_AND_LIMITATIONS.txt'))

cat('\nCOHORT AND DESIGN:\n'); print(meta,row.names=FALSE)
cat('\nPRIMARY POOLED MIXTURES (all 3 hypotheses, BH-adjusted):\n')
print(primary[,c('exposure_group','n','residual_df','percent_change','lower95_pct',
                 'upper95_pct','p_value','p_value_BH_3_tests')],row.names=FALSE)
cat('\nPRIMARY SINGLE-EXPOSURE EPIDEMIOLOGY:\n')
print(single[,c('chemical','percent_change_per_doubling','p_value','p_value_BH_FDR')],row.names=FALSE)
cat('\nSENSITIVITY <=20 mg/L:\n')
print(sensitivity[,c('exposure_group','n','residual_df','percent_change',
                     'lower95_pct','upper95_pct','p_value')],row.names=FALSE)
cat('\n[saved] ',outdir,'\n',sep='')
