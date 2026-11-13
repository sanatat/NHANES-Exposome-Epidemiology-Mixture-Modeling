#!/usr/bin/env Rscript
# NHANES 2017-18 CUE-CAN pilot: detection-limit QC + exploratory survey sensitivity.
# No causal claims, no cancer outcomes. Needs the pre-existing data/processed/analysis_joint.csv.
command <- grep('^--file=', commandArgs(trailingOnly=FALSE), value=TRUE)
if (!length(command)) stop('Run using Rscript 04_sensitivity_qc.R')
root <- dirname(normalizePath(sub('^--file=', '', command[1])))
if (!requireNamespace('survey', quietly=TRUE)) stop('Install survey: sudo apt install r-cran-survey')
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust')
file <- file.path(root, 'data', 'processed', 'analysis_joint.csv')
if (!file.exists(file)) stop('Missing analysis_joint.csv. First run 01_download_build.py')
outdir <- file.path(root, 'results')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
d <- read.csv(file, check.names=FALSE)
exposures <- c('LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB')
flags <- c('LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL','URDMBPLC','URDMEPLC','URDMIBLC')
needed <- c('SEQN','survey_weight','SDMVSTRA','SDMVPSU','LBXHSCRP','RIDAGEYR','RIAGENDR',
            'RIDRETH3','BMXBMI','INDFMPIR','URXUCR',exposures,flags)
if (!all(needed %in% names(d))) stop('Missing columns: ',paste(setdiff(needed,names(d)),collapse=', '))
d <- d[complete.cases(d[,needed]),,drop=FALSE]
d <- d[is.finite(d$LBXHSCRP) & d$LBXHSCRP>0 & is.finite(d$URXUCR) & d$URXUCR>0 &
       is.finite(d$survey_weight) & d$survey_weight>0,,drop=FALSE]
for (x in exposures) d <- d[is.finite(d[[x]]) & d[[x]]>0,,drop=FALSE]
if (nrow(d)<100) stop('Unexpectedly small analytic cohort.')
if (anyDuplicated(d$SEQN)) stop('Duplicate SEQN in joint analysis.')

d$log_hscrp <- log(d$LBXHSCRP)
d$sex <- factor(d$RIAGENDR)
d$race <- factor(d$RIDRETH3)
d$age <- d$RIDAGEYR
d$bmi <- d$BMXBMI
d$pir <- d$INDFMPIR
d$log_urine_cr <- log(d$URXUCR)
confounders <- c('age','sex','race','bmi','pir','log_urine_cr')
make_design <- function(data) svydesign(ids=~SDMVPSU, strata=~SDMVSTRA,
                                       weights=~survey_weight, data=data, nest=TRUE)
base_design <- make_design(d)

# LOD information from CDC lab quality indicators; handle sample weights.
qc <- lapply(seq_along(exposures), function(i) {
  f <- flags[i]
  x <- d[[f]]
  ok <- !is.na(x) & x %in% c(0,1)
  data.frame(chemical=exposures[i], lod_flag=f,
             n_complete=nrow(d),n_with_flag=sum(ok),
             n_below_LOD=sum(x[ok]==1),
             pct_below_LOD_unweighted=if(sum(ok)>0) 100*mean(x[ok]==1) else NA_real_,
             pct_below_LOD_weighted=if(sum(ok)>0) 100*weighted.mean(x[ok]==1,d$survey_weight[ok]) else NA_real_)
})
qc <- do.call(rbind, qc)
write.csv(qc, file.path(outdir,'sensitivity_lod_qc.csv'),row.names=FALSE)

# Quantile boundaries are fixed from the original cohort so excluding high CRP
# does NOT redefine the units of the exposure shifts.
qcols <- paste0('q_',tolower(exposures))
for (i in seq_along(exposures)) {
  z <- d[[exposures[i]]]
  breaks <- unique(as.numeric(quantile(z,probs=seq(0,1,0.25),na.rm=TRUE)))
  if(length(breaks)<4) stop('Insufficient unique quantile boundaries: ',exposures[i])
  d[[qcols[i]]] <- as.integer(cut(z, breaks=breaks, include.lowest=TRUE, labels=FALSE)) - 1L
  if(anyNA(d[[qcols[i]]])) stop('Quantile assignment failure for ',exposures[i])
}

run_mix <- function(data, chems, label, scenario) {
  des <- make_design(data)
  qterms <- paste0('q_', tolower(chems))
  form <- reformulate(c(qterms,confounders),response='log_hscrp')
  fit <- svyglm(form,design=des,family=gaussian())
  b <- coef(fit)[qterms]
  if(anyNA(b)) stop('Aliased mixture exposure in ', label, ' / ',scenario)
  covariance <- vcov(fit)[qterms,qterms,drop=FALSE]
  psi <- sum(b)
  se <- sqrt(sum(covariance))
  df_design <- degf(des)
  df_residual <- df_design - (fit$rank - 1L)
  critical <- qt(0.975,df_design)
  delta <- 100*(exp(psi)-1)
  data.frame(scenario=scenario, mixture=label, n=nrow(data),
             n_exposures=length(chems), design_df=df_design,
             model_rank=fit$rank, residual_df=df_residual,
             psi=psi, SE=se,
             pct_difference=delta,
             lower95_pct=100*(exp(psi-critical*se)-1),
             upper95_pct=100*(exp(psi+critical*se)-1),
             p_exploratory=2*pt(-abs(psi/se),df=df_design),
             inference='Exploratory design-df Wald; not causal')
}

models <- list(
  list(label='All 7 chemicals', chems=exposures),
  list(label='PFAS only, same people', chems=exposures[1:4]),
  list(label='Phthalates only, same people', chems=exposures[5:7])
)
scenarios <- list('Full complete-case cohort'=d,'hs-CRP <= 10 mg/L'=d[d$LBXHSCRP<=10,,drop=FALSE])
results <- list()
for (sc in names(scenarios)) {
  ds <- scenarios[[sc]]
  if(nrow(ds)<100) stop('Too few individuals for sensitivity ',sc)
  for (m in models) {
    results[[length(results)+1]] <- run_mix(ds,m$chems,m$label,sc)
  }
}
result <- do.call(rbind,results)
write.csv(result,file.path(outdir,'sensitivity_mixture_comparison.csv'),row.names=FALSE)

summary <- c(
  paste('Date:',Sys.Date()),
  paste('Participants:',nrow(d)),
  paste('High hs-CRP (>10 mg/L):',sum(d$LBXHSCRP>10)),
  paste('Weighted high hs-CRP fraction:',round(weighted.mean(d$LBXHSCRP>10,d$survey_weight)*100,2),'%'),
  paste('Base design PSU minus strata:',degf(base_design)),
  'Subsample B weight; observations are people sampled using a complex survey.',
  'Quantile thresholds fixed using full joint analytic cohort.',
  'Serum PFAS and urinary phthalates are separate matrices; urinary creatinine adjusted.',
  'Same covariates and same base people across mixtures for a fair exploratory comparison.',
  'Adjusted for age, sex, race/ethnicity, BMI, PIR and log urine creatinine.',
  'All p-values in sensitivity file use design df Wald approximation: exploratory.',
  'No causal interpretation; high-CRP restriction is a sensitivity analysis, not default.'
)
writeLines(summary,file.path(outdir,'sensitivity_notes.txt'))
cat('\nQUALITY CONTROL: percentage below limit of detection\n')
print(qc[,c('chemical','n_below_LOD','pct_below_LOD_weighted')],row.names=FALSE)
cat('\nSENSITIVITY COMPARISON\n')
print(result[,c('scenario','mixture','n','residual_df','pct_difference','lower95_pct','upper95_pct','p_exploratory')],row.names=FALSE)
cat('\nSaved in',outdir,'\n')
