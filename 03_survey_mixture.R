#!/usr/bin/env Rscript
# Survey-design-aware epidemiology + linear-additive quantile g-computation.
# This is an association analysis, not proof of effects/causality.
args <- commandArgs(trailingOnly = TRUE)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
root <- if (length(args) && args[1] != "") normalizePath(args[1]) else dirname(normalizePath(script_path))
if (!requireNamespace("survey", quietly=TRUE)) {
  stop("R package 'survey' is missing. On Linux Mint: sudo apt install r-cran-survey")
}
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu="adjust")

proc <- file.path(root, 'data', 'processed')
results <- file.path(root, 'results')
dir.create(results, recursive=TRUE, showWarnings=FALSE)
joint_path <- file.path(proc, 'analysis_joint.csv')
use_joint <- file.exists(joint_path)
input <- if (use_joint) joint_path else file.path(proc, 'analysis_pfas.csv')
d <- read.csv(input, check.names=FALSE)
exposures <- c('LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS')
if (use_joint) exposures <- c(exposures, 'URXMBP','URXMEP','URXMIB')
vars_needed <- c('survey_weight', 'SDMVPSU','SDMVSTRA', 'LBXHSCRP',
                 'RIDAGEYR','RIAGENDR','RIDRETH3','BMXBMI','INDFMPIR', exposures)
if (use_joint) vars_needed <- c(vars_needed, 'URXUCR')
if (!all(vars_needed %in% colnames(d))) stop('Missing columns: ', paste(setdiff(vars_needed, names(d)), collapse=', '))
d <- d[complete.cases(d[,vars_needed]),,drop=FALSE]
d <- d[d$LBXHSCRP>0 & d$survey_weight>0,,drop=FALSE]
for (variable in exposures) d <- d[d[[variable]] > 0,,drop=FALSE]
if (use_joint) d <- d[d$URXUCR>0,,drop=FALSE]
if (nrow(d) < 100) stop('Fewer than 100 complete participants; review QC.')

d$log_hscrp <- log(d$LBXHSCRP)
d$sex <- factor(d$RIAGENDR)
d$race <- factor(d$RIDRETH3)
d$age <- d$RIDAGEYR
d$bmi <- d$BMXBMI
d$pir <- d$INDFMPIR
# Include log urine creatinine as a covariate in models with urinary analytes.
confounders <- c('age','sex','race','bmi','pir')
if (use_joint) {
  d$log_urine_cr <- log(d$URXUCR)
  confounders <- c(confounders, 'log_urine_cr')
}

# Survey PSUs/strata and subsample-B weight, never the full-MEC weight.
des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~survey_weight,
                 data=d, nest=TRUE)
if (degf(des) < 4) stop('Too few survey design degrees of freedom.')

# Epidemiology: adjusted associations per doubling in chemical concentration.
continuous <- paste0('log2_', tolower(exposures))
for (i in seq_along(exposures)) d[[continuous[i]]] <- log2(d[[exposures[i]]])
des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~survey_weight,
                 data=d, nest=TRUE)
form <- reformulate(c(continuous, confounders), response='log_hscrp')
fit <- svyglm(form, design=des, family=gaussian())
# A simultaneous 7-exposure model with many adjustment parameters can consume
# more than the NHANES 15 design df, producing NaN p-values under the default
# residual-df rule. The explicitly chosen *design df* is an exploratory Wald
# approximation documented by survey::summary.svyglm; coefficients/SEs do not
# change. For defensible individual-exposure inference, see the parsimonious
# one-exposure-at-a-time analyses below (default residual df).
model_default_df <- degf(des) - (fit$rank - 1L)
if (model_default_df <= 0) {
  warning('Full joint exposure model has ', fit$rank,
          ' fitted coefficients vs ', degf(des),
          ' design df; full-model p-values are exploratory (using design df).')
}
co <- summary(fit, df.resid=degf(des))$coefficients
epi <- data.frame(term=rownames(co), beta=co[,1], SE=co[,2], statistic=co[,3],
                  p_value=co[,4], row.names=NULL)
epi$percent_change_estimate <- 100 * (exp(epi$beta)-1)
epi$percent_change_lower95 <- 100 * (exp(epi$beta - qt(.975,degf(des))*epi$SE)-1)
epi$percent_change_upper95 <- 100 * (exp(epi$beta + qt(.975,degf(des))*epi$SE)-1)
write.csv(epi, file.path(results,'epidemiology_survey_svyglm.csv'), row.names=FALSE)
# Stable alias clarifies that the joint 7-exposure model uses design-df Wald p-values.
write.csv(epi, file.path(results,'epidemiology_joint_exposure_exploratory.csv'), row.names=FALSE)

# Primary epidemiology analysis: separate regression for each chemical, each
# controlling for the SAME pre-specified potential confounders. This avoids a
# high-dimensional, simultaneous model exhausting the design residual df.
# These coefficients are NOT conditioned on the other six exposures.
single_results <- vector('list', length(exposures))
for (i in seq_along(exposures)) {
  term <- continuous[i]
  f1 <- svyglm(reformulate(c(term, confounders), response='log_hscrp'),
               design=des, family=gaussian())
  denom_df <- degf(des) - (f1$rank - 1L)
  if (!is.finite(denom_df) || denom_df <= 0) {
    stop('Insufficient residual design degrees of freedom for single exposure ',
         exposures[i], ' (df=', denom_df, '). Simplify the model or add cycles.')
  }
  c1 <- summary(f1, df.resid=denom_df)$coefficients
  b <- unname(c1[term, 1]); e <- unname(c1[term, 2])
  tval <- unname(c1[term, 3]); pv <- unname(c1[term, 4])
  crit <- qt(.975, df=denom_df)
  single_results[[i]] <- data.frame(
    chemical=exposures[i], term=term, n=nrow(d), residual_design_df=denom_df,
    beta=b, SE=e, statistic=tval, p_value=pv,
    percent_change_per_doubling=100 * (exp(b)-1),
    percent_change_lower95=100 * (exp(b-crit*e)-1),
    percent_change_upper95=100 * (exp(b+crit*e)-1)
  )
}
singles <- do.call(rbind, single_results)
singles$p_value_FDR_BH <- p.adjust(singles$p_value, method='BH')
write.csv(singles, file.path(results,'epidemiology_single_exposure_adjusted.csv'), row.names=FALSE)

# Mixture modeling: quartile-g-computation, additive linear model (no interactions).
# One-unit simultaneous shift = all exposures increase by ONE quartile bin.
# Survey-based SE of the sum computed from full covariance matrix of svyglm.
qvars <- paste0('q_', tolower(exposures))
for (i in seq_along(exposures)) {
  x <- d[[exposures[i]]]
  breaks <- unique(as.numeric(quantile(x, probs=seq(0,1,0.25), na.rm=TRUE)))
  if (length(breaks) < 4) stop('Too few distinct exposure quantiles: ', exposures[i])
  # cut() handles true tied concentrations without arbitrarily splitting ties.
  d[[qvars[i]]] <- as.integer(cut(x, breaks=breaks, include.lowest=TRUE, labels=FALSE)) - 1L
}
des_q <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~survey_weight,
                   data=d, nest=TRUE)
fit_q <- svyglm(reformulate(c(qvars, confounders), response='log_hscrp'),
                design=des_q, family=gaussian())
beta <- coef(fit_q)[qvars]
if (anyNA(beta)) stop('Exposure quantile coefficient missing; investigate collinearity.')
v <- vcov(fit_q)[qvars,qvars,drop=FALSE]
psi <- sum(beta)
se <- sqrt(sum(v))
tdf <- degf(des_q)
critical <- qt(.975, tdf)
mixture <- data.frame(
  analysis=if (use_joint) 'PFAS + urine phthalates' else 'PFAS only',
  n_participants=nrow(d), n_exposures=length(exposures),
  survey_degrees_freedom=tdf,
  psi_log_hscrp=psi, standard_error=se,
  lower95=psi-critical*se, upper95=psi+critical*se,
  p_value=2*pt(-abs(psi/se), df=tdf),
  percent_difference=100*(exp(psi)-1),
  percent_difference_lower95=100*(exp(psi-critical*se)-1),
  percent_difference_upper95=100*(exp(psi+critical*se)-1)
)
write.csv(mixture, file.path(results,'quantile_gcomputation_joint_effect.csv'), row.names=FALSE)
positive <- pmax(beta, 0)
negative <- pmax(-beta, 0)
weights <- data.frame(chemical=exposures, coefficient=as.numeric(beta),
  positive_direction_weight=if (sum(positive)>0) as.numeric(positive/sum(positive)) else NA_real_,
  negative_direction_weight=if (sum(negative)>0) as.numeric(negative/sum(negative)) else NA_real_)
write.csv(weights, file.path(results,'mixture_component_weights.csv'), row.names=FALSE)
writeLines(c(
  paste('Input:', basename(input)),
  paste('Complete cases:', nrow(d)),
  paste('Exposures:', paste(exposures, collapse=', ')),
  paste('Design df:', tdf),
  'Design: strata SDMVSTRA, PSU SDMVPSU, subsample B weight WTSB2YR',
  'Outcome: log hs-CRP, not cancer',
  'Covariates: age, sex, race/ethnicity, BMI, income/poverty ratio,',
  'log urine creatinine (joint model only)',
  'Model: additive quantile g-computation via survey-weighted svyglm',
  'Quantile bins are derived from the unweighted analytic sample, with tied values kept together.',
  'Mixture psi is the sum of quartile exposure coefficients; SE from covariance.',
  'Conditional associations, not causal effects. No nonlinear interactions modeled.',
  'Missingness handled by complete cases; evaluate nonresponse before population inference.',
  paste('Simultaneous model residual design df (default):', model_default_df),
  'Joint 7-exposure Wald p-values use survey::summary(..., df.resid=degf(des)) as an exploratory approximation.',
  'Primary epidemiology file uses one chemical at a time, with same confounders, default residual df.',
  'Multiple single-exposure comparisons are adjusted using Benjamini-Hochberg FDR.',
  'Mixture p-value uses design t df; mixture inference exploratory given only 15 PSU-stratum df.'
), file.path(results,'statistical_analysis_notes.txt'))
cat('\nEPIDEMIOLOGY (survey weighted):\n')
print(epi[epi$term %in% continuous,,drop=FALSE])
cat('\nSINGLE-CHEMICAL EPIDEMIOLOGY (covariate-adjusted; residual survey df):\n')
print(singles[,c('chemical','residual_design_df','percent_change_per_doubling','p_value','p_value_FDR_BH'),drop=FALSE])
cat('\nMIXTURE (additive quantile g-computation, survey weighted):\n')
print(mixture)
cat('\nOutputs:', results, '\n')
