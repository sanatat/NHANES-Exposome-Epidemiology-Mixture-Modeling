#!/usr/bin/env Rscript
# NHANES 2015-2018 / v8: post-hoc creatinine adjustment, robustness, and
# restricted cubic spline dose-response for urinary mono-n-butyl phthalate (MBP).
# Does NOT replace the v6 primary models or provide external validation.
# Do not select specifications by p-value; report ALL models below.

cmd <- grep('^--file=', commandArgs(trailingOnly=FALSE), value=TRUE)
if (length(cmd) != 1L) stop('Run: Rscript 12_creatinine_robustness_splines.R')
root <- dirname(normalizePath(sub('^--file=', '', cmd[[1L]])))
if (!requireNamespace('survey', quietly=TRUE)) stop('R survey package missing; sudo apt install r-cran-survey')
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust')

input <- file.path(root,'data','processed','analysis_joint_pooled_2015_2018.csv')
output <- file.path(root,'results','pooled_2015_2018')
dir.create(output,recursive=TRUE,showWarnings=FALSE)
if (!file.exists(input)) stop('Missing pooled input. First run: .venv/bin/python 08_build_pooled_2015_2018.py')

d0 <- read.csv(input, check.names=FALSE)
all_exposures <- c('LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB')
all_flags <- c('LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL',
               'URDMBPLC','URDMEPLC','URDMIBLC')
needed <- c('SEQN','cycle_start','survey_weight','pooled_strata','pooled_psu',
  'crp_harmonized_mg_L','RIDAGEYR','RIAGENDR','RIDRETH3',
  'INDFMPIR','BMXBMI','URXUCR',all_exposures,all_flags)
if (!all(needed %in% names(d0))) stop('Missing required columns: ',paste(setdiff(needed,names(d0)),collapse=', '))
if (anyDuplicated(d0$SEQN)) stop('Duplicate SEQN in pooled file')
if (!setequal(unique(d0$cycle_start),c(2015,2017))) stop('Unexpected pooled cycles')
if (any(!is.finite(d0$survey_weight) | d0$survey_weight <= 0)) stop('Invalid sampling weights')
for (flag in all_flags) {
  if (any(!is.na(d0[[flag]]) & !(d0[[flag]] %in% c(0,1)))) stop('Invalid LOD flag: ',flag)
}

# Preserve the EXACT v6 analytic sample definition (avoid changing sample
# composition across creatinine or covariate specifications).
base <- d0[complete.cases(d0[,needed,drop=FALSE]) &
  d0$RIDAGEYR >= 20 & d0$URXUCR > 0 & d0$URXMBP > 0 &
  d0$crp_harmonized_mg_L > 0 & d0$crp_harmonized_mg_L <= 20,,drop=FALSE]
main <- base[base$crp_harmonized_mg_L <= 10,,drop=FALSE]
if (nrow(main)<1500 || nrow(base)<1500) stop('Unexpectedly small pooled analytic cohort')
if (!all(c(2015,2017) %in% main$cycle_start)) stop('Both cycles required')

prepare <- function(d) {
  d$log_crp <- log(d$crp_harmonized_mg_L)
  d$log2_mbp <- log2(d$URXMBP)
  d$age <- d$RIDAGEYR
  d$bmi <- d$BMXBMI
  d$pir <- d$INDFMPIR
  d$sex <- factor(d$RIAGENDR)
  d$race_ethnicity <- factor(d$RIDRETH3)
  d$cycle <- factor(d$cycle_start)
  d$log_urine_creatinine <- log(d$URXUCR)
  # URXMBP in ng/mL, URXUCR in mg/dL:
  # ng/g creatinine = (ng/mL)*100000/(mg/dL).
  d$mbp_ng_per_g_creatinine <- 100000 * d$URXMBP / d$URXUCR
  d$log2_mbp_per_g_creatinine <- log2(d$mbp_ng_per_g_creatinine)
  if (any(!is.finite(d$log2_mbp_per_g_creatinine))) stop('Bad creatinine-standardized MBP')
  d
}
base <- prepare(base)
main <- prepare(main)

wquant <- function(x,w,ps) {
  if (length(x)!=length(w) || anyNA(x) || anyNA(w) || any(w <= 0)) stop('Bad weighted-quantile inputs')
  j <- order(x,seq_along(x))
  xs <- x[j]; ws <- w[j]; cs <- cumsum(ws)/sum(ws)
  vapply(ps,function(p) xs[which(cs>=p)[[1L]]],numeric(1))
}

# Fixed in primary cohort and reused for sensitivity. Exposure-only choices;
# never use outcome associations or p-values to choose winsor limits/knots.
mbp_winsor <- wquant(main$URXMBP,main$survey_weight,c(.01,.99))
knots <- wquant(main$log2_mbp,main$survey_weight,c(.05,.35,.65,.95))
if (length(unique(knots)) != 4L) stop('Restricted cubic spline knots must be distinct')
if (mbp_winsor[[1L]] >= mbp_winsor[[2L]]) stop('Degenerate winsor limits')

# Four-knot restricted cubic spline (natural cubic spline):
# x + two truncated-cubic terms. The extra two terms represent NONLINEARITY.
# Scaled by (k4-k1)^2 for numerical conditioning.
rcs4 <- function(x,k) {
  if (length(k)!=4 || any(diff(k)<=0)) stop('Invalid RCS knots')
  tp <- function(v) pmax(v,0)^3
  span <- (k[[4L]]-k[[1L]])^2
  make <- function(j) (
    tp(x-k[[j]]) -
    tp(x-k[[3L]])*(k[[4L]]-k[[j]])/(k[[4L]]-k[[3L]]) +
    tp(x-k[[4L]])*(k[[3L]]-k[[j]])/(k[[4L]]-k[[3L]])
    ) / span
  data.frame(mbp_rcs1=make(1L),mbp_rcs2=make(2L))
}
add_locked_terms <- function(d) {
  d$log2_mbp_winsor <- log2(pmax(mbp_winsor[[1L]],pmin(d$URXMBP,mbp_winsor[[2L]])))
  b <- rcs4(d$log2_mbp,knots)
  d$mbp_rcs1 <- b$mbp_rcs1
  d$mbp_rcs2 <- b$mbp_rcs2
  d
}
main <- add_locked_terms(main)
base <- add_locked_terms(base)

write.csv(data.frame(parameter=c('mbp_winsor_q01_ng_mL','mbp_winsor_q99_ng_mL',
      'spline_knot_05_log2_ng_mL','spline_knot_35_log2_ng_mL',
      'spline_knot_65_log2_ng_mL','spline_knot_95_log2_ng_mL'),
      value=c(mbp_winsor,knots)),
      file.path(output,'v8_locked_exposure_cutoffs.csv'),row.names=FALSE)

conf_full <- c('age','sex','race_ethnicity','pir','bmi','cycle','log_urine_creatinine')
conf_no_cr <- setdiff(conf_full,'log_urine_creatinine')
conf_no_bmi <- setdiff(conf_full,'bmi')
conf_no_pir <- setdiff(conf_full,'pir')
conf_core <- c('age','sex','race_ethnicity','cycle','log_urine_creatinine')
make_design <- function(d) svydesign(ids=~pooled_psu,strata=~pooled_strata,
                                     weights=~survey_weight,data=d,nest=TRUE)

fit_effect <- function(d,exposure,adjustment,scenario,kind) {
  des <- make_design(d)
  form <- reformulate(c(exposure,adjustment),response='log_crp')
  fit <- svyglm(form,design=des,family=gaussian())
  cf <- coef(fit)
  if (anyNA(cf) || !(exposure %in% names(cf))) stop('Aliased model: ',scenario)
  rdf <- degf(des) - (fit$rank-1L)
  if (rdf <= 0) stop('Nonpositive design residual df: ',scenario,' / ',rdf)
  beta <- as.numeric(cf[[exposure]])
  se <- as.numeric(sqrt(vcov(fit)[exposure,exposure]))
  if (!is.finite(se) || se<=0) stop('Invalid uncertainty in: ',scenario)
  crit <- qt(.975,rdf)
  data.frame(scenario=scenario,model_kind=kind,n=nrow(d),design_df=degf(des),
    residual_df=rdf,exposure_variable=exposure,
    urine_dilution_handling=if ('log_urine_creatinine' %in% adjustment)
      'log_creatinine_covariate' else if (exposure=='log2_mbp_per_g_creatinine')
      'creatinine_standardized_ratio' else 'none',
    adjustment=paste(adjustment,collapse='+'),
    beta_log_crp=beta,se_log_crp=se,percent_change_per_doubling=100*expm1(beta),
    lower95_pct=100*expm1(beta-crit*se),
    upper95_pct=100*expm1(beta+crit*se),
    p_value_exploratory=2*pt(-abs(beta/se),df=rdf),stringsAsFactors=FALSE)
}

specs <- list(
  list(id='Reference v6: raw MBP + log urinary creatinine',x='log2_mbp',adj=conf_full),
  list(id='MBP per g creatinine (no extra creatinine covariate)',x='log2_mbp_per_g_creatinine',adj=conf_no_cr),
  list(id='Raw MBP without creatinine dilution correction',x='log2_mbp',adj=conf_no_cr),
  list(id='Reference minus BMI',x='log2_mbp',adj=conf_no_bmi),
  list(id='Reference minus income PIR',x='log2_mbp',adj=conf_no_pir),
  list(id='Demographic core + urine creatinine',x='log2_mbp',adj=conf_core),
  list(id='Winsorized MBP at locked weighted 1-99 pct',x='log2_mbp_winsor',adj=conf_full)
)
results <- list()
for (s in specs) {
  results[[length(results)+1L]] <- fit_effect(main,s$x,s$adj,
    paste0('Primary <=10: ',s$id),'primary sensitivity')
}
# Re-run the same seven specifications on an alternative outcome threshold;
# no refitting cutoffs or sample exclusions by p-value.
for (s in specs) {
  results[[length(results)+1L]] <- fit_effect(base,s$x,s$adj,
    paste0('Sensitivity <=20: ',s$id),'alternative hs-CRP threshold')
}
robust <- do.call(rbind,results)
write.csv(robust,file.path(output,'v8_MBP_robustness_all_specifications.csv'),row.names=FALSE)

# Audit reference-model agreement with v6 (do not silently overwrite v6).
v6fn <- file.path(output,'PRIMARY_pooled_single_chemical_epidemiology.csv')
if (file.exists(v6fn)) {
  v6 <- read.csv(v6fn)
  prev <- v6[v6$chemical=='URXMBP','percent_change_per_doubling']
  if (length(prev) != 1L) stop('Cannot compare to v6 MBP row')
  delta <- abs(prev - robust$percent_change_per_doubling[[1L]])
  if (delta > 1e-5) stop(sprintf('v8 reference mismatch vs v6: %.6f percentage points. Investigate!',delta))
  cat(sprintf('[QC] Reference MBP estimate matches v6: %.6f%%\n',prev))
} else {
  warning('v6 result CSV absent: reference model replication cannot be checked')
}

# RCS: compare full adjustment reference linear model to nested linear + 2
# nonlinear components (4 weighted knots locked from PRIMARY cohort).
design <- make_design(main)
linear_form <- reformulate(c('log2_mbp',conf_full),response='log_crp')
rcs_form <- reformulate(c('log2_mbp','mbp_rcs1','mbp_rcs2',conf_full),response='log_crp')
linear_fit <- svyglm(linear_form,design=design,family=gaussian())
spline_fit <- svyglm(rcs_form,design=design,family=gaussian())
if (anyNA(coef(spline_fit)) || anyNA(coef(linear_fit))) stop('Aliased spline terms')
rdf <- degf(design)-(spline_fit$rank-1L)
if (rdf <= 0) stop('Insufficient residual design df for RCS nonlinearity test')
nonlinear_test <- regTermTest(spline_fit,~mbp_rcs1+mbp_rcs2,df=rdf,method='Wald')
p_nl <- as.numeric(nonlinear_test$p)
if (length(p_nl)!=1L || !is.finite(p_nl)) stop('Nonlinearity p-value was not computed')

# Descriptive adjusted dose-response contrasts relative to the weighted
# median MBP. These are NOT direct predictions of absolute hs-CRP levels.
xgrid <- seq(wquant(main$log2_mbp,main$survey_weight,.01),
             wquant(main$log2_mbp,main$survey_weight,.99),length.out=120L)
xref <- wquant(main$log2_mbp,main$survey_weight,.50)[[1L]]
z <- rcs4(xgrid,knots)
z0 <- rcs4(xref,knots)
contrast <- cbind(log2_mbp=xgrid-xref,
  mbp_rcs1=z$mbp_rcs1-z0$mbp_rcs1,
  mbp_rcs2=z$mbp_rcs2-z0$mbp_rcs2)
terms <- colnames(contrast)
if (!all(terms %in% names(coef(spline_fit)))) stop('Spline terms missing from fit')
cov3 <- vcov(spline_fit)[terms,terms,drop=FALSE]
beta3 <- coef(spline_fit)[terms]
delta <- as.numeric(as.matrix(contrast) %*% beta3)
se <- sqrt(pmax(0,rowSums((as.matrix(contrast) %*% cov3)*as.matrix(contrast))))
crit <- qt(.975,rdf)
curve <- data.frame(mbp_ng_mL=2^xgrid,log2_mbp=xgrid,
  reference_mbp_ng_mL=2^xref,
  relative_hscrp_pct=100*expm1(delta),
  lower95_pct=100*expm1(delta-crit*se),
  upper95_pct=100*expm1(delta+crit*se))
write.csv(curve,file.path(output,'v8_MBP_RCS_curve_exploratory.csv'),row.names=FALSE)
write.csv(data.frame(n=nrow(main),design_df=degf(design),residual_df=rdf,
  spline_knots=4,nonlinear_parameters=2,
  p_nonlinearity_exploratory=p_nl,
  mbp_reference_ng_mL=2^xref,
  interpretation='POST-HOC: Wald test of additional nonlinear spline terms; NOT independent confirmation'),
  file.path(output,'v8_MBP_RCS_nonlinearity_test.csv'),row.names=FALSE)

# Safe plotting; the tabular CSV still exists even if plotting device fails.
png(file.path(output,'v8_MBP_exposure_response_spline.png'),width=1050,height=680,res=135)
par(mar=c(4.6,5,3.2,1.5))
xx <- curve$mbp_ng_mL
ylim <- range(curve$lower95_pct,curve$upper95_pct,finite=TRUE)
plot(xx,curve$relative_hscrp_pct,type='n',log='x',ylim=ylim,
     xlab='Urinary MBP (ng/mL; logarithmic scale)',
     ylab='Relative difference in hs-CRP (%)',
     main='Exploratory survey-weighted MBP spline; reference = weighted median')
polygon(c(xx,rev(xx)),c(curve$lower95_pct,rev(curve$upper95_pct)),
        col=grDevices::adjustcolor('gray65',alpha.f=.45),border=NA)
lines(xx,curve$relative_hscrp_pct,lwd=2,col='steelblue4')
abline(h=0,lty=2,col='gray40')
abline(v=2^xref,lty=3,col='gray40')
legend('topleft',legend=c('Adjusted contrast','95% design-based CI'),
       col=c('steelblue4','gray65'),lty=c(1,1),lwd=c(2,8),bty='n')
dev.off()

writeLines(c(
  'NHANES v8: POST-HOC/EXPLORATORY robustness, dilution adjustment, spline analysis.',
  'v6 is still the primary model; v7 cycle analysis is post-hoc, not independent validation.',
  sprintf('Participants primary hs-CRP <=10: %d; sensitivity <=20: %d.',nrow(main),nrow(base)),
  'Outcome = harmonized log hs-CRP. No cancer outcome; no causal claim.',
  'Two 2015-2018 NHANES cycles; preserves pooled public strata/PSU and WTSB2YR/2.',
  'Single-chemical models all use same complete-case cohort within each threshold.',
  'Urine dilution alternatives are DIFFERENT estimands: concentration+creatinine covariate,',
  'creatinine-normalized exposure with NO second creatinine covariate, and unadjusted concentration.',
  'Urinary MBP original units ng/mL; urine creatinine mg/dL; MBP ratio ng/g = 100000*MBP/creatinine.',
  'Covariates: age, sex, population group, PIR, BMI, cycle (and log creatinine where specified).',
  'Additional covariate variants omit BMI, PIR, or both as sensitivity analyses, not preferred models.',
  sprintf('Winsorization thresholds prelocked on primary weighted exposure: %.4g and %.4g ng/mL.',
          mbp_winsor[[1L]],mbp_winsor[[2L]]),
  'Restricted cubic spline has four primary exposure-weighted knots (5/35/65/95 percentiles).',
  'Test of nonlinearity = design-based Wald joint test of TWO extra nonlinear spline terms.',
  'Spline curves are % contrasts relative to weighted median MBP, not absolute predicted hs-CRP.',
  'Do NOT select a significant result among 14 alternative model specifications.',
  'All p-values in v8 are exploratory and unadjusted for post-hoc search; report effect sizes and CIs.',
  'Urine creatinine itself can vary by sex, BMI, age and kidney health; ratio standardization may bias.',
  'Cross-sectional spot-urine biomarkers and hs-CRP cannot establish temporal order or causation.',
  'Changing thresholds and winsorization can cause selection effects; interpret cautiously.',
  'CDC official: https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx',
  'CDC urine creatinine units: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/ALB_CR_I.htm',
  'R survey Wald term test: https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/regTermTest.html'
),file.path(output,'v8_methods_and_limits.txt'))

cat('\nV8 SAMPLE QC:\n')
print(data.frame(scenario=c('primary <=10','sensitivity <=20'),n=c(nrow(main),nrow(base)),
  design_df=c(degf(make_design(main)),degf(make_design(base)))),row.names=FALSE)
cat('\nV8 MBP ROBUSTNESS (all post-hoc; no significance selection):\n')
print(robust[,c('scenario','n','residual_df','percent_change_per_doubling',
      'lower95_pct','upper95_pct','p_value_exploratory')],row.names=FALSE)
cat('\nV8 RESTRICTED CUBIC SPLINE:\n')
print(data.frame(n=nrow(main),residual_df=rdf,
  p_nonlinearity_exploratory=p_nl,reference_mbp_ng_mL=2^xref),row.names=FALSE)
cat('\n[saved] ',output,'\n',sep='')
