#!/usr/bin/env Rscript
# NHANES 2015-2018: exploratory within-cycle consistency of MBP association.
# This is NOT independent replication: the MBP signal was selected using pooled data.
# The established v6 primary three mixture tests and seven single-exposure tests
# remain PRIMARY. This script cannot be used to redefine them or cherry-pick p-values.
args <- grep('^--file=', commandArgs(trailingOnly=FALSE), value=TRUE)
if (length(args) != 1L) stop('Run: Rscript 11_cycle_consistency.R')
root <- dirname(normalizePath(sub('^--file=', '', args[[1L]])))
if (!requireNamespace('survey',quietly=TRUE)) stop('Install R package survey: sudo apt install r-cran-survey')
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust')
fn <- file.path(root,'data','processed','analysis_joint_pooled_2015_2018.csv')
if (!file.exists(fn)) stop('Pooled CSV absent. First run: .venv/bin/python 08_build_pooled_2015_2018.py')
out <- file.path(root,'results','pooled_2015_2018')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(fn,check.names=FALSE)
chem <- c('LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB')
need <- c('SEQN','cycle_start','survey_weight','pooled_strata','pooled_psu','crp_harmonized_mg_L',
  'RIDAGEYR','RIAGENDR','RIDRETH3','INDFMPIR','BMXBMI','URXUCR',chem)
if (!all(need %in% names(d))) stop('Missing columns: ',paste(setdiff(need,names(d)),collapse=', '))
if (anyDuplicated(d$SEQN)) stop('Duplicate participants in pooled dataset')
if (!setequal(unique(d$cycle_start), c(2015,2017))) stop('Expected 2015 and 2017 cycle starts')
if (any(!is.finite(d$survey_weight) | d$survey_weight <= 0)) stop('Invalid survey weights')
if (any(!is.finite(d$URXUCR) | d$URXUCR <= 0)) stop('Invalid urine creatinine')
# Match the EXACT v6 primary cohort and confounder set.
d <- d[complete.cases(d[,need,drop=FALSE]) & d$RIDAGEYR>=20 &
  d$crp_harmonized_mg_L>0 & d$crp_harmonized_mg_L<=10,,drop=FALSE]
if (nrow(d) < 500) stop('Unexpectedly small primary cohort')
d$log_crp <- log(d$crp_harmonized_mg_L)
d$age <- d$RIDAGEYR
d$bmi <- d$BMXBMI
d$pir <- d$INDFMPIR
d$sex <- factor(d$RIAGENDR)
d$race_ethnicity <- factor(d$RIDRETH3)
d$cycle <- factor(d$cycle_start)
d$log_urine_creatinine <- log(d$URXUCR)
for (nm in chem) d[[paste0('log2_',tolower(nm))]] <- log2(d[[nm]])
conf_base <- c('age','sex','race_ethnicity','pir','bmi','log_urine_creatinine')
design <- function(x) svydesign(ids=~pooled_psu,strata=~pooled_strata,
                                 weights=~survey_weight,data=x,nest=TRUE)
get_effect <- function(fit, des, term, lab, n) {
  cf <- coef(fit)
  if (!(term %in% names(cf)) || !is.finite(cf[[term]])) stop('Coefficient missing/aliased: ', lab)
  df <- degf(des) - (fit$rank - 1L)
  if (df<=0) stop('No residual design df for ',lab,' (df=',df,')')
  beta <- as.numeric(cf[[term]])
  se <- as.numeric(sqrt(vcov(fit)[term,term]))
  if (!is.finite(se) || se<=0) stop('Invalid standard error: ',lab)
  crit <- qt(.975,df)
  data.frame(cycle=lab,n=n,design_df=degf(des),residual_df=df,
    percent_change_per_doubling=100*expm1(beta),
    lower95_pct=100*expm1(beta-crit*se),
    upper95_pct=100*expm1(beta+crit*se),
    beta=beta,se=se,p_value=2*pt(-abs(beta/se),df=df))
}
# Within-cycle models use identical covariates except cycle, which is constant.
# Fit ALL seven chemicals to avoid hiding less favorable within-cycle results.
cycle_results <- list()
for (cy in c(2015,2017)) {
  sub <- d[d$cycle_start==cy,,drop=FALSE]
  des <- design(sub)
  lab <- paste0(cy,'-',cy+1)
  for (nm in chem) {
    trm <- paste0('log2_',tolower(nm))
    f <- svyglm(reformulate(c(trm,conf_base),response='log_crp'),
                design=des,family=gaussian())
    row <- get_effect(f,des,trm,lab,nrow(sub))
    row$chemical <- nm
    cycle_results[[length(cycle_results)+1L]] <- row
  }
}
cyc <- do.call(rbind,cycle_results)
cyc$p_value_BH_within_cycle <- ave(cyc$p_value,cyc$cycle,
                                     FUN=function(x) p.adjust(x,method='BH'))
cyc <- cyc[,c('cycle','chemical','n','design_df','residual_df',
  'percent_change_per_doubling','lower95_pct','upper95_pct',
  'beta','se','p_value','p_value_BH_within_cycle')]
write.csv(cyc,file.path(out,'v7_exploratory_cycle_exposures.csv'),row.names=FALSE)

# Test for an MBP-by-cycle difference in the pooled cohort; this is an
# exploratory heterogeneity test selected AFTER observing the pooled MBP result.
des <- design(d)
f <- svyglm(log_crp ~ log2_urxmbp * cycle + age + sex + race_ethnicity +
              pir + bmi + log_urine_creatinine,
            design=des,family=gaussian())
term <- grep('log2_urxmbp:cycle',names(coef(f)),value=TRUE,fixed=TRUE)
if (length(term)!=1L) stop('Expected a single MBP-by-cycle interaction coefficient')
int <- get_effect(f,des,term,'2017-vs-2015 MBP slope interaction',nrow(d))
int$interpretation <- 'Exploratory coefficient contrast; not a new primary finding'
write.csv(int,file.path(out,'v7_exploratory_MBP_cycle_interaction.csv'),row.names=FALSE)

writeLines(c(
 'Exploratory, POST-HOC consistency analysis of the already-seen pooled MBP association.',
 'The two cycles are subsets of discovery data, NOT independent external validation.',
 'Main result remains v6 PRIMARY_pooled_mixtures.csv and PRIMARY_pooled_single_chemical_epidemiology.csv.',
 'No selection of outcome, exposures, covariates, or significance thresholds based on these new p-values.',
 'hs-CRP outcome: harmonized, >0 and <=10 mg/L, log-transformed.',
 'All 7 chemicals tested separately for each cycle with BH correction within each cycle.',
 'Confounders: age, sex, race/ethnicity, PIR, BMI, log urine creatinine.',
 'Within-cycle models omit cycle covariate because it is constant.',
 'Survey strata and PSUs and four-year scaled weights preserved from v6.',
 'Caution: 2015 hs-CRP calibration outside validated range was excluded upstream.',
 'An interaction p-value >0.05 is not proof of equivalence across cycles.',
 'Cross-sectional observational data do not establish causality.'
),file.path(out,'v7_exploratory_replication_notes.txt'))
cat('\nEXPLORATORY WITHIN-CYCLE MBP (selection occurred in pooled data; not independent validation):\n')
print(cyc[cyc$chemical=='URXMBP',c('cycle','n','residual_df','percent_change_per_doubling',
          'lower95_pct','upper95_pct','p_value','p_value_BH_within_cycle')],row.names=FALSE)
cat('\nEXPLORATORY MBP-by-cycle interaction:\n')
print(int[,c('n','residual_df','beta','se','p_value')],row.names=FALSE)
cat('\n[saved]',out,'\n')
