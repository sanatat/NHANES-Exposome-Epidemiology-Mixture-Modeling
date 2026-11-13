#!/usr/bin/env Rscript
# NHANES 2017-18 exploratory exposome: independent flag audit plus
# simplified, survey-weighted additive quantile mixture sensitivity models.
# Previously generated files are NOT changed. NO CAUSAL INFERENCE.
args <- grep('^--file=', commandArgs(trailingOnly=FALSE), value=TRUE)
if (length(args) != 1) stop('Run: Rscript 05_lod_and_parsimonious_mixture.R')
root <- dirname(normalizePath(sub('^--file=', '', args[1])))
if (!requireNamespace('survey', quietly=TRUE)) stop('Install: sudo apt install r-cran-survey')
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust')
input <- file.path(root, 'data/processed/analysis_joint.csv')
if (!file.exists(input)) stop('Missing ',input,' (run the existing download script first).')
outdir <- file.path(root, 'results')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

d <- read.csv(input, check.names=FALSE)
chemicals <- c('LBXPFHS','LBXPFNA','LBXNFOA','LBXNFOS','URXMBP','URXMEP','URXMIB')
flags <- c('LBDPFHSL','LBDPFNAL','LBDNFOAL','LBDNFOSL','URDMBPLC','URDMEPLC','URDMIBLC')
needed <- c('SEQN','survey_weight','SDMVSTRA','SDMVPSU','LBXHSCRP',
            'RIDAGEYR','RIAGENDR','BMXBMI','URXUCR',chemicals,flags)
if (!all(needed %in% names(d))) stop('Missing: ',paste(setdiff(needed,names(d)),collapse=', '))
if (anyDuplicated(d$SEQN)) stop('Duplicate SEQN values')
# Match the original complete-case analytic cohort, excluding only invalid records.
keep <- complete.cases(d[,needed,drop=FALSE]) &
        is.finite(d$survey_weight) & d$survey_weight>0 &
        is.finite(d$LBXHSCRP) & d$LBXHSCRP>0 &
        is.finite(d$URXUCR) & d$URXUCR>0 & d$RIDAGEYR>=20
for (chem in chemicals) keep <- keep & is.finite(d[[chem]]) & d[[chem]]>0
d <- d[keep,,drop=FALSE]
if (nrow(d)<100) stop('Unexpectedly small analysis cohort (',nrow(d),').')

# The actual denominator must include both flag==0 (detected) and flag==1 (non-detect).
# Output raw counts separately before any regression is run.
qc <- lapply(seq_along(chemicals),function(i) {
  x <- suppressWarnings(as.numeric(as.character(d[[flags[i]]])))
  w <- d$survey_weight
  eligible <- !is.na(x) & x %in% c(0,1)
  n0 <- sum(eligible & x==0)
  n1 <- sum(eligible & x==1)
  nonstandard <- sum(!is.na(x) & !(x %in% c(0,1)))
  n_missing <- sum(is.na(x))
  denom <- sum(w[eligible])
  pct_weighted <- if (denom>0) 100*sum(w[eligible & x==1])/denom else NA_real_
  pct_unweighted <- if(sum(eligible)>0) 100*n1/sum(eligible) else NA_real_
  data.frame(chemical=chemicals[i],flag=flags[i],n_total=nrow(d),
             n_detected_0=n0,n_below_LOD_1=n1,n_missing=n_missing,
             n_other=nonstandard,n_valid_flags=sum(eligible),
             pct_below_LOD_unweighted=pct_unweighted,
             pct_below_LOD_weighted=pct_weighted)
})
qc <- do.call(rbind,qc)
write.csv(qc,file.path(outdir,'v4_lod_qc.csv'),row.names=FALSE)
cat('\nLOD AUDIT: BOTH 0 AND 1 ARE COUNTED IN THE DENOMINATOR\n')
print(qc[,c('chemical','n_detected_0','n_below_LOD_1','n_missing','n_other',
           'pct_below_LOD_unweighted','pct_below_LOD_weighted')],row.names=FALSE)
if(any(qc$n_valid_flags != nrow(d))) {
  stop('Unexpected missing/invalid LOD flags; inspect v4_lod_qc.csv and the downloaded CDC files. ',
       'No model results were calculated.')
}

# Mixture quartiles are set on the full complete-case cohort, reused unchanged
# in the high-CRP sensitivity subset.
qnames <- paste0('q_',tolower(chemicals))
for(i in seq_along(chemicals)) {
  breaks <- unique(as.numeric(quantile(d[[chemicals[i]]],probs=seq(0,1,.25),na.rm=TRUE)))
  if(length(breaks) < 4) stop('Too few unique quartile cutoffs for ',chemicals[i])
  d[[qnames[i]]] <- as.integer(cut(d[[chemicals[i]]],breaks=breaks,
                                      include.lowest=TRUE,labels=FALSE))-1L
  if(anyNA(d[[qnames[i]]])) stop('Quartile assignment failed: ', chemicals[i])
}
d$log_hscrp <- log(d$LBXHSCRP)
d$age <- d$RIDAGEYR
d$sex <- factor(d$RIAGENDR)
d$bmi <- d$BMXBMI
d$log_urine_cr <- log(d$URXUCR)
# PRE-SPECIFIED REDUCED COVARIATES: lower parameter count allows positive
# survey residual df for the 7-chemical mixture. Omitting PIR/race leaves
# potential residual confounding: this is a sensitivity exercise, not the
# main fully adjusted or causal model.
reduced_covars <- c('age','sex','bmi','log_urine_cr')
run_model <- function(ds,selected,label,scenario) {
  design <- svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~survey_weight,
                      data=ds,nest=TRUE)
  q <- paste0('q_',tolower(selected))
  fit <- svyglm(reformulate(c(q,reduced_covars),response='log_hscrp'),
                design=design,family=gaussian())
  beta <- coef(fit)[q]
  if(anyNA(beta)) stop('Aliased exposure term in ',label,', ',scenario)
  v <- vcov(fit)[q,q,drop=FALSE]
  psi <- sum(beta)
  se <- sqrt(sum(v))
  df_design <- degf(design)
  df_res <- df_design - (fit$rank - 1L)
  crit <- if(df_res>0) qt(.975,df_res) else NA_real_
  p <- if(df_res>0 && se>0) 2*pt(-abs(psi/se),df_res) else NA_real_
  data.frame(scenario=scenario,mixture=label,n=nrow(ds),n_exposures=length(selected),
             design_df=df_design,model_rank=fit$rank,residual_df=df_res,
             psi=psi,SE=se,percent_change=100*expm1(psi),
             lower95_pct=100*expm1(psi-crit*se),
             upper95_pct=100*expm1(psi+crit*se),p_value_residual_df=p,
             interpretation='Exploratory association; reduced adjustment, not causal')
}
models <- list(list('All 7 chemicals',chemicals),
               list('PFAS only',chemicals[1:4]),
               list('Phthalates only',chemicals[5:7]))
scenarios <- list('Full cohort'=d,'hs-CRP <= 10 mg/L'=d[d$LBXHSCRP<=10,,drop=FALSE])
res <- list()
for(s in names(scenarios)) {
  ds <- scenarios[[s]]
  if(nrow(ds)<100) stop('Too few participants: ',s)
  for(m in models) res[[length(res)+1]] <- run_model(ds,m[[2]],m[[1]],s)
}
res <- do.call(rbind,res)
write.csv(res,file.path(outdir,'v4_reduced_covariate_mixture.csv'),row.names=FALSE)
cat('\nREDUCED-COVARIATE MIXTURE SENSITIVITY\n')
print(res[,c('scenario','mixture','n','residual_df','percent_change',
             'lower95_pct','upper95_pct','p_value_residual_df')],row.names=FALSE)
writeLines(c(
  'NHANES 2017-2018: survey-weighted, additive QUARTILE mixture pilot',
  paste('Complete-case adult cohort:',nrow(d)),
  paste('High-CRP restriction removed:',sum(d$LBXHSCRP>10),'participants'),
  'LOD percentage denominators: flags 0 or 1 only, using full cohort weights.',
  'Quartile cutoffs: full cohort and held fixed for both sensitivity scenarios.',
  paste('REDUCED adjustment covariates:',paste(reduced_covars,collapse=', ')),
  'Important limitation: race/ethnicity and income omitted to conserve residual survey df.',
  'Therefore estimates may be residually confounded. These models are sensitivity analyses, not a substitute for adequately adjusted epidemiology.',
  'Survey df: # PSU - # strata; model df for tests: design_df - (fit rank - 1).',
  'If residual df <= 0, this script intentionally reports NA CIs and p-values.',
  'Exposure coefficients summed for a simultaneous one-quantile category increase.',
  'Serum PFAS and urine phthalate metabolite measurements are not the same biological matrix.',
  'This is cross-sectional NHANES; associations are not causal and do not establish health effects.',
  'The earlier 7-exposure FULL-covariate model has negative residual df and must not be used for confirmatory inference.'
),file.path(outdir,'v4_methodological_notes.txt'))
cat('\nNew outputs in: ',outdir,'\n',sep='')
