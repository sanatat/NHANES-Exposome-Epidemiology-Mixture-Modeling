# Results recorded from local Linux Mint console outputs

> **Provenance.** These numeric results were transcribed from terminal outputs supplied by the project operator in the conversation. The underlying NHANES and ToxCast exports are intentionally **not bundled**, and **we have not independently re-executed the scientific analysis on the raw datasets** in this distribution. Re-run the scripts to reproduce and verify them before citing or publishing a scientific manuscript.

## NHANES 2015–2018 pooled cohort (v6)

- Complete-case pooled participants: **2,784** (2015–16: 1,430; 2017–18: 1,354).
- Primary analysis, hs-CRP ≤10 mg/L: **2,568** (2015–16: 1,319; 2017–18: 1,249).
- hs-CRP ≤20 mg/L sensitivity: **2,719**.
- Survey design degrees of freedom: **30**; effective sample size reported by the script: **1,206.903** in primary cohort.
- Serum PFAS: PFHxS, PFNA, PFOA, PFOS; urinary phthalate metabolites: MBP, MEP, MiBP.

### Primary survey-weighted additive quantile mixture models

| Exposure family | Estimated hs-CRP difference | 95% CI | p | BH p (3 tests) |
|---|---:|---|---:|---:|
| All seven exposures | +3.69% | −5.66% to +13.97% | 0.4199 | 0.6298 |
| Four PFAS | −0.24% | −6.61% to +6.57% | 0.9401 | 0.9401 |
| Three phthalate metabolites | +3.92% | −4.75% to +13.38% | 0.3629 | 0.6298 |

No mixture association met a two-sided 0.05 threshold. These estimates are *associations*, not causal effects.

### Selected single-exposure signal (v6–v8)

- Urinary **MBP**: +5.25% hs-CRP per doubling (reference model); unadjusted for multiple testing p=0.0181; **FDR-adjusted p=0.1267** over seven chemicals.
- Cycle 2015–16: +0.99% (95% CI −5.49%, +7.92%); cycle 2017–18: +9.95% (+0.93%, +19.79%); cycle interaction p=0.1088. This is not independent replication because the candidate was selected after inspecting pooled data.
- v8: 14 *post-hoc correlated* sensitivity specifications. 13/14 nominal p-values <0.05; this **does not** represent 13 independent confirmations or supersede FDR-adjusted primary results. No evidence of nonlinearity by restricted cubic spline (exploratory p=0.4980).

### Cross-cycle prediction (v10)

Additive chemicals did **not** improve weighted out-of-cycle MSE in any of six matched-algorithm comparisons (Ridge, Elastic Net, Random Forest; each direction). Positive deltas below mean **increased** MSE (worse):

| Model | 2015→2017: ΔMSE (chemical minus base) | 2017→2015: ΔMSE |
|---|---:|---:|
| Ridge | +0.0054 | +0.0060 |
| Elastic Net | +0.0035 | +0.0052 |
| Random Forest | +0.0010 | +0.0132 |

### Toxicology bridge / screening evidence (v11–v19)

- Exact CAS cytotoxicity table matches: **6/7** (MiBP lacked an exact CAS match in the queried table).
- Annotation candidates: **172 distinct endpoints**, of which **44** met a targeted-gene shortlist criterion (target mapping is not proof of direct binding).
- Third-party ToxCast mirror query: **170/308** chemical × endpoint pairs observed; **138/308** not recorded in queried source. Missing is **not** negative.
- Across 44 endpoints, **52** chemical–endpoint pairs had at least one positive screening hit; **17** in core 12. There were **86** positive source records. Do not regard distinct sample records as independent biological replication.
- v18–v19 initially marked 79 records as containing unresolved numeric MC6 flag IDs and seven with missing information. Flag IDs require correct decoding and individual curve review. v20 code includes a decoder but its results were **not provided/re-run on real data** when this document was assembled.
- PFNA/PFOS/PFOA NRF2/ARE positives deserve quality review; AC50 values are in vitro assay estimates, **not environmental human risk thresholds**. `fitc=100`/LOEC should not be reported as a fitted AC50.

## Essential qualifications

This is an exploratory, cross-sectional analysis with post-hoc model development. The outcome is hs-CRP, **not cancer**. Chemical assays in NHANES were **targeted HPLC–MS/MS**, not raw untargeted LC–HRMS/GC–HRMS processing. PFAS and phthalate metabolites use different biological matrices and timescales. Exposure mixtures may be correlated, and serum/urine concentrations are not direct measures of causal dose. The public ToxCast mirror should be checked against EPA's official release, concentration–response QC, and original assay annotations before any mechanistic assertion. No risk scoring or human causality is established here.
