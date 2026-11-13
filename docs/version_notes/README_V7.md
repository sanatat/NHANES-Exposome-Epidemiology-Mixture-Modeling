# v7 — Exploratory consistency of MBP across NHANES cycles

This is an **add-on to v6**, not a replacement for its primary analyses.

- **Input:** `data/processed/analysis_joint_pooled_2015_2018.csv` generated in v6.
- **Run:** `Rscript 11_cycle_consistency.R`.
- **Outputs:** `results/pooled_2015_2018/v7_exploratory_cycle_exposures.csv`, `v7_exploratory_MBP_cycle_interaction.csv`, and `v7_exploratory_replication_notes.txt`.

The pooled v6 MBP analysis had unadjusted `p = 0.0181` and BH FDR `p = 0.1267`. Therefore, **MBP is not statistically significant under the pre-declared 7-test FDR analysis**.

This add-on tests all seven exposures separately within each cycle and an MBP-by-cycle slope interaction, with survey weights and the same confounder set. These are explicitly **post-hoc consistency checks** (not external replication, not new primary hypotheses); report all results even if null or contradictory. Do not select methods based on the p-value. No cancer or causal inference is claimed.

Neither the pooled input dataset nor earlier results are modified.
