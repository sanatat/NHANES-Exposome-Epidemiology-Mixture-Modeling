# v12 — Mechanistic Exposomics (EPA ToxCast invitrodb v4.3)

**Goal.** Build an auditable bridge between *chemical-specific NHANES biomarkers* and published EPA *in vitro assay endpoints*. The analysis does **not** estimate human toxicity, disease causation, or exposure-adjusted potency. The ToxCast results are **chemical-level** rather than person-level observations.

## Requirements

Use the existing `~/nhanes_exposome_project/.venv`, containing pandas, requests, openpyxl, matplotlib. No changes to previously generated v6–v11 results.

```bash
cd ~/nhanes_exposome_project
bash RUN_V12.sh
```

Default run downloads and SHA-256-verifies **only two small official EPA workbooks**:
- Assay annotations: https://clowder.edap-cluster.com/files/68af6bd3e4b02565fc7c3aa8
- Target mappings: https://clowder.edap-cluster.com/files/68af6bd3e4b02565fc7c3aa0

Results: `results/v12_mechanistic_exposomics/v12_candidate_domain_inventory.csv`, `v12_annotated_assay_inventory.csv`, `v12_annotation_inventory.png`, `v12_provenance_audit.json`. Endpoint/domain matches are **keyword candidates**, not validated pathway annotations. With no full assay results, **no chemical hitcounts are computed.**

## Full dose–response analysis (optional, large download)

The official `INVITRODB_SUMMARY.zip` is **approximately 7.5 GB**. It is never auto-downloaded.

- Source page: https://clowder.edap-cluster.com/files/68af6b70e4b02565fc7c3a98
- Verified official ZIP SHA-256 (published on file page): `c61d46cf685d5d84e72d0d4f6f46746ab1714d6f924e43480d390925ae827c59`
- Extract just `mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv` (the archive member has `AUG2024` in its filename, though published in the v4.3 release in August 2025).

One option to unpack only the CSV after downloading the ZIP manually:

```bash
unzip -j ~/Downloads/INVITRODB_SUMMARY.zip '*mc5-6_winning_model_fits-flags*.csv' -d ~/nhanes_exposome_project/data/toxicology/

cd ~/nhanes_exposome_project
.venv/bin/python 23_v12_mechanistic_toxcast.py \
  --mc56 data/toxicology/mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv
```

Full output includes a seven-chemical × three-domain coverage/hit-count table, assay-level exact-CAS records and a heatmap. Activity criterion is a **continuous hitcall ≥0.9** at each tested endpoint; counts are **NOT** risk scores, and do not distinguish assay quality or potency without an expert review. Multiple assays can represent the same biological pathway. No tests performed should not be coded as negative tests. Unannotated endpoints are not coded as inactive. No father/mother parent chemical results are substituted for phthalate metabolites.

## Quality interpretation

- The three domain labels (inflammation, oxidative stress, immune response) use **exploratory keyword filters** on assay annotations and gene mappings. These are not confirmatory pathway assessments. Manually review the matched assay and gene mappings before citing specific mechanism claims.
- `casrn` is a chemical identifier and does not guarantee identical salt/protonation/formulation or assay-relevant exposure.
- Cytotoxicity-burst summary from v11 is NOT equivalent to the total number of biological endpoints tested.
- There is **no connection at the individual participant level** between ToxCast endpoints and NHANES records.
- A finding of no activity in a particular assay panel should not be interpreted as evidence that a chemical is safe.
- The full CSV path is a local file; downloading it can require significant disk space (both ZIP and extracted CSV).

**Sources:** EPA (2025), *ToxCast invitrodb v4.3*, DOI: https://doi.org/10.23645/epacomptox.6062623.v14 ; https://www.epa.gov/comptox-tools/exploring-toxcast-data
