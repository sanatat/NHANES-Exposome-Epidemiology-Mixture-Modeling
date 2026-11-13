# V11 — EPA ToxCast cytotoxicity evidence bridge (NHANES 2015–2018)

**Purpose:** connect the seven existing NHANES exposure biomarkers to official **US EPA invitrodb v4.3 (released August 2025)** cytotoxicity *burst* summaries. This patch does not alter v5–v10 outputs, retrain models, claim a cancer result, or chase low p-values. The source data are **chemical-level in vitro screens**, not data measured in NHANES participants.

## Install and run (Linux Mint)

This patch goes into the existing `~/nhanes_exposome_project` directory. `v10` is optional for this analysis (the pooled v6 CSV is required for human exposure summaries).

```bash
cd ~/nhanes_exposome_project
.venv/bin/python -m pip install 'openpyxl>=3.1,<4'
.venv/bin/python 22_v11_offline_tests.py
.venv/bin/python 20_v11_epa_toxcast_bridge.py
```

The last command tries to download **only a ~943 kB official cytotoxicity Excel table**, validates the exact SHA-256 against EPA's published file digest, and caches it at `data/toxicology/cytotox_invitrodb_v4_3_AUG2025.xlsx`. No API key required for the public download. Network access to Clowder can sometimes be restricted. In that event, **do not download the ~7.5 GB full release automatically**. Open the official public file page in a browser and click **Download**:

https://clowder.edap-cluster.com/files/68af6bd3e4b02565fc7c3aa4

Save the file to `~/nhanes_exposome_project/data/toxicology/cytotox_invitrodb_v4_3_AUG2025.xlsx` and re-run the last command. The cached file is not overwritten. The script checks SHA-256:

`3bc444f658c80e4ac9a1189f083c0120365b504e10d0aa4e0f07feda294c08dd`

If you already have a corresponding file in a custom location, run `--cytotox /path/to/file.xlsx`; the script flags any checksum mismatch instead of falsely certifying it as official.

## What this adds

1. **Chemical identity crosswalk**: NHANES column, biomarker, CAS registry number, family and specimen. Structural and salt/isomer comparability still requires evaluation for biological interpretation.
2. **ToxCast toxicity-screening coverage**: exact-CAS matches only; reported `ntested`, `nhit`, cytotoxicity median/lower bound **if supplied** in the official spreadsheet. Missing data remain NA, not zero or safe.
3. **Parent/metabolite audit**: DBP, DEP and DiBP data are separated into **CONTEXT ONLY**; they are *not* substituted for MBP, MEP or MiBP biomarker data.
4. **Human exposure statistics from NHANES**: weighted geometric mean and percentage below limit of detection, with primary sample hs-CRP > 0 and <=10 mg/L, from the existing v6 pooled CSV.
5. **Two-panel figure**: NHANES LOD frequencies and EPA fraction of active cytotoxicity *burst* assay endpoints. The panels are separate because unlike endpoints, specimen types and units cannot be directly compared as a common risk scale.
6. **Traceability**: source timestamp, SHA-256, exact-match coverage, workbook sheet and original columns recorded in a JSON audit.

## Outputs

`results/v11_toxicology_bridge/`

- `v11_verified_chemical_crosswalk.csv`
- `v11_exact_metabolites_and_PFAS.csv`
- `v11_parent_plasticizers_CONTEXT_ONLY.csv`
- `v11_two_evidence_sources.png`
- `v11_provenance_and_qc.json`

**Optional deeper layer:** if you separately obtain the full official invitrodb v4.3 multi-concentration `mc5-6_winning_model_fits-flags_...csv`, run:

```bash
.venv/bin/python 20_v11_epa_toxcast_bridge.py \
  --mc56 /path/to/mc5-6_winning_model_fits-flags_invitrodbv4_3_AUG2024.csv
```

This streams the large CSV in chunks, filters exact-CAS rows, and writes additional assay-level evidence and coverage. **No large download is triggered automatically.** A hit call >=0.9 denotes ToxCast activity evidence *for that endpoint*, not clinical hazard. Results are not normalized by assay selection or cytotoxicity interference and should not be naively compared as risk percentages.

EPA release and documentation: https://www.epa.gov/comptox-tools/exploring-toxcast-data

EPA release DOI: https://doi.org/10.23645/epacomptox.6062623.v14

CDC source names: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PFAS_J.htm and https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/PHTHTE_J.htm

## Critical scientific limitations

- These seven assays are **NHANES targeted chemical biomarkers**, not untargeted LC-HRMS/GC-HRMS instrument data. This step adds orthogonal toxicology context, not mass-spectrometry preprocessing.
- Human urinary metabolites and parent plasticizers are chemically distinct. **Do not assign a parent's in vitro potency to its urinary metabolite.**
- ToxCast screening results vary by assay coverage, chemical sample quality and test concentrations. `cytotox_median_um = 1000` can represent a default ceiling and must not be interpreted uncritically as measured potency.
- ToxCast effects cannot be converted to population health risk, causal mechanisms or individual-level outcomes just by joining on chemical identity. The bridge is descriptive, not inferential.
- It is not legitimate to report any of this study as a confirmed health or mechanistic effect for MBP, nor to use the integration to obtain significant NHANES p-values.
- **All actual EPA values are downloaded at runtime on the user's machine**; our automated testing covers schema handling and exact-CAS matching on synthetic data only. No fake EPA bioactivity metrics are bundled.

### Expected GitHub figure (available after executing v11)

![Separate human exposure detection and EPA in vitro screening coverage](results/v11_toxicology_bridge/v11_two_evidence_sources.png)

If GitHub shows a broken image before you run the pipeline, that means the figure has not yet been generated and committed. Do not replace it with a simulated or invented figure.
