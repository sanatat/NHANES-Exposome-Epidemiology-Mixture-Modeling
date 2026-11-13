# V19 — EPA ToxCast Level 6 caution flag taxonomy (exploratory)

## What it does
Reads existing V18 positive-record quality CSV. Splits records with reported QC flags into recognized EPA tcpl caution flag types and indicated-but-unresolved types. Never labels positive hits as verified.

## Run
```
cd ~/nhanes_exposome_project
.venv/bin/python 38_v19_tests.py
.venv/bin/python 37_v19_flag_taxonomy.py
```

## Output
- `results/v19_flag_taxonomy/v19_flag_frequency.csv`
- `results/v19_flag_taxonomy/v19_positive_pairs_review.csv`
- `results/v19_flag_taxonomy/v19_positive_records_flag_taxonomy.csv`
- `results/v19_flag_taxonomy/v19_audit.json`

## Scientific limitations
Source is a third-party mirror. Flag parsing is textual and cannot validate the curve fit. If a record has only an indication of flags but no recognized EPA flag names, it is labeled unresolved. Missing flag text does **not** imply clean or validated data. No assay potency or risk ranking is produced.

EPA tcpl documentation: https://usepa.r-universe.dev/tcpl/doc/Introduction_Appendices.html
