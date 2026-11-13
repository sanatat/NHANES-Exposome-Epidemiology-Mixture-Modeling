# V18 — Screening-hit quality metadata audit

Reads existing V15 raw ToxCast records and V17 selected hit pairs. No downloads; original results v1-v17 unchanged.

Run:

```bash
cd ~/nhanes_exposome_project
.venv/bin/python 36_v18_tests.py
.venv/bin/python 35_v18_quality_evidence.py
```

Outputs in `results/v18_curve_quality_audit/`: record-level evidence with source QC fields; per-pair candidate audit (all 44 and core 12); JSON schema availability audit.

**Limits:** Fields from the third-party mirror are not official validation. If fit flags or concentration ranges are missing, those checks remain unperformed. An absence of flags is not proof of quality. The review queue is not a chemical risk ranking or evidence of causal pathways. Validate the original EPA invitrodb v4.3 mc5-6 and relevant mc6 flags/curve plots before publication.
