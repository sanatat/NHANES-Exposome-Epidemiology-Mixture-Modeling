# V13 — EPA Assay Candidate Evidence Audit (annotation-only)

Requires outputs from `23_v12_mechanistic_toxcast.py`. This patch does not download new data and does not make bioactivity/hit or disease causality claims.

Run in `~/nhanes_exposome_project`:

```bash
.venv/bin/python 26_v13_synthetic_tests.py
.venv/bin/python 25_v13_audit_assay_candidates.py
```

Outputs at `results/v13_mechanism_evidence_audit/`:

- `v13_disjoint_domain_overlap.csv`: seven mutually exclusive intersections, so summed counts are not misread as unique assays.
- `v13_target_gene_supported_shortlist.csv`: keyword-tagged endpoints with exact official target gene symbols matching **exploratory** predefined panels.
- `v13_keyword_only_review.csv`: annotations tagged by keywords without exact gene support.
- `v13_all_tagged_endpoints.csv`, `v13_domain_evidence_summary.csv`, `v13_audit.json`, `v13_annotation_evidence.png`.

**Strict interpretation:** These are metadata candidate lists. A target gene annotation does **not** demonstrate that any of seven chemicals has been tested or has activated the endpoint. For actual activity, acquire official per-chemical invitrodb v4.3 mc5-6 results and rerun v12 with `--mc56`. EPA's CTX Bioactivity API is an alternative but requires an individual API key; request it via `ccte_api@epa.gov`. Do not replace missing toxicity with zero. The panel is for review, not a validated pathway enrichment analysis.
