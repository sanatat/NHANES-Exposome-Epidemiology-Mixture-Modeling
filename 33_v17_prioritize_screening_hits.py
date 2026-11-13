#!/usr/bin/env python3
"""V17: transparent review queue for ToxCast hits in V16.

Offline, observational only. A high hitc is not a verified biological mechanism.
No evidence of independent replicate agreement is inferred from multiple spid values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DEFAULT_V16 = ROOT / "results/v16_bioactivity_quality_audit/v16_per_chemical_endpoint_evidence.csv"
DEFAULT_RAW = ROOT / "results/v15_measured_bioactivity/v15_matched_raw_assays.csv"
DEFAULT_OUT = ROOT / "results/v17_positive_hit_review"

NEEDED = {
    "chemical", "casrn", "aeid", "tier", "assay_component_endpoint_name",
    "mapped_gene_symbols", "n_source_records", "n_valid_hitc",
    "n_high_positive", "n_high_opposite", "n_below_positive_threshold",
    "n_missing_hitc", "n_unique_sample_ids", "n_repeated_sample_id_records",
    "min_ac50_uM_from_verified_units", "n_positive_records_with_verified_ac50",
    "evidence_status",
}


def select_hits(frame: pd.DataFrame) -> pd.DataFrame:
    missing = NEEDED - set(frame.columns)
    if missing:
        raise ValueError("Missing V16 fields: " + repr(sorted(missing)))
    d = frame.copy()
    d["aeid"] = pd.to_numeric(d["aeid"], errors="raise").astype(int)
    if d.duplicated(["casrn", "aeid"]).any():
        raise ValueError("Duplicate exact chemical × endpoint pairs")
    numeric = ["n_source_records", "n_valid_hitc", "n_high_positive", "n_high_opposite",
               "n_below_positive_threshold", "n_missing_hitc", "n_repeated_sample_id_records",
               "n_positive_records_with_verified_ac50"]
    for col in numeric:
        d[col] = pd.to_numeric(d[col], errors="raise")
    if (d["n_high_positive"] > d["n_valid_hitc"]).any():
        raise ValueError("Positive calls cannot exceed valid hitc records")
    d = d[d["n_high_positive"] > 0].copy()
    if d.empty:
        return d
    d["positive_fraction_of_valid_records"] = d["n_high_positive"] / d["n_valid_hitc"]
    d["positive_ac50_coverage_fraction"] = d["n_positive_records_with_verified_ac50"] / d["n_high_positive"]
    d["review_notes"] = d.apply(review_notes, axis=1)
    d["flag_review_status"] = "UNKNOWN__MC6_FIT_FLAGS_NOT_RETRIEVED"
    d["cytotoxicity_interference_status"] = "UNKNOWN__ASSAY_SPECIFIC_CONTEXT_NOT_VALIDATED"
    d["official_epa_crosscheck_status"] = "PENDING__THIRD_PARTY_MIRROR"
    # Sort for inspection; ordering is NOT chemical risk rank or biological replication.
    tier_order = {"CORE_12": 0}
    d["_tier"] = d["tier"].map(tier_order).fillna(1)
    d = d.sort_values(["_tier", "chemical", "aeid"], kind="stable").drop(columns="_tier")
    return d


def review_notes(r):
    notes = []
    if r.n_high_positive < r.n_valid_hitc:
        notes.append("mixed_positive_and_other_records")
    if r.n_source_records == 1:
        notes.append("single_source_record")
    if r.n_repeated_sample_id_records > 0:
        notes.append("repeated_sample_id_check")
    if r.n_positive_records_with_verified_ac50 < r.n_high_positive:
        notes.append("positive_ac50_missing_or_unit_unverified")
    if r.n_high_opposite > 0:
        notes.append("opposite_orientation_needs_assay_review")
    notes.append("mc6_flags_and_concentration_range_not_checked")
    return "; ".join(notes)


def details(raw: pd.DataFrame, selected: pd.DataFrame) -> pd.DataFrame:
    if selected.empty:
        return pd.DataFrame()
    required = {"casrn", "aeid", "hitc_numeric"}
    if not required.issubset(raw):
        raise ValueError("V15 raw records missing: " + repr(sorted(required - set(raw))))
    r = raw.copy()
    r["aeid"] = pd.to_numeric(r["aeid"], errors="raise").astype(int)
    r["hitc_numeric"] = pd.to_numeric(r["hitc_numeric"], errors="coerce")
    r = r.merge(selected[["chemical", "casrn", "aeid", "tier"]], on=["casrn", "aeid"], validate="many_to_one")
    r = r[r["hitc_numeric"] >= .9].copy()
    if "ac50_uM_if_verified" not in r:
        r["ac50_uM_if_verified"] = np.nan
    cols = ["chemical", "casrn", "aeid", "tier", "hitc_numeric", "ac50_uM_if_verified"]
    for extra in ["spid", "sample_id", "tx_sample_id", "ac50_original", "ac50_unit_verified", "conc_unit", "fitc", "flags", "mc6_flags", "max_conc", "max_conc_um"]:
        if extra in r.columns:
            cols.append(extra)
    return r[cols].sort_values(["chemical", "aeid"]).reset_index(drop=True)


def run(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--v16", type=Path, default=DEFAULT_V16)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    for f in (args.v16, args.raw):
        if not f.is_file():
            raise FileNotFoundError(f"Missing data file: {f}")
    v16 = pd.read_csv(args.v16, low_memory=False)
    raw = pd.read_csv(args.raw, low_memory=False)
    candidates = select_hits(v16)
    records = details(raw, candidates)
    if len(records) != int(candidates["n_high_positive"].sum()):
        raise ValueError("Mismatch between V16 counts and positive V15 raw records; refuse report")
    core = candidates[candidates["tier"] == "CORE_12"].copy()
    args.out.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(args.out / "v17_all44_positive_endpoint_pairs.csv", index=False)
    core.to_csv(args.out / "v17_core12_positive_review.csv", index=False)
    records.to_csv(args.out / "v17_positive_source_records.csv", index=False)
    audit = {
        "source": "V15 third-party ToxCast mirror / V16 counts",
        "total_pairs_input": len(v16),
        "positive_pairs_all44": len(candidates),
        "positive_pairs_core12": len(core),
        "positive_source_records_all44": len(records),
        "core12_positive_pairs_by_chemical": {str(k): int(v) for k, v in core.groupby("chemical").size().items()},
        "limitations": [
            "hitc >= 0.9 is a screening hit, not verified pathway activation or human risk",
            "the sample records are not necessarily independent biological confirmations",
            "fit quality flags, curve plots and concentration ranges require original EPA mc6 or tcpl verification",
            "ac50 figures are exploratory, only when the source unit was explicitly verified as µM",
            "cytotoxicity assay-specific interference and physiological relevance not established",
            "this is not a chemical hazard ranking and does not change NHANES epidemiological results",
        ],
    }
    (args.out / "v17_evidence_gaps.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(f"[V17] review candidate pairs: {len(candidates)} (core12: {len(core)})")
    print(f"[V17] positive raw source records: {len(records)}")
    if len(core):
        display_cols = ["chemical", "aeid", "mapped_gene_symbols", "n_high_positive", "n_valid_hitc", "min_ac50_uM_from_verified_units", "review_notes"]
        print("\n[V17] CORE12 PAIRS (AC50 in µM only when verified):")
        print(core[display_cols].to_string(index=False, max_colwidth=66))
    print(f"\n[V17] results saved: {args.out}")
    print("[V17] WARNING: mc6 flags, concentration ranges, and original EPA export remain unverified")
    return audit


if __name__ == "__main__":
    run()
