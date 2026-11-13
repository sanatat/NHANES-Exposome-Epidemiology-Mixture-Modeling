#!/usr/bin/env python3
"""V16: reproducible observational audit of ToxCast invitrodb v4.3 mirror output.

This DOES NOT infer effects in humans. It keeps observations at sample/endpoint
level and avoids interpreting missing tests as inactive compounds. No web access.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "results/v15_measured_bioactivity/v15_matched_raw_assays.csv"
MATRIX = ROOT / "results/v15_measured_bioactivity/v15_chemical_x_aeid_matrix.csv"
PANEL = ROOT / "results/v14_mechanistic_panel/v14_chemical_x_endpoint_evidence.csv"
AUDIT = ROOT / "results/v15_measured_bioactivity/v15_audit.json"
OUT = ROOT / "results/v16_bioactivity_quality_audit"
CHEM_ORDER = ["PFHxS", "PFNA", "PFOA", "PFOS", "MBP", "MEP", "MiBP"]


def status_of(n_records, valid, positive, opposite):
    if not n_records:
        return "NOT_RECORDED_IN_SOURCE"
    if not valid:
        return "RECORDED_NO_VALID_HITC"
    if opposite:
        return "OPPOSITE_DIRECTION_REVIEW"
    if positive == valid:
        return "ALL_VALID_HIGH_POSITIVE"
    if positive:
        return "MIXED_HIGH_AND_BELOW"
    return "NO_HIGH_POSITIVE_HITC"


def calculate(matrix, raw, panel):
    required_m = {"chemical", "casrn", "aeid", "tier", "retrieval_status"}
    required_r = {"casrn", "aeid", "hitc_numeric"}
    if not required_m.issubset(matrix):
        raise ValueError(f"V15 matrix missing columns: {required_m-set(matrix)}")
    if not required_r.issubset(raw):
        raise ValueError(f"V15 raw records missing columns: {required_r-set(raw)}")
    matrix = matrix.copy()
    raw = raw.copy()
    matrix["aeid"] = pd.to_numeric(matrix["aeid"], errors="raise").astype("int64")
    raw["aeid"] = pd.to_numeric(raw["aeid"], errors="raise").astype("int64")
    matrix["casrn"] = matrix["casrn"].astype(str).str.strip()
    raw["casrn"] = raw["casrn"].astype(str).str.strip()
    if matrix.duplicated(["casrn", "aeid"]).any():
        raise ValueError("Duplicate chemical/AEID keys in V15 matrix")
    if raw[["casrn", "aeid"]].merge(matrix[["casrn", "aeid"]].drop_duplicates(), on=["casrn", "aeid"], how="left", indicator=True)["_merge"].eq("left_only").any():
        raise ValueError("Raw bioactivity includes CAS/AEID not present in the curated V15 matrix")
    raw["hitc_numeric"] = pd.to_numeric(raw["hitc_numeric"], errors="coerce")
    valid_hitc = raw["hitc_numeric"].dropna()
    if ((valid_hitc < -1) | (valid_hitc > 1)).any():
        raise ValueError("Hitc outside documented [-1, 1] interval")

    recs=[]
    for (cas, aeid), g in raw.groupby(["casrn", "aeid"], sort=False):
        h=g["hitc_numeric"]
        v=h.notna()
        pos = h.ge(.9)
        neg = h.le(-.9)
        sample_col = next((c for c in ("spid", "sample_id", "tx_sample_id") if c in g.columns), None)
        ids = g[sample_col].dropna().astype(str) if sample_col else pd.Series([], dtype=str)
        ac50 = pd.to_numeric(g.get("ac50_uM_if_verified", pd.Series(np.nan,index=g.index)), errors="coerce")
        eligible=ac50[pos & ac50.gt(0)]
        recs.append({
            "casrn":cas,"aeid":aeid,"n_source_records":len(g),
            "n_valid_hitc":int(v.sum()),"n_high_positive":int(pos.sum()),
            "n_high_opposite":int(neg.sum()),"n_below_positive_threshold":int((v & ~pos & ~neg).sum()),
            "n_missing_hitc":int((~v).sum()),
            "sample_identifier_column": sample_col or "UNAVAILABLE",
            "n_unique_sample_ids":int(ids.nunique()) if sample_col else np.nan,
            "n_repeated_sample_id_records":int(len(ids)-ids.nunique()) if sample_col else np.nan,
            "min_ac50_uM_from_verified_units":float(eligible.min()) if len(eligible) else np.nan,
            "n_positive_records_with_verified_ac50":int(eligible.count()),
        })
    agg=pd.DataFrame(recs)
    base=matrix[["chemical", "casrn", "aeid", "tier"]].copy()
    details=[]
    for c in ("assay_component_endpoint_name", "mapped_gene_symbols", "priority_topic", "interpretation_note"):
        if c in panel.columns: details.append(c)
    if details:
        meta=panel[["casrn", "aeid", *details]].drop_duplicates(["casrn","aeid"])
        base=base.merge(meta, how="left", on=["casrn", "aeid"], validate="one_to_one")
    if len(agg):
        base=base.merge(agg, how="left", on=["casrn", "aeid"], validate="one_to_one")
    for c in ("n_source_records","n_valid_hitc","n_high_positive","n_high_opposite","n_below_positive_threshold","n_missing_hitc","n_positive_records_with_verified_ac50"):
        if c not in base: base[c]=0
        base[c]=base[c].fillna(0).astype(int)
    base["sample_identifier_column"]=base.get("sample_identifier_column", pd.Series("UNAVAILABLE", index=base.index)).fillna("UNAVAILABLE")
    base["evidence_status"]=[status_of(r.n_source_records,r.n_valid_hitc,r.n_high_positive,r.n_high_opposite) for r in base.itertuples(index=False)]
    expected=matrix[["casrn", "aeid", "retrieval_status"]].copy()
    x=base[["casrn", "aeid", "n_source_records"]].merge(expected,on=["casrn","aeid"],validate="one_to_one")
    if ((x["n_source_records"]>0) != x["retrieval_status"].eq("OBSERVED_IN_SOURCE")).any():
        raise ValueError("V15 matrix and raw source disagree on observed vs missing pairs")
    order={c:i for i,c in enumerate(CHEM_ORDER)}
    base["_order"]=base["chemical"].map(order).fillna(99)
    base=base.sort_values(["_order","tier","aeid"]).drop(columns=["_order"]).reset_index(drop=True)
    return base


def summarize(pairs):
    grouped=[]
    for chemical,g in pairs.groupby("chemical",sort=False):
        core=g[g.tier.eq("CORE_12")]
        grouped.append({
            "chemical":chemical,
            "n_endpoint_candidates":len(g),
            "n_endpoint_observed":int(g.n_source_records.gt(0).sum()),
            "n_endpoint_unobserved":int(g.n_source_records.eq(0).sum()),
            "n_endpoint_with_positive_hitc":int(g.n_high_positive.gt(0).sum()),
            "n_endpoint_mixed_positive_nonpositive":int(g.evidence_status.eq("MIXED_HIGH_AND_BELOW").sum()),
            "n_endpoint_opposite_direction_review":int(g.n_high_opposite.gt(0).sum()),
            "n_core12_observed":int(core.n_source_records.gt(0).sum()),
            "n_core12_with_positive_hitc":int(core.n_high_positive.gt(0).sum()),
            "n_core12_source_records":int(core.n_source_records.sum()),
            "n_core12_high_positive_records":int(core.n_high_positive.sum()),
        })
    return pd.DataFrame(grouped)


def draw_core(pairs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    core=pairs[pairs.tier.eq("CORE_12")].copy()
    aeids=sorted(core.aeid.unique())
    cmap=ListedColormap(["#e9eaee", "#d7e7f4", "#f5d7a6", "#d1dfc0", "#b4a3db"])
    mapping={"NOT_RECORDED_IN_SOURCE":0,"RECORDED_NO_VALID_HITC":0,
             "NO_HIGH_POSITIVE_HITC":1,"MIXED_HIGH_AND_BELOW":2,
             "ALL_VALID_HIGH_POSITIVE":3,"OPPOSITE_DIRECTION_REVIEW":4}
    tab=core.set_index(["chemical","aeid"])
    arr=np.zeros((len(CHEM_ORDER),len(aeids)),dtype=int)
    labels={}
    for i,chem in enumerate(CHEM_ORDER):
        for j,aeid in enumerate(aeids):
            if (chem,aeid) not in tab.index:
                continue
            r=tab.loc[(chem,aeid)]
            arr[i,j]=mapping[r.evidence_status]
            if r.n_source_records:
                labels[(i,j)]=f"{int(r.n_high_positive)}/{int(r.n_valid_hitc)}"
            else:
                labels[(i,j)]="—"
    fig,ax=plt.subplots(figsize=(14,5))
    ax.imshow(arr,aspect="auto",cmap=cmap,norm=BoundaryNorm(np.arange(-.5,5.5,1),5))
    ax.set_xticks(np.arange(len(aeids)),[str(a) for a in aeids],rotation=0)
    ax.set_yticks(np.arange(len(CHEM_ORDER)),CHEM_ORDER)
    ax.set_xlabel("ToxCast assay endpoint ID (AEID)")
    ax.set_title("V16: Exploratory ToxCast core 12 | cells = high hitc / valid records")
    for (i,j),txt in labels.items(): ax.text(j,i,txt,ha="center",va="center",fontsize=8)
    fig.text(.02,.02,"Gray = no source record; blue = no high-positive; orange = mixed; green = all high-positive; purple = opposite-direction review",fontsize=8)
    fig.subplots_adjust(left=.08,right=.99,top=.88,bottom=.18)
    fig.savefig(path,dpi=175)
    plt.close(fig)


def run(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument("--raw", type=Path,default=RAW)
    p.add_argument("--matrix", type=Path,default=MATRIX)
    p.add_argument("--panel", type=Path,default=PANEL)
    p.add_argument("--out",type=Path,default=OUT)
    a=p.parse_args(argv)
    for f in (a.raw,a.matrix,a.panel):
        if not f.is_file():raise FileNotFoundError(f"Missing input: {f}")
    matrix=pd.read_csv(a.matrix,low_memory=False)
    raw=pd.read_csv(a.raw,low_memory=False)
    panel=pd.read_csv(a.panel,low_memory=False)
    pairs=calculate(matrix,raw,panel)
    summ=summarize(pairs)
    core=pairs[pairs.tier.eq("CORE_12")]
    a.out.mkdir(parents=True,exist_ok=True)
    pairs.to_csv(a.out/"v16_per_chemical_endpoint_evidence.csv",index=False)
    summ.to_csv(a.out/"v16_chemical_coverage_and_hits.csv",index=False)
    core.to_csv(a.out/"v16_core12_status.csv",index=False)
    draw_core(pairs,a.out/"v16_core12_hitcall_heatmap.png")
    out={"n_possible_pairs":len(pairs),"n_observed_pairs":int(pairs.n_source_records.gt(0).sum()),
         "n_missing_pairs":int(pairs.n_source_records.eq(0).sum()),
         "n_core_pairs":len(core),"n_core_observed":int(core.n_source_records.gt(0).sum()),
         "n_core_pairs_with_high_hitc":int(core.n_high_positive.gt(0).sum()),
         "n_core_high_hitc_source_records":int(core.n_high_positive.sum()),
         "source":"Secondary public ToxCast invitrodb v4.3 mirror queried by V15; NOT validated against EPA primary download",
         "caveats":["Multiple records may reflect sample IDs or study replicates, not independent biological confirmations.",
                    "No-record is not inactive.","Positive screening hit is not proof of human toxicity or inflammation.",
                    "Active hitcall can be influenced by cytotoxicity and assay interference.",
                    "Check tcpl mc6 flags, representative-sample designation and official EPA export before publication.",
                    "Potency is retained only if original unit was verified as micromolar."]}
    if AUDIT.is_file():
        out["v15_audit_source"]=json.loads(AUDIT.read_text()).get("source","UNKNOWN")
    (a.out/"v16_audit.json").write_text(json.dumps(out,indent=2))
    print("[V16] observed",out["n_observed_pairs"],"/",out["n_possible_pairs"],"pairs")
    print("[V16] core12 observed",out["n_core_observed"],"/",out["n_core_pairs"],"pairs; pairs with high hitc:",out["n_core_pairs_with_high_hitc"])
    print("[V16] Per-chemical core12 audit:\n",summ[["chemical","n_core12_observed","n_core12_with_positive_hitc","n_core12_source_records","n_core12_high_positive_records"]].to_string(index=False))
    print("[V16] saved:",a.out)
    return out

if __name__=="__main__":run()
