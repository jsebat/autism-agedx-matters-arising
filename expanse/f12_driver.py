#!/usr/bin/env python3
"""Driver analysis (driver_pgs.py, unchanged logic) with the PGS of the two published autism factors of Zhang, Grove et al. 2025:
F1 (earlier-diagnosed factor), F2 (later-diagnosed factor) and F2 - F1 (the contrast their model says tracks age at diagnosis).
Scores: SBayesRC weights on GSA v1 / GSA v2 / Twist EUR (pgs_warrier/scores), standardized to non-autistic relatives of the same
platform; pTDT in complete trios (as agedx_strata.py).
BIAS NOTE: F1 contains SPARK cases diagnosed < 6 y and F2 SPARK cases diagnosed > 10 y, so the age-at-diagnosis split is favoured
by sample overlap; any other split that matches or beats it does so despite that advantage.
Output (aggregate): ~/rvh2/agedx/f12_driver_*.tsv"""
import os
import numpy as np, pandas as pd
B = "/expanse/projects/sebat1/s3/data/sebat"; SC = "/expanse/projects/sebat1/jsebat/pgs_warrier/scores"; OUT = os.path.expanduser("~/rvh2/agedx")
R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift"
d = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype={"sample_id": str, "IID": str}, low_memory=False,
                usecols=["source", "IID", "sample_id", "case_status", "predicted_ancestry"])
d = d[d.source.isin(["spark_gsa_v1", "spark_gsa_v2", "spark_twist"]) & d.predicted_ancestry.eq("EUR")].drop_duplicates(["source", "sample_id"])
parts = []
for src in ["spark_gsa_v1", "spark_gsa_v2", "spark_twist"]:
    g = d[d.source == src].copy()
    for f in ["F1", "F2"]:
        s = pd.read_csv(f"{SC}/{src}.warrier_{f}.sscore", sep="\t", dtype={"#IID": str, "IID": str})
        ic = "#IID" if "#IID" in s.columns else "IID"; sc = [c for c in s.columns if c.endswith("_SUM") or c.endswith("_AVG")][-1]
        g[f] = g.IID.map(s.set_index(ic)[sc])
    parts.append(g)
d = pd.concat(parts)
for p in ["F1", "F2"]:
    ref = d[~d.case_status.eq("case")].groupby("source")[p].agg(["mean", "std"]); d[p] = (d[p] - d.source.map(ref["mean"])) / d.source.map(ref["std"])
d["F2mF1"] = d.F2 - d.F1
ref = d[~d.case_status.eq("case")].groupby("source")["F2mF1"].agg(["mean", "std"]); d["F2mF1"] = (d.F2mF1 - d.source.map(ref["mean"])) / d.source.map(ref["std"])
ped = pd.read_csv(f"{R}/pheno/spark_pedigree_ids.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
out = []
for src, g in d.groupby("source"):
    ids = set(g.sample_id); g = g.copy(); by = g.set_index("sample_id")
    mo, fa = g.sample_id.map(ped.biomother_sp_id), g.sample_id.map(ped.biofather_sp_id); ok = mo.isin(ids) & fa.isin(ids)
    for p in ["F1", "F2", "F2mF1"]:
        mid = (mo.map(by[p]) + fa.map(by[p])) / 2; g[p + "_ptdt"] = np.where(ok, (g[p] - mid) / np.nanstd(mid[ok]), np.nan)
    out.append(g)
d = pd.concat(out)
c = pd.read_parquet(f"{OUT}/cases.parquet")
c = c.merge(d[["source", "sample_id", "F1", "F2", "F2mF1", "F1_ptdt", "F2_ptdt", "F2mF1_ptdt"]], on=["source", "sample_id"], how="inner")
print("cases with factor PGS:", len(c), flush=True)
c.to_parquet(f"{OUT}/cases_f12.parquet", index=False)  # participant-level; stays on Expanse
src = open(os.path.expanduser("~/rvh2/driver_pgs.py")).read()
src = src.replace('PGS = ["asd_2026", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]', 'PGS = ["F1", "F2", "F2mF1"]')
src = src.replace('c = pd.read_parquet(f"{OUT}/cases.parquet")', 'c = pd.read_parquet(f"{OUT}/cases_f12.parquet")')
src = src.replace('f"{OUT}/driver_', 'f"{OUT}/f12_driver_')
assert 'PGS = ["F1", "F2", "F2mF1"]' in src and "cases_f12" in src and "f12_driver_summary" in src
exec(compile(src, "driver_pgs_f12", "exec"))
D = pd.read_csv(f"{OUT}/f12_driver_profiles.tsv", sep="\t")
pd.set_option("display.width", 250)
for meas in ["population", "within-family pTDT"]:
    x = D[(D.measure == meas) & (D.pgs == "F2mF1")].copy(); x["chi2"] = x.z ** 2; x = x.sort_values("chi2", ascending=False); x["rank"] = range(1, len(x) + 1)
    print("\nF2 - F1 contrast,", meas); print(x[["split", "n_high", "n_low", "diff", "se", "z", "rank"]].round(3).to_string(index=False))
