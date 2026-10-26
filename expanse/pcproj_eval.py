#!/usr/bin/env python3
"""Projection artefact check (SPARK Twist, EUR). PCs computed in a KING-unrelated set and projected onto everyone (PC-AiR logic).
  1. Who is in the unrelated (PCA) set: by role, by number of genotyped parents, by age-at-diagnosis stratum.
  2. Shrinkage: SD of each projected PC in samples outside vs inside the PCA set (ratio < 1 = shrinkage toward 0).
  3. Does the artefact separate the GWAS groups? Mean and SD of PCs by case group (0/1/2 parents; before 6 / after 10) and
     in non-autistic relatives; and the PC-score difference after 10 vs before 6 among cases (in SD units of the PCA set).
Output: cnv_drift/pcproj/eval_*.tsv (aggregate)."""
import numpy as np, pandas as pd, os
W = "/expanse/projects/sebat1/jsebat/cnv_drift/pcproj"
s = pd.read_csv(f"{W}/proj_all.sscore", sep="\t"); idc = [c for c in s.columns if c.endswith("IID")][0]
pcs = [c for c in s.columns if c.endswith("_SUM") or c.startswith("PC")][:10]; s = s.rename(columns={idc: "IID"})
print("score columns:", pcs)
inset = pd.read_csv(f"{W}/king.king.cutoff.in.id", sep="\t", dtype=str); iic = [c for c in inset.columns if "IID" in c][0]
s["in_pca"] = s.IID.astype(str).isin(set(inset[iic]))
d = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "sample_id", "case_status"])
d = d[d.source == "spark_twist"].drop_duplicates("IID")
s = s.merge(d, on="IID", how="left")
c = pd.read_parquet(os.path.expanduser("~/rvh2/agedx/cases.parquet")); c = c[c.source == "spark_twist"][["IID", "n_par", "stratum"]]
s = s.merge(c, on="IID", how="left")
s["group"] = np.where(s.case_status.eq("case"), "case, " + s.stratum.fillna("no dx age") + ", " + s.n_par.astype("Int64").astype(str) + " parents", "non-autistic relative")
ref = s[s.in_pca][pcs].std()
for p in pcs: s[p] = s[p] / ref[p]          # PC units: SD in the PCA (unrelated) set
print("\n1. share in the PCA set:")
print(s.groupby("group").in_pca.agg(["mean", "size"]).round(3).to_string())
r2 = pd.DataFrame({"SD outside / SD inside PCA set": s[~s.in_pca][pcs].std() / s[s.in_pca][pcs].std()}).T.round(3)
print("\n2. shrinkage of projected PCs (1 = none):"); print(r2.to_string())
g = s.groupby("group")[pcs]; M = g.mean().round(3); SD = g.std().round(3)
print("\n3a. mean PC (in PCA-set SD units) by group:"); print(M.to_string()); print("\n3b. SD of PC by group:"); print(SD.to_string())
cc = s[s.case_status.eq("case") & s.stratum.isin(["before 6", "after 10"])]
diff = (cc[cc.stratum == "after 10"][pcs].mean() - cc[cc.stratum == "before 6"][pcs].mean()).round(3)
sdr = (cc[cc.stratum == "after 10"][pcs].std() / cc[cc.stratum == "before 6"][pcs].std()).round(3)
print("\n3c. cases after 10 minus before 6: mean PC difference"); print(diff.to_string()); print("SD ratio after10/before6"); print(sdr.to_string())
r2.to_csv(f"{W}/eval_shrinkage.tsv", sep="\t"); M.to_csv(f"{W}/eval_pc_means.tsv", sep="\t"); SD.to_csv(f"{W}/eval_pc_sd.tsv", sep="\t")
s.groupby("group").in_pca.agg(["mean", "size"]).to_csv(f"{W}/eval_inset.tsv", sep="\t")
