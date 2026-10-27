#!/usr/bin/env python3
"""(All-data version: also SRC = spark_wgs, spark_wgs_2026_08, one row per person across strata.)
SNP heritability of the case characteristics we split on, estimated the way Zhang, Grove et al. estimated it for age at
diagnosis (GCTA-GREML among autistic participants; qcovar 10 PCs, covar sex), per platform (SRC). Unrelated autistic SPARK EUR
participants (KING set from cnv_drift/mgwas; GRM cut-off 0.05). Traits:
  dx_age        age at diagnosis in years, standardized (their trait); dx_age_u22 = same, registrants < 22 y (their restriction)
  adi_pct       area deprivation index national percentile (standardized)
  adult_reg     registered as an adult (>= 18 y)            [binary; liability scale with K = sample prevalence]
  adhd          parent/self-reported ADHD                    [binary]
  late_dx       diagnosed > 10 y vs < 6 y (their strata)     [binary]
Output: warrier_rep/<src>/greml_pheno.txt, greml_names.txt, greml_qcovar.txt, greml_covar.txt, greml_prev.txt"""
import os
import numpy as np, pandas as pd
SRC = os.environ["SRC"]; W = f"/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/{SRC}"; M = f"/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/{SRC}"
R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"; B = "/expanse/projects/sebat1/s3/data/sebat"
fam = pd.read_csv(f"{W}/hm3w.fam", sep=r"\s+", header=None, dtype=str).iloc[:, :2]; fam.columns = ["FID", "IID"]
d = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "sample_id"])
arrays = set(d[d.source.isin(["spark_gsa_v1", "spark_gsa_v2", "spark_twist"])].sample_id); iwgs = set(d[d.source.eq("spark_wgs")].sample_id)
if SRC == "spark_wgs_2026_08":   # SDSM barcodes -> SP IDs; one row per person across strata (exclude anyone on an array or iWGS v1.1)
    meta = pd.read_csv("/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv", sep="	", dtype=str).set_index("sample_id")
    fam["sp"] = fam.IID.map(meta.ind_id); fam.loc[fam.sp.isin(arrays | iwgs), "sp"] = np.nan
else:
    d = d[d.source == SRC].drop_duplicates("IID").set_index("IID"); fam["sp"] = fam.IID.map(d.sample_id)
    if SRC == "spark_wgs": fam.loc[fam.sp.isin(arrays), "sp"] = np.nan
reg = pd.read_csv(f"{R}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
ho = pd.read_csv(f"{R}/spark_heldout.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
env = pd.read_csv(f"{R}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
sv = pd.read_csv(f"{R}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
unrel = set(pd.read_csv(f"{M}/king.king.cutoff.in.id", sep="\t", dtype=str).iloc[:, 0])
f = fam.copy(); f["asd"] = f.sp.map(reg.asd).astype(str).str.upper().eq("TRUE")
f = f[f.asd & f.IID.isin(unrel)].copy()
dxm = pd.to_numeric(f.sp.map(reg.diagnosis_age), errors="coerce"); areg = pd.to_numeric(f.sp.map(reg.age_at_registration_years), errors="coerce")
z = lambda v: (v - v.mean()) / v.std()
T = {"dx_age": z(dxm / 12), "dx_age_u22": z((dxm / 12).where(areg < 22)), "adi_pct": z(pd.to_numeric(f.sp.map(env.adi_pct), errors="coerce")),
     "adult_reg": (areg >= 18).astype(float).where(areg.notna()), "adhd": pd.to_numeric(f.sp.map(ho.behav_adhd), errors="coerce"),
     "late_dx": pd.Series(np.where(dxm > 120, 1.0, np.where(dxm < 72, 0.0, np.nan)), index=f.index),
     # 2026-10-26: the SES stratifiers named in the text, defined as their split GWAS (spark_survey codes: s09 low vs high income third, s10 mother without degree)
     "mat_age": z(pd.to_numeric(f.sp.map(env.maternal_age), errors="coerce")), "pat_age": z(pd.to_numeric(f.sp.map(env.paternal_age), errors="coerce")),
     "low_income": pd.to_numeric(f.sp.map(sv.s09), errors="coerce"), "mother_nodeg": pd.to_numeric(f.sp.map(sv.s10), errors="coerce")}
P = f[["FID", "IID"]].copy()
for k, v in T.items(): P[k] = v.values
P.fillna("NA").to_csv(f"{W}/greml_pheno.txt", sep=" ", index=False, header=False)
open(f"{W}/greml_names.txt", "w").write("\n".join(T) + "\n")
prev = {k: float(np.nanmean(v)) for k, v in T.items() if k in ("adult_reg", "adhd", "late_dx", "low_income", "mother_nodeg")}
pd.Series(prev).to_csv(f"{W}/greml_prev.txt", sep="\t", header=False)
f[["FID", "IID"]].to_csv(f"{W}/greml_keep.txt", sep=" ", index=False, header=False)
pc = pd.read_csv(f"{W}/pcair.sscore", sep="\t"); ic = [c for c in pc.columns if "IID" in c][0]; pcs = [c for c in pc.columns if c.endswith("_AVG")][:10]
f[["FID", "IID"]].merge(pc[[ic] + pcs].rename(columns={ic: "IID"}), on="IID", how="left").fillna("NA").to_csv(f"{W}/greml_qcovar.txt", sep=" ", index=False, header=False)
c = f[["FID", "IID"]].copy(); c["sex"] = f.sp.map(reg.sex).values; c.fillna("NA").to_csv(f"{W}/greml_covar.txt", sep=" ", index=False, header=False)
print(SRC, "unrelated cases", len(f), {k: int(v.notna().sum()) for k, v in T.items()}, "prev", {k: round(v, 3) for k, v in prev.items()})
