#!/usr/bin/env python3
"""Matched-design GWAS, step 2 (per platform, SRC=<source>): phenotype/covariate files for 4 splits x 2 halves x 2 designs.
Splits of autistic cases: age at diagnosis (after 10 / before 6), area deprivation (top / bottom third, thresholds over all
SPARK EUR cases), adult at registration (yes / no), ADHD comorbidity (yes / no).
  between-family: cases of the half that are in the KING-unrelated set vs non-autistic members of the unrelated set
                  (same controls for both halves); covariates PC1-10 (computed in the unrelated set) + sex.
  trio (TDT):     affected offspring of the half with both parents genotyped on this platform; fam file with pedigree links.
Output: mgwas/<src>/pheno_<split>_<half>.txt (between-family), tdt_<split>_<half>.pheno, covar.txt, trio.fam"""
import os
import numpy as np, pandas as pd
SRC = os.environ["SRC"]; W = f"/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/{SRC}"; R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
B = "/expanse/projects/sebat1/s3/data/sebat"
d = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "sample_id", "case_status"])
d = d[d.source == SRC].drop_duplicates("IID")
c = pd.read_parquet(os.path.expanduser("~/rvh2/agedx/cases.parquet"))
allc = c.copy()
ho = pd.read_csv(f"{R}/spark_heldout.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
env = pd.read_csv(f"{R}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
allc["adhd"] = pd.to_numeric(allc.sample_id.map(ho.behav_adhd), errors="coerce"); allc["adi"] = pd.to_numeric(allc.sample_id.map(env.adi_pct), errors="coerce")
lo, hi = allc.adi.quantile(1 / 3), allc.adi.quantile(2 / 3)
c = allc[allc.source == SRC].copy()
SPL = {"agedx": {"late": c.stratum.eq("after 10"), "early": c.stratum.eq("before 6")},
       "adi": {"high": c.adi >= hi, "low": c.adi <= lo},
       "adult": {"adult": c.adult_at_reg.eq(1), "child": c.adult_at_reg.eq(0)},
       "adhd": {"yes": c.adhd.eq(1), "no": c.adhd.eq(0)}}
unrel = set(pd.read_csv(f"{W}/king.king.cutoff.in.id", sep="\t", dtype=str).iloc[:, 0])
pcs = pd.read_csv(f"{W}/pca_unrel.eigenvec", sep="\t", dtype={"#IID": str}).rename(columns={"#IID": "IID"})
ped = pd.read_csv(f"{B}/SPARK_iWES_v3/SPARK_iWES_v3.ped", sep="\t", header=None, dtype=str)
ped2 = pd.read_csv(f"{B}/SPARK_iWGS_v1.1/SPARK_iWGS_v1.1.ped", sep="\t", header=None, dtype=str)
ped = pd.concat([ped, ped2]).drop_duplicates(1).set_index(1)
d["sex"] = d.sample_id.map(ped[4]).map({"1": 1, "2": 2})
cov = pcs.merge(d[["IID", "sex"]], on="IID", how="left"); cov["sex"] = cov.sex.fillna(0).astype(int)
cov[["IID", "sex"] + [f"PC{i}" for i in range(1, 11)]].to_csv(f"{W}/covar.txt", sep="\t", index=False)
controls = set(d[(d.case_status != "case") & d.IID.isin(unrel)].IID)
# trio fam: FID = SF family, parents = genotyped parents' IIDs on this platform
fam = pd.read_csv(f"{W}/hm3_bed.fam", sep=r"\s+", header=None, dtype=str); fam.columns = ["FID", "IID", "PAT", "MAT", "SEX", "PHE"]
s2i = d.set_index("sample_id").IID; i2s = d.set_index("IID").sample_id
sp = fam.IID.map(i2s)
fam["FID"] = sp.map(ped[0]).fillna(fam.IID)
fam["PAT"] = sp.map(ped[2]).map(s2i).fillna("0"); fam["MAT"] = sp.map(ped[3]).map(s2i).fillna("0")
fam["SEX"] = sp.map(ped[4]).fillna("0"); fam["PHE"] = "-9"
fam.to_csv(f"{W}/trio.fam", sep=" ", index=False, header=False)
trio_kids = set(fam[(fam.PAT != "0") & (fam.MAT != "0")].IID)
rows = []
for sp_, halves in SPL.items():
    for h, m in halves.items():
        ids = set(c[m.fillna(False)].IID)
        cas = ids & unrel
        ph = pd.DataFrame({"IID": sorted(cas | controls)}); ph["Y"] = np.where(ph.IID.isin(cas), 2, 1)
        ph.to_csv(f"{W}/pheno_{sp_}_{h}.txt", sep="\t", index=False)
        tk = ids & trio_kids
        tp = fam[["FID", "IID"]].copy(); tp["Y"] = np.where(tp.IID.isin(tk), 2, 1)
        tp.to_csv(f"{W}/tdt_{sp_}_{h}.pheno", sep=" ", index=False, header=False)
        rows.append(dict(source=SRC, split=sp_, half=h, cases_in_half=len(ids), between_cases=len(cas), between_controls=len(controls), trio_offspring=len(tk)))
S = pd.DataFrame(rows); S.to_csv(f"{W}/design_counts.tsv", sep="\t", index=False); print(S.to_string(index=False))
