#!/usr/bin/env python3
"""Matched-design GWAS, round 2 (per platform, SRC=<source>): splits made ORTHOGONAL to age at diagnosis / adult registration, to
test whether a 'second factor' needs late diagnosis at all. Between-family design only (same unrelated controls and covariates
as round 1):
  adiE   area deprivation top vs bottom third, among cases diagnosed before 6
  adhdE  ADHD yes vs no, among cases diagnosed before 6
  sexE   female vs male, among cases diagnosed before 6
  agedxC diagnosed after 10 vs before 6, among cases registered as children (< 18 y)
Also (SRC=spark_gsa_v1 only) phi correlations between all split indicators across SPARK EUR cases (aggregate).
Output: mgwas/<src>/pheno2_<split>_<half>.txt; rows appended to design_counts.tsv (round-1 rows kept)."""
import os
import numpy as np, pandas as pd
SRC = os.environ["SRC"]; W = f"/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/{SRC}"; R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
B = "/expanse/projects/sebat1/s3/data/sebat"
d = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "sample_id", "case_status"])
d = d[d.source == SRC].drop_duplicates("IID")
allc = pd.read_parquet(os.path.expanduser("~/rvh2/agedx/cases.parquet"))
ho = pd.read_csv(f"{R}/spark_heldout.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
env = pd.read_csv(f"{R}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
allc["adhd"] = pd.to_numeric(allc.sample_id.map(ho.behav_adhd), errors="coerce"); allc["adi"] = pd.to_numeric(allc.sample_id.map(env.adi_pct), errors="coerce")
lo, hi = allc.adi.quantile(1 / 3), allc.adi.quantile(2 / 3)


def ind(c):
    early = c.stratum.eq("before 6"); child = c.adult_at_reg.eq(0)
    return {"adiE": {"high": early & (c.adi >= hi), "low": early & (c.adi <= lo)},
            "adhdE": {"yes": early & c.adhd.eq(1), "no": early & c.adhd.eq(0)},
            "sexE": {"female": early & c.female.eq(1), "male": early & c.female.eq(0)},
            "agedxC": {"late": child & c.stratum.eq("after 10"), "early": child & c.stratum.eq("before 6")}}


if SRC == "spark_gsa_v1":
    x = pd.DataFrame({"late dx (vs before 6)": np.where(allc.stratum.eq("after 10"), 1.0, np.where(allc.stratum.eq("before 6"), 0.0, np.nan)),
                      "adult at registration": allc.adult_at_reg.astype(float), "high ADI (vs low third)": np.where(allc.adi >= hi, 1.0, np.where(allc.adi <= lo, 0.0, np.nan)),
                      "ADHD": allc.adhd, "female": allc.female.astype(float), "no parents genotyped (vs 2)": np.where(allc.n_par.eq(0), 1.0, np.where(allc.n_par.eq(2), 0.0, np.nan))})
    x.corr(min_periods=200).round(3).to_csv(os.path.expanduser("~/rvh2/agedx/split_phi.tsv"), sep="\t"); print(x.corr(min_periods=200).round(2).to_string())
    e = allc[allc.stratum.eq("after 10")]; print("late-dx cases: % adult reg", round(e.adult_at_reg.mean() * 100, 1), "% high-ADI third", round((e.adi >= hi).mean() * 100, 1), "% ADHD", round(e.adhd.mean() * 100, 1))
c = allc[allc.source == SRC].copy()
unrel = set(pd.read_csv(f"{W}/king.king.cutoff.in.id", sep="\t", dtype=str).iloc[:, 0])
controls = set(d[(d.case_status != "case") & d.IID.isin(unrel)].IID)
rows = []
for sp_, halves in ind(c).items():
    for h, m in halves.items():
        ids = set(c[m.fillna(False)].IID); cas = ids & unrel
        ph = pd.DataFrame({"IID": sorted(cas | controls)}); ph["Y"] = np.where(ph.IID.isin(cas), 2, 1)
        ph.to_csv(f"{W}/pheno2_{sp_}_{h}.txt", sep="\t", index=False)
        rows.append(dict(source=SRC, split=sp_, half=h, cases_in_half=len(ids), between_cases=len(cas), between_controls=len(controls), trio_offspring=0))
S = pd.DataFrame(rows)
old = pd.read_csv(f"{W}/design_counts.tsv", sep="\t"); old = old[~old.split.isin(S.split.unique())]
pd.concat([old, S]).to_csv(f"{W}/design_counts.tsv", sep="\t", index=False); print(S.to_string(index=False))
