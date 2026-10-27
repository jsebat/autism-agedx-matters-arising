#!/usr/bin/env python3
"""Supplementary table: within-family associations of the genetic predictors with autism (direct effects; indirect genetic effects,
population stratification and participation cancel within families).
(a) Polygenic scores: pTDT (child minus mean of both parents' scores, in s.d. of the mid-parent score within genotype set) for each of the nine
    scores, in autistic children and in their unaffected siblings (both parents genotyped, European ancestry, SPARK and SSC); difference
    autistic minus unaffected from OLS pTDT ~ autistic with family-clustered standard errors.
(b) Rare variants: autistic minus unaffected sibling in the same family, for each burden class (sibling-family cohorts: SPARK WES IDT / Twist,
    iWGS v1.1, SPARK-WGS-2026_08, SSC; out_dnm_v3), classes weighted and standardized to siblings within cohort as in rare_full.py;
    mean of within-family differences, family bootstrap 95% CI (2,000) and paired t-test.
Output (aggregate): ~/rvh2/agedx/within_family_table_{pgs,rare}.tsv"""
import os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import ttest_1samp
OUT = os.path.expanduser("~/rvh2/agedx"); H = os.path.expanduser("~/rvh2")
A = pd.read_parquet(f"{OUT}/full_probands.parquet").drop_duplicates("sample_id")
PG = ["asd", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
# SES polygenic score (not one of the nine CCA scores): occupational status (CAMSIS; Akimova et al. 2025), scored by score_occstatus.sbatch on the same
# genotype sets; each person takes the score from their own source and is standardized within source to non-autistic participants, as for the nine.
OCC = "/expanse/projects/sebat1/jsebat/pgs_occstatus"
if os.path.exists(f"{OCC}/spark_gsa_v1.occstatus.sscore"):
    ids = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, usecols=["source", "IID", "sample_id"]).drop_duplicates(["source", "IID"])
    meta = pd.read_csv("/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv", sep="\t", dtype=str).set_index("sample_id").ind_id
    parts = []
    for src in A.source.unique():
        f = f"{OCC}/{src}.occstatus.sscore"
        if not os.path.exists(f): continue
        sc = pd.read_csv(f, sep="\t", dtype=str); ic = [c for c in sc.columns if c.lstrip("#") == "IID"][0]
        sc["v"] = pd.to_numeric(sc["SCORE1_AVG"], errors="coerce")
        m = sc[ic].map(meta) if src == "spark_wgs_2026_08" else sc[ic].map(ids[ids.source.eq(src)].set_index("IID").sample_id)
        parts.append(pd.DataFrame(dict(source=src, sample_id=m, occ=sc.v)).dropna())
    O = pd.concat(parts).drop_duplicates(["source", "sample_id"])
    A = A.merge(O, on=["source", "sample_id"], how="left")
    ref = A[~A.case.astype(bool)].groupby("source").occ.agg(["mean", "std"])
    A["occstatus"] = (A.occ - A.source.map(ref["mean"])) / A.source.map(ref["std"]); PG = PG + ["occstatus"]
    print("occupational-status score coverage by source:", A.occstatus.notna().groupby(A.source).mean().round(3).to_dict(), flush=True)
allA = A.set_index("sample_id")
c = A[A.predicted_ancestry.eq("EUR") & A.n_parents_genotyped.eq(2) & A.mother.notna() & A.father.notna()].copy()
c["autistic"] = c.case.astype(float); c["fam"] = pd.factorize(c.family.fillna(c.sample_id))[0]
print("children with both parents genotyped:", c.groupby(["cohort", "autistic"]).size().to_dict(), flush=True)
rows = []
for p in PG:
    mid = (c.mother.map(allA[p]) + c.father.map(allA[p])) / 2; sd = mid.groupby(c.source).transform("std")
    c["y"] = (c[p] - mid) / sd; x = c.dropna(subset=["y"])
    for who, d in [("autistic", x[x.autistic == 1]), ("unaffected siblings", x[x.autistic == 0])]:
        m = smf.ols("y ~ 1", data=d).fit(cov_type="cluster", cov_kwds={"groups": d.fam.to_numpy()})
        rows.append(dict(pgs=p, group=who, n=len(d), ptdt=m.params.iloc[0], se=m.bse.iloc[0], p=m.pvalues.iloc[0]))
    m = smf.ols("y ~ autistic", data=x).fit(cov_type="cluster", cov_kwds={"groups": x.fam.to_numpy()})
    rows.append(dict(pgs=p, group="autistic minus unaffected", n=len(x), ptdt=m.params["autistic"], se=m.bse["autistic"], p=m.pvalues["autistic"]))
T = pd.DataFrame(rows); T["lo95"] = T.ptdt - 1.96 * T.se; T["hi95"] = T.ptdt + 1.96 * T.se
T.to_csv(f"{OUT}/within_family_table_pgs.tsv", sep="\t", index=False); pd.set_option("display.width", 200); print(T.round(4).to_string(index=False), flush=True)

# ---- (b) rare variants, within sibling families
TIERS = ["lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4"]
cell = lambda modes, sets: [f"{t}__{m}__{g}" for t in TIERS for m in modes for g in sets]
CAT = {"de novo, NDD genes": cell(["dnm"], ["ndd"]), "de novo, other genes": cell(["dnm"], ["nonndd"]),
       "inherited, NDD genes": cell(["inh"], ["ndd"]), "inherited, other genes": cell(["inh"], ["nonndd"]),
       "inherited LoF tier 1": [f"lof_t1__inh__{g}" for g in ("ndd", "nonndd")]}
W = pd.read_csv("/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift_full/results/collider_weights.tsv", sep="\t", index_col=0).iloc[:, 0]
D = pd.concat([pd.read_parquet(f"{H}/out_dnm_v3/{k}.parquet").assign(cohort=k) for k in ["spark_wgs", "spark_wes_twist", "spark_wes_idt", "ssc", "spark_wgs_new"]], ignore_index=True)
for k, ts in CAT.items():
    ts = [t for t in ts if t in W.index]; D[k] = (D[ts].fillna(0).to_numpy(float) * W[ts].to_numpy()).sum(1)
    ref = D[D.case.eq(0)].groupby("cohort")[k].agg(["mean", "std"]); D[k] = (D[k] - D.cohort.map(ref["mean"])) / D.cohort.map(ref["std"])
D["fk"] = D.cohort + "|" + D.family_key.astype(str)
g = D.groupby(["fk", "case"])[list(CAT)].mean().unstack("case").dropna()
rng = np.random.default_rng(7); rows = []
for k in CAT:
    d = (g[(k, 1)] - g[(k, 0)]).to_numpy(); bs = [rng.choice(d, len(d)).mean() for _ in range(2000)]
    rows.append(dict(cls=k, n_families=len(d), diff=d.mean(), lo95=np.percentile(bs, 2.5), hi95=np.percentile(bs, 97.5), p=ttest_1samp(d, 0).pvalue))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/within_family_table_rare.tsv", sep="\t", index=False); print(R.round(4).to_string(index=False))
print("families by cohort:", pd.Series([f.split("|")[0] for f in g.index]).value_counts().to_dict())
