#!/usr/bin/env python3
"""Between-family rare-variant signature on the EXPANDED data : rare_between.py redone with
  probands = expanded trio probands (rare_full_cases.parquet: legacy SPARK trios without sibling requirement, SPARK-WGS-2026_08 and SSC
             sibling-family trios), registered before 22 (SSC: all) -- main; all ages as sensitivity,
  siblings = unaffected siblings of the sibling-family cohorts (~/rvh2/out_dnm_v3, incl. SSC), classes with the same weights and the same
             within-cohort sibling standardization.
Splits (SPARK survey codes, family-level value of the autistic child) for probands and siblings; parental age (per 10 years) for SPARK
(spark_env) and SSC (harmonised SFARI_phenotypes table). SNV/indel only; de novo = SNV/indel. Output (aggregate):
~/rvh2/agedx/rare_between_full{_all}.tsv, rare_between_full_parentage{_all}.tsv"""
import os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
H = os.path.expanduser("~/rvh2"); OUT = H + "/agedx"; P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
U22 = os.environ.get("U22", "1") == "1"; TAG = "" if U22 else "_all"
TIERS = ["lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4"]
cell = lambda modes, sets: [f"{t}__{m}__{g}" for t in TIERS for m in modes for g in sets]
CAT = {"de novo, NDD genes": cell(["dnm"], ["ndd"]), "de novo, other genes": cell(["dnm"], ["nonndd"]),
       "inherited, NDD genes": cell(["inh"], ["ndd"]), "inherited, other genes": cell(["inh"], ["nonndd"])}
W = pd.read_csv("/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift_full/results/collider_weights.tsv", sep="\t", index_col=0).iloc[:, 0]
pro = pd.read_parquet(f"{OUT}/rare_full_cases.parquet")
if U22: pro = pro[pro.u22]
sibs = pd.concat([pd.read_parquet(f"{H}/out_dnm_v3/{c}.parquet").query("case == 0").assign(cohort=c) for c in ["spark_wgs", "spark_wes_twist", "spark_wes_idt", "ssc", "spark_wgs_new"]], ignore_index=True)
meta = pd.read_csv("/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv", sep="\t", dtype=str).set_index("sample_id")
sibs["sp"] = np.where(sibs.cohort.eq("spark_wgs_new"), sibs["sample"].map(meta.ind_id), sibs["sample"])
for k, ts in CAT.items():
    ts = [t for t in ts if t in W.index]
    sibs[k] = (sibs[ts].fillna(0).to_numpy(float) * W[ts].to_numpy()).sum(1)
    ref = sibs.groupby("cohort")[k].agg(["mean", "std"]); sibs[k] = (sibs[k] - sibs.cohort.map(ref["mean"])) / sibs.cohort.map(ref["std"])
pro["who"] = "probands"; sibs["who"] = "unaffected siblings"
cols = ["sp", "cohort", "family_key", "who"] + list(CAT)
df = pd.concat([pro[cols], sibs[cols]], ignore_index=True)
df["fam"] = df.cohort + "|" + df.family_key.astype(str)
sv = pd.read_csv(f"{P}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
env = pd.read_csv(f"{P}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
hz = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/SFARI_phenotypes/spark_ssc_phenotypes_child.csv", low_memory=False)
hz = hz[hz.cohort.str.upper().eq("SSC")]
ssc_ids = pd.read_csv(f"{P}/ssc_ids_sanders.csv", dtype=str, encoding="utf-8-sig").drop_duplicates("SS_RUID")
idmap = ssc_ids.set_index("SS_RUID").individual   # SS_RUID -> <family>.p1
def hzcol(c):
    s = hz.set_index("individual_id")[c] if "individual_id" in hz.columns else None
    return s
SPL = {"age at diagnosis (>10 vs <6 y)": "s00", "household income (low vs high third)": "s09", "mother no degree vs degree": "s10", "maternal age (older vs younger third)": "s11",
       "paternal age (older vs younger third)": "s12", "deprivation (high vs low third)": "s08", "ADHD (yes vs no)": "s27", "intellectual disability": "s13"}
rows, prow = [], []
fi = pd.factorize(df.fam)[0]; nu = fi.max() + 1; rng = np.random.default_rng(11); BW = [np.bincount(rng.integers(0, nu, nu), minlength=nu)[fi] for _ in range(500)]
for nm, code in SPL.items():
    cs = df[df.who.eq("probands")].assign(g=lambda d: pd.to_numeric(d.sp.map(sv[code]), errors="coerce"))
    fg = cs.dropna(subset=["g"]).groupby("fam").g.first(); g = df.fam.map(fg).to_numpy(float)
    for who in ["probands", "unaffected siblings"]:
        m = df.who.eq(who).to_numpy()
        for col in CAT:
            v = df[col].to_numpy(float); ok = m & np.isfinite(v) & np.isfinite(g); h, l = ok & (g == 1), ok & (g == 0)
            if h.sum() < 30 or l.sum() < 30: continue
            d0 = v[h].mean() - v[l].mean(); bs = [np.average(v[h], weights=w[h]) - np.average(v[l], weights=w[l]) for w in BW if w[h].sum() and w[l].sum()]
            se = float(np.std(bs)); rows.append(dict(split=nm, who=who, category=col, n_high=int(h.sum()), n_low=int(l.sum()), diff=d0, se=se, z=d0 / se))
# parental age per 10 years: SPARK env; SSC harmonised (keyed by individual_id = Sanders 'individual')
ma = pd.to_numeric(df.sp.map(env.maternal_age), errors="coerce"); pa = pd.to_numeric(df.sp.map(env.paternal_age), errors="coerce")
hzi = hz.dropna(subset=["repository_id"]).drop_duplicates("repository_id").set_index("repository_id")   # repository_id = SS_RUID, as in our SSC tables
ma = ma.fillna(pd.to_numeric(df.sp.map(hzi.mother_age_at_birth), errors="coerce")); pa = pa.fillna(pd.to_numeric(df.sp.map(hzi.father_age_at_birth), errors="coerce"))
df["maternal_age"] = ma / 10; df["paternal_age"] = pa / 10
for par in ["maternal_age", "paternal_age"]:
    for who in ["probands", "unaffected siblings"]:
        for col in CAT:
            x = df[df.who.eq(who)].dropna(subset=[par, col])
            if len(x) < 200: continue
            f = smf.ols(f"Q('{col}') ~ {par} + C(cohort)", data=x).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(x.fam)[0]})
            prow.append(dict(parent=par, who=who, category=col, n=len(x), n_ssc=int(x.cohort.eq("ssc").sum()), beta_per_10y=f.params[par], se=f.bse[par], z=f.tvalues[par], p=f.pvalues[par]))
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/rare_between_full{TAG}.tsv", sep="\t", index=False)
Q = pd.DataFrame(prow); Q.to_csv(f"{OUT}/rare_between_full_parentage{TAG}.tsv", sep="\t", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_rows", 500)
print("probands", int(df.who.eq("probands").sum()), "| siblings", int(df.who.eq("unaffected siblings").sum()), df.groupby(["who", "cohort"]).size().to_dict())
for who in ["probands", "unaffected siblings"]:
    print(f"\n== between-family z, {who}"); print(D[D.who == who].pivot_table(index="split", columns="category", values="z").round(1).to_string())
print("\n== parental age (per 10 y)"); print(Q.round(4).to_string(index=False))
