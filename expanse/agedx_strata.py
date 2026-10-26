#!/usr/bin/env python3
"""SPARK autistic individuals split as in Zhang, Grove et al. Nature 2025 (age at diagnosis: before 6 / 6-10 / after 10;
'before 6' and 'after 10' define their Factors 1 and 2). For each stratum:
  composition: genotyping platform, number of biological parents genotyped on the same platform (0/1/2), adult at
  registration (>= 18 y), diagnosis era, sex, ID, minimal language, ADHD and mood/anxiety comorbidity, age at diagnosis;
  PGS: (a) population-style mean of each PGS in cases, standardized to non-autistic parents + siblings of the same platform
  (the comparison their stratified GWAS make: cases vs unscreened relatives); (b) within-family pTDT in complete trios
  (child minus mid-parent, immune to PCs, family structure and population stratification).
PGS (pgs_results_v2): ASD 2026 leave-one-out, ADHD 2022, MDD 2025, bipolar 2024, SCZ 2022, EA 2022, IQ 2018, AFB 2018,
height (control). European ancestry only. Output (aggregate): ~/rvh2/agedx/*.tsv"""
import os
import numpy as np, pandas as pd
R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift"; OUT = os.path.expanduser("~/rvh2/agedx"); os.makedirs(OUT, exist_ok=True)
PGS = ["asd_2026", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
d = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype={"sample_id": str, "FID": str, "IID": str}, low_memory=False)
d = d[d.source.str.startswith("spark")].copy()
d = d.merge(pd.read_csv("/expanse/projects/sebat1/jsebat/pgs_asd2026_loo/asd2026_loo_scores.tsv", sep="\t", dtype={"IID": str})[["source", "IID", "asd_2026"]], on=["source", "IID"], how="left")
d = d[d.predicted_ancestry.eq("EUR")].drop_duplicates(["source", "sample_id"])
ped = pd.read_csv(f"{R}/pheno/spark_pedigree_ids.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
dx = pd.read_csv(f"{R}/pheno/spark_dx_derived.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
sev = pd.read_csv(f"{R}/pheno/spark_severity.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
ho = pd.read_csv(f"{R}/pheno/spark_heldout.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
print("dx columns:", list(dx.columns)[:40]); print("heldout columns:", [c for c in ho.columns][:40], flush=True)
# standardize PGS within platform to NON-autistic individuals (parents + siblings), like a case vs relative-control contrast
for p in PGS:
    ref = d[~d.case_status.eq("case")].groupby("source")[p].agg(["mean", "std"])
    d[p] = (d[p] - d.source.map(ref["mean"])) / d.source.map(ref["std"])
parts = []
for src, g in d.groupby("source"):
    ids = set(g.sample_id); g = g.copy()
    g["mo"] = g.sample_id.map(ped.biomother_sp_id); g["fa"] = g.sample_id.map(ped.biofather_sp_id)
    g["n_par"] = g.mo.isin(ids).astype(int) + g.fa.isin(ids).astype(int)
    by = g.set_index("sample_id")
    for p in PGS:
        ok = g.mo.isin(ids) & g.fa.isin(ids); mid = (g.mo.map(by[p]) + g.fa.map(by[p])) / 2
        g[p + "_ptdt"] = np.where(ok, (g[p] - mid) / np.nanstd(mid[ok]), np.nan)
    parts.append(g)
d = pd.concat(parts)
c = d[d.case_status.eq("case")].copy()
c["dx_age"] = c.sample_id.map(dx.dx_age_y); c["dx_year"] = c.sample_id.map(dx.dx_year)
for col in ["ID", "sex", "validity_flag"]:
    if col in dx.columns: c[col] = c.sample_id.map(dx[col])
for col in ["minimal_language", "cog_impair", "fsiq"]:
    if col in sev.columns: c[col] = c.sample_id.map(sev[col])
for col in ["behav_adhd_ageadj", "mood_or_anx_ageadj"]:
    if col in ho.columns: c[col] = c.sample_id.map(ho[col])
if "age_at_registration" in dx.columns: c["adult_at_reg"] = (c.sample_id.map(dx.age_at_registration) >= 18).astype(float)
elif "reg_age_y" in dx.columns: c["adult_at_reg"] = (c.sample_id.map(dx.reg_age_y) >= 18).astype(float)
c = c[c.dx_age.notna()]
c["stratum"] = pd.cut(c.dx_age, [-1, 6, 11, 200], labels=["before 6", "6-10", "after 10"], right=False).astype(str)
c["era"] = pd.cut(c.dx_year, [1950, 2010, 2015, 2019, 2027], labels=["<2010", "2010-14", "2015-18", "2019-26"], right=False).astype(str)
c["female"] = (c.get("sex", pd.Series(index=c.index, dtype=object)) == "Female").astype(float)
ST = ["before 6", "6-10", "after 10"]
c.to_parquet(f"{OUT}/cases.parquet", index=False)   # participant-level; stays on Expanse
comp = []
for st in ST:
    x = c[c.stratum == st]
    r = dict(stratum=st, n=len(x), median_dx_age=x.dx_age.median(), pct_0_parents=(x.n_par == 0).mean() * 100, pct_1_parent=(x.n_par == 1).mean() * 100, pct_2_parents=(x.n_par == 2).mean() * 100,
             pct_adult_at_reg=x.adult_at_reg.mean() * 100 if "adult_at_reg" in x else np.nan, female_pct=x.female.mean() * 100)
    for e in ["<2010", "2010-14", "2015-18", "2019-26"]: r[f"pct_dx_{e}"] = (x.era == e).mean() * 100
    for s in ["spark_gsa_v1", "spark_gsa_v2", "spark_twist", "spark_wgs"]: r[f"pct_{s}"] = (x.source == s).mean() * 100
    for col in ["ID", "minimal_language", "cog_impair", "validity_flag"]:
        if col in x: r[f"pct_{col}"] = pd.to_numeric(x[col], errors="coerce").mean() * 100
    for col in ["behav_adhd_ageadj", "mood_or_anx_ageadj", "fsiq"]:
        if col in x: r[f"mean_{col}"] = pd.to_numeric(x[col], errors="coerce").mean()
    comp.append(r)
C = pd.DataFrame(comp).set_index("stratum").T; C.to_csv(f"{OUT}/composition.tsv", sep="\t")
pd.set_option("display.width", 220); print(C.round(2).to_string(), flush=True)
def bm(v, fam, B=400, seed=1):
    ok = np.isfinite(v); v, fam = v[ok], fam[ok]
    if len(v) < 30: return np.nan, np.nan, len(v)
    fi, uf = pd.factorize(fam); rng = np.random.default_rng(seed); bs = []
    for _ in range(B):
        w = np.bincount(rng.integers(0, len(uf), len(uf)), minlength=len(uf))[fi]; bs.append(np.average(v, weights=w))
    return float(v.mean()), float(np.std(bs)), int(len(v))
rows = []
for st in ST:
    for grp, x in [("all cases", c[c.stratum == st]), ("cases with 0 parents genotyped", c[(c.stratum == st) & (c.n_par == 0)]), ("cases with 2 parents genotyped", c[(c.stratum == st) & (c.n_par == 2)])]:
        for p in PGS:
            m, se, n = bm(x[p].to_numpy(float), x.FID.astype(str).to_numpy()); rows.append(dict(stratum=st, group=grp, measure="population mean (SD vs non-autistic relatives)", pgs=p, mean=m, se=se, n=n))
    x = c[(c.stratum == st) & (c.n_par == 2)]
    for p in PGS:
        m, se, n = bm(x[p + "_ptdt"].to_numpy(float), x.FID.astype(str).to_numpy()); rows.append(dict(stratum=st, group="complete trios", measure="pTDT (within family)", pgs=p, mean=m, se=se, n=n))
P = pd.DataFrame(rows); P.to_csv(f"{OUT}/pgs_by_stratum.tsv", sep="\t", index=False)
for meas, grp in [("population mean (SD vs non-autistic relatives)", "all cases"), ("population mean (SD vs non-autistic relatives)", "cases with 0 parents genotyped"),
                  ("population mean (SD vs non-autistic relatives)", "cases with 2 parents genotyped"), ("pTDT (within family)", "complete trios")]:
    q = P[(P.measure == meas) & (P.group == grp)]
    print(f"\n{meas} | {grp}"); print(q.pivot(index="pgs", columns="stratum", values="mean").reindex(columns=ST).round(3).to_string())
    print("n:", q.groupby("stratum").n.max().reindex(ST).to_dict())
