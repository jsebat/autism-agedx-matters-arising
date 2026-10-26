#!/usr/bin/env python3
"""'What holds' (pTDT) on the FULL SPARK data (arrays + iWGS v1.1 + SPARK-WGS-2026_08). Zhang et al. (Suppl. Table 24): ADHD (and
depression) polygenic scores are over-transmitted to later-diagnosed autistic people. Test whether diagnosis after 10 vs before 6 years predicts
pTDT once the child's ADHD (and adult registration, sex, diagnosis era) are in the model.
pTDT = child score minus mean of both parents' scores, in SD of the mid-parent score within platform. Probands of European ancestry with both
parents genotyped. Samples: registered before 22 (Zhang et al.'s criterion for within-case analyses; main) and all ages. OLS, family-clustered SEs.
SSC has no age at diagnosis and is not included. Output (aggregate): ~/rvh2/agedx/ptdt_full.tsv"""
import os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
A = pd.read_parquet(f"{OUT}/full_probands.parquet").drop_duplicates("sample_id")
allA = A.set_index("sample_id")
c = A[A.cohort.eq("SPARK") & A.case & A.predicted_ancestry.eq("EUR") & A.n_parents_genotyped.eq(2)].copy()
PG = ["adhd_2022", "mdd_2025", "edu_2022", "asd"]
for p in PG:
    mid = (c.mother.map(allA[p]) + c.father.map(allA[p])) / 2
    sd = mid.groupby(c.source).transform("std")
    c[p + "_ptdt"] = (c[p] - mid) / sd
rd = lambda f: pd.read_csv(f"{P}/{f}", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
reg, ho, env, dx = rd("spark_registration.tsv"), rd("spark_heldout.tsv"), rd("spark_env.tsv"), rd("spark_dx_derived.tsv")
s = c.sample_id; num = lambda ser: pd.to_numeric(s.map(ser), errors="coerce")
dxm = num(reg.diagnosis_age); c["late"] = np.where(dxm > 120, 1.0, np.where(dxm < 72, 0.0, np.nan))
c["adhd"] = num(ho.behav_adhd); c["regage"] = num(reg.age_at_registration_years); c["adult_at_reg"] = c.regage.ge(18).astype(float)
c["female"] = num(dx.sex.map({"Female": 1, "Male": 0})); c["dx_recent"] = (num(env.birth_year) + dxm / 12).ge(2019).astype(float)
c["fam"] = pd.factorize(c.family.fillna(c.sample_id))[0]
rows = []
for samp, d in [("registered before 22", c[c.regage.lt(22)]), ("all ages", c)]:
    for p in PG:
        y = p + "_ptdt"
        for lab, f in [("age only", f"{y} ~ late"), ("+ ADHD", f"{y} ~ late + adhd"), ("+ ADHD, adult registration, sex, era", f"{y} ~ late + adhd + adult_at_reg + female + dx_recent")]:
            if samp == "registered before 22": f = f.replace(" + adult_at_reg", "")
            x = d.dropna(subset=[y, "late", "adhd", "female", "dx_recent"])
            m = smf.ols(f, data=x).fit(cov_type="cluster", cov_kwds={"groups": x.fam.to_numpy()})
            for term in ["late", "adhd"]:
                if term in m.params: rows.append(dict(sample=samp, pgs=p, model=lab, term=term, beta=m.params[term], se=m.bse[term], p=m.pvalues[term], n=int(m.nobs), n_late=int(x.late.sum())))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/ptdt_full.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); print(T[T.term.eq("late")].round(4).to_string(index=False))
print(c.groupby("source").size().to_dict())
