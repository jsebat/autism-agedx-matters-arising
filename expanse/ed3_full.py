#!/usr/bin/env python3
"""Extended Data Fig. 3a on the FULL SPARK data (registered before 22): polygenic-score difference between the halves of three
splits (maternal age, household income, age at diagnosis >10 vs <6 y), between families (autistic participants, PGS standardized within platform;
platform fixed effects) and within families (pTDT, both parents genotyped). OLS with family-clustered SEs. Output: ~/rvh2/agedx/ed3a_full.tsv"""
import os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
A = pd.read_parquet(f"{OUT}/full_probands.parquet").drop_duplicates("sample_id"); allA = A.set_index("sample_id")
reg = pd.read_csv(f"{P}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
sv = pd.read_csv(f"{P}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
c = A[A.cohort.eq("SPARK") & A.case & A.predicted_ancestry.eq("EUR")].copy()
c = c[pd.to_numeric(c.sample_id.map(reg.age_at_registration_years), errors="coerce").lt(22)].copy()
PG = ["afb_2018", "edu_2022", "adhd_2022", "asd"]
for p in PG:
    mid = (c.mother.map(allA[p]) + c.father.map(allA[p])) / 2; c[p + "_ptdt"] = (c[p] - mid) / mid.groupby(c.source).transform("std")
c["fam"] = pd.factorize(c.family.fillna(c.sample_id))[0]
SPL = {"maternal age (older vs younger third)": "s11", "household income (low vs high third)": "s09", "age at diagnosis (>10 vs <6 y)": "s00"}
rows = []
for nm, code in SPL.items():
    c["g"] = pd.to_numeric(c.sample_id.map(sv[code]), errors="coerce")
    for p in PG:
        for meas, y, f in [("population", p, f"{p} ~ g + C(source)"), ("within-family pTDT", p + "_ptdt", f"{p}_ptdt ~ g")]:
            x = c.dropna(subset=[y, "g"])
            m = smf.ols(f, data=x).fit(cov_type="cluster", cov_kwds={"groups": x.fam.to_numpy()})
            rows.append(dict(split=nm, pgs="asd_2026" if p == "asd" else p, measure=meas, diff=m.params["g"], se=m.bse["g"], p=m.pvalues["g"], n=int(m.nobs)))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/ed3a_full.tsv", sep="\t", index=False); print(T.round(3).to_string(index=False))
