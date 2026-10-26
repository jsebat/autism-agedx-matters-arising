#!/usr/bin/env python3
"""Their within-family evidence (Zhang, Grove et al. 2025, Discussion / Supplementary Table 24): ADHD PGS is over-transmitted mainly
to later-diagnosed autistic individuals. Test whether age at diagnosis predicts the ADHD-PGS pTDT once ADHD comorbidity (and adult
registration, sex, diagnosis era) are in the model. SPARK EUR complete trios (cases.parquet; pTDT = child minus mid-parent, SD of
mid-parent). OLS with family-clustered SEs. Also MDD and EA PGS for comparison.
Output (aggregate): ~/rvh2/agedx/adhd_ptdt_conditional.tsv"""
import os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
OUT = os.path.expanduser("~/rvh2/agedx"); R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
c = pd.read_parquet(f"{OUT}/cases.parquet")
ho = pd.read_csv(f"{R}/spark_heldout.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
c["adhd"] = pd.to_numeric(c.sample_id.map(ho.behav_adhd), errors="coerce")
c["late"] = np.where(c.stratum.eq("after 10"), 1.0, np.where(c.stratum.eq("before 6"), 0.0, np.nan))
c["dx_recent"] = pd.to_numeric(c.dx_year, errors="coerce").ge(2019).astype(float)
rows = []
for p in ["adhd_2022", "mdd_2025", "edu_2022", "asd_2026"]:
    y = p + "_ptdt"
    for lab, f in [("age only", f"{y} ~ late"), ("+ ADHD comorbidity", f"{y} ~ late + adhd"),
                   ("+ ADHD, adult registration, sex, era", f"{y} ~ late + adhd + adult_at_reg + female + dx_recent")]:
        x = c.dropna(subset=[y, "late", "adhd", "adult_at_reg", "female"])
        m = smf.ols(f, data=x).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(x.FID.astype(str))[0]})
        for term in ["late", "adhd"]:
            if term in m.params: rows.append(dict(pgs=p, model=lab, term=term, beta=m.params[term], se=m.bse[term], z=m.tvalues[term], p=m.pvalues[term], n=int(m.nobs)))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/adhd_ptdt_conditional.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(T.round(4).to_string(index=False))
