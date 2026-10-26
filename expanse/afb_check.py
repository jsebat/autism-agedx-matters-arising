#!/usr/bin/env python3
"""Does the maternal-age split capture the genetics of reproductive timing? Among SPARK EUR autistic cases (cases.parquet; PGS in SD of
non-autistic relatives, pTDT in complete trios): PGS difference between cases of older (top third maternal age at birth) vs younger
(bottom third) mothers, and low vs high household income, for AFB 2018, EA 2022, IQ 2018, ADHD 2022, ASD 2026. Family-bootstrap SE.
Output (aggregate): ~/rvh2/agedx/afb_check.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
c = pd.read_parquet(f"{OUT}/cases.parquet"); sv = pd.read_csv(f"{R}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
SPL = {"maternal age (older vs younger third)": "s11", "household income (low vs high third)": "s09", "age at diagnosis (>10 vs <6 y)": "s00"}
fi = pd.factorize(c.FID.astype(str))[0]; nu = fi.max() + 1; rng = np.random.default_rng(1); BW = [np.bincount(rng.integers(0, nu, nu), minlength=nu)[fi] for _ in range(200)]
rows = []
for nm, code in SPL.items():
    g = pd.to_numeric(c.sample_id.map(sv[code]), errors="coerce").to_numpy(float)
    for p in ["afb_2018", "edu_2022", "iq_2018", "adhd_2022", "asd_2026"]:
        for meas, col in [("population", p), ("within-family pTDT", p + "_ptdt")]:
            v = c[col].to_numpy(float); ok = np.isfinite(v) & np.isfinite(g); h, l = ok & (g == 1), ok & (g == 0)
            d0 = v[h].mean() - v[l].mean(); bs = [np.average(v[h], weights=w[h]) - np.average(v[l], weights=w[l]) for w in BW if w[h].sum() and w[l].sum()]
            rows.append(dict(split=nm, pgs=p, measure=meas, n_high=int(h.sum()), n_low=int(l.sum()), diff=d0, se=float(np.std(bs)), z=d0 / float(np.std(bs))))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/afb_check.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(T.round(3).to_string(index=False))
