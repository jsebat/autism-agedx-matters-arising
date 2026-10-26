#!/usr/bin/env python3
"""Overlap-free factor-PGS test: SPARK cases diagnosed at 6-10 y were in NEITHER factor's discovery GWAS (F1 used SPARK < 6 y, F2
SPARK > 10 y; Zhang, Grove et al. 2025, Box 1), so among them case-vs-case contrasts of the F2 - F1 score are not inflated by sample
overlap. Within '< 6 y' (all in F1) and '> 10 y' (all in F2) the inflation is shared by every case of the stratum, so
case-vs-case contrasts there are overlap-balanced too (ST env var). Same splits as driver_pgs.py (age at diagnosis excluded), population comparison among 6-10 y cases only.
Output (aggregate): ~/rvh2/agedx/f12_mid_driver.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx")
ST = os.environ.get("ST", "6-10"); c = pd.read_parquet(f"{OUT}/cases_f12.parquet"); c = c[c.stratum.eq(ST)].copy(); print(ST, "cases:", len(c))
src = open(os.path.expanduser("~/rvh2/driver_pgs.py")).read()
i = src.index("fam = c.FID"); pre = src[src.index("sev = pd.read_csv"):i]
R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
exec(pre)
SPL.pop("age at diagnosis (after 10 vs before 6)", None)
fi = pd.factorize(c.FID.astype(str))[0]; nu = fi.max() + 1; rng = np.random.default_rng(1); BW = [np.bincount(rng.integers(0, nu, nu), minlength=nu)[fi] for _ in range(200)]
rows = []
for nm, g in SPL.items():
    g = np.asarray(g, float)
    for p in ["F1", "F2", "F2mF1"]:
        v = c[p].to_numpy(float); ok = np.isfinite(v) & np.isfinite(g); h, l = ok & (g == 1), ok & (g == 0)
        if h.sum() < 30 or l.sum() < 30: continue
        d0 = v[h].mean() - v[l].mean(); bs = [np.average(v[h], weights=w[h]) - np.average(v[l], weights=w[l]) for w in BW if w[h].sum() and w[l].sum()]
        rows.append(dict(split=nm, pgs=p, n_high=int(h.sum()), n_low=int(l.sum()), diff=d0, se=float(np.std(bs)), z=d0 / float(np.std(bs))))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/f12_mid_driver_{ST.replace(' ', '')}.tsv", sep="\t", index=False)
x = T[T.pgs == "F2mF1"].assign(a=lambda d: d.z.abs()).sort_values("a", ascending=False)
pd.set_option("display.width", 200); print(x.drop(columns="a").round(3).to_string(index=False))
