#!/usr/bin/env python3
"""Oaxaca-style decomposition of the 'after 10' vs 'before 6' gap in case PGS (SD vs non-autistic relatives of the same
platform; the contrast made by the Zhang, Grove et al. 2025 SPARK stratified GWAS). Pooled OLS on cases in the two strata:
  PGS ~ after10 + blocks; blocks = family structure (0/1/2 parents genotyped), adult at registration, female, diagnosis era,
  platform. Explained share of block k = sum over its terms of (mean_after10 - mean_before6) x coefficient (pooled model);
  'residual' = coefficient on after10 in the full model. Family-clustered bootstrap SEs (B = 200). Also reports the gap
  sequentially adjusted block by block, and the within-family pTDT difference for comparison.
Output: ~/rvh2/agedx/decomposition.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx")
c = pd.read_parquet(f"{OUT}/cases.parquet"); c = c[c.stratum.isin(["before 6", "after 10"])].copy()
c["after10"] = (c.stratum == "after 10").astype(float)
c["female"] = c.female.astype(float); c["adult"] = c.adult_at_reg.fillna(0).astype(float)
c["era"] = c.era.replace("nan", "unknown")
B = {"family structure": pd.get_dummies(c.n_par.astype(int), prefix="par", drop_first=True),
     "adult at registration": c[["adult"]], "sex": c[["female"]],
     "diagnosis era": pd.get_dummies(c.era, prefix="era", drop_first=True),
     "platform": pd.get_dummies(c.source, prefix="plat", drop_first=True)}
B = {k: v.astype(float) for k, v in B.items()}
PGS = ["mdd_2025", "adhd_2022", "edu_2022", "iq_2018", "asd_2026", "bip_2024", "scz_2022"]
fam = c.FID.astype(str).to_numpy(); fi, uf = pd.factorize(fam); rng = np.random.default_rng(3)
def fit(y, X):
    A = np.column_stack([np.ones(len(y)), X]); return np.linalg.lstsq(A, y, rcond=None)[0]
rows = []
for p in PGS:
    ok = c[p].notna().to_numpy(); y = c[p].to_numpy(float)[ok]; a10 = c.after10.to_numpy()[ok]
    Xb = {k: v.to_numpy()[ok] for k, v in B.items()}
    def decomp(idx):
        yy, aa = y[idx], a10[idx]; XX = np.column_stack([aa] + [Xb[k][idx] for k in B])
        b = fit(yy, XX); gap = yy[aa == 1].mean() - yy[aa == 0].mean(); out = {"total gap": gap, "residual (stratum)": b[1]}; j = 2
        for k in B:
            m = Xb[k][idx]; d_ = m[aa == 1].mean(0) - m[aa == 0].mean(0); kk = m.shape[1]
            out[k] = float(d_ @ b[j:j + kk]); j += kk
        seq = {}; cols = [aa]
        for k in B:
            cols.append(Xb[k][idx]); seq["after10 | + " + k] = fit(yy, np.column_stack(cols))[1]
        return out, seq
    base, seq = decomp(np.arange(len(y)))
    fio = fi[ok]; nu = fio.max() + 1; bs = []
    for _ in range(200):
        w = np.bincount(rng.integers(0, nu, nu), minlength=nu)[fio]; idx = np.repeat(np.arange(len(y)), w); bs.append(decomp(idx)[0])
    bs = pd.DataFrame(bs)
    for k, v in base.items(): rows.append(dict(pgs=p, component=k, estimate=v, se=bs[k].std(), kind="decomposition"))
    for k, v in seq.items(): rows.append(dict(pgs=p, component=k, estimate=v, se=np.nan, kind="sequential"))
    print(p, {k: round(v, 3) for k, v in base.items()}, flush=True)
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/decomposition.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(D[D.kind == "decomposition"].assign(v=lambda x: x.estimate.round(3).astype(str) + " (" + x.se.round(3).astype(str) + ")").pivot(index="component", columns="pgs", values="v").to_string())
print(D[D.kind == "sequential"].pivot(index="component", columns="pgs", values="estimate").round(3).to_string())
print("\nexplanatory variables: mean in before 6 vs after 10")
for k, v in B.items():
    for col in v.columns: print(f"  {k:24s} {col:22s} {v[col][c.after10 == 0].mean():.3f} {v[col][c.after10 == 1].mean():.3f}")
