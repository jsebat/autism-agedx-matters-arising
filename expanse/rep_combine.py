#!/usr/bin/env python3
"""Replication step 4: inverse-variance meta of the per-platform fastGWA-mlm-binary results (their meta: SPARK v3 + v4 releases;
here GSA v1 + GSA v2 + Twist), per-SNP effective N = sum over platforms of 4 / (1/cases + 1/controls).
Output: warrier_rep/combined/rep_<gwas>.sumstats.txt (SNP A1 A2 N Z P), counts_all.tsv, lambda per platform and combined."""
import glob, os
import numpy as np, pandas as pd
from scipy.stats import norm, chi2
B = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep"; OUT = f"{B}/combined"; os.makedirs(OUT, exist_ok=True)
SRCS = ["spark_gsa_v1", "spark_gsa_v2", "spark_twist"]
CF = os.environ.get("COUNTS", "counts.tsv"); TAG = "" if CF == "counts.tsv" else "_" + CF.replace("counts_", "").replace(".tsv", "")   # COUNTS=counts_null.tsv for the random-split null
C = pd.concat([pd.read_csv(f"{B}/{s}/{CF}", sep="\t") for s in SRCS]); C.to_csv(f"{OUT}/counts{TAG}_all.tsv", sep="\t", index=False)
cc = C.set_index(["source", "gwas"]); lam = []
for g in C.gwas.unique():
    parts = []
    for s in SRCS:
        fs = glob.glob(f"{B}/{s}/fg_{g}.fastGWA")
        if not fs: continue
        x = pd.read_csv(fs[0], sep="\t", low_memory=False)
        b, se = ("BETA", "SE") if "BETA" in x.columns else ("T", "SE_T")
        x = x[np.isfinite(x[b]) & (x[se] > 0)][["SNP", "A1", "A2", b, se]].rename(columns={b: "b", se: "se"})
        r = cc.loc[(s, g)]; x["neff"] = 4 / (1 / max(r.cases, 1) + 1 / max(r.controls, 1))
        lam.append(dict(gwas=g, source=s, lambda_gc=float(np.median((x.b / x.se) ** 2) / chi2.ppf(0.5, 1)), n_snp=len(x)))
        parts.append(x)
    if not parts: continue
    a = pd.concat(parts); ref = a.drop_duplicates("SNP").set_index("SNP")[["A1", "A2"]]
    a = a.join(ref, on="SNP", rsuffix="_ref"); a.loc[a.A1 != a.A1_ref, "b"] *= -1; a["w"] = 1 / a.se ** 2; a["wb"] = a.w * a.b
    m = a.groupby("SNP").agg(sw=("w", "sum"), swb=("wb", "sum"), N=("neff", "sum")).join(ref)
    m["Z"] = m.swb / np.sqrt(m.sw); m["P"] = 2 * norm.sf(m.Z.abs())
    m.reset_index()[["SNP", "A1", "A2", "N", "Z", "P"]].to_csv(f"{OUT}/rep_{g}.sumstats.txt", sep="\t", index=False)
    lam.append(dict(gwas=g, source="combined", lambda_gc=float(np.median(m.Z ** 2) / chi2.ppf(0.5, 1)), n_snp=len(m)))
    print(g, "done", flush=True)
L = pd.DataFrame(lam); L.to_csv(f"{OUT}/lambda{TAG}.tsv", sep="\t", index=False)
print(L.pivot(index="gwas", columns="source", values="lambda_gc").round(3).to_string())
print(C.groupby("gwas")[["cases", "controls"]].sum().to_string())
