#!/usr/bin/env python3
"""ALL-DATA split GWAS: inverse-variance meta of the per-stratum fastGWA-mlm-binary results over SIX strata
(SPARK GSA v1, GSA v2, Twist, iWGS v1.1, SPARK-WGS-2026_08, SSC), one row per person across strata. Same method as rep_combine.py
(arrays-only results in warrier_rep/combined stay as the Zhang-era replication). A GWAS absent from a stratum (e.g. SSC has no age strata)
is combined over the strata that have it. Per-SNP effective N = sum over strata of 4 / (1/cases + 1/controls).
Output: warrier_rep/combined_all/rep_<gwas>.sumstats.txt, counts_all.tsv, lambda.tsv"""
import glob, os
import numpy as np, pandas as pd
from scipy.stats import norm, chi2
B = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep"; OUT = f"{B}/combined_all"; os.makedirs(OUT, exist_ok=True)
CF = {"spark_gsa_v1": ["counts.tsv", "counts_survey.tsv", "counts_null.tsv"], "spark_gsa_v2": ["counts.tsv", "counts_survey.tsv", "counts_null.tsv"],
      "spark_twist": ["counts.tsv", "counts_survey.tsv", "counts_null.tsv"], "spark_wgs": ["counts_all.tsv", "counts_allnull.tsv"],
      "spark_wgs_2026_08": ["counts_all.tsv", "counts_allnull.tsv"], "ssc": ["counts_all.tsv"]}
C = pd.concat([pd.read_csv(f"{B}/{s}/{f}", sep="\t") for s, fs in CF.items() for f in fs if os.path.exists(f"{B}/{s}/{f}")]).drop_duplicates(["source", "gwas"])
C.to_csv(f"{OUT}/counts_all.tsv", sep="\t", index=False)
cc = C.set_index(["source", "gwas"]); lam = []
for g in sorted(C.gwas.unique()):
    parts = []
    for s in CF:
        fs = glob.glob(f"{B}/{s}/fg_{g}.fastGWA")
        if not fs or (s, g) not in cc.index: continue
        x = pd.read_csv(fs[0], sep="\t", low_memory=False)
        b, se = ("BETA", "SE") if "BETA" in x.columns else ("T", "SE_T")
        x = x[np.isfinite(x[b]) & (x[se] > 0)][["SNP", "A1", "A2", b, se]].rename(columns={b: "b", se: "se"})
        r = cc.loc[(s, g)]; x["neff"] = 4 / (1 / max(r.cases, 1) + 1 / max(r.controls, 1))
        lam.append(dict(gwas=g, source=s, lambda_gc=float(np.median((x.b / x.se) ** 2) / chi2.ppf(0.5, 1)), n_snp=len(x)))
        parts.append(x)
    if not parts: continue
    a = pd.concat(parts); ref = a.drop_duplicates("SNP").set_index("SNP")[["A1", "A2"]]
    a = a.join(ref, on="SNP", rsuffix="_ref"); a.loc[a.A1 != a.A1_ref, "b"] *= -1; a["w"] = 1 / a.se ** 2; a["wb"] = a.w * a.b
    m = a.groupby("SNP").agg(sw=("w", "sum"), swb=("wb", "sum"), N=("neff", "sum"), k=("w", "size")).join(ref)
    m["Z"] = m.swb / np.sqrt(m.sw); m["P"] = 2 * norm.sf(m.Z.abs())
    m.reset_index()[["SNP", "A1", "A2", "N", "Z", "P"]].to_csv(f"{OUT}/rep_{g}.sumstats.txt", sep="\t", index=False)
    lam.append(dict(gwas=g, source="combined", lambda_gc=float(np.median(m.Z ** 2) / chi2.ppf(0.5, 1)), n_snp=len(m), n_strata=len(parts)))
    print(g, len(parts), "strata", flush=True)
L = pd.DataFrame(lam); L.to_csv(f"{OUT}/lambda.tsv", sep="\t", index=False)
print(C.groupby("gwas")[["cases", "controls"]].sum().loc[lambda d: d.index.str.startswith("W_")].to_string())
