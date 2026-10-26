#!/usr/bin/env python3
"""Matched-design GWAS, step 4: combine platforms and write LDSC-ready sumstats (SNP A1 A2 N Z) per split x half x design.
  between-family: inverse-variance meta of log(OR) across GSA v1 / GSA v2 / Twist (alleles aligned on A1); per-SNP effective
                  N = sum over contributing platforms of 4 / (1/cases + 1/controls).
  trio (TDT):     transmitted (T) / untransmitted (U) counts of A1 summed across platforms; Z = (T - U) / sqrt(T + U);
                  N = 2 x affected offspring (1:1 case / pseudo-control) on contributing platforms.
Output: cnv_drift/mgwas/combined/<design>_<split>_<half>.sumstats.txt and design_counts_all.tsv"""
import glob, os
import numpy as np, pandas as pd
from scipy.stats import norm
M = "/expanse/projects/sebat1/jsebat/cnv_drift/mgwas"; OUT = f"{M}/combined"; os.makedirs(OUT, exist_ok=True)
SRCS = ["spark_gsa_v1", "spark_gsa_v2", "spark_twist"]
DC = pd.concat([pd.read_csv(f"{M}/{s}/design_counts.tsv", sep="\t") for s in SRCS]); DC.to_csv(f"{OUT}/design_counts_all.tsv", sep="\t", index=False)
dc = DC.set_index(["source", "split", "half"])
for split, half in DC[["split", "half"]].drop_duplicates().itertuples(index=False):
    # between-family
    parts = []
    for s in SRCS:
        fs = glob.glob(f"{M}/{s}/bf_{split}_{half}.Y.glm.logistic*")
        if not fs: continue
        x = pd.read_csv(fs[0], sep="\t", low_memory=False)
        x = x[x.P.notna() & x["LOG(OR)_SE"].notna()].copy()
        x["A2"] = np.where(x.A1 == x.ALT, x.REF, x.ALT); x["b"] = np.log(x.OR); x["se"] = x["LOG(OR)_SE"]
        r = dc.loc[(s, split, half)]; x["neff"] = 4 / (1 / max(r.between_cases, 1) + 1 / max(r.between_controls, 1))
        parts.append(x[["ID", "A1", "A2", "b", "se", "neff"]].assign(src=s))
    if parts:
        a = pd.concat(parts); ref = a.drop_duplicates("ID").set_index("ID")[["A1", "A2"]]
        a = a.join(ref, on="ID", rsuffix="_ref"); flip = a.A1 != a.A1_ref; a.loc[flip, "b"] *= -1
        a["w"] = 1 / a.se ** 2
        g = a.groupby("ID").agg(sw=("w", "sum"), swb=("b", lambda v: 0), N=("neff", "sum"))
        g["swb"] = (a.w * a.b).groupby(a.ID).sum(); g["Z"] = (g.swb / g.sw) / np.sqrt(1 / g.sw)
        g = g.join(ref); g.reset_index().rename(columns={"ID": "SNP"}).assign(P=lambda q: 2 * norm.sf(q.Z.abs()))[["SNP", "A1", "A2", "N", "Z", "P"]].to_csv(f"{OUT}/between_{split}_{half}.sumstats.txt", sep="\t", index=False)
    # trio TDT
    parts = []
    for s in SRCS:
        f = f"{M}/{s}/tdt_{split}_{half}.tdt"
        if not os.path.exists(f): continue
        x = pd.read_csv(f, sep=r"\s+", low_memory=False); x = x[(x["T"] + x["U"]) > 0]
        r = dc.loc[(s, split, half)]
        parts.append(x[["SNP", "A1", "A2", "T", "U"]].assign(n=2 * r.trio_offspring))
    if parts:
        a = pd.concat(parts); ref = a.drop_duplicates("SNP").set_index("SNP")[["A1", "A2"]]
        a = a.join(ref, on="SNP", rsuffix="_ref"); flip = a.A1 != a.A1_ref
        a.loc[flip, ["T", "U"]] = a.loc[flip, ["U", "T"]].values
        g = a.groupby("SNP").agg(T=("T", "sum"), U=("U", "sum"), N=("n", "sum")).join(ref)
        g["Z"] = (g["T"] - g["U"]) / np.sqrt(g["T"] + g["U"])
        # sibling-null version: per-SNP null transmission rate p0 from the same platforms' non-autistic trio offspring
        sn = []
        for s in SRCS:
            f = f"{M}/{s}/tdt_sibnull.tdt"
            if os.path.exists(f): sn.append(pd.read_csv(f, sep=r"\s+", low_memory=False)[["SNP", "A1", "T", "U"]])
        if sn:
            q = pd.concat(sn).join(ref, on="SNP", rsuffix="_ref"); fl = q.A1 != q.A1_ref
            q.loc[fl, ["T", "U"]] = q.loc[fl, ["U", "T"]].values; q = q.groupby("SNP")[["T", "U"]].sum()
            gg = g.join(q, rsuffix="_sib", how="inner"); n = gg["T"] + gg["U"]; ns = gg["T_sib"] + gg["U_sib"]
            gg = gg[(n > 0) & (ns > 0)]; n = gg["T"] + gg["U"]; ns = gg["T_sib"] + gg["U_sib"]; p0 = (gg["T_sib"] + 0.5) / (ns + 1)
            gg["Z"] = (gg["T"] - n * p0) / np.sqrt(n * p0 * (1 - p0) + n ** 2 * p0 * (1 - p0) / ns)
            gg.reset_index().assign(P=lambda q_: 2 * norm.sf(q_.Z.abs()))[["SNP", "A1", "A2", "N", "Z", "P"]].to_csv(f"{OUT}/sibtrio_{split}_{half}.sumstats.txt", sep="	", index=False)
        g.reset_index().assign(P=lambda q: 2 * norm.sf(q.Z.abs()))[["SNP", "A1", "A2", "N", "Z", "P"]].to_csv(f"{OUT}/trio_{split}_{half}.sumstats.txt", sep="\t", index=False)
    print(split, half, "combined", flush=True)
print(DC.groupby(["split", "half"])[["between_cases", "between_controls", "trio_offspring"]].sum().to_string())
