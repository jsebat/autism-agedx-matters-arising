#!/usr/bin/env python3
"""Inverse-variance combination of the per-stratum GCTA-GREML SNP heritabilities over ALL SPARK strata (GSA v1, GSA v2, Twist, iWGS v1.1,
SPARK-WGS-2026_08; one row per person across strata). Binary traits use the liability-scale estimate (V(G)/Vp_L), except late_dx (observed scale, matching greml_ivw). Arrays-only rows are also
reported (Zhang-era). Output (aggregate): ~/rvh2/agedx/greml_all.tsv"""
import os, re
import numpy as np, pandas as pd
B = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep"; OUT = os.path.expanduser("~/rvh2/agedx")
SRCS = ["spark_gsa_v1", "spark_gsa_v2", "spark_twist", "spark_wgs", "spark_wgs_2026_08"]
TR = ["dx_age", "dx_age_u22", "adi_pct", "adult_reg", "adhd", "late_dx", "mat_age", "pat_age", "low_income", "mother_nodeg"]
def read(f):
    t = open(f).read(); est = {}
    for line in t.splitlines():
        p = line.split()
        if p and p[0] in ("V(G)/Vp", "V(G)/Vp_L") and len(p) >= 3: est[p[0]] = (float(p[1]), float(p[2]))
        if p and p[0] == "n": est["n"] = int(p[1])
    return est
rows = []
for s in SRCS:
    for t in TR:
        f = f"{B}/{s}/greml_{t}.hsq"
        if not os.path.exists(f): continue
        e = read(f); k = "V(G)/Vp_L" if ("V(G)/Vp_L" in e and t != "late_dx") else "V(G)/Vp"   # late_dx on the observed scale, as in the arrays-only analysis (greml_ivw)
        if k in e: rows.append(dict(source=s, trait=t, h2=e[k][0], se=e[k][1], n=e.get("n", np.nan), scale="liability" if k.endswith("_L") else "observed"))
D = pd.DataFrame(rows); out = []
for t, g in D.groupby("trait"):
    for lab, gg in [("all SPARK strata", g), ("arrays only (Zhang-era)", g[g.source.isin(SRCS[:3])])]:
        gg = gg[(gg.se > 0) & (gg.n >= 500)]   # strata with < 500 unrelated cases (iWGS v1.1 after removing array-genotyped people) are uninformative
        if not len(gg): continue
        w = 1 / gg.se ** 2; h = float((w * gg.h2).sum() / w.sum()); se = float(np.sqrt(1 / w.sum()))
        out.append(dict(trait=t, sample=lab, h2=h, se=se, z=h / se, n=int(gg.n.sum()), strata=len(gg)))
T = pd.DataFrame(out); T.to_csv(f"{OUT}/greml_all.tsv", sep="\t", index=False); D.to_csv(f"{OUT}/greml_by_stratum.tsv", sep="\t", index=False)
print(T.round(3).to_string(index=False)); print(D.round(3).to_string(index=False))
