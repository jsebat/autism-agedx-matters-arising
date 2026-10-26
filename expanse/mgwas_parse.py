#!/usr/bin/env python3
"""Parse LDSC --rg logs of the matched-design GWAS into tables: all pairwise rg, within-split pairs by design, rg with external
traits, and h2 (observed scale) of each GWAS. Output: cnv_drift/mgwas/combined/rg_all.tsv, rg_pairs.tsv, rg_external.tsv"""
import glob, os, re, io
import pandas as pd
C = "/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/combined"
rows = []
for f in sorted(glob.glob(f"{C}/rg/*.log")):
    t = open(f).read()
    if "Summary of Genetic Correlation Results" not in t: print("no summary:", os.path.basename(f)); continue
    blk = t.split("Summary of Genetic Correlation Results")[1].strip().split("\n\n")[0]
    rows.append(pd.read_csv(io.StringIO(blk), sep=r"\s+"))
A = pd.concat(rows, ignore_index=True)
nm = lambda s: os.path.basename(str(s)).replace(".sumstats.gz", "")
A["p1"] = A.p1.map(nm); A["p2"] = A.p2.map(nm); A = A[A.p1 != A.p2]
A.to_csv(f"{C}/rg_all.tsv", sep="\t", index=False)
SPL = {"agedx": ("early", "late"), "adi": ("low", "high"), "adult": ("child", "adult"), "adhd": ("no", "yes")}
pr = []
for sp, (a, b) in SPL.items():
    for d1 in ("between", "sibtrio"):
        for d2 in ("between", "sibtrio"):
            for h1, h2 in [(a, b), (a, a), (b, b)]:
                if (h1, d1) == (h2, d2): continue
                x = A[(A.p1 == f"{d1}_{sp}_{h1}") & (A.p2 == f"{d2}_{sp}_{h2}")]
                if len(x): r = x.iloc[0]; pr.append(dict(split=sp, a=f"{d1} {h1}", b=f"{d2} {h2}", rg=r.rg, se=r.se, p=r.p))
P = pd.DataFrame(pr); P.to_csv(f"{C}/rg_pairs.tsv", sep="\t", index=False)
E = A[A.p2.str.startswith("ext_") & ~A.p1.str.startswith("ext_")].copy(); E["trait"] = E.p2.str.replace("ext_", "", regex=False)
E.to_csv(f"{C}/rg_external.tsv", sep="\t", index=False)
H = A.drop_duplicates("p2")[["p2", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se"]].rename(columns={"p2": "gwas"})
H.to_csv(f"{C}/h2.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(H[~H.gwas.str.startswith("ext_")].round(4).to_string(index=False))
print(P.round(3).to_string(index=False))
print(E.pivot_table(index="p1", columns="trait", values="rg").round(2).to_string())
print(E.pivot_table(index="p1", columns="trait", values="se").round(2).to_string())
