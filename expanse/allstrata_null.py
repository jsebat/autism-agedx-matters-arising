#!/usr/bin/env python3
"""Random-split null for the new SPARK strata (spark_wgs, spark_wgs_2026_08), with the same names as null_pheno.py on the arrays so the null
can be meta-analysed over all strata: ageN<r> (sizes of over10 / under6 in this stratum; r = 1..8, seeds 2001-2008) and adiN<r> (sizes of the ADI
thirds s08; r = 1..4, seeds 3001-3004). Same case/control definitions, exclusions (one row per person) and C1/C2 as allstrata_pheno.py.
Output: warrier_rep/<src>/allnull_pheno.txt, allnull_names.txt, counts_allnull.tsv"""
import os
import numpy as np, pandas as pd
SRC = os.environ["SRC"]; W = f"/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/{SRC}"
src = open(os.path.expanduser("~/rvh2/allstrata_pheno.py")).read()
exec(src[src.index("import os, hashlib"):src.index("for c, v in codes.items():")])   # reuse f, case, ctrl, c1, c2, S, codes (SPARK strata only)
cnt = {k: int(m.fillna(False).astype(bool).sum()) for k, (m, cc) in S.items()}
s08 = codes["s08"]; n_adi_hi, n_adi_lo = int((case & s08.eq(1)).sum()), int((case & s08.eq(0)).sum())
pool = np.where(case)[0]; P = f[["FID", "IID"]].copy(); rows = []
def add(k, idx, cc):
    m = np.zeros(len(f), bool); m[idx] = True; cc = cc.fillna(False).astype(bool).to_numpy()
    P[k] = np.where(m, 1, np.where(cc, 0, np.nan)); rows.append(dict(source=SRC, gwas=k, cases=int(m.sum()), controls=int(cc.sum())))
for pre, (na, nb), reps, seed0 in [("ageN", (cnt["W_over10"], cnt["W_under6"]), 8, 2000), ("adiN", (n_adi_hi, n_adi_lo), 4, 3000)]:
    for r in range(1, reps + 1):
        perm = np.random.default_rng(seed0 + r).permutation(pool)
        add(f"{pre}{r}_A_c2", perm[:na], c2); add(f"{pre}{r}_B_c1", perm[na:na + nb], c1)
P.fillna("NA").to_csv(f"{W}/allnull_pheno.txt", sep=" ", index=False, header=False)
open(f"{W}/allnull_names.txt", "w").write("\n".join(P.columns[2:]) + "\n")
pd.DataFrame(rows).to_csv(f"{W}/counts_allnull.tsv", sep="\t", index=False); print(SRC, len(rows), "null GWAS", rows[:2])
