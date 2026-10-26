#!/usr/bin/env python3
"""Random-split null for the replication (their pipeline), per platform (SRC): cases split at random into two disjoint sets whose
sizes match a real split on this platform; each set tested against its own control half (C1/C2, the family-level split used in
rep_pheno.py), so the true rg between the two GWAS is 1 and any shortfall is estimation noise.
  ageN<r>  sizes of over10 / under6 (r = 1..16; replicates 1-4 of rep_pheno.py are seeds 1001-1004, these are 2001-2016)
  adiN<r>  sizes of ADI high / low thirds (r = 1..8)
Output: warrier_rep/<src>/null_pheno.txt, null_names.txt, counts_null.tsv"""
import os
import numpy as np, pandas as pd
SRC = os.environ["SRC"]; W = f"/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/{SRC}"; R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
B = "/expanse/projects/sebat1/s3/data/sebat"
fam = pd.read_csv(f"{W}/hm3w.fam", sep=r"\s+", header=None, dtype=str).iloc[:, :2]; fam.columns = ["FID", "IID"]
d = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "sample_id"])
d = d[d.source == SRC].drop_duplicates("IID").set_index("IID"); fam["sp"] = fam.IID.map(d.sample_id)
reg = pd.read_csv(f"{R}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
f = fam.copy(); f["asd"] = f.sp.map(reg.asd).astype(str).str.upper().eq("TRUE"); f["known"] = f.sp.isin(reg.index)
case, ctrl = f.asd & f.known, (~f.asd) & f.known
fam_id = f.sp.map(reg.family_sf_id).fillna(f.IID).astype(str)
uf = np.sort(fam_id.unique()); rng = np.random.default_rng(20261023); c1fam = set(uf[rng.random(len(uf)) < 0.5])   # same C1/C2 as rep_pheno.py
c1, c2 = ctrl & fam_id.isin(c1fam), ctrl & ~fam_id.isin(c1fam)
cnt = pd.read_csv(f"{W}/counts.tsv", sep="\t").set_index("gwas").cases
pool = np.where(case)[0]; P = f[["FID", "IID"]].copy(); rows = []
def add(k, idx, cc):
    m = np.zeros(len(f), bool); m[idx] = True; P[k] = np.where(m, 1, np.where(cc, 0, np.nan)); rows.append(dict(source=SRC, gwas=k, cases=int(m.sum()), controls=int(cc.sum())))
for pre, (a, b), reps, seed0 in [("ageN", ("W_over10", "W_under6"), 16, 2000), ("adiN", ("adi_high", "adi_low"), 8, 3000)]:
    na, nb = int(cnt[a]), int(cnt[b])
    for r in range(1, reps + 1):
        perm = np.random.default_rng(seed0 + r).permutation(pool)
        add(f"{pre}{r}_A_c2", perm[:na], c2); add(f"{pre}{r}_B_c1", perm[na:na + nb], c1)
P.fillna("NA").to_csv(f"{W}/null_pheno.txt", sep=" ", index=False, header=False)
open(f"{W}/null_names.txt", "w").write("\n".join(P.columns[2:]) + "\n")
pd.DataFrame(rows).to_csv(f"{W}/counts_null.tsv", sep="\t", index=False); print(SRC, len(rows), "null GWAS")
