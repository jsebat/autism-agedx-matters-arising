#!/usr/bin/env python3
"""Replication of Zhang, Grove et al. 2025 SPARK stratified GWAS, step 2 per platform (SRC=<source>): phenotype and covariate files
exactly as in SPARK_replicationGWAS.R, then the same design applied to alternative splits of the autistic participants.
  cases: asd == TRUE in SPARK individuals_registration (release 2026-09-28); controls: ALL genotyped non-autistic participants
  (parents and siblings), coded 0 in every analysis; cases outside the split half are NA.
  their strata: over10 = diagnosis_age > 120 months, under6 = diagnosis_age < 72 months (and all cases)
  alternative splits (same controls): area deprivation top/bottom third, adult/child at registration, ADHD yes/no; and splits
  orthogonal to age at diagnosis: ADI, ADHD and sex WITHIN cases diagnosed < 72 months; over10 vs under6 WITHIN child registrants.
  qcovar: 10 PC-AiR-style PCs + age_at_registration_years; covar: sex (their 'batch' = imputation batch; here one platform per run).
Output: warrier_rep/<src>/{pheno.txt, pheno_names.txt, qcovar.txt, covar.txt, counts.tsv}"""
import os
import numpy as np, pandas as pd
SRC = os.environ["SRC"]; W = f"/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/{SRC}"; R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
B = "/expanse/projects/sebat1/s3/data/sebat"
fam = pd.read_csv(f"{W}/hm3w.fam", sep=r"\s+", header=None, dtype=str).iloc[:, :2]; fam.columns = ["FID", "IID"]
d = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "sample_id"])
d = d[d.source == SRC].drop_duplicates("IID").set_index("IID")
fam["sp"] = fam.IID.map(d.sample_id)
reg = pd.read_csv(f"{R}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
ho = pd.read_csv(f"{R}/spark_heldout.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
env = pd.read_csv(f"{R}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
allc = pd.read_parquet(os.path.expanduser("~/rvh2/agedx/cases.parquet"))
adi_all = pd.to_numeric(allc.sample_id.map(env.adi_pct), errors="coerce"); lo, hi = adi_all.quantile(1 / 3), adi_all.quantile(2 / 3)
f = fam.copy()
f["asd"] = f.sp.map(reg.asd).astype(str).str.upper().eq("TRUE")
f["known"] = f.sp.isin(reg.index)
f["dxm"] = pd.to_numeric(f.sp.map(reg.diagnosis_age), errors="coerce")
f["areg"] = pd.to_numeric(f.sp.map(reg.age_at_registration_years), errors="coerce")
f["sex"] = f.sp.map(reg.sex).map({"Male": 1, "Female": 2})
f["adi"] = pd.to_numeric(f.sp.map(env.adi_pct), errors="coerce"); f["adhd"] = pd.to_numeric(f.sp.map(ho.behav_adhd), errors="coerce")
case, ctrl = f.asd & f.known, (~f.asd) & f.known
u6, o10, child = f.dxm < 72, f.dxm > 120, f.areg < 18
S = {"W_all": case, "W_under6": case & u6, "W_over10": case & o10,
     "adi_high": case & (f.adi >= hi), "adi_low": case & (f.adi <= lo), "adult_adult": case & (f.areg >= 18), "adult_child": case & child,
     "adhd_yes": case & f.adhd.eq(1), "adhd_no": case & f.adhd.eq(0),
     "adiE_high": case & u6 & (f.adi >= hi), "adiE_low": case & u6 & (f.adi <= lo), "adhdE_yes": case & u6 & f.adhd.eq(1), "adhdE_no": case & u6 & f.adhd.eq(0),
     "sexE_female": case & u6 & f.sex.eq(2), "sexE_male": case & u6 & f.sex.eq(1), "agedxC_late": case & child & o10, "agedxC_early": case & child & u6}
P = f[["FID", "IID"]].copy(); rows = []
for k, m in S.items():
    P[k] = np.where(m, 1, np.where(ctrl, 0, np.nan))
    rows.append(dict(source=SRC, gwas=k, cases=int(m.sum()), controls=int(ctrl.sum())))
# Robustness 1 - disjoint controls: non-autistic participants split at random BY FAMILY into C1 / C2; the first half of each pair is
# tested against C1 and the second against C2, so no pair of GWAS shares a control (over10 also against C1 for cross-split rg).
fam_id = f.sp.map(reg.family_sf_id).fillna(f.IID).astype(str)
uf = np.sort(fam_id.unique()); rng = np.random.default_rng(20261023); c1fam = set(uf[rng.random(len(uf)) < 0.5])
c1, c2 = ctrl & fam_id.isin(c1fam), ctrl & ~fam_id.isin(c1fam)
PAIRS = [("W_under6", "W_over10"), ("adi_low", "adi_high"), ("adult_child", "adult_adult"), ("adhd_no", "adhd_yes"), ("adiE_low", "adiE_high"),
         ("adhdE_no", "adhdE_yes"), ("sexE_male", "sexE_female"), ("agedxC_early", "agedxC_late")]
def add(k, m, cc):
    P[k] = np.where(m, 1, np.where(cc, 0, np.nan)); rows.append(dict(source=SRC, gwas=k, cases=int(m.sum()), controls=int(cc.sum())))
for a, b in PAIRS:
    add(f"{a}_c1", S[a], c1); add(f"{b}_c2", S[b], c2)
add("W_over10_c1", S["W_over10"], c1)
# Robustness 2 - random-split null: cases with a known diagnosis age split at random into two disjoint sets of the SAME sizes as
# over10 / under6 on this platform (true rg = 1), disjoint controls as above; 4 replicates.
pool = np.where(case & f.dxm.notna())[0]; n_o, n_u = int(S["W_over10"].sum()), int(S["W_under6"].sum())
for r in range(1, 5):
    perm = np.random.default_rng(1000 + r).permutation(pool); ia, ib = perm[:n_o], perm[n_o:n_o + n_u]
    ma = np.zeros(len(f), bool); ma[ia] = True; mb = np.zeros(len(f), bool); mb[ib] = True
    add(f"rand{r}_A_c2", pd.Series(ma, index=f.index), c2); add(f"rand{r}_B_c1", pd.Series(mb, index=f.index), c1)
S = {k: None for k in P.columns[2:]}
P.fillna("NA").to_csv(f"{W}/pheno.txt", sep=" ", index=False, header=False)
open(f"{W}/pheno_names.txt", "w").write("\n".join(S) + "\n")
pc = pd.read_csv(f"{W}/pcair.sscore", sep="\t"); ic = [c for c in pc.columns if "IID" in c][0]; pcs = [c for c in pc.columns if c.endswith("_AVG")][:10]
q = f[["FID", "IID"]].merge(pc[[ic] + pcs].rename(columns={ic: "IID"}), on="IID", how="left"); q["areg"] = f.areg.values
q.fillna("NA").to_csv(f"{W}/qcovar.txt", sep=" ", index=False, header=False)
c = f[["FID", "IID"]].copy(); c["sex"] = f.sex.map({1: "M", 2: "F"}); c.fillna("NA").to_csv(f"{W}/covar.txt", sep=" ", index=False, header=False)
C = pd.DataFrame(rows); C.to_csv(f"{W}/counts.tsv", sep="\t", index=False); print(C.to_string(index=False))
print("samples", len(f), "with registration", int(f.known.sum()), "PCs missing", int(q[pcs[0]].isna().sum()), "age_reg missing", int(f.areg.isna().sum()))
