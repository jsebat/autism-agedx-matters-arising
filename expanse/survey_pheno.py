#!/usr/bin/env python3
"""Broad survey, per platform (SRC): for every split in pheno/spark_survey.tsv (built by survey_build.py from the SPARK release) two
GWAS with their pipeline - high/yes half vs control half C2, low/no half vs control half C1 (the family-level control split of
rep_pheno.py, so the two GWAS of a pair share no control). Plus one platform-specific design split: no parent vs both parents
genotyped on this platform (s99).
Output: warrier_rep/<src>/survey_pheno.txt, survey_names.txt, counts_survey.tsv"""
import os
import numpy as np, pandas as pd
SRC = os.environ["SRC"]; W = f"/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/{SRC}"; R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
B = "/expanse/projects/sebat1/s3/data/sebat"
fam = pd.read_csv(f"{W}/hm3w.fam", sep=r"\s+", header=None, dtype=str).iloc[:, :2]; fam.columns = ["FID", "IID"]
d = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "sample_id"])
d = d[d.source == SRC].drop_duplicates("IID").set_index("IID"); fam["sp"] = fam.IID.map(d.sample_id)
reg = pd.read_csv(f"{R}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
ped = pd.read_csv(f"{R}/spark_pedigree_ids.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
sv = pd.read_csv(f"{R}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
f = fam.copy(); f["asd"] = f.sp.map(reg.asd).astype(str).str.upper().eq("TRUE"); f["known"] = f.sp.isin(reg.index)
case, ctrl = f.asd & f.known, (~f.asd) & f.known
fam_id = f.sp.map(reg.family_sf_id).fillna(f.IID).astype(str)
uf = np.sort(fam_id.unique()); rng = np.random.default_rng(20261023); c1fam = set(uf[rng.random(len(uf)) < 0.5])   # same C1/C2 as rep_pheno.py
c1, c2 = ctrl & fam_id.isin(c1fam), ctrl & ~fam_id.isin(c1fam)
geno = set(f.sp.dropna()); npar = f.sp.map(ped.biomother_sp_id).isin(geno).astype(int) + f.sp.map(ped.biofather_sp_id).isin(geno).astype(int)
codes = {c: pd.to_numeric(f.sp.map(sv[c]), errors="coerce") for c in sv.columns}
codes["s99"] = pd.Series(np.where(npar == 0, 1.0, np.where(npar == 2, 0.0, np.nan)), index=f.index)
P = f[["FID", "IID"]].copy(); rows = []
for c, v in codes.items():
    hi, lo = case & v.eq(1), case & v.eq(0)
    for k, m, cc in [(f"{c}_hi_c2", hi, c2), (f"{c}_lo_c1", lo, c1)]:
        P[k] = np.where(m, 1, np.where(cc, 0, np.nan)); rows.append(dict(source=SRC, gwas=k, cases=int(m.sum()), controls=int(cc.sum())))
keep = [r["gwas"] for r in rows if r["cases"] >= 150]   # too few cases on a platform -> not run there
P[["FID", "IID"] + keep].fillna("NA").to_csv(f"{W}/survey_pheno.txt", sep=" ", index=False, header=False)
open(f"{W}/survey_names.txt", "w").write("\n".join(keep) + "\n")
C = pd.DataFrame(rows); C = C[C.gwas.isin(keep)]; C.to_csv(f"{W}/counts_survey.tsv", sep="\t", index=False)
print(SRC, len(keep), "survey GWAS; dropped", len(rows) - len(keep))
