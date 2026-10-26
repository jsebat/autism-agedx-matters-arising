#!/usr/bin/env python3
"""Split-GWAS phenotypes for the strata beyond the three array platforms (use all data): SRC = spark_wgs (iWGS v1.1),
spark_wgs_2026_08, ssc. Same design as rep_pheno.py + survey_pheno.py (their pipeline: cases vs ALL non-autistic relatives; disjoint control
halves C1/C2 split by family), so that each new stratum enters the per-stratum inverse-variance meta-analysis next to GSA v1, GSA v2, Twist.
One row per person across ALL strata: iWGS v1.1 keeps only people without array genotypes; SPARK-WGS-2026_08 only people in neither; SSC is
separate. Others are set to NA (genotype files, GRM and PCs unchanged).
SPARK strata: W_all, W_under6, W_over10, W_under6_c1 / W_over10_c2 / W_over10_c1, 4 random-split null pairs, and every survey split (spark_survey.tsv)
as <code>_hi_c2 / <code>_lo_c1. SSC: W_all and the survey splits SSC can define (thirds within SSC): s09-s13, s16-s18, s23-s26, s31-s33, s40-s41.
qcovar: 10 PC-AiR-style PCs (+ age at registration for SPARK); covar: sex.
Output: warrier_rep/<src>/{allpheno.txt, allpheno_names.txt, qcovar.txt, covar.txt, counts_all.tsv}"""
import os, hashlib
import numpy as np, pandas as pd
SRC = os.environ["SRC"]; W = f"/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/{SRC}"; R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
B = "/expanse/projects/sebat1/s3/data/sebat"; J = "/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08"
fam = pd.read_csv(f"{W}/hm3w.fam", sep=r"\s+", header=None, dtype=str).iloc[:, :2]; fam.columns = ["FID", "IID"]
old = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "sample_id"])
arrays = set(old[old.source.isin(["spark_gsa_v1", "spark_gsa_v2", "spark_twist"])].sample_id); iwgs = set(old[old.source.eq("spark_wgs")].sample_id)
half = lambda famid: famid.map(lambda x: int(hashlib.md5(str(x).encode()).hexdigest(), 16) % 2 == 0)   # deterministic family-level C1/C2
f = fam.copy(); rows = []
tert = lambda v: pd.Series(np.where(v >= v.quantile(2 / 3), 1.0, np.where(v <= v.quantile(1 / 3), 0.0, np.nan)), index=v.index)
if SRC.startswith("spark"):
    if SRC == "spark_wgs_2026_08":
        meta = pd.read_csv(f"{J}/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv", sep="\t", dtype=str).set_index("sample_id")
        f["sp"] = f.IID.map(meta.ind_id); excl = arrays | iwgs
    else:
        d = old[old.source.eq(SRC)].drop_duplicates("IID").set_index("IID"); f["sp"] = f.IID.map(d.sample_id); excl = arrays
    reg = pd.read_csv(f"{R}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
    ped = pd.read_csv(f"{R}/spark_pedigree_ids.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
    sv = pd.read_csv(f"{R}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
    new = ~f.sp.isin(excl) & f.sp.notna()
    f["asd"] = f.sp.map(reg.asd).astype(str).str.upper().eq("TRUE"); known = f.sp.isin(reg.index) & new
    case, ctrl = f.asd & known, (~f.asd) & known
    dxm = pd.to_numeric(f.sp.map(reg.diagnosis_age), errors="coerce")
    famid = f.sp.map(reg.family_sf_id).fillna(f.IID); c1 = ctrl & half(famid); c2 = ctrl & ~half(famid)
    S = {"W_all": (case, ctrl), "W_under6": (case & (dxm < 72), ctrl), "W_over10": (case & (dxm > 120), ctrl),
         "W_under6_c1": (case & (dxm < 72), c1), "W_over10_c2": (case & (dxm > 120), c2), "W_over10_c1": (case & (dxm > 120), c1)}
    pool = np.where(case & dxm.notna())[0]; n_o, n_u = int((case & (dxm > 120)).sum()), int((case & (dxm < 72)).sum())
    for r in range(1, 5):
        perm = np.random.default_rng(1000 + r).permutation(pool); ma = np.zeros(len(f), bool); ma[perm[:n_o]] = True; mb = np.zeros(len(f), bool); mb[perm[n_o:n_o + n_u]] = True
        S[f"rand{r}_A_c2"] = (pd.Series(ma, index=f.index), c2); S[f"rand{r}_B_c1"] = (pd.Series(mb, index=f.index), c1)
    geno = set(f.sp[new].dropna()); npar = f.sp.map(ped.biomother_sp_id).isin(geno).astype(int) + f.sp.map(ped.biofather_sp_id).isin(geno).astype(int)
    codes = {c: pd.to_numeric(f.sp.map(sv[c]), errors="coerce") for c in sv.columns}
    codes["s99"] = pd.Series(np.where(npar == 0, 1.0, np.where(npar == 2, 0.0, np.nan)), index=f.index)
    areg = pd.to_numeric(f.sp.map(reg.age_at_registration_years), errors="coerce"); sex = f.sp.map(reg.sex).map({"Male": "M", "Female": "F"})
else:   # SSC
    ids = pd.read_csv(f"{R}/ssc_ids_sanders.csv", dtype=str, encoding="utf-8-sig").drop_duplicates("SS_RUID").set_index("SS_RUID")
    rel = f.IID.map(ids.relation); known = rel.notna()
    case = known & rel.astype(str).str.match(r"^p\d"); ctrl = known & ~case
    famid = f.IID.map(ids.family).fillna(f.IID); c1 = ctrl & half(famid); c2 = ctrl & ~half(famid)
    S = {"W_all": (case, ctrl)}
    X = pd.read_csv(f"{R}/ssc_cca_chars.tsv", sep="\t", index_col=0).drop_duplicates("SS_RUID").set_index("SS_RUID")
    hz = pd.read_csv(f"{B}/SFARI_phenotypes/spark_ssc_phenotypes_child.csv", low_memory=False)
    hz = hz[hz.cohort.str.upper().eq("SSC")].dropna(subset=["repository_id"]).drop_duplicates("repository_id").set_index("repository_id")
    senv = pd.read_csv(f"{R}/ssc_env.tsv", sep="\t").drop_duplicates("id").set_index("id"); sder = pd.read_csv(f"{R}/ssc_derived.tsv", sep="\t").drop_duplicates("id").set_index("id")
    ind = f.IID.map(ids.individual); cs = case
    def onc(v):   # thirds / binary defined among SSC cases only
        v = v.where(cs); return v
    num = lambda s: pd.to_numeric(s, errors="coerce")
    me = f.IID.map(hz.mother_education).astype(str).str.lower(); deg = me.isin(["baccalaureate", "graduate"])
    codes = {"s09": 1 - tert(onc(num(f.IID.map(hz.household_income)))), "s10": pd.Series(np.where(f.IID.map(hz.mother_education).isna(), np.nan, (~deg).astype(float)), index=f.index).where(cs),
             "s11": tert(onc(num(f.IID.map(hz.mother_age_at_birth)))), "s12": tert(onc(num(f.IID.map(hz.father_age_at_birth)))),
             "s13": onc(num(f.IID.map(X["intellectual disability"]))), "s16": tert(onc(num(f.IID.map(X["age walked"])))), "s17": tert(onc(num(f.IID.map(X["age first words"])))),
             "s18": onc(num(f.IID.map(X["language regression"]))), "s23": tert(onc(num(ind.map(sder.scq)))), "s24": tert(onc(num(f.IID.map(X["RBS-R"])))),
             "s25": tert(onc(num(f.IID.map(X["SRS"])))), "s26": 1 - tert(onc(num(f.IID.map(X["Vineland ABC"])))), "s31": tert(onc(num(f.IID.map(X["CBCL internalizing"])))),
             "s32": tert(onc(num(f.IID.map(X["CBCL externalizing"])))), "s33": onc(num(f.IID.map(X["seizures"]))),
             "s40": tert(onc(num(ind.map(senv.prenatal_pm25_mean)))), "s41": tert(onc(num(ind.map(senv.prenatal_no2_mean))))}
    areg = None; sex = f.IID.map(ids.sex).astype(str).str.lower().map({"male": "M", "female": "F"})
for c, v in codes.items():
    S[f"{c}_hi_c2"] = (case & v.eq(1), c2); S[f"{c}_lo_c1"] = (case & v.eq(0), c1)
P = f[["FID", "IID"]].copy(); keep = []
for k, (m, cc) in S.items():
    m = m.fillna(False).astype(bool); cc = cc.fillna(False).astype(bool)
    n = int(m.sum()); rows.append(dict(source=SRC, gwas=k, cases=n, controls=int(cc.sum())))
    if n >= 150: P[k] = np.where(m, 1, np.where(cc, 0, np.nan)); keep.append(k)
P[["FID", "IID"] + keep].fillna("NA").to_csv(f"{W}/allpheno.txt", sep=" ", index=False, header=False)
open(f"{W}/allpheno_names.txt", "w").write("\n".join(keep) + "\n")
pc = pd.read_csv(f"{W}/pcair.sscore", sep="\t"); ic = [c for c in pc.columns if "IID" in c][0]; pcs = [c for c in pc.columns if c.endswith("_AVG")][:10]
q = f[["FID", "IID"]].merge(pc[[ic] + pcs].rename(columns={ic: "IID"}), on="IID", how="left")
if areg is not None: q["areg"] = areg.values
q.fillna("NA").to_csv(f"{W}/qcovar.txt", sep=" ", index=False, header=False)
cv = f[["FID", "IID"]].copy(); cv["sex"] = sex.values; cv.fillna("NA").to_csv(f"{W}/covar.txt", sep=" ", index=False, header=False)
C = pd.DataFrame(rows); C = C[C.gwas.isin(keep)]; C.to_csv(f"{W}/counts_all.tsv", sep="\t", index=False)
print(SRC, "samples", len(f), "| cases", int(case.sum()), "controls", int(ctrl.sum()), "| GWAS kept", len(keep), "of", len(S))
print(C.head(12).to_string(index=False))
