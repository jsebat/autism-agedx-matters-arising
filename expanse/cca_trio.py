#!/usr/bin/env python3
"""Combined CCA (polygenic scores + rare variants) on the EXPANDED trio probands, main sample = registered before 22 (Zhang et al.'s criterion).
Probands: SPARK trio probands from rare_full_cases.parquet (legacy expanded trios + new-WGS sibling-family trios) with PGS in full_probands.parquet.
SSC is not included: its characteristics differ from SPARK's 34. Blocks: YSET = all | pgs | rare on the SAME probands.
Adjustment: genotyping platform, sequencing cohort and family structure (0/1/2 parents genotyped) regressed out of both blocks.
Output (aggregate): ~/rvh2/agedx/cca_trio_{YSET}{_all}_{cv,structure}.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
YSET = os.environ.get("YSET", "all"); U22 = os.environ.get("U22", "1") == "1"
R = pd.read_parquet(f"{OUT}/rare_full_cases.parquet"); R = R[~R.cohort.eq("ssc")]
if U22: R = R[R.u22]
A = pd.read_parquet(f"{OUT}/full_probands.parquet")
A = A[A.cohort.eq("SPARK") & A.predicted_ancestry.eq("EUR")].drop_duplicates("sample_id")
RC = ["de novo, NDD genes", "de novo, other genes", "inherited, NDD genes", "inherited, other genes", "inherited LoF tier 1"]
pro = R[["sp", "cohort", "family_key"] + RC].rename(columns={"cohort": "seq_cohort"}).merge(A, left_on="sp", right_on="sample_id", how="inner")
PG = ["asd", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
pro = pro.rename(columns={"asd": "asd_2026"}); PG = ["asd_2026"] + PG[1:]
pro["sample"] = pro.sample_id; pro = pro.reset_index(drop=True)
SUB = YSET + ("" if U22 else "_all")
print("probands with rare + PGS:", len(pro), pro.seq_cohort.value_counts().to_dict(), flush=True)
rd = lambda f: pd.read_csv(f"{P}/{f}", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
dx, sev, ho, env, reg, sv = rd("spark_dx_derived.tsv"), rd("spark_severity.tsv"), rd("spark_heldout.tsv"), rd("spark_env.tsv"), rd("spark_registration.tsv"), rd("spark_survey.tsv")
s = pro["sample"]; num = lambda ser: pd.to_numeric(s.map(ser), errors="coerce")
XC = {"age at diagnosis": num(dx.dx_age_y), "age at registration": num(reg.age_at_registration_years), "registration year": num(reg.registration_year),
      "birth year": num(env.birth_year), "maternal age": num(env.maternal_age), "paternal age": num(env.paternal_age), "deprivation (ADI)": num(env.adi_pct),
      "FSIQ": num(sev.fsiq), "SCQ": num(sev.scq_res), "RBS-R": num(sev.rbsr_res), "SRS": num(sev.srs_res), "Vineland ABC": num(sev.vabc_res),
      "CBCL internalizing": num(sev.cbcl_int_res), "CBCL externalizing": num(sev.cbcl_ext_res), "DCDQ": num(sev.dcdq_res),
      "age walked": num(sev.walked_age_mos), "age first words": num(sev.used_words_age_mos), "prenatal PM2.5": num(env.prenatal_pm25)}
XB = {"female": num(dx.sex.map({"Female": 1, "Male": 0})), "intellectual disability": num(dx.ID), "minimal language": num(dx.minimal_language),
      "language regression": num(env.regression), "ADHD": num(ho.behav_adhd), "mood/anxiety": num(ho.mood_or_anx), "seizures": num(ho.seizures),
      "sleep disorder": num(ho.sleep_dx), "feeding disorder": num(ho.feeding_dx), "multiplex family": num(ho.multiplex), "specialist diagnosis": num(ho.dx_specialist),
      "validity flag": num(dx.validity_flag), "low household income (third)": num(sv.s09),
      "mother without degree": num(sv.s10), "ODD/conduct": num(sv.s29), "macrocephaly": num(sv.s37)}
X = pd.DataFrame({**XC, **XB}, index=pro.index)
miss = X.isna().mean(); Xi = X.copy()
for c in X.columns:
    Xi[c] = X[c].fillna(X[c].median() if c in XC else 0.5)
    if miss[c] > 0.2: Xi[c + " [missing]"] = X[c].isna().astype(float)
Xi = Xi.loc[:, Xi.std() > 0]
fam = pd.factorize(pro.family_key)[0]
def zs(M): M = np.asarray(M, float); return (M - M.mean(0)) / M.std(0)
def fit(Xm, Ym, lam, k=6):
    k = min(k, Ym.shape[1], Xm.shape[1]); n = len(Xm); Sxx = Xm.T @ Xm / n + lam * np.eye(Xm.shape[1]); Syy = Ym.T @ Ym / n + lam * np.eye(Ym.shape[1]); Sxy = Xm.T @ Ym / n
    ix = np.linalg.inv(np.linalg.cholesky(Sxx)); iy = np.linalg.inv(np.linalg.cholesky(Syy))
    U, d, Vt = np.linalg.svd(ix @ Sxy @ iy.T); A_ = ix.T @ U[:, :k]; B = iy.T @ Vt.T[:, :k]; return A_, B, d[:k]
def cv(Xm, Ym, famk, lam, k=6, seed=0):
    k = min(k, Ym.shape[1], Xm.shape[1])
    rng = np.random.default_rng(seed); uf = np.unique(famk); fold = dict(zip(uf, rng.integers(0, 5, len(uf)))); f = np.array([fold[x] for x in famk])
    U, V = np.zeros((len(Xm), k)), np.zeros((len(Xm), k))
    for q in range(5):
        tr, te = f != q, f == q; mx, sx = Xm[tr].mean(0), Xm[tr].std(0) + 1e-9; my, sy = Ym[tr].mean(0), Ym[tr].std(0) + 1e-9
        A_, B, _ = fit((Xm[tr] - mx) / sx, (Ym[tr] - my) / sy, lam, k); U[te] = ((Xm[te] - mx) / sx) @ A_; V[te] = ((Ym[te] - my) / sy) @ B
    return np.array([abs(np.corrcoef(U[:, j], V[:, j])[0, 1]) for j in range(k)])
Yd = pro[{"all": PG + RC, "pgs": PG, "rare": RC}[YSET]]; okr = Yd.notna().all(axis=1).to_numpy()
Xm, Ym, famk = Xi.to_numpy(float)[okr], Yd.to_numpy(float)[okr], fam[okr]
Z = pd.get_dummies(pro.loc[okr, ["source", "seq_cohort", "n_parents_genotyped"]].astype(str), drop_first=True).to_numpy(float); Z = np.column_stack([np.ones(len(Z)), Z])
res = lambda M: M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]
Xm, Ym = res(Xm), res(Ym); keep = Xm.std(0) > 1e-8; Xi = Xi.loc[:, keep]; Xm = Xm[:, keep]
lams = [0.01, 0.05, 0.1, 0.3, 1.0]
if os.environ.get("SINGLE") == "1":
    # Strength of age at diagnosis alone: the same ridge CCA with age at diagnosis as the whole characteristics block (out-of-sample genetic prediction),
    # observed values only, same adjustment, penalty choice and 5-fold family CV; 200-permutation null. Output cca_trio_<SUB>_single.tsv
    import sys
    Xraw = X.loc[okr].reset_index(drop=True); rows = []
    for c in ["age at diagnosis"]:
        m = Xraw[c].notna().to_numpy(); x = Xraw[c].to_numpy(float)[m][:, None]; zz = Z[m]
        xr = x - zz @ np.linalg.lstsq(zz, x, rcond=None)[0]; yr = Ym[m]; fk = famk[m]
        sc = {l: cv(xr, yr, fk, l)[0] for l in lams}; l1 = max(sc, key=sc.get); r = cv(xr, yr, fk, l1)[0]
        rg_ = np.random.default_rng(7); nl = np.array([cv(xr, yr[rg_.permutation(len(yr))], fk, l1, seed=1)[0] for _ in range(200)])
        rows.append(dict(characteristic=c, n=int(m.sum()), n_families=len(np.unique(fk)), cv_r=r, shared_variance=r * r, null_95=np.percentile(nl, 95),
                         p_perm=(1 + (nl >= r).sum()) / 201, **{"lambda": l1}))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/cca_trio_{SUB}_single.tsv", sep="\t", index=False); sys.exit(0)
scores = {l: cv(Xm, Ym, famk, l)[0] for l in lams}; lam = max(scores, key=scores.get)
obs = cv(Xm, Ym, famk, lam); rng = np.random.default_rng(99)
null = np.array([cv(Xm, Ym[rng.permutation(len(Ym))], famk, lam, seed=1) for _ in range(200)])
pv = [(1 + (null[:, j] >= obs[j]).sum()) / (1 + len(null)) for j in range(len(obs))]
CV = pd.DataFrame(dict(component=range(1, len(obs) + 1), cv_canonical_r=obs, null_mean=null.mean(0), null_95=np.percentile(null, 95, axis=0), p_perm=pv))
CV["lambda"] = lam; CV["n"] = len(Xm); CV["n_families"] = len(np.unique(famk)); CV.to_csv(f"{OUT}/cca_trio_{SUB}_cv.tsv", sep="\t", index=False)
A_, B, d = fit(zs(Xm), zs(Ym), lam); U, V = zs(Xm) @ A_, zs(Ym) @ B
st = []
for j in range(min(4, A_.shape[1])):
    for i, c in enumerate(Xi.columns): st.append(dict(component=j + 1, block="characteristics", variable=c, r=float(np.corrcoef(U[:, j], zs(Xm)[:, i])[0, 1])))
    for i, c in enumerate(Yd.columns): st.append(dict(component=j + 1, block="genetics", variable=c, r=float(np.corrcoef(V[:, j], zs(Ym)[:, i])[0, 1])))
ST = pd.DataFrame(st); ST.to_csv(f"{OUT}/cca_trio_{SUB}_structure.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(f"===== SUBSET {SUB}: n = {len(Xm)}, families = {len(np.unique(famk))}, lambda = {lam}"); print(CV.round(3).to_string(index=False))
for j in range(1, 5):
    x = ST[ST.component == j]; print(f"-- component {j}")
    for blk in ["genetics", "characteristics"]:
        y = x[(x.block == blk) & ~x.variable.str.contains(r"\[missing\]")].assign(a=lambda q: q.r.abs()).sort_values("a", ascending=False)
        print(f"   {blk}: " + "; ".join(f"{v} {r:+.2f}" for v, r in zip(y.variable.head(7), y.r.head(7))))
        if blk == "characteristics": print(f"   age at diagnosis r = {y.set_index('variable').r['age at diagnosis']:+.2f}, rank {list(y.variable).index('age at diagnosis') + 1}/{len(y)}")
