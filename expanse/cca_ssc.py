#!/usr/bin/env python3
"""SSC replication of the polygenic (and rare) CCA axes, and SPARK on the same harmonised characteristics (include SSC).
Harmonised characteristics (SSC and SPARK): female, FSIQ, intellectual disability, SRS, RBS-R, Vineland ABC, CBCL internalizing/externalizing,
language regression, seizures, age walked, age first words, maternal age, paternal age, mother without degree, low household income (bottom third
within cohort). SSC has no age at diagnosis and no ADI, so SSC tests the EDU/SES and DEV axes, not AADx.
COHORT (env): ssc | spark (SPARK, registered before 22, platform + family structure regressed out). YSET: pgs | all (pgs + rare; SSC sibling-family
trio probands). Same ridge CCA, 5-fold family CV and 200-permutation null as cca_full.py. Output: ~/rvh2/agedx/cca_harm_{COHORT}_{YSET}_{cv,structure}.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
COH = os.environ.get("COHORT", "ssc"); YSET = os.environ.get("YSET", "pgs")
A = pd.read_parquet(f"{OUT}/full_probands.parquet").drop_duplicates("sample_id")
PG = ["asd", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
RC = ["de novo, NDD genes", "de novo, other genes", "inherited, NDD genes", "inherited, other genes", "inherited LoF tier 1"]
num = lambda s: pd.to_numeric(s, errors="coerce")
if COH == "ssc":
    pro = A[A.cohort.eq("SSC") & A.case & A.predicted_ancestry.eq("EUR")].copy().reset_index(drop=True)
    S = pd.read_csv(f"{P}/ssc_cca_chars.tsv", sep="\t", index_col=0).drop_duplicates("SS_RUID").set_index("SS_RUID")
    hz = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/SFARI_phenotypes/spark_ssc_phenotypes_child.csv", low_memory=False)
    hz = hz[hz.cohort.str.upper().eq("SSC")].dropna(subset=["repository_id"]).drop_duplicates("repository_id").set_index("repository_id")
    print("SSC mother_education:", hz.mother_education.value_counts().to_dict()); print("SSC household_income:", hz.household_income.value_counts().sort_index().to_dict())
    s = pro.sample_id
    X = pd.DataFrame({c: num(s.map(S[c])) for c in S.columns}, index=pro.index)
    X["maternal age"] = num(s.map(hz.mother_age_at_birth)); X["paternal age"] = num(s.map(hz.father_age_at_birth))
    me = s.map(hz.mother_education).astype(str).str.lower()
    X["mother without degree"] = (~me.isin(["baccalaureate", "graduate", "post-graduate", "masters", "doctorate", "professional", "graduate-degree"])).astype(float).where(s.map(hz.mother_education).notna())
    inc = num(s.map(hz.household_income)); X["low household income (third)"] = (inc <= inc.quantile(1 / 3)).astype(float).where(inc.notna())
    fam = pd.factorize(pro.family.fillna(pro.sample_id))[0]; Z = np.ones((len(pro), 1))
    if YSET == "all":
        R = pd.read_parquet(f"{OUT}/rare_full_cases.parquet"); R = R[R.cohort.eq("ssc")].set_index("sp")
        for c in RC: pro[c] = pro.sample_id.map(R[c])
else:
    reg = pd.read_csv(f"{P}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
    pro = A[A.cohort.eq("SPARK") & A.case & A.predicted_ancestry.eq("EUR")].copy()
    pro = pro[num(pro.sample_id.map(reg.age_at_registration_years)).lt(22)].reset_index(drop=True)
    rd = lambda f: pd.read_csv(f"{P}/{f}", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
    dx, sev, ho, env, sv = rd("spark_dx_derived.tsv"), rd("spark_severity.tsv"), rd("spark_heldout.tsv"), rd("spark_env.tsv"), rd("spark_survey.tsv")
    s = pro.sample_id; m_ = lambda ser: num(s.map(ser))
    X = pd.DataFrame({"female": m_(dx.sex.map({"Female": 1, "Male": 0})), "FSIQ": m_(sev.fsiq), "intellectual disability": m_(dx.ID), "SRS": m_(sev.srs_res), "RBS-R": m_(sev.rbsr_res),
                      "Vineland ABC": m_(sev.vabc_res), "CBCL internalizing": m_(sev.cbcl_int_res), "CBCL externalizing": m_(sev.cbcl_ext_res), "language regression": m_(env.regression),
                      "seizures": m_(ho.seizures), "age walked": m_(sev.walked_age_mos), "age first words": m_(sev.used_words_age_mos), "maternal age": m_(env.maternal_age),
                      "paternal age": m_(env.paternal_age), "mother without degree": m_(sv.s10), "low household income (third)": m_(sv.s09)}, index=pro.index)
    fam = pd.factorize(pro.family.fillna(pro.sample_id))[0]
    Z = pd.get_dummies(pro[["source", "n_parents_genotyped"]].astype(str), drop_first=True).to_numpy(float); Z = np.column_stack([np.ones(len(Z)), Z])
print(COH, "probands with PGS:", len(pro))
miss = X.isna().mean(); Xi = X.copy()
for c in X.columns:
    Xi[c] = X[c].fillna(X[c].median())
    if miss[c] > 0.2: Xi[c + " [missing]"] = X[c].isna().astype(float)
Xi = Xi.loc[:, Xi.std() > 0]
Yd = pro[PG + (RC if YSET == "all" else [])].rename(columns={"asd": "asd_2026"}); okr = Yd.notna().all(axis=1).to_numpy()
Xm, Ym, famk = Xi.to_numpy(float)[okr], Yd.to_numpy(float)[okr], fam[okr]; Zk = Z[okr]
res = lambda M: M - Zk @ np.linalg.lstsq(Zk, M, rcond=None)[0]
Xm, Ym = res(Xm), res(Ym); keep = Xm.std(0) > 1e-8; Xi = Xi.loc[:, keep]; Xm = Xm[:, keep]
def zs(M): M = np.asarray(M, float); return (M - M.mean(0)) / M.std(0)
def fit(Xm, Ym, lam, k=6):
    k = min(k, Ym.shape[1], Xm.shape[1]); n = len(Xm); Sxx = Xm.T @ Xm / n + lam * np.eye(Xm.shape[1]); Syy = Ym.T @ Ym / n + lam * np.eye(Ym.shape[1]); Sxy = Xm.T @ Ym / n
    ix = np.linalg.inv(np.linalg.cholesky(Sxx)); iy = np.linalg.inv(np.linalg.cholesky(Syy))
    U, d, Vt = np.linalg.svd(ix @ Sxy @ iy.T); return ix.T @ U[:, :k], iy.T @ Vt.T[:, :k], d[:k]
def cv(Xm, Ym, famk, lam, k=6, seed=0):
    k = min(k, Ym.shape[1], Xm.shape[1]); rng = np.random.default_rng(seed); uf = np.unique(famk); fold = dict(zip(uf, rng.integers(0, 5, len(uf)))); f = np.array([fold[x] for x in famk])
    U, V = np.zeros((len(Xm), k)), np.zeros((len(Xm), k))
    for q in range(5):
        tr, te = f != q, f == q; mx, sx = Xm[tr].mean(0), Xm[tr].std(0) + 1e-9; my, sy = Ym[tr].mean(0), Ym[tr].std(0) + 1e-9
        A_, B, _ = fit((Xm[tr] - mx) / sx, (Ym[tr] - my) / sy, lam, k); U[te] = ((Xm[te] - mx) / sx) @ A_; V[te] = ((Ym[te] - my) / sy) @ B
    return np.array([abs(np.corrcoef(U[:, j], V[:, j])[0, 1]) for j in range(k)])
lams = [0.01, 0.05, 0.1, 0.3, 1.0]; sc = {l: cv(Xm, Ym, famk, l)[0] for l in lams}; lam = max(sc, key=sc.get)
obs = cv(Xm, Ym, famk, lam); rng = np.random.default_rng(99); null = np.array([cv(Xm, Ym[rng.permutation(len(Ym))], famk, lam, seed=1) for _ in range(200)])
pv = [(1 + (null[:, j] >= obs[j]).sum()) / (1 + len(null)) for j in range(len(obs))]
CV = pd.DataFrame(dict(component=range(1, len(obs) + 1), cv_canonical_r=obs, null_95=np.percentile(null, 95, axis=0), p_perm=pv)); CV["n"] = len(Xm)
TAGF = f"cca_harm_{COH}_{YSET}"; CV.to_csv(f"{OUT}/{TAGF}_cv.tsv", sep="\t", index=False)
A_, B, d = fit(zs(Xm), zs(Ym), lam); U, V = zs(Xm) @ A_, zs(Ym) @ B; st = []
for j in range(min(4, A_.shape[1])):
    for i, c in enumerate(Xi.columns): st.append(dict(component=j + 1, block="characteristics", variable=c, r=float(np.corrcoef(U[:, j], zs(Xm)[:, i])[0, 1])))
    for i, c in enumerate(Yd.columns): st.append(dict(component=j + 1, block="genetics", variable=c, r=float(np.corrcoef(V[:, j], zs(Ym)[:, i])[0, 1])))
ST = pd.DataFrame(st); ST.to_csv(f"{OUT}/{TAGF}_structure.tsv", sep="\t", index=False)
print(f"===== {TAGF}: n = {len(Xm)}, lambda = {lam}"); print(CV.round(3).head(4).to_string(index=False))
for j in range(1, 4):
    x = ST[ST.component == j]
    for blk in ["genetics", "characteristics"]:
        y = x[(x.block == blk) & ~x.variable.str.contains(r"\[missing\]")].assign(a=lambda q: q.r.abs()).sort_values("a", ascending=False)
        print(f"-- {j} {blk}: " + "; ".join(f"{v} {r:+.2f}" for v, r in zip(y.variable.head(7), y.r.head(7))))
