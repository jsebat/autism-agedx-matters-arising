#!/usr/bin/env python3
"""Polygenic-score CCA on the FULL SPARK sample (use all data; Zhang-era SPARK and family structure as sensitivity analyses).
Same method as cca_axes.py (ridge CCA, 5-fold family CV, 200-permutation null, structure coefficients), between families, polygenic scores only.
Sample: SPARK autistic participants of European ancestry with PGS from full_probands.parquet (arrays + iWGS v1.1 + SPARK-WGS-2026_08),
PGS standardised within source to non-autistic relatives. Characteristics: the 34 used before (genetic diagnosis excluded).
SUBSET (env): all | zhang (participants available to Zhang et al.: arrays + iWGS v1.1) | new (SPARK-WGS-2026_08 only) | p0 | p1 | p2 (parents genotyped)
Output (aggregate): ~/rvh2/agedx/cca_full_{SUBSET}_{cv,structure}.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
SUB = os.environ.get("SUBSET", "all")
A = pd.read_parquet(f"{OUT}/full_probands.parquet")
pro = A[A.cohort.eq("SPARK") & A.case & A.predicted_ancestry.eq("EUR")].copy()
regage = pd.to_numeric(pro.sample_id.map(pd.read_csv(f"{P}/spark_registration.tsv", sep="	").drop_duplicates("subject_sp_id").set_index("subject_sp_id").age_at_registration_years), errors="coerce")
u22 = regage.lt(22)   # Zhang et al. restrict their within-case age-at-diagnosis analyses to age_at_registration_years < 22
pro = {"all": pro, "zhang": pro[pro.zhang_era], "new": pro[~pro.zhang_era], "p0": pro[pro.n_parents_genotyped.eq(0)],
       "p1": pro[pro.n_parents_genotyped.eq(1)], "p2": pro[pro.n_parents_genotyped.eq(2)], "u22": pro[u22], "p0u22": pro[u22 & pro.n_parents_genotyped.eq(0)],
       "zhangu22": pro[u22 & pro.zhang_era], **{f"src_{k}": pro[pro.source.eq(k)] for k in ["spark_gsa_v1", "spark_gsa_v2", "spark_twist", "spark_wgs", "spark_wgs_2026_08"]}}[SUB].reset_index(drop=True)
PG = ["asd", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
WITHIN = os.environ.get("WITHIN", "0") == "1"   # within families: pTDT = child score minus mean of both parents' scores (both parents genotyped)
if WITHIN:
    allA = pd.read_parquet(f"{OUT}/full_probands.parquet").drop_duplicates("sample_id").set_index("sample_id")
    pro = pro[pro.n_parents_genotyped.eq(2)].copy()
    for p_ in PG:
        mid = (pro.mother.map(allA[p_]) + pro.father.map(allA[p_])) / 2
        pro[p_] = pro[p_] - mid
    pro = pro.reset_index(drop=True); SUB = SUB + "_within"
pro = pro.rename(columns={"asd": "asd_2026"}); PG = ["asd_2026"] + PG[1:]
pro["sample"] = pro.sample_id; pro["family_key"] = pro.family.fillna(pro.sample_id)
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
Yd = pro[PG]; okr = Yd.notna().all(axis=1).to_numpy()
Xm, Ym, famk = Xi.to_numpy(float)[okr], Yd.to_numpy(float)[okr], fam[okr]
PARTIAL = os.environ.get("PARTIAL", "")   # "source": regress genotyping platform out of both blocks (platform is confounded with registration year and ASD-score version)
if PARTIAL in ("source", "source_family"):
    # "source": genotyping platform; "source_family": platform + family structure (0/1/2 parents genotyped), as indicator variables
    zc = ["source"] + (["n_parents_genotyped"] if PARTIAL == "source_family" else [])
    Z = pd.get_dummies(pro.loc[okr, zc].astype(str), drop_first=True).to_numpy(float); Z = np.column_stack([np.ones(len(Z)), Z])
    res = lambda M: M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]
    Xm, Ym = res(Xm), res(Ym)
    keep = Xm.std(0) > 1e-8; Xi = Xi.loc[:, keep]; Xm = Xm[:, keep]
    SUB = SUB + ("_partial" if PARTIAL == "source" else "_partialfam")
lams = [0.01, 0.05, 0.1, 0.3, 1.0]
if os.environ.get("SINGLE") == "1":
    # Strength of single characteristics: the same ridge CCA with one characteristic as the whole characteristics block (= out-of-sample polygenic
    # prediction of that characteristic), observed values only, same adjustment, penalty choice, 5-fold family CV; 200-permutation null.
    import sys
    Xraw = X.loc[okr].reset_index(drop=True); Zs = Z if PARTIAL in ("source", "source_family") else np.ones((okr.sum(), 1)); rows = []
    for c in ["age at diagnosis", "deprivation (ADI)", "mother without degree", "low household income (third)", "maternal age", "minimal language", "ADHD", "birth year"]:
        m = Xraw[c].notna().to_numpy(); x = Xraw[c].to_numpy(float)[m][:, None]; zz = Zs[m]
        xr = x - zz @ np.linalg.lstsq(zz, x, rcond=None)[0]; yr = Ym[m]; fk = famk[m]
        sc = {l: cv(xr, yr, fk, l)[0] for l in lams}; l1 = max(sc, key=sc.get); r = cv(xr, yr, fk, l1)[0]
        rg_ = np.random.default_rng(7); nl = np.array([cv(xr, yr[rg_.permutation(len(yr))], fk, l1, seed=1)[0] for _ in range(200)])
        rows.append(dict(characteristic=c, n=int(m.sum()), n_families=len(np.unique(fk)), cv_r=r, shared_variance=r * r, null_95=np.percentile(nl, 95),
                         p_perm=(1 + (nl >= r).sum()) / 201, **{"lambda": l1}))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/cca_full_{SUB}_single.tsv", sep="\t", index=False); sys.exit(0)
scores = {l: cv(Xm, Ym, famk, l)[0] for l in lams}; lam = max(scores, key=scores.get)
obs = cv(Xm, Ym, famk, lam); rng = np.random.default_rng(99)
null = np.array([cv(Xm, Ym[rng.permutation(len(Ym))], famk, lam, seed=1) for _ in range(200)])
pv = [(1 + (null[:, j] >= obs[j]).sum()) / (1 + len(null)) for j in range(len(obs))]
CV = pd.DataFrame(dict(component=range(1, len(obs) + 1), cv_canonical_r=obs, null_mean=null.mean(0), null_95=np.percentile(null, 95, axis=0), p_perm=pv))
CV["lambda"] = lam; CV["n"] = len(Xm); CV["n_families"] = len(np.unique(famk)); CV.to_csv(f"{OUT}/cca_full_{SUB}_cv.tsv", sep="\t", index=False)
A_, B, d = fit(zs(Xm), zs(Ym), lam); U, V = zs(Xm) @ A_, zs(Ym) @ B
st = []
for j in range(min(4, A_.shape[1])):
    for i, c in enumerate(Xi.columns): st.append(dict(component=j + 1, block="characteristics", variable=c, r=float(np.corrcoef(U[:, j], zs(Xm)[:, i])[0, 1])))
    for i, c in enumerate(Yd.columns): st.append(dict(component=j + 1, block="genetics", variable=c, r=float(np.corrcoef(V[:, j], zs(Ym)[:, i])[0, 1])))
ST = pd.DataFrame(st); ST.to_csv(f"{OUT}/cca_full_{SUB}_structure.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(f"===== SUBSET {SUB}: n = {len(Xm)}, families = {len(np.unique(famk))}, lambda = {lam}"); print(CV.round(3).to_string(index=False))
for j in range(1, 5):
    x = ST[ST.component == j]; print(f"-- component {j}")
    for blk in ["genetics", "characteristics"]:
        y = x[(x.block == blk) & ~x.variable.str.contains(r"\[missing\]")].assign(a=lambda q: q.r.abs()).sort_values("a", ascending=False)
        print(f"   {blk}: " + "; ".join(f"{v} {r:+.2f}" for v, r in zip(y.variable.head(7), y.r.head(7))))
        if blk == "characteristics": print(f"   age at diagnosis r = {y.set_index('variable').r['age at diagnosis']:+.2f}, rank {list(y.variable).index('age at diagnosis') + 1}/{len(y)}")
# ---- power analysis (env POWER=1): each characteristic's correlation with the GENETIC variate of axes 1-4, observed values only
if os.environ.get("POWER") == "1":
    from scipy.stats import norm
    Vg = zs(Ym) @ B   # genetic-side canonical variates, full-sample fit
    Xraw = X.loc[okr].reset_index(drop=True)
    Zp = Z if PARTIAL in ("source", "source_family") else np.ones((okr.sum(), 1))
    zc = norm.ppf(1 - 0.05 / 2 / Xraw.shape[1]) + norm.ppf(0.80)   # Bonferroni over characteristics, 80% power
    rows = []
    for c in Xraw.columns:
        m = Xraw[c].notna().to_numpy(); n = int(m.sum())
        if n < 50: continue
        x = Xraw[c].to_numpy(float)[m]; zz = Zp[m]
        xr = x - zz @ np.linalg.lstsq(zz, x, rcond=None)[0]
        for j in range(min(4, Vg.shape[1])):
            v = Vg[m, j]; vr = v - zz @ np.linalg.lstsq(zz, v, rcond=None)[0]
            r = float(np.corrcoef(xr, vr)[0, 1]); se = 1 / np.sqrt(n - 3); fz = np.arctanh(r)
            rows.append(dict(characteristic=c, axis=j + 1, n=n, r=r, lo95=float(np.tanh(fz - 1.96 * se)), hi95=float(np.tanh(fz + 1.96 * se)),
                             p=float(2 * norm.sf(abs(fz) / se)), mde_r_80pct_bonf=float(np.tanh(zc * se))))
    PW = pd.DataFrame(rows); PW.to_csv(f"{OUT}/cca_full_{SUB}_power.tsv", sep="\t", index=False)
    print(PW[PW.axis == 1].sort_values("n").round(3).to_string(index=False))
    print(PW[(PW.characteristic == "age at diagnosis")].round(3).to_string(index=False))
# ---- missingness by family structure (env MISSBIAS=1): is recording non-random with respect to trios / single parent / case-only, and does it bias associations?
if os.environ.get("MISSBIAS") == "1":
    from scipy.stats import norm, chi2 as chi2d
    Vg = zs(Ym) @ B
    Xraw = X.loc[okr].reset_index(drop=True); fs = pro.loc[okr, "n_parents_genotyped"].to_numpy(); src = pro.loc[okr, "source"].astype(str).to_numpy()
    Zs = pd.get_dummies(pd.Series(src), drop_first=True).to_numpy(float); Zs = np.column_stack([np.ones(len(Zs)), Zs])
    resid = lambda y, z: y - z @ np.linalg.lstsq(z, y, rcond=None)[0]
    rows = []
    for c in Xraw.columns:
        rec = Xraw[c].notna().to_numpy()
        row = dict(characteristic=c, pct_recorded=100 * rec.mean(), **{f"pct_recorded_{k}par": 100 * rec[fs == k].mean() for k in (0, 1, 2)})
        # missingness vs genetic variates (platform regressed out)
        for j in (0, 1):
            row[f"r_missing_axis{j + 1}"] = float(np.corrcoef(resid((~rec).astype(float), Zs), resid(Vg[:, j], Zs))[0, 1])
        # association with the EDU/SES genetic variate within each family-structure stratum (platform regressed out) + heterogeneity
        est = []
        for k in (0, 1, 2):
            m = rec & (fs == k); n = int(m.sum())
            if n < 100: row[f"r_axis1_{k}par"] = np.nan; continue
            x = resid(Xraw[c].to_numpy(float)[m], Zs[m]); v = resid(Vg[m, 0], Zs[m]); r = float(np.corrcoef(x, v)[0, 1])
            row[f"r_axis1_{k}par"] = r; row[f"n_{k}par"] = n; est.append((np.arctanh(r), 1 / (n - 3)))
        if len(est) >= 2:
            z = np.array([e[0] for e in est]); w = 1 / np.array([e[1] for e in est]); zb = (w * z).sum() / w.sum(); Q = float((w * (z - zb) ** 2).sum())
            row["Q_het"] = Q; row["p_het"] = float(chi2d.sf(Q, len(est) - 1))
        rows.append(row)
    MB = pd.DataFrame(rows); MB.to_csv(f"{OUT}/cca_full_{SUB}_missbias.tsv", sep="\t", index=False)
    pd.set_option("display.width", 260); print(MB.sort_values("pct_recorded").round(3).to_string(index=False))
