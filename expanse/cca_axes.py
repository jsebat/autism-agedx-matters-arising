#!/usr/bin/env python3
"""(YSET=all|pgs|rare selects the genetics block.) Strongest axes linking characteristics of autistic people to their genetics, common and rare combined (Matters Arising on Zhang,
Grove et al. 2025). Regularized canonical correlation analysis (CCA) between
  X  characteristics of SPARK autistic probands (age at diagnosis, registration, sex, parental ages, SES, deprivation, cognition,
     development, severity, psychiatric and medical comorbidity, design variables)
  Y  genetics: polygenic scores (ASD 2026, ADHD, MDD, bipolar, SCZ, EA, IQ, age at first birth, height) + rare-variant classes
     (de novo NDD / other genes, inherited NDD / other genes, inherited LoF tier 1)
Two versions: BETWEEN families (absolute PGS standardized to non-autistic relatives; absolute rare burden standardized to siblings) and
WITHIN families (pTDT for PGS; proband minus sibling for rare variants).
Sample: SPARK EUR probands in complete trios with an unaffected sibling and both genotype and sequence data.
Missing X: median imputation (continuous) or 0.5 (binary / tertile codes) plus a missingness indicator for variables > 20% missing.
CCA: ridge-regularized (lambda chosen by 5-fold family CV of the first canonical correlation); out-of-sample canonical correlations
by 5-fold CV for the first 6 components; permutation null (Y rows shuffled, 200 permutations) for the CV correlations.
Structure coefficients (correlation of each variable with the canonical variate, full-sample fit) define the axes.
Output (aggregate): ~/rvh2/agedx/cca_{between,within}_{cv,structure}.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
# rare-variant tier definitions (SNV/indel), as in the lab's exploratory rare_survey.py
TI = ["lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4"]
cell = lambda tiers, modes, sets: [f"{t}__{m}__{g}" for t in tiers for m in modes for g in sets]
CAT = {"de novo, NDD genes": cell(TI, ["dnm"], ["ndd"]), "de novo, other genes": cell(TI, ["dnm"], ["nonndd"]),
       "de novo LoF": cell(["lof_t1", "lof_t2"], ["dnm"], ["ndd", "nonndd"]), "de novo missense": cell(TI[2:], ["dnm"], ["ndd", "nonndd"]),
       "inherited, NDD genes": cell(TI, ["inh"], ["ndd"]), "inherited, other genes": cell(TI, ["inh"], ["nonndd"]),
       "inherited LoF tier 1": cell(["lof_t1"], ["inh"], ["ndd", "nonndd"]), "inherited missense": cell(TI[2:], ["inh"], ["ndd", "nonndd"]),
       "all rare SNV/indel": cell(TI, ["dnm", "inh"], ["ndd", "nonndd"])}
def load(root, cats):
    a = pd.read_parquet(f"/expanse/projects/sebat1/jsebat/{root}/assembled_children.parquet"); a = a[a.cohort.str.startswith("spark")].copy()
    W = pd.read_csv(f"/expanse/projects/sebat1/jsebat/{root}/results/collider_weights.tsv", sep="\t", index_col=0).iloc[:, 0]
    for k, ts in cats.items():
        ts = [t for t in ts if t in a.columns and t in W.index]
        a[k] = (a[ts].fillna(0).to_numpy(float) * W[ts].to_numpy()).sum(1)
        ref = a[a.case.eq(0)].groupby("cohort")[k].agg(["mean", "std"]); a[k] = (a[k] - a.cohort.map(ref["mean"])) / a.cohort.map(ref["std"])
    return a
RC = ["de novo, NDD genes", "de novo, other genes", "inherited, NDD genes", "inherited, other genes", "inherited LoF tier 1"]
a = load("pgs_rare_h2_drift_full", {k: CAT[k] for k in RC})
sibm = a[a.case.eq(0)].groupby(["cohort", "family_key"])[RC].mean()
pro = a[a.case.eq(1)].copy(); m = sibm.reindex(pd.MultiIndex.from_frame(pro[["cohort", "family_key"]])); ok = m.notna().all(axis=1).to_numpy()
pro = pro[ok].copy()
for c in RC: pro[c + " (excess)"] = pro[c].to_numpy() - m[ok][c].to_numpy()
cs = pd.read_parquet(f"{OUT}/cases.parquet").drop_duplicates("sample_id").set_index("sample_id")
PG = ["asd_2026", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
for p in PG: pro[p] = pro["sample"].map(cs[p]); pro[p + "_ptdt"] = pro["sample"].map(cs[p + "_ptdt"])
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
      "validity flag": num(dx.validity_flag), "genetic diagnosis": num(ho.genetic_dx_reported), "low household income (third)": num(sv.s09),
      "mother without degree": num(sv.s10), "ODD/conduct": num(sv.s29), "macrocephaly": num(sv.s37)}
X = pd.DataFrame({**XC, **XB}, index=pro.index)
DROP = [v for v in os.environ.get("DROP", "").split(",") if v]   # e.g. DROP="genetic diagnosis": reported from the same sequencing, circular for rare variants
X = X.drop(columns=DROP)
miss = X.isna().mean(); Xi = X.copy()
for c in X.columns:
    Xi[c] = X[c].fillna(X[c].median() if c in XC else 0.5)
    if miss[c] > 0.2: Xi[c + " [missing]"] = X[c].isna().astype(float)
Xi = Xi.loc[:, Xi.std() > 0]
YSET = os.environ.get("YSET", "all")   # all | pgs | rare
GB = {"all": PG + RC, "pgs": PG, "rare": RC}[YSET]; GW = {"all": [p + "_ptdt" for p in PG] + [c + " (excess)" for c in RC], "pgs": [p + "_ptdt" for p in PG], "rare": [c + " (excess)" for c in RC]}[YSET]
YB = pro[GB]; YW = pro[GW]; SUF = ("" if YSET == "all" else "_" + YSET) + ("_common" if os.environ.get("COMMON") == "1" else "") + ("_nogdx" if "genetic diagnosis" in DROP else "")
fam =pd.factorize(pro.cohort + "|" + pro.family_key)[0]
def zs(M): M = np.asarray(M, float); return (M - M.mean(0)) / M.std(0)
def fit(Xm, Ym, lam, k=6):
    k = min(k, Ym.shape[1], Xm.shape[1]); n = len(Xm); Sxx = Xm.T @ Xm / n + lam * np.eye(Xm.shape[1]); Syy = Ym.T @ Ym / n + lam * np.eye(Ym.shape[1]); Sxy = Xm.T @ Ym / n
    ix = np.linalg.inv(np.linalg.cholesky(Sxx)); iy = np.linalg.inv(np.linalg.cholesky(Syy))
    U, d, Vt = np.linalg.svd(ix @ Sxy @ iy.T); A = ix.T @ U[:, :k]; B = iy.T @ Vt.T[:, :k]; return A, B, d[:k]
def cv(Xm, Ym, lam, k=6, seed=0):
    k = min(k, Ym.shape[1], Xm.shape[1])
    rng = np.random.default_rng(seed); uf = np.unique(fam); fold = dict(zip(uf, rng.integers(0, 5, len(uf)))); f = np.array([fold[x] for x in fam])
    U, V = np.zeros((len(Xm), k)), np.zeros((len(Xm), k))
    for q in range(5):
        tr, te = f != q, f == q; mx, sx = Xm[tr].mean(0), Xm[tr].std(0) + 1e-9; my, sy = Ym[tr].mean(0), Ym[tr].std(0) + 1e-9
        A, B, _ = fit((Xm[tr] - mx) / sx, (Ym[tr] - my) / sy, lam, k); U[te] = ((Xm[te] - mx) / sx) @ A; V[te] = ((Ym[te] - my) / sy) @ B
    return np.array([abs(np.corrcoef(U[:, j], V[:, j])[0, 1]) for j in range(k)])
for tag, Yd in [("between", YB), ("within", YW)]:
    okr = Yd.notna().all(axis=1).to_numpy().copy()
    if os.environ.get("COMMON") == "1":   # restrict to probands with polygenic scores, so rare-only CCA uses the same sample as PGS-only
        okr &= pro[[p if tag == "between" else p + "_ptdt" for p in PG]].notna().all(axis=1).to_numpy()
    Xm, Ym = Xi.to_numpy(float)[okr], Yd.to_numpy(float)[okr]; famk = fam[okr]
    fam_save = fam; fam = famk
    lams = [0.01, 0.05, 0.1, 0.3, 1.0]; scores = {l: cv(Xm, Ym, l)[0] for l in lams}; lam = max(scores, key=scores.get)
    obs = cv(Xm, Ym, lam); rng = np.random.default_rng(99)
    null = np.array([cv(Xm, Ym[rng.permutation(len(Ym))], lam, seed=1) for _ in range(200)])
    pv = [(1 + (null[:, j] >= obs[j]).sum()) / (1 + len(null)) for j in range(len(obs))]
    CV = pd.DataFrame(dict(component=range(1, len(obs) + 1), cv_canonical_r=obs, null_mean=null.mean(0), null_95=np.percentile(null, 95, axis=0), p_perm=pv))
    CV["lambda"] = lam; CV["n"] = len(Xm); CV.to_csv(f"{OUT}/cca_{tag}{SUF}_cv.tsv", sep="\t", index=False)
    A, B, d = fit(zs(Xm), zs(Ym), lam); U, V = zs(Xm) @ A, zs(Ym) @ B
    st = []
    for j in range(min(4, A.shape[1])):
        sgn = np.sign(np.corrcoef(V[:, j], zs(Ym)[:, list(Yd.columns).index(Yd.columns[0])])[0, 1]) or 1
        for i, c in enumerate(Xi.columns): st.append(dict(component=j + 1, block="characteristics", variable=c, r=float(np.corrcoef(U[:, j], zs(Xm)[:, i])[0, 1])))
        for i, c in enumerate(Yd.columns): st.append(dict(component=j + 1, block="genetics", variable=c, r=float(np.corrcoef(V[:, j], zs(Ym)[:, i])[0, 1])))
    ST = pd.DataFrame(st); ST.to_csv(f"{OUT}/cca_{tag}{SUF}_structure.tsv", sep="\t", index=False)
    fam = fam_save
    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 400)
    print(f"\n===== {tag.upper()} families: n = {len(Xm)}, lambda = {lam}"); print(CV.round(3).to_string(index=False))
    for j in range(1, 5):
        x = ST[ST.component == j]; print(f"\n-- component {j}")
        for blk in ["genetics", "characteristics"]:
            y = x[x.block == blk].assign(a=lambda q: q.r.abs()).sort_values("a", ascending=False).head(8)
            print(f"   {blk}: " + "; ".join(f"{v} {r:+.2f}" for v, r in zip(y.variable, y.r)))
        ad = x[x.variable == "age at diagnosis"].r.iloc[0]; print(f"   age at diagnosis structure r = {ad:+.2f}")
