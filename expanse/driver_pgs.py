#!/usr/bin/env python3
"""Which characteristic of SPARK autistic cases drives polygenic heterogeneity? Age at diagnosis (the key stratum of Zhang,
Grove et al. 2025) vs every other characteristic we can measure.
Cases: EUR SPARK autistic individuals with PGS (cases.parquet from agedx_strata.py; PGS in SD of non-autistic relatives on the
same platform; pTDT in complete trios). Candidate splits (high vs low group; continuous measures split top vs bottom third):
  age at diagnosis (after 10 vs before 6, as in the paper), family structure (0 vs 2 parents genotyped), adult at registration,
  sex, diagnosis era (2019-26 vs <2010), platform (Twist vs GSA), ID, minimal language, FSIQ, SCQ, Vineland, RBS-R, CBCL total,
  ADHD, mood/anxiety, seizures, regression, multiplex family, diagnosis by specialist, validity flag, area deprivation.
1. Profile per split: PGS difference (high minus low) for 9 PGS, population-style and within-family (pTDT), family-bootstrap SE;
   divergence = root-mean-square standardized difference (sum z^2 / k), and cosine similarity of the population profile to the
   age-at-diagnosis profile.
2. Joint model (cases, population-style): each PGS ~ all characteristics (median-imputed + missing indicators, standardized);
   per characteristic: partial R^2 summed over the 9 PGS, and the age-at-diagnosis coefficient with and without the others.
Output (aggregate): ~/rvh2/agedx/driver_*.tsv"""
import os
import numpy as np, pandas as pd
R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"; OUT = os.path.expanduser("~/rvh2/agedx")
PGS = ["asd_2026", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
c = pd.read_parquet(f"{OUT}/cases.parquet")
sev = pd.read_csv(f"{R}/spark_severity.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
ho = pd.read_csv(f"{R}/spark_heldout.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
env = pd.read_csv(f"{R}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
dx = pd.read_csv(f"{R}/spark_dx_derived.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
for col in ["scq_res", "vabc_res", "rbsr_res", "cbcl_tot_res", "fsiq"]: c[col] = c.sample_id.map(sev[col])
for col in ["behav_adhd", "mood_or_anx", "seizures", "multiplex", "dx_specialist"]: c[col] = c.sample_id.map(ho[col])
c["regression"] = c.sample_id.map(env["regression"]); c["adi_pct"] = c.sample_id.map(env["adi_pct"])
c["validity_flag"] = c.sample_id.map(dx["validity_flag"]); c["ID"] = c.sample_id.map(dx["ID"]); c["minimal_language"] = c.sample_id.map(dx["minimal_language"])
c["platform_twist"] = np.where(c.source.eq("spark_twist"), 1.0, np.where(c.source.isin(["spark_gsa_v1", "spark_gsa_v2"]), 0.0, np.nan))
def b(x): return pd.to_numeric(x, errors="coerce")
def tert(x):
    x = b(x); lo, hi = x.quantile(1 / 3), x.quantile(2 / 3); return np.where(x >= hi, 1.0, np.where(x <= lo, 0.0, np.nan))
SPL = {
    "age at diagnosis (after 10 vs before 6)": np.where(c.stratum.eq("after 10"), 1.0, np.where(c.stratum.eq("before 6"), 0.0, np.nan)),
    "family structure (0 vs 2 parents genotyped)": np.where(c.n_par.eq(0), 1.0, np.where(c.n_par.eq(2), 0.0, np.nan)),
    "adult at registration": b(c.adult_at_reg),
    "female": b(c.female),
    "diagnosed 2019-26 vs before 2010": np.where(c.era.eq("2019-26"), 1.0, np.where(c.era.eq("<2010"), 0.0, np.nan)),
    "platform Twist vs GSA": c.platform_twist,
    "intellectual disability": b(c.ID), "minimal language": b(c.minimal_language), "FSIQ (top vs bottom third)": tert(c.fsiq),
    "SCQ (top vs bottom third)": tert(c.scq_res), "Vineland ABC (top vs bottom third)": tert(c.vabc_res), "RBS-R (top vs bottom third)": tert(c.rbsr_res),
    "CBCL total (top vs bottom third)": tert(c.cbcl_tot_res), "ADHD": b(c.behav_adhd), "mood or anxiety disorder": b(c.mood_or_anx), "seizures": b(c.seizures),
    "language regression": b(c.regression), "multiplex family": b(c.multiplex), "diagnosed by specialist": b(c.dx_specialist),
    "ASD validity flag": b(c.validity_flag), "area deprivation (top vs bottom third)": tert(c.adi_pct)}
fam = c.FID.astype(str).to_numpy(); fi, uf = pd.factorize(fam)
def diff(v, g, B=200, seed=1):
    ok = np.isfinite(v) & np.isfinite(g); v, g, f = v[ok], g[ok], fi[ok]
    if (g == 1).sum() < 30 or (g == 0).sum() < 30: return np.nan, np.nan, int((g == 1).sum()), int((g == 0).sum())
    d0 = v[g == 1].mean() - v[g == 0].mean(); rng = np.random.default_rng(seed); nu = fi.max() + 1; bs = []
    for _ in range(B):
        w = np.bincount(rng.integers(0, nu, nu), minlength=nu)[f]
        m1, m0 = w[g == 1].sum(), w[g == 0].sum()
        if m1 and m0: bs.append((w[g == 1] * v[g == 1]).sum() / m1 - (w[g == 0] * v[g == 0]).sum() / m0)
    return d0, float(np.std(bs)), int((g == 1).sum()), int((g == 0).sum())
rows = []
for nm, g in SPL.items():
    g = np.asarray(g, float)
    for meas, suf in [("population", ""), ("within-family pTDT", "_ptdt")]:
        for p in PGS:
            d0, se, n1, n0 = diff(c[p + suf].to_numpy(float), g); rows.append(dict(split=nm, measure=meas, pgs=p, diff=d0, se=se, z=d0 / se if se else np.nan, n_high=n1, n_low=n0))
    print(nm, "done", flush=True)
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/driver_profiles.tsv", sep="\t", index=False)
ref = D[(D.split.str.startswith("age at diagnosis")) & (D.measure == "population")].set_index("pgs")["diff"].reindex(PGS)
summ = []
for (nm, meas), x in D.groupby(["split", "measure"], sort=False):
    x = x.set_index("pgs").reindex(PGS); z = x.z.dropna()
    v = x["diff"].to_numpy(float); cos = np.nansum(v * ref.to_numpy()) / np.sqrt(np.nansum(v ** 2) * np.nansum(ref.to_numpy() ** 2)) if meas == "population" else np.nan
    summ.append(dict(split=nm, measure=meas, n_high=int(x.n_high.max()), n_low=int(x.n_low.max()), rms_diff=float(np.sqrt(np.nanmean(v ** 2))), chi2_profile=float((z ** 2).sum()),
                     max_abs_z=float(z.abs().max()) if len(z) else np.nan, cosine_vs_age_at_dx=cos))
S = pd.DataFrame(summ); S.to_csv(f"{OUT}/driver_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(S[S.measure == "population"].sort_values("chi2_profile", ascending=False).round(3).to_string(index=False))
print(S[S.measure != "population"].sort_values("chi2_profile", ascending=False).round(3).to_string(index=False))
print(D[D.measure == "population"].pivot(index="split", columns="pgs", values="z").reindex(columns=PGS).round(1).to_string())
print(D[D.measure != "population"].pivot(index="split", columns="pgs", values="z").reindex(columns=PGS).round(1).to_string())
# ---- joint model (population-style)
X = pd.DataFrame({nm: np.asarray(g, float) for nm, g in SPL.items()}, index=c.index)
Xc = pd.DataFrame(index=c.index)
for col in X.columns:
    v = X[col]; Xc[col] = v.fillna(v.median()) if v.notna().any() else 0.0
    if v.isna().any(): Xc[col + " (missing)"] = v.isna().astype(float)
Xc = (Xc - Xc.mean()) / Xc.std().replace(0, 1)
def ols_r2(y, A):
    A1 = np.column_stack([np.ones(len(y)), A]); bb = np.linalg.lstsq(A1, y, rcond=None)[0]; r = y - A1 @ bb; return 1 - r.var() / y.var(), bb
jrows = []
for p in PGS:
    y = c[p].to_numpy(float); ok = np.isfinite(y); yy = y[ok]; A = Xc.to_numpy()[ok]
    r2full, bfull = ols_r2(yy, A)
    for nm in SPL:
        keep = [i for i, cc in enumerate(Xc.columns) if not cc.startswith(nm)]
        r2red, _ = ols_r2(yy, A[:, keep]); jrows.append(dict(pgs=p, characteristic=nm, partial_r2=r2full - r2red, beta_joint=bfull[1 + list(Xc.columns).index(nm)]))
J = pd.DataFrame(jrows); J.to_csv(f"{OUT}/driver_joint.tsv", sep="\t", index=False)
print("\njoint model: partial R2 (x1000) by characteristic, summed over 9 PGS")
print((J.groupby("characteristic").partial_r2.sum() * 1000).sort_values(ascending=False).round(2).to_string())
print("\njoint betas (SD of PGS per SD of characteristic):"); print(J.pivot(index="characteristic", columns="pgs", values="beta_joint").reindex(columns=PGS).round(3).to_string())
