#!/usr/bin/env python3
"""Formal comparison: does a split stratify rare-variant burden more strongly than age at diagnosis?
Every candidate is defined the same way where possible (top vs bottom third over all SPARK autistic participants):
  age at diagnosis (thirds of diagnosis age; and the paper's > 10 vs < 6 y), maternal age, paternal age, household income
  (low vs high third), deprivation (ADI, high vs low third), mother's education (no degree vs degree; binary), ADHD (binary).
Outcome vector per split (4 classes, SD units of unaffected siblings): de novo NDD genes, de novo other genes, inherited NDD genes,
inherited other genes. Three outcome sets: (A) probands' absolute burden (between families), (B) unaffected siblings' absolute burden
(family-level signature), (C) within-family excess (proband minus sibling mean).
Statistics per split and outcome set: d = high minus low half (4-vector); norm ||d|| (overall stratification strength) and the Wald
chi-square d' S^-1 d (4 df; S = bootstrap covariance). Comparison with age at diagnosis: Delta = ||d_split|| - ||d_age||, paired family
bootstrap (the same resampled families for all splits, 1,000 replicates): 95% CI and two-sided bootstrap P.
Output (aggregate): ~/rvh2/agedx/rare_compare.tsv"""
import os
import numpy as np, pandas as pd
from scipy.stats import chi2
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
CATS = ["de novo, NDD genes", "de novo, other genes", "inherited, NDD genes", "inherited, other genes"]
a = load("pgs_rare_h2_drift_full", {k: CAT[k] for k in CATS})
sv = pd.read_csv(f"{P}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
reg = pd.read_csv(f"{P}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
dxm = pd.to_numeric(reg.diagnosis_age, errors="coerce").where(reg.asd.astype(str).str.upper().eq("TRUE"))
lo, hi = dxm.quantile(1 / 3), dxm.quantile(2 / 3); dx3 = pd.Series(np.where(dxm >= hi, 1.0, np.where(dxm <= lo, 0.0, np.nan)), index=dxm.index)
SPL = {"age at diagnosis (top vs bottom third)": dx3, "age at diagnosis (>10 vs <6 y; paper)": sv.s00, "maternal age (older vs younger third)": sv.s11,
       "paternal age (older vs younger third)": sv.s12, "household income (low vs high third)": sv.s09, "deprivation (high vs low third)": sv.s08,
       "mother no degree vs degree": sv.s10, "ADHD (yes vs no)": sv.s27,
       # developmental / autistic-trait splits (survey codes; continuous measures top vs bottom third)
       "late walking (top vs bottom third)": sv.s16, "late first words (top vs bottom third)": sv.s17, "motor coordination problems (DCDQ)": sv.s20,
       "low Vineland ABC (bottom vs top third)": sv.s26, "SCQ (top vs bottom third)": sv.s23, "RBS-R (top vs bottom third)": sv.s24,
       "SRS (top vs bottom third)": sv.s25, "intellectual disability": sv.s13, "minimal language": sv.s15, "language regression": sv.s18,
       "early onset concern": sv.s19, "motor delay": sv.s21}
fam = (a.cohort + "|" + a.family_key); fi = pd.factorize(fam)[0]; nu = fi.max() + 1
case = a.case.eq(1).to_numpy(); Y = a[CATS].to_numpy(float)
sibm = a[a.case.eq(0)].groupby(["cohort", "family_key"])[CATS].mean()
EX = Y - sibm.reindex(pd.MultiIndex.from_frame(a[["cohort", "family_key"]])).to_numpy(float)   # within-family excess (meaningful for cases)
G = {}
for nm, s in SPL.items():   # family-level split value from the family's autistic child
    cs = a[case].assign(g=pd.to_numeric(a.loc[case, "sample"].map(s), errors="coerce").to_numpy())
    G[nm] = fam.map(cs.dropna(subset=["g"]).groupby((cs.cohort + "|" + cs.family_key)[cs.g.notna()]).g.first()).to_numpy(float)
SETS = {"A probands, absolute": (case, Y), "B unaffected siblings, absolute": (~case, Y), "C within-family excess (proband - sibling)": (case, EX)}
def dvec(w, m, V, g):
    ok = m & np.isfinite(g) & np.isfinite(V).all(1); h, l = ok & (g == 1), ok & (g == 0)
    return (V[h] * w[h, None]).sum(0) / w[h].sum() - (V[l] * w[l, None]).sum(0) / w[l].sum()
rng = np.random.default_rng(2026); B = 1000; ones = np.ones(len(a))
W = [np.bincount(rng.integers(0, nu, nu), minlength=nu)[fi].astype(float) for _ in range(B)]
rows = []
for sname, (m, V) in SETS.items():
    est = {nm: dvec(ones, m, V, g) for nm, g in G.items()}
    boot = {nm: np.array([dvec(w, m, V, g) for w in W]) for nm, g in G.items()}
    ref = "age at diagnosis (top vs bottom third)"
    nref = np.linalg.norm(boot[ref], axis=1)
    for nm in G:
        d = est[nm]; S = np.cov(boot[nm].T); wald = float(d @ np.linalg.solve(S, d)); nb = np.linalg.norm(boot[nm], axis=1)
        wd = float(d[:2] @ np.linalg.solve(S[:2, :2], d[:2])); wi = float(d[2:] @ np.linalg.solve(S[2:, 2:], d[2:]))   # de novo classes only / inherited only (2 df each)
        delta = np.linalg.norm(d) - np.linalg.norm(est[ref]); db = nb - nref
        p = 2 * min((db <= 0).mean(), (db >= 0).mean()) if nm != ref else np.nan
        h_ = m & (G[nm] == 1); l_ = m & (G[nm] == 0)
        rows.append(dict(outcome=sname, split=nm, n_high=int(h_.sum()), n_low=int(l_.sum()), wald_denovo_2df=wd, p_denovo=float(chi2.sf(wd, 2)), wald_inherited_2df=wi, p_inherited=float(chi2.sf(wi, 2)), **{f"d_{c}": v for c, v in zip(CATS, d)},
                         norm=float(np.linalg.norm(d)), wald_chi2_4df=wald, p_wald=float(chi2.sf(wald, 4)),
                         delta_vs_age=delta, delta_lo95=float(np.percentile(db, 2.5)), delta_hi95=float(np.percentile(db, 97.5)), p_delta=p))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/rare_compare.tsv", sep="\t", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 30)
print(T[["outcome", "split", "n_high", "n_low", "norm", "wald_chi2_4df", "p_wald", "delta_vs_age", "delta_lo95", "delta_hi95", "p_delta"]].round(4).to_string(index=False))
