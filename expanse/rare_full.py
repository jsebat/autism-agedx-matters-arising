#!/usr/bin/env python3
"""Rare-variant comparison on the EXPANDED trio set (Matters Arising, 2026-10-25), main sample = Zhang et al.'s criterion for their
age-at-diagnosis analyses (registered before age 22; env U22=0 for all ages).
Probands: legacy SPARK trio probands from the lab's expanded collection (~/rvh2/out_trio: WES IDT, WES Twist, iWGS v1.1; sibling not required)
  + sibling-family trio probands of SSC and SPARK-WGS-2026_08 (~/rvh2/out_dnm_v3; new-WGS v3 strand filter). One row per person
  (priority iWGS v1.1 > Twist > IDT > new WGS). SNV/indel only; de novo = SNV/indel de novo calls (no CNV).
Classes (as before): De novo - NDD genes, De novo - other genes, Inherited - NDD genes, Inherited - other genes; weighted sums of the 24 cells
  with the within-family log-OR weights (pgs_rare_h2_drift_full/results/collider_weights.tsv), standardised to unaffected siblings of the
  same cohort (out_dnm_v3 siblings).
Splits: SPARK survey codes (spark_survey.tsv) + age at diagnosis thirds; SSC has no survey codes, so SSC probands enter
  only the CCA. Per split: 4-vector d (high minus low), Wald chi2 (4 df; 2 df de novo; 2 df inherited) from a 1,000-replicate family
  bootstrap, and the paired bootstrap difference in norm vs age at diagnosis.
Outputs: ~/rvh2/agedx/rare_full_cases.parquet (participant-level, Expanse only), rare_compare_full{_all}.tsv, rare_full_counts.txt (aggregate)"""
import os, io
import numpy as np, pandas as pd
from scipy.stats import chi2
H = os.path.expanduser("~/rvh2"); OUT = H + "/agedx"; P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
J = "/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08"
U22 = os.environ.get("U22", "1") == "1"; TAG = "" if U22 else "_all"
buf = io.StringIO(); pr = lambda *x: print(*x, file=buf, flush=True)
TIERS = ["lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4"]
cell = lambda modes, sets: [f"{t}__{m}__{g}" for t in TIERS for m in modes for g in sets]
CAT = {"de novo, NDD genes": cell(["dnm"], ["ndd"]), "de novo, other genes": cell(["dnm"], ["nonndd"]),
       "inherited, NDD genes": cell(["inh"], ["ndd"]), "inherited, other genes": cell(["inh"], ["nonndd"]),
       "inherited LoF tier 1": [f"lof_t1__inh__{g}" for g in ("ndd", "nonndd")]}
CATS = list(CAT)[:4]
W = pd.read_csv("/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift_full/results/collider_weights.tsv", sep="\t", index_col=0).iloc[:, 0]
ALLC = sorted({c for v in CAT.values() for c in v})
# ---- probands and siblings
parts = []
for c in ["spark_wgs", "spark_wes_twist", "spark_wes_idt"]:
    x = pd.read_parquet(f"{H}/out_trio/{c}.parquet"); x["family_key"] = x.family_id.astype(str); parts.append(x[["sample", "cohort", "case", "family_key"] + ALLC])
for c in ["ssc", "spark_wgs_new"]:
    x = pd.read_parquet(f"{H}/out_dnm_v3/{c}.parquet"); parts.append(x[x.case.eq(1)][["sample", "cohort", "case", "family_key"] + ALLC])
pro = pd.concat(parts, ignore_index=True)
meta = pd.read_csv(f"{J}/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv", sep="\t", dtype=str).set_index("sample_id")
pro["sp"] = np.where(pro.cohort.eq("spark_wgs_new"), pro["sample"].map(meta.ind_id), pro["sample"])
pri = {"spark_wgs": 0, "spark_wes_twist": 1, "spark_wes_idt": 2, "spark_wgs_new": 3, "ssc": 0}
pro["pri"] = pro.cohort.map(pri); n0 = len(pro); pro = pro.sort_values("pri").drop_duplicates("sp").drop(columns="pri")
pr("probands before/after de-duplication:", n0, len(pro))
sibs = pd.concat([pd.read_parquet(f"{H}/out_dnm_v3/{c}.parquet").query("case == 0").assign(cohort=c) for c in ["spark_wgs", "spark_wes_twist", "spark_wes_idt", "ssc", "spark_wgs_new"]], ignore_index=True)
for k, ts in CAT.items():
    ts = [t for t in ts if t in W.index]
    pro[k] = (pro[ts].fillna(0).to_numpy(float) * W[ts].to_numpy()).sum(1); sibs[k] = (sibs[ts].fillna(0).to_numpy(float) * W[ts].to_numpy()).sum(1)
    ref = sibs.groupby("cohort")[k].agg(["mean", "std"]); pro[k] = (pro[k] - pro.cohort.map(ref["mean"])) / pro.cohort.map(ref["std"])
# ---- registration age, restriction
reg = pd.read_csv(f"{P}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
pro["reg_age"] = pd.to_numeric(pro.sp.map(reg.age_at_registration_years), errors="coerce")
pro["u22"] = pro.reg_age.lt(22) | pro.cohort.eq("ssc")
pr("probands by cohort (all / registered before 22):"); pr(pro.groupby("cohort").agg(all=("sp", "size"), u22=("u22", "sum")).to_string())
if not os.environ.get("ADJ"): pro.to_parquet(f"{OUT}/rare_full_cases.parquet", index=False)   # sensitivity runs never rewrite the shared table
a = pro[pro.u22].copy() if U22 else pro.copy()
a = a[~a.cohort.eq("ssc")].reset_index(drop=True)   # splits need SPARK survey codes
# ---- splits, age-at-diagnosis thirds defined among SPARK autistic registrants of the same age criterion
sv = pd.read_csv(f"{P}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
asd = reg.asd.astype(str).str.upper().eq("TRUE") & (pd.to_numeric(reg.age_at_registration_years, errors="coerce").lt(22) if U22 else True)
dxm = pd.to_numeric(reg.diagnosis_age, errors="coerce").where(asd)
lo, hi = dxm.quantile(1 / 3), dxm.quantile(2 / 3); dx3 = pd.Series(np.where(dxm >= hi, 1.0, np.where(dxm <= lo, 0.0, np.nan)), index=dxm.index)
SPL = {"age at diagnosis (top vs bottom third)": dx3, "age at diagnosis (>10 vs <6 y; paper)": sv.s00, "maternal age (older vs younger third)": sv.s11,
       "paternal age (older vs younger third)": sv.s12, "household income (low vs high third)": sv.s09, "deprivation (high vs low third)": sv.s08,
       "mother no degree vs degree": sv.s10, "ADHD (yes vs no)": sv.s27, "late walking (top vs bottom third)": sv.s16, "late first words (top vs bottom third)": sv.s17,
       "motor coordination problems (DCDQ)": sv.s20, "low Vineland ABC (bottom vs top third)": sv.s26, "SCQ (top vs bottom third)": sv.s23, "RBS-R (top vs bottom third)": sv.s24,
       "SRS (top vs bottom third)": sv.s25, "intellectual disability": sv.s13, "minimal language": sv.s15, "language regression": sv.s18,
       "early onset concern": sv.s19, "motor delay": sv.s21}
fi = pd.factorize(a.cohort + "|" + a.family_key)[0]; nu = fi.max() + 1; Y = a[CATS].to_numpy(float).copy()
G = {nm: pd.to_numeric(a.sp.map(s), errors="coerce").to_numpy(float) for nm, s in SPL.items()}
ADJ = os.environ.get("ADJ", ""); ADJ = "" if ADJ == "pooled" else ADJ   # "pooled" = original analysis without rewriting the shared table; review fix: "" pooled (original) | "cohort" within-cohort contrasts pooled | "cohort_parage" + burdens residualised on parental ages within cohort
coh = a.cohort.astype(str).to_numpy(); COH = np.unique(coh)
if ADJ == "cohort_parage":
    env_ = pd.read_csv(f"{P}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
    pa = np.column_stack([pd.to_numeric(a.sp.map(env_[c]), errors="coerce").to_numpy(float) for c in ("maternal_age", "paternal_age")])
    keep_pa = np.isfinite(pa).all(1)
    for c_ in COH:
        m = (coh == c_) & keep_pa & np.isfinite(Y).all(1); Zc = np.column_stack([np.ones(m.sum()), pa[m]])
        Y[m] = Y[m] - Zc @ np.linalg.lstsq(Zc, Y[m], rcond=None)[0] + Y[m].mean(0)
    Y[~keep_pa] = np.nan
    TAG += "_parage"
if ADJ: TAG += "_strat"
def dvec(w, g):
    ok = np.isfinite(g) & np.isfinite(Y).all(1)
    if not ADJ:
        h, l = ok & (g == 1), ok & (g == 0)
        return (Y[h] * w[h, None]).sum(0) / w[h].sum() - (Y[l] * w[l, None]).sum(0) / w[l].sum()
    num, den = 0.0, 0.0   # within-cohort high-minus-low, weighted by the harmonic weight of the two (bootstrap-weighted) group sizes
    for c_ in COH:
        h, l = ok & (g == 1) & (coh == c_), ok & (g == 0) & (coh == c_); nh, nl = w[h].sum(), w[l].sum()
        if nh < 5 or nl < 5: continue
        d_ = (Y[h] * w[h, None]).sum(0) / nh - (Y[l] * w[l, None]).sum(0) / nl; wc = nh * nl / (nh + nl); num = num + wc * d_; den += wc
    return num / den
rng = np.random.default_rng(2026); B = 1000; ones = np.ones(len(a))
WB = [np.bincount(rng.integers(0, nu, nu), minlength=nu)[fi].astype(float) for _ in range(B)]
est = {nm: dvec(ones, g) for nm, g in G.items()}; boot = {nm: np.array([dvec(w, g) for w in WB]) for nm, g in G.items()}
ref = "age at diagnosis (top vs bottom third)"; nref = np.linalg.norm(boot[ref], axis=1); rows = []
for nm, g in G.items():
    d = est[nm]; S = np.cov(boot[nm].T); wald = float(d @ np.linalg.solve(S, d)); nb = np.linalg.norm(boot[nm], axis=1)
    wd = float(d[:2] @ np.linalg.solve(S[:2, :2], d[:2])); wi = float(d[2:] @ np.linalg.solve(S[2:, 2:], d[2:])); db = nb - nref
    ok = np.isfinite(g); rows.append(dict(outcome="A probands, absolute", split=nm, n_high=int((ok & (g == 1)).sum()), n_low=int((ok & (g == 0)).sum()),
        **{f"d_{c}": v for c, v in zip(CATS, d)}, norm=float(np.linalg.norm(d)), wald_chi2_4df=wald, p_wald=float(chi2.sf(wald, 4)),
        wald_denovo_2df=wd, p_denovo=float(chi2.sf(wd, 2)), wald_inherited_2df=wi, p_inherited=float(chi2.sf(wi, 2)),
        delta_vs_age=float(np.linalg.norm(d) - np.linalg.norm(est[ref])), delta_lo95=float(np.percentile(db, 2.5)), delta_hi95=float(np.percentile(db, 97.5)),
        p_delta=float(max(2 * min((db <= 0).mean(), (db >= 0).mean()), 1 / B)) if nm != ref else np.nan, p_delta_is_floor=bool(nm != ref and min((db <= 0).mean(), (db >= 0).mean()) == 0),
        **{f"d_lo95_{c}": float(np.percentile(boot[nm][:, i], 2.5)) for i, c in enumerate(CATS)}, **{f"d_hi95_{c}": float(np.percentile(boot[nm][:, i], 97.5)) for i, c in enumerate(CATS)}))
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/rare_compare_full{TAG}.tsv", sep="\t", index=False)
pr(f"\nsplit comparison ({'registered before 22' if U22 else 'all ages'}), SPARK trio probands n = {len(a)}")
pd.set_option("display.width", 250); pr(T[["split", "n_high", "n_low", "wald_chi2_4df", "p_wald", "wald_inherited_2df", "wald_denovo_2df", "delta_vs_age", "p_delta"]].round(4).sort_values("wald_chi2_4df", ascending=False).to_string(index=False))
open(f"{OUT}/rare_full_counts{TAG}.txt", "w").write(buf.getvalue()); print(buf.getvalue())
