#!/usr/bin/env python3
"""Between-family rare-variant signature of the stratifying axes (Antaki et al. Nat Genet 2022, PMID 35654974: parental-age effects on
autism risk reflect de novo mutations in the paternal germline AND inherited risk that influences parental behaviour).
Unlike rare_survey.py (proband minus sibling), this compares ABSOLUTE burden between families: probands and, separately, their
unaffected siblings, in the high vs low half of each family-level split (split value of the family's autistic child). If unaffected
siblings show the same gradient, the axis is a family-level genetic signature rather than autism-specific transmission.
Scores: same rare-variant categories and weights as rare_survey.py, standardized to all unaffected siblings within cohort.
Also continuous regressions on maternal and paternal age at the child's birth (per 10 years), family-clustered SE.
Output (aggregate): ~/rvh2/agedx/rare_between.tsv, rare_between_parentage.tsv"""
import os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
src = open(os.path.expanduser("~/rvh2/rare_survey.py")).read()
exec(src[src.index("TI = "):src.index("def excess")])                  # CAT and load()
sv = pd.read_csv(f"{P}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
env = pd.read_csv(f"{P}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
a = load("pgs_rare_h2_drift_full", CAT)
a2 = load("pgs_rare_h2_drift_fullcnv", {"CNV (all classes)": ["cnv_del__pgc18__GSA", "cnv_dup__pgc18__GSA", "cnv_del__ndd__GSA", "cnv_del__pgc18__XHMM", "cnv_dup__pgc18__XHMM", "cnv_del__ndd__XHMM"]})
SPL = ["s00", "s09", "s10", "s11", "s12", "s08", "s27", "s13", "s04"]
rows, prow = [], []
for df, cols, tag in [(a, list(CAT), "SNV/indel"), (a2, ["CNV (all classes)"], "CNV")]:
    fam = df.cohort + "|" + df.family_key
    for code in SPL:   # family-level split value: that of the family's autistic child(ren) (first non-missing)
        cs = df[df.case.eq(1)].assign(g=lambda d: pd.to_numeric(d["sample"].map(sv[code]), errors="coerce"))
        fg = cs.dropna(subset=["g"]).groupby(cs.cohort + "|" + cs.family_key).g.first()
        g = fam.map(fg).to_numpy(float)
        fi = pd.factorize(fam)[0]; nu = fi.max() + 1; rng = np.random.default_rng(11); BW = [np.bincount(rng.integers(0, nu, nu), minlength=nu)[fi] for _ in range(300)]
        for who, m in [("probands", df.case.eq(1).to_numpy()), ("unaffected siblings", df.case.eq(0).to_numpy())]:
            for col in cols:
                v = df[col].to_numpy(float); ok = m & np.isfinite(v) & np.isfinite(g); h, l = ok & (g == 1), ok & (g == 0)
                if h.sum() < 30 or l.sum() < 30: continue
                d0 = v[h].mean() - v[l].mean(); bs = [np.average(v[h], weights=w[h]) - np.average(v[l], weights=w[l]) for w in BW if w[h].sum() and w[l].sum()]
                se = float(np.std(bs)); rows.append(dict(set=tag, code=code, who=who, category=col, n_high=int(h.sum()), n_low=int(l.sum()), diff=d0, se=se, z=d0 / se))
    for par in ["maternal_age", "paternal_age"]:
        df = df.assign(x=pd.to_numeric(df["sample"].map(env[par]), errors="coerce") / 10, fam=fam)
        for who, sel in [("probands", df.case.eq(1)), ("unaffected siblings", df.case.eq(0))]:
            for col in cols:
                x = df[sel].dropna(subset=["x", col])
                if len(x) < 200: continue
                f = smf.ols(f"Q('{col}') ~ x", data=x).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(x.fam)[0]})
                prow.append(dict(set=tag, parent=par, who=who, category=col, n=len(x), beta_per_10y=f.params["x"], se=f.bse["x"], z=f.tvalues["x"], p=f.pvalues["x"]))
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/rare_between.tsv", sep="\t", index=False)
Q = pd.DataFrame(prow); Q.to_csv(f"{OUT}/rare_between_parentage.tsv", sep="\t", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_rows", 500)
for who in ["probands", "unaffected siblings"]:
    print(f"\n== between-family z, {who} (high minus low half)"); print(D[D.who == who].pivot_table(index="code", columns="category", values="z").round(1).to_string())
print("\n== parental age (per 10 y)"); print(Q.pivot_table(index=["parent", "who"], columns="category", values="z").round(1).to_string())
print(Q[Q.category.isin(["de novo, NDD genes", "de novo, other genes", "inherited, NDD genes", "inherited, other genes", "inherited LoF tier 1"])].round(3).to_string(index=False))
