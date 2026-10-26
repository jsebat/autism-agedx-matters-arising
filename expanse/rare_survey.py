#!/usr/bin/env python3
"""Do the survey splits (Matters Arising on Zhang, Grove et al. 2025) stratify rare-variant burden? Their own trio analysis found no
association of age at diagnosis with de novo or inherited rare variants (6,206 SPARK trios, their Supplementary Table 11).
Families: SPARK case / unaffected-sibling discordant complete trios with the 24-cell rare SNV/indel model (pgs_rare_h2_drift_full:
6 tiers x de novo / inherited x NDD / other genes, filter v3) and, separately, rare CNVs (pgs_rare_h2_drift_fullcnv; GSA v1 + exome).
Category scores = sum of pooled within-family log-OR weights x counts (results/collider_weights.tsv) over the category's cells,
standardized to unaffected siblings within cohort:
  de novo NDD genes | de novo other genes | de novo LoF (tiers 1-2) | de novo missense (tiers 1-4)
  inherited NDD genes | inherited other genes | inherited LoF tier 1 | inherited missense (tiers 1-4) | all SNV/indel | CNV (all classes)
Within-family excess per autistic child = child score - mean score of the unaffected siblings. Per split (spark_survey.tsv codes):
excess in high/yes half - excess in low/no half, family-bootstrap SE (500), z.
Output (aggregate): ~/rvh2/agedx/rare_survey.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
sv = pd.read_csv(f"{P}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
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
def excess(a, cols):
    sib = a[a.case.eq(0)].groupby(["cohort", "family_key"])[cols].mean()
    c = a[a.case.eq(1)].copy(); m = sib.reindex(pd.MultiIndex.from_frame(c[["cohort", "family_key"]])); ok = m.notna().all(axis=1).to_numpy()
    c = c[ok].copy(); m = m[ok]
    for col in cols: c[col] = c[col].to_numpy() - m[col].to_numpy()
    return c
def run(c, cols, tag, B=500):
    fi = pd.factorize(c.cohort + "|" + c.family_key)[0]; nu = fi.max() + 1; rng = np.random.default_rng(7)
    BW = [np.bincount(rng.integers(0, nu, nu), minlength=nu)[fi] for _ in range(B)]; rows = []
    for code in sv.columns:
        g = pd.to_numeric(c["sample"].map(sv[code]), errors="coerce").to_numpy(float)
        for col in cols:
            v = c[col].to_numpy(float); ok = np.isfinite(v) & np.isfinite(g); h, l = ok & (g == 1), ok & (g == 0)
            if h.sum() < 30 or l.sum() < 30: continue
            d0 = v[h].mean() - v[l].mean(); bs = [np.average(v[h], weights=w[h]) - np.average(v[l], weights=w[l]) for w in BW if w[h].sum() and w[l].sum()]
            se = float(np.std(bs)); rows.append(dict(set=tag, code=code, category=col, n_high=int(h.sum()), n_low=int(l.sum()), excess_high=v[h].mean(), excess_low=v[l].mean(), diff=d0, se=se, z=d0 / se))
    return pd.DataFrame(rows)
a = load("pgs_rare_h2_drift_full", CAT); c = excess(a, list(CAT)); print("SNV/indel trio probands:", len(c), flush=True)
D1 = run(c, list(CAT), "SNV/indel")
a2 = load("pgs_rare_h2_drift_fullcnv", {"CNV (all classes)": ["cnv_del__pgc18__GSA", "cnv_dup__pgc18__GSA", "cnv_del__ndd__GSA", "cnv_del__pgc18__XHMM", "cnv_dup__pgc18__XHMM", "cnv_del__ndd__XHMM"]})
c2 = excess(a2, ["CNV (all classes)"]); print("CNV probands:", len(c2), flush=True)
D = pd.concat([D1, run(c2, ["CNV (all classes)"], "CNV")]); D.to_csv(f"{OUT}/rare_survey.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
Z = D.pivot_table(index="code", columns="category", values="z").round(1)
print(Z.to_string())
