#!/usr/bin/env python3
"""Rare-variant analogue of driver_pgs.py: which characteristic of SPARK autistic cases separates rare-variant burden, and where
does age at diagnosis (after 10 vs before 6; Zhang, Grove et al. 2025) rank?
Families: case/unaffected-sibling discordant SPARK complete trios with the 24-term rare SNV/indel model (pgs_rare_h2_drift_full;
6 tiers x de novo/inherited x NDD/other genes, filter v3). Score = sum of pooled within-family log-OR weights x counts
(results/collider_weights.tsv), standardized to unaffected siblings within cohort; also de novo and inherited parts.
Within-family excess per case = case score - mean score of unaffected siblings in the same family (immune to ancestry, platform,
family structure). Same-family PGS sibling differences (ASD, ADHD, EA, IQ, SCZ) for comparison.
Per split: excess in high group - excess in low group, family-bootstrap SE (200), z; rank of age at diagnosis among splits.
CNV (secondary; pgs_rare_h2_drift_fullcnv, GSA v1 + exome XHMM, GSA v2 missing): CNV part of that root's weighted score.
Output (aggregate): ~/rvh2/agedx/rare_driver_*.tsv"""
import os
import numpy as np, pandas as pd
R = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift"; OUT = os.path.expanduser("~/rvh2/agedx")
P = f"{R}/pheno"
sev = pd.read_csv(f"{P}/spark_severity.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
ho = pd.read_csv(f"{P}/spark_heldout.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
env = pd.read_csv(f"{P}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
dx = pd.read_csv(f"{P}/spark_dx_derived.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
print("dx columns:", list(dx.columns), flush=True)


def b(x): return pd.to_numeric(x, errors="coerce")


def tert(x):
    x = b(x); lo, hi = x.quantile(1 / 3), x.quantile(2 / 3); return np.where(x >= hi, 1.0, np.where(x <= lo, 0.0, np.nan))


def load(root):
    a = pd.read_parquet(f"/expanse/projects/sebat1/jsebat/{root}/assembled_children.parquet")
    a = a[a.cohort.str.startswith("spark")].copy()
    W = pd.read_csv(f"/expanse/projects/sebat1/jsebat/{root}/results/collider_weights.tsv", sep="\t", index_col=0).iloc[:, 0]
    W = W[[t for t in W.index if t in a.columns]]
    parts = {"rare SNV/indel (all)": [t for t in W.index if not t.startswith("cnv")], "de novo": [t for t in W.index if "__dnm__" in t],
             "inherited": [t for t in W.index if "__inh__" in t], "CNV": [t for t in W.index if t.startswith("cnv")]}
    for k, ts in parts.items():
        if not ts: continue
        a[k] = (a[ts].fillna(0).to_numpy(float) * W[ts].to_numpy()).sum(1)
        ref = a[a.case.eq(0)].groupby("cohort")[k].agg(["mean", "std"]); a[k] = (a[k] - a.cohort.map(ref["mean"])) / a.cohort.map(ref["std"])
    print(root, "terms", len(W), {k: len(v) for k, v in parts.items()}, flush=True)
    return a, [k for k, v in parts.items() if v]


def excess(a, cols):
    sib = a[a.case.eq(0)].groupby(["cohort", "family_key"])[cols].mean()
    c = a[a.case.eq(1)].copy(); k = pd.MultiIndex.from_frame(c[["cohort", "family_key"]])
    m = sib.reindex(k); ok = m.notna().all(axis=1).to_numpy(); c = c[ok].copy(); m = m[ok]
    for col in cols: c[col] = c[col].to_numpy() - m[col].to_numpy()
    return c


def splits(c):
    s = c["sample"]
    da = b(s.map(dx.dx_age_y)); dy = b(s.map(dx.dx_year))
    regcol = "age_at_registration" if "age_at_registration" in dx.columns else "reg_age_y"
    reg = b(s.map(dx[regcol]))
    return {"age at diagnosis (after 10 vs before 6)": np.where(da >= 11, 1.0, np.where(da < 6, 0.0, np.nan)),
            "adult at registration": np.where(reg.isna(), np.nan, (reg >= 18).astype(float)),
            "female": np.where(c.male.isna(), np.nan, 1 - c.male.astype(float)),
            "diagnosed 2019-26 vs before 2010": np.where(dy >= 2019, 1.0, np.where(dy < 2010, 0.0, np.nan)),
            "WGS vs exome": np.where(c.cohort.str.contains("wgs"), 1.0, 0.0),
            "intellectual disability": b(s.map(dx.ID)), "minimal language": b(s.map(dx.minimal_language)), "FSIQ (top vs bottom third)": tert(s.map(sev.fsiq)),
            "SCQ (top vs bottom third)": tert(s.map(sev.scq_res)), "Vineland ABC (top vs bottom third)": tert(s.map(sev.vabc_res)), "RBS-R (top vs bottom third)": tert(s.map(sev.rbsr_res)),
            "CBCL total (top vs bottom third)": tert(s.map(sev.cbcl_tot_res)), "ADHD": b(s.map(ho.behav_adhd)), "mood or anxiety disorder": b(s.map(ho.mood_or_anx)), "seizures": b(s.map(ho.seizures)),
            "language regression": b(s.map(env.regression)), "multiplex family": b(s.map(ho.multiplex)), "diagnosed by specialist": b(s.map(ho.dx_specialist)),
            "ASD validity flag": b(s.map(dx.validity_flag)), "area deprivation (top vs bottom third)": tert(s.map(env.adi_pct))}


def run(c, cols, tag, B=200):
    fi = pd.factorize(c.cohort + "|" + c.family_key)[0]; nu = fi.max() + 1; rng = np.random.default_rng(1)
    BW = [np.bincount(rng.integers(0, nu, nu), minlength=nu)[fi] for _ in range(B)]
    rows = []
    for nm, g in splits(c).items():
        g = np.asarray(g, float)
        for col in cols:
            v = c[col].to_numpy(float); ok = np.isfinite(v) & np.isfinite(g); h, l = ok & (g == 1), ok & (g == 0)
            if h.sum() < 30 or l.sum() < 30: continue
            d0 = v[h].mean() - v[l].mean(); bs = []
            for w in BW:
                if w[h].sum() and w[l].sum(): bs.append(np.average(v[h], weights=w[h]) - np.average(v[l], weights=w[l]))
            se = float(np.std(bs))
            rows.append(dict(set=tag, split=nm, score=col, n_high=int(h.sum()), n_low=int(l.sum()), excess_high=v[h].mean(), excess_low=v[l].mean(), diff=d0, se=se, z=d0 / se))
    print(tag, "splits done", flush=True)
    return pd.DataFrame(rows)


a, cols = load("pgs_rare_h2_drift_full")
pg = ["asd_loSSC_z", "adhd_2022_z", "edu_2022_z", "iq_2018_z", "scz_2022_z"]
for p in pg:  # standardize PGS to unaffected siblings within cohort
    ref = a[a.case.eq(0)].groupby("cohort")[p].agg(["mean", "std"]); a[p] = (a[p] - a.cohort.map(ref["mean"])) / a.cohort.map(ref["std"])
c = excess(a, cols + pg); print("SNV families: cases", len(c), flush=True)
D = run(c, cols + pg, "SNV/indel trios")
a2, cols2 = load("pgs_rare_h2_drift_fullcnv"); c2 = excess(a2, ["CNV"]); print("CNV families: cases", len(c2), flush=True)
D = pd.concat([D, run(c2, ["CNV"], "CNV (GSA v1 + exome XHMM)")]); D.to_csv(f"{OUT}/rare_driver_splits.tsv", sep="\t", index=False)
ov = pd.DataFrame([dict(set="SNV/indel trios", score=k, mean_excess=c[k].mean(), se=c[k].std() / np.sqrt(len(c)), n=len(c)) for k in cols + pg] +
                  [dict(set="CNV", score="CNV", mean_excess=c2.CNV.mean(), se=c2.CNV.std() / np.sqrt(len(c2)), n=len(c2))])
ov.to_csv(f"{OUT}/rare_driver_overall.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(ov.round(4).to_string(index=False))
for sc in cols + pg + ["CNV"]:
    x = D[D.score == sc].copy(); x["chi2"] = x.z ** 2; x = x.sort_values("chi2", ascending=False); x["rank"] = range(1, len(x) + 1)
    print("\n==", sc); print(x[["split", "n_high", "n_low", "excess_high", "excess_low", "diff", "se", "z", "rank"]].round(3).to_string(index=False))
