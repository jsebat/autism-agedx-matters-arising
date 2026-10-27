#!/usr/bin/env python3
"""Split vectors for Fig. 2a. Every split is drawn as an arrow from its low/no half to its high/yes half in a two-dimensional
external-genetics space: x = psychiatric (mean difference in rg with ADHD and MDD), y = reproductive timing / education (mean difference
in rg with age at first birth and educational attainment); Zhang F2 minus F1 for comparison. Arrows oriented so x >= 0.
Inputs: results/full_data/all_rg_ext_all.tsv (results/replication/rg_ext_all.tsv with DATA=arrays), survey_splits.tsv.
Output: results/axis_vectors.tsv (read by fig2_main.py)."""
import os
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results")
E = pd.read_csv(os.path.join(R, "full_data", "all_rg_ext_all.tsv") if os.environ.get("DATA", "all") == "all" else os.path.join(R, "replication", "rg_ext_all.tsv"), sep="\t")   # all SPARK strata by default
rg = E.pivot_table(index="gwas", columns="trait", values="rg"); se = E.pivot_table(index="gwas", columns="trait", values="se")
sv = pd.read_csv(os.path.join(R, "survey_splits.tsv"), sep="\t").set_index("code")
def idx(g):
    psy = (rg.loc[g, "adhd_2022"] + rg.loc[g, "mdd_2025"]) / 2; rep = (rg.loc[g, "afb_2018"] + rg.loc[g, "edu_2022"]) / 2
    return psy, rep
rows = []
for code in sorted({g.split("_")[1] for g in rg.index if g.startswith("rep_s")}):
    hi, lo = f"rep_{code}_hi_c2", f"rep_{code}_lo_c1"
    need = ["adhd_2022", "mdd_2025", "afb_2018", "edu_2022"]
    if any(g not in rg.index or g not in se.index or rg.loc[g, need].isna().any() for g in (hi, lo)): continue
    (p1, r1), (p0, r0) = idx(hi), idx(lo); dp, dr = p1 - p0, r1 - r0
    s = np.sqrt(se.loc[hi, ["adhd_2022", "mdd_2025"]].pow(2).sum() + se.loc[lo, ["adhd_2022", "mdd_2025"]].pow(2).sum()) / 2
    sr = np.sqrt(se.loc[hi, ["afb_2018", "edu_2022"]].pow(2).sum() + se.loc[lo, ["afb_2018", "edu_2022"]].pow(2).sum()) / 2
    rows.append(dict(code=code, split=sv.split.get(code, "design: no parent vs both parents genotyped"), d_psych=dp, d_repro=dr, se_psych=s, se_repro=sr))
V = pd.DataFrame(rows)
(p1, r1), (p0, r0) = idx("W_F2"), idx("W_F1")
theirs = dict(code="F2-F1", split="Zhang F2 minus Zhang F1", d_psych=p1 - p0, d_repro=r1 - r0, se_psych=np.nan, se_repro=np.nan)
V = pd.concat([V, pd.DataFrame([theirs])], ignore_index=True)
flip = V.d_psych < 0; V.loc[flip, ["d_psych", "d_repro"]] *= -1; V["flipped"] = flip
V["angle_from_psych_axis_deg"] = np.degrees(np.arctan2(V.d_repro.abs(), V.d_psych))
V.to_csv(os.path.join(R, "axis_vectors.tsv"), sep="\t", index=False)
print(V[V.code.isin(["F2-F1", "s00", "s09", "s10", "s11", "s12", "s08", "s27", "s28", "s29", "s13", "s01"])][["code", "split", "d_psych", "d_repro", "se_psych", "se_repro", "angle_from_psych_axis_deg"]].round(2).to_string(index=False))
