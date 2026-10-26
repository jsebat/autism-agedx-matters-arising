#!/usr/bin/env python3
"""The missing EDU/SES axis (not captured by Zhang F1/F2). (a) Every split drawn as an arrow from its low/no half to its high/yes half in a two-dimensional
external-genetics space: x = psychiatric (mean difference in rg with ADHD and MDD), y = reproductive timing / education (mean difference
in rg with age at first birth and educational attainment); their factor contrast (F2 minus F1) for comparison. Arrows oriented so x >= 0.
(b) Rare-variant stratification of the same splits (paired family bootstrap; rare_compare.tsv): Wald chi-square (4 df, probands'
absolute burden) and siblings' result. Inputs: results/replication/rg_ext_all.tsv, survey_splits.tsv, rare_compare.tsv.
Output: Manuscripts/Nature_LTE/figures/L_AxisMissing.(png|pdf), results/axis_vectors.tsv"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results"); MS = os.path.join(H, "..", "..", "..", "Manuscripts", "Nature_LTE")
BLUE, ORANGE, PURPLE, GREY, DGREY, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#6a3d9a", "#c4c4c4", "#7a7a7a", "#222222", "#6b6b6b", "#dddddd"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False})
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

LAB = {"F2-F1": ("Zhang F2 - F1\n(later - earlier)", PURPLE), "s00": ("age at diagnosis\n(>10 vs <6 y)", ORANGE), "s09": ("household income", BLUE),
       "s10": ("mother's education", BLUE), "s11": ("maternal age", BLUE), "s12": ("paternal age", BLUE), "s08": ("ADI", BLUE),
       "s27": ("ADHD", DGREY), "s28": ("mood/anxiety", DGREY), "s29": ("ODD/conduct", DGREY)}
OFF = {"F2-F1": (4, 6), "s00": (5, -4), "s09": (4, 2), "s10": (4, -6), "s11": (-62, -2), "s12": (4, 0), "s08": (-50, -10), "s27": (4, 2), "s28": (-8, 8), "s29": (4, 0)}
fig = plt.figure(figsize=(7.2, 3.7)); gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.55)
a = fig.add_subplot(gs[0])
for r in V.itertuples():
    if r.code in LAB: continue
    a.annotate("", xy=(r.d_psych, r.d_repro), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.7, mutation_scale=6))
for r in V.itertuples():
    if r.code not in LAB: continue
    lab, col = LAB[r.code]; lw = 2.6 if r.code in ("F2-F1", "s00") else 1.6
    a.annotate("", xy=(r.d_psych, r.d_repro), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, mutation_scale=9))
    if np.isfinite(r.se_psych): a.errorbar(r.d_psych, r.d_repro, xerr=1.96 * r.se_psych, yerr=1.96 * r.se_repro, fmt="none", ecolor=col, elinewidth=0.5, alpha=0.5)
    a.annotate(lab, (r.d_psych, r.d_repro), xytext=OFF[r.code], textcoords="offset points", fontsize=6, color=col, fontweight="bold" if r.code in ("F2-F1", "s00") else None)
a.axhline(0, color=MUTED, lw=0.5); a.axvline(0, color=MUTED, lw=0.5)
a.set_xlabel("difference in psychiatric genetic signal\n(rg with ADHD and MDD, high minus low half)")
a.set_ylabel("missing EDU/SES axis: difference in reproductive-timing /\neducation signal (rg with age at first birth and educational attainment)")
a.set_xlim(-0.05, 0.95); a.set_ylim(-1.2, 0.45)
a.set_title("a  Zhang F1/F2 and the age split run along the psychiatric axis;\n    socioeconomic splits reveal the missing EDU/SES axis", loc="left", fontweight="bold", fontsize=7.2)
a.grid(color=GRID, lw=0.4); a.set_axisbelow(True)
b = fig.add_subplot(gs[1]); C = pd.read_csv(os.path.join(R, "rare_compare.tsv"), sep="\t")
order = ["age at diagnosis (top vs bottom third)", "age at diagnosis (>10 vs <6 y; paper)", "ADHD (yes vs no)", "paternal age (older vs younger third)",
         "mother no degree vs degree", "deprivation (high vs low third)", "maternal age (older vs younger third)", "household income (low vs high third)"]
for k, (oc, dy, mk) in enumerate([("A probands, absolute", 0.18, "o"), ("B unaffected siblings, absolute", 0.0, "s"), ("C within-family excess (proband - sibling)", -0.18, "^")]):
    x = C[C.outcome == oc].set_index("split").reindex(order)
    cols = [ORANGE if s.startswith("age at") else (DGREY if s.startswith("ADHD") else BLUE) for s in order]
    b.scatter(-np.log10(x.p_wald), np.arange(len(order)) + dy, c=cols, marker=mk, s=18, zorder=3, edgecolor="white", linewidth=0.4)
    b.scatter([], [], marker=mk, color=MUTED, s=18, label=oc.split(" ", 1)[1].replace(", absolute", ""))
b.axvline(-np.log10(0.05), color=MUTED, lw=0.5, ls=":")
b.set_yticks(range(len(order))); b.set_yticklabels([s.replace(" (top vs bottom third)", " (thirds)").replace(" (older vs younger third)", "").replace(" (low vs high third)", "").replace(" (high vs low third)", "") for s in order], fontsize=6)
b.set_xlabel("-log10 P, rare-variant stratification\n(Wald, 4 classes: de novo / inherited x NDD / other genes)")
b.legend(frameon=False, fontsize=5.6, loc="lower right"); b.grid(axis="x", color=GRID, lw=0.4); b.set_axisbelow(True)
b.set_title("b  Rare variants follow the missing EDU/SES axis,\n    not age at diagnosis", loc="left", fontweight="bold", fontsize=7.2)
fig.savefig(os.path.join(MS, "figures", "L_AxisMissing.png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(MS, "figures", "L_AxisMissing.pdf"), bbox_inches="tight")
print("figure written")
