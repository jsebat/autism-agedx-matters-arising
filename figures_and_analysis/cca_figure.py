#!/usr/bin/env python3
"""Extended Data figure: the strongest axes linking characteristics of autistic people to their genetics (regularized CCA), FULL data
, main sample = registered before 22 (Zhang et al.'s criterion), genotyping platform (+ sequencing cohort) and family structure
regressed out. (a) Out-of-sample canonical correlation per axis vs permutation null: polygenic scores between families (cca2.py, ANALYSIS=full), within
families (pTDT), and polygenic + rare variants in trio probands (cca2.py, ANALYSIS=trio). (b) Structure coefficients of the four axes of the combined
analysis for selected genetic and phenotypic variables; age at diagnosis highlighted.
Inputs: results/full_data/cca_full_u22_partialfam_*, cca_full_u22_within_partial_*, cca_trio_all_*. Output: Manuscripts/Nature_LTE/figures/L_CCA.(png|pdf)"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results", "full_data"); MS = os.path.join(H, "..", "..", "..", "Manuscripts", "Nature_LTE")
BLUE, ORANGE, AQUA, GREY, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#252525", "#222222", "#6b6b6b", "#dddddd"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(7.2, 4.4)); gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1.6], wspace=0.55)
a = fig.add_subplot(gs[0])
SETS = [("cca_full_u22_partialfam", BLUE, -0.26, "polygenic scores, between families"), ("cca_full_u22_within_partial", AQUA, 0.0, "polygenic scores, within families (pTDT)"),
        ("cca_trio_all", GREY, 0.26, "polygenic + rare, trio probands")]
for tag, col, dx, lab in SETS:
    cv = pd.read_csv(os.path.join(R, f"{tag}_cv.tsv"), sep="\t").head(5)
    a.bar(cv.component + dx, cv.cv_canonical_r, width=0.25, color=col, label=f"{lab} (n={int(cv.n.iloc[0]):,})")
    a.scatter(cv.component + dx, cv.null_95, marker="_", color=ORANGE, s=30, zorder=3)
a.scatter([], [], marker="_", color=ORANGE, s=30, label="permutation null, 95th percentile")
a.set_ylim(0, 0.6); a.set_xticks(range(1, 6)); a.set_xlabel("canonical axis"); a.set_ylabel("out-of-sample canonical correlation (5-fold CV)")
a.legend(frameon=False, fontsize=5.2, loc="upper right"); a.set_title("a  Strength of each axis", loc="left", fontweight="bold", fontsize=7.2)
b = fig.add_subplot(gs[1]); st = pd.read_csv(os.path.join(R, "cca_trio_all_structure.tsv"), sep="\t")
G = ["edu_2022", "afb_2018", "iq_2018", "adhd_2022", "mdd_2025", "scz_2022", "asd_2026", "de novo, NDD genes", "inherited, NDD genes", "inherited, other genes"]
GL = ["EA PGS", "age at first birth PGS", "IQ PGS", "ADHD PGS", "MDD PGS", "SCZ PGS", "autism PGS", "De novo – NDD genes", "Inherited – NDD genes", "Inherited – other genes"]
X = ["age at diagnosis", "deprivation (ADI)", "mother without degree", "low household income (third)", "maternal age", "paternal age", "intellectual disability",
     "age walked", "minimal language", "ADHD", "CBCL externalizing"]
AXES = [(1, "EDU/SES"), (2, "DEV (developmental)"), (4, "PSYCH (psychiatric)"), (3, "LANG (language)")]   # combined-analysis components
rows = [("genetics", v, l) for v, l in zip(G, GL)] + [("characteristics", v, v) for v in X]
M = np.array([[st[(st.component == k) & (st.variable == v)].r.iloc[0] for k, _ in AXES] for _, v, _ in rows])
for j in range(M.shape[1]):   # orient each axis so its largest genetic loading is positive
    g = M[:len(G), j]
    if g[np.argmax(np.abs(g))] < 0: M[:, j] *= -1
im = b.imshow(M.T, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
for (i, j), v in np.ndenumerate(M.T): b.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5.0, color="white" if abs(v) > 0.6 else INK)
b.set_xticks(range(len(rows))); b.set_xticklabels([l for _, _, l in rows], rotation=60, ha="right", fontsize=5.5)
for k, (blk, v, _) in enumerate(rows):
    if v == "age at diagnosis": b.get_xticklabels()[k].set_color(ORANGE); b.get_xticklabels()[k].set_fontweight("bold")
b.axvline(len(G) - 0.5, color=INK, lw=0.8)
b.set_xlabel("genetics (polygenic + rare)                                   characteristics of autistic participants", fontsize=6.2)
b.set_yticks(range(len(AXES))); b.set_yticklabels([f"axis {k}: {n}" for k, n in AXES], fontsize=6)
plt.colorbar(im, ax=b, fraction=0.025, pad=0.01, label="structure coefficient")
b.set_title("b  What defines the axes (polygenic + rare, trio probands)", loc="left", fontweight="bold", fontsize=7.2)
fig.savefig(os.path.join(MS, "figures", "L_CCA.png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(MS, "figures", "L_CCA.pdf"), bbox_inches="tight")
print("figure written")
