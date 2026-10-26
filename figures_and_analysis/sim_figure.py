#!/usr/bin/env python3
"""Supplementary Figure 5: two factors from one autism liability (sim_twofactor.py, sim_twofactor2.py; summary-level simulation, no participant data).
(a) rg between the SPARK-like early- and late-diagnosed GWAS against the SNP heritability of the stratifier, by how strongly the stratifier drives age at
diagnosis (a); observed values for reference. (b) Share of replicates in which two correlated factors fit better than one (AIC), and the fitted factor
correlation, for selected scenarios. Inputs: results/full_data/sim_twofactor.tsv, sim_twofactor2.tsv. Output: Manuscripts/Nature_LTE/figures/L_Sim.(png|pdf)"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results", "full_data"); F = os.path.join(H, "..", "..", "..", "Manuscripts", "Nature_LTE", "figures")
INK, MUTED, GRID, ORANGE, BLUE, TEAL, GREY = "#222222", "#6b6b6b", "#dddddd", "#eb6834", "#2a78d6", "#1b8a8f", "#9a9a9a"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK, "text.color": INK})
s1 = pd.read_csv(os.path.join(R, "sim_twofactor.tsv"), sep="\t"); s2 = pd.read_csv(os.path.join(R, "sim_twofactor2.tsv"), sep="\t")
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1, 1.25], wspace=0.75))
# (a) rg vs stratifier heritability (rho = 0), one line per link strength a; a = 1 (split directly on the stratifier) from the extended grid
x0 = s1[s1.rho.eq(0)]
for av, col in [(0.0, GREY), (0.3, BLUE), (0.6, TEAL)]:
    d = x0[x0.a.eq(av)].sort_values("h2S")
    a.errorbar(d.h2S, d.rg_spark_median, yerr=[d.rg_spark_median - d.rg_spark_q05, d.rg_spark_q95 - d.rg_spark_median], fmt="o-", color=col, ms=3.5, lw=1.2, elinewidth=0.7, label=f"a = {av}")
dd = s2[s2.name.isin(["A direct, h2S 0.3", "A direct, h2S 0.5"])].assign(h2S=[0.3, 0.5])
d1 = pd.concat([pd.DataFrame(dict(h2S=[0.0], rg_spark_median=[1.0], rg_spark_q05=[1.0], rg_spark_q95=[1.0])), dd])
a.errorbar(d1.h2S, d1.rg_spark_median, yerr=[d1.rg_spark_median - d1.rg_spark_q05, d1.rg_spark_q95 - d1.rg_spark_median], fmt="s-", color=ORANGE, ms=3.5, lw=1.2, elinewidth=0.7,
           label="a = 1 (split directly\non the stratifier)")
a.axhline(0.11, color=INK, lw=0.7, ls="--"); a.text(0.01, 0.13, "observed in SPARK, 0.11", fontsize=5.8, color=INK)
a.axhline(1, color=MUTED, lw=0.5, ls=":")
a.set_xlabel("SNP heritability of the stratifier"); a.set_ylabel("rg, early- vs late-diagnosed GWAS (SPARK-like)")
a.set_xticks([0, 0.1, 0.3, 0.5]); a.set_ylim(0, 1.08); a.set_xlim(-0.03, 0.55)
a.set_title("a  One liability: rg falls with stratifier genetics", loc="left", fontweight="bold", fontsize=7.2)
a.legend(frameon=False, fontsize=5.8, loc="lower left", bbox_to_anchor=(0.0, 0.17), title="age at diagnosis\ndriven by stratifier", title_fontsize=5.8)
a.grid(axis="y", color=GRID, lw=0.4); a.set_axisbelow(True)
# (b) two-factor preference and factor correlation for selected scenarios
rows = [("Non-heritable stratifier", s1[(s1.h2S == 0) & (s1.rho == 0) & (s1.a == 0.6)].iloc[0]),
        ("h² 0.1, a 0.3", s1[(s1.h2S == 0.1) & (s1.rho == 0) & (s1.a == 0.3)].iloc[0]),
        ("h² 0.3, a 0.6", s1[(s1.h2S == 0.3) & (s1.rho == 0) & (s1.a == 0.6)].iloc[0]),
        ("Split on stratifier, h² 0.3", s2[s2.name.eq("A direct, h2S 0.3")].iloc[0]),
        ("Split on stratifier, h² 0.5", s2[s2.name.eq("A direct, h2S 0.5")].iloc[0]),
        ("Three stratifiers, h² 0.3", s2[s2.name.eq("B three stratifiers, h2S 0.3 each")].iloc[0]),
        ("Late-only ascertainment only", s2[s2.name.eq("D late-only asc on S2, gamma 1.0, a 0")].iloc[0]),
        ("Three stratifiers h² 0.5\n+ late-only ascertainment", s2[s2.name.str.startswith("E combined")].iloc[0])]
y = np.arange(len(rows))[::-1]
for yi, (lab, r) in zip(y, rows):
    b.barh(yi, r.frac_aic_2f, color=TEAL if r.frac_aic_2f > 0.5 else GREY, height=0.6)
    b.text(1.03, yi, f"rg {r.rg_spark_median:.2f}; factor r {r.phi_2f_median:.2f}", va="center", fontsize=5.6, color=MUTED)
b.set_yticks(y); b.set_yticklabels([l for l, _ in rows], fontsize=5.8); b.set_xlim(0, 1); b.axvline(0.05, color=MUTED, lw=0.5, ls=":")
b.set_xlabel("replicates preferring two correlated factors (AIC)")
b.set_title("b  Two factors without two aetiologies", loc="left", fontweight="bold", fontsize=7.2)
b.grid(axis="x", color=GRID, lw=0.4); b.set_axisbelow(True)
fig.savefig(os.path.join(F, "L_Sim.png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(F, "L_Sim.pdf"), bbox_inches="tight"); plt.close(fig)
print("figure written")
