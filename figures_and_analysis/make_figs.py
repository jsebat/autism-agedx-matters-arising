#!/usr/bin/env python3
"""Figures for the Matters Arising on Zhang, Grove et al. Nature 2025 (age at diagnosis and two autism polygenic factors).
Inputs: aggregate tables in ../results (copied from Expanse ~/rvh2/agedx, cnv_drift/mgwas/combined, cnv_drift/pcproj).
Outputs: ../figures/Fig1.png, Fig2.png, ED1.png, ED2.png, ED3.png (+ .pdf)."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2, norm
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results"); F = os.path.join(H, "..", "figures"); os.makedirs(F, exist_ok=True)
BLUE, ORANGE, AQUA, YELLOW, MAG, GREEN = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"
INK, MUTED, GRID = "#222222", "#6b6b6b", "#dddddd"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
AGE = "age at diagnosis (after 10 vs before 6)"
SHORT = {AGE: "Age at diagnosis (>10 vs <6 y)", "family structure (0 vs 2 parents genotyped)": "No parent genotyped (vs both)", "adult at registration": "Registered as adult",
         "female": "Female", "diagnosed 2019-26 vs before 2010": "Diagnosed 2019-26 (vs <2010)", "platform Twist vs GSA": "Platform (Twist vs GSA)", "WGS vs exome": "Platform (WGS vs exome)",
         "intellectual disability": "Intellectual disability", "minimal language": "Minimal language", "FSIQ (top vs bottom third)": "FSIQ", "SCQ (top vs bottom third)": "SCQ",
         "Vineland ABC (top vs bottom third)": "Vineland ABC", "RBS-R (top vs bottom third)": "RBS-R", "CBCL total (top vs bottom third)": "CBCL total", "ADHD": "ADHD",
         "mood or anxiety disorder": "Mood or anxiety disorder", "seizures": "Seizures", "language regression": "Language regression", "multiplex family": "Multiplex family",
         "diagnosed by specialist": "Diagnosed by specialist", "ASD validity flag": "SPARK validity flag", "area deprivation (top vs bottom third)": "Area deprivation (ADI)"}


def lp_chi(x, df): return -np.log10(np.maximum(chi2.sf(x, df), 1e-300))
def lp_z(z): return -np.log10(np.maximum(2 * norm.sf(np.abs(z)), 1e-300))


def save(fig, name):
    fig.savefig(os.path.join(F, name + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(F, name + ".pdf"), bbox_inches="tight"); plt.close(fig)


# ---------------- Figure 1 ----------------
comp = pd.read_csv(os.path.join(R, "composition.tsv"), sep="\t", index_col=0)
S = pd.read_csv(os.path.join(R, "driver_summary.tsv"), sep="\t")
RD = pd.read_csv(os.path.join(R, "rare_driver_splits.tsv"), sep="\t")
fig = plt.figure(figsize=(7.4, 4.0))
gs = fig.add_gridspec(1, 2, width_ratios=[0.75, 2.75], wspace=0.95)
ax = fig.add_subplot(gs[0])
rows = [("pct_0_parents", "No parent genotyped"), ("pct_adult_at_reg", "Registered as adult"), ("female_pct", "Female"), ("pct_dx_2019-26", "Diagnosed 2019-26"),
        ("pct_spark_twist", "Genotyped on Twist"), ("pct_ID", "Intellectual disability")]
y = np.arange(len(rows))[::-1]
for (k, lab), yy in zip(rows, y):
    e, l = comp.loc[k, "before 6"], comp.loc[k, "after 10"]
    ax.plot([e, l], [yy, yy], color=GRID, lw=1.5, zorder=1)
    ax.scatter([e], [yy], s=22, color=BLUE, zorder=2, edgecolor="white", linewidth=0.6); ax.scatter([l], [yy], s=22, color=ORANGE, zorder=2, edgecolor="white", linewidth=0.6)
ax.set_yticks(y); ax.set_yticklabels([r[1] for r in rows]); ax.set_xlim(0, 80); ax.set_xlabel("% of autistic participants")
ax.scatter([], [], color=BLUE, s=22, label=f"Diagnosed <6 y (n={int(comp.loc['n','before 6']):,})"); ax.scatter([], [], color=ORANGE, s=22, label=f"Diagnosed >10 y (n={int(comp.loc['n','after 10']):,})")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.13), fontsize=6, handletextpad=0.2)
ax.set_title("a  SPARK strata differ in design", loc="left", fontweight="bold", fontsize=7.5)
ax.grid(axis="x", color=GRID, lw=0.4); ax.set_axisbelow(True)
# b: heterogeneity by split, three measures
pop = S[S.measure == "population"].set_index("split"); wf = S[S.measure != "population"].set_index("split")
rr = RD[(RD.score == "rare SNV/indel (all)")].set_index("split"); cn = RD[RD.score == "CNV"].set_index("split")
splits = [s for s in pop.sort_values("chi2_profile", ascending=False).index]
for s in rr.index:
    if s not in splits: splits.append(s)
cols = [("PGS profile\npopulation\n(cases vs relatives)", pop.chi2_profile, 9), ("PGS profile\nwithin family\n(pTDT)", wf.chi2_profile, 9),
        ("Rare SNV/indel\nwithin family\n(case vs sibling)", rr.z, None), ("Rare CNV\nwithin family\n(case vs sibling)", cn.z, None)]
sub = gs[1].subgridspec(1, 4, wspace=0.12)
yy = np.arange(len(splits))[::-1]
for j, (tt, ser, df) in enumerate(cols):
    a = fig.add_subplot(sub[j])
    v = np.array([(lp_chi(ser[s], df) if df else lp_z(ser[s])) if (s in ser.index and np.isfinite(ser[s]) and (df is None or ser[s] > 0)) else np.nan for s in splits])
    colr = [ORANGE if s == AGE else BLUE for s in splits]
    a.barh(yy, np.nan_to_num(v), color=colr, height=0.62)
    for k, s in enumerate(splits):
        if not np.isfinite(v[k]): a.text(0.2, yy[k], "n/a", va="center", fontsize=5.5, color=MUTED)
    rk = pd.Series(v, index=splits).rank(ascending=False)
    a.text(0.98, 0.02, f"rank {int(rk[AGE])}/{int(np.isfinite(v).sum())}", transform=a.transAxes, ha="right", va="bottom", fontsize=6, color=ORANGE)
    a.set_xscale("symlog", linthresh=10) if j == 0 else None
    a.axvline(-np.log10(0.05 / len(splits)), color=MUTED, lw=0.5, ls=":")
    a.set_yticks(yy); a.set_yticklabels([SHORT.get(s, s) for s in splits] if j == 0 else []); a.tick_params(axis="y", length=0)
    a.set_title(tt, fontsize=6); a.set_xlabel("-log10 P")
    a.grid(axis="x", color=GRID, lw=0.4); a.set_axisbelow(True)
    if j == 0:
        a.text(-0.02, 1.2, "b  Which case characteristic separates autism genetics?", transform=a.transAxes, ha="right" if False else "left", fontweight="bold", fontsize=7.5)
for t in fig.axes[1].get_yticklabels():
    if t.get_text() == SHORT[AGE]: t.set_color(ORANGE); t.set_fontweight("bold")
save(fig, "Fig1")

# ---------------- Figure 2 (matched-design LDSC) ----------------
A = pd.read_csv(os.path.join(R, "rg_all.tsv"), sep="\t")
def rg(a, b):
    x = A[((A.p1 == a) & (A.p2 == b)) | ((A.p1 == b) & (A.p2 == a))]
    return (x.rg.iloc[0], x.se.iloc[0]) if len(x) else (np.nan, np.nan)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.55))
# a: within-pair rg per split
a = axs[0]
PAIRS = [("agedx", "early", "late", "Age at diagnosis\n(<6 vs >10 y)"), ("adi", "low", "high", "Area deprivation\n(low vs high third)"), ("adhd", "no", "yes", "ADHD\n(no vs yes)"),
         ("adiE", "low", "high", "ADI, within cases\ndiagnosed <6 y"), ("adhdE", "no", "yes", "ADHD, within cases\ndiagnosed <6 y"), ("sexE", "male", "female", "Sex, within cases\ndiagnosed <6 y"),
         ("agedxC", "early", "late", "Age at diagnosis, within\nchild registrants")]
lab, est = [], []
for sp, h1, h2, l in PAIRS:
    r, s = rg(f"between_{sp}_{h1}", f"between_{sp}_{h2}")
    if np.isfinite(r): lab.append(l); est.append((r, s, sp))
yy = np.arange(len(lab))[::-1]
for k, (r, s, sp) in enumerate(est):
    c = ORANGE if sp.startswith("agedx") else BLUE
    a.errorbar(np.clip(r, -1.5, 2), yy[k], xerr=1.96 * s, fmt="o", color=c, ms=4, elinewidth=0.9, capsize=0)
a.axvline(1, color=MUTED, lw=0.5, ls="--"); a.axvline(0, color=MUTED, lw=0.5)
a.set_yticks(yy); a.set_yticklabels(lab, fontsize=6); a.set_xlim(-2, 2.5); a.set_xlabel("Genetic correlation between the two halves (rg, 95% CI)")
a.set_title("a  A low rg between halves is not specific\n    to age at diagnosis (between-family GWAS)", loc="left", fontweight="bold", fontsize=7.5)
a.grid(axis="x", color=GRID, lw=0.4); a.set_axisbelow(True)
# b: correlation with late-diagnosis GWAS vs within-family autism signal
a = axs[1]
GW = [("between_adult_adult", "Registered as adult"), ("between_adi_high", "High ADI"), ("between_adhd_yes", "ADHD"), ("between_adiE_high", "High ADI, dx <6 y"),
      ("between_adhdE_yes", "ADHD, dx <6 y"), ("between_sexE_female", "Female, dx <6 y"), ("between_agedxC_late", "Dx >10 y, child registrants"), ("between_agedx_early", "Dx <6 y")]
GW = [g for g in GW if np.isfinite(rg(g[0], "between_agedx_late")[0]) or g[0] == "between_agedx_early"]
yy = np.arange(len(GW))[::-1]
for k, (g, l) in enumerate(GW):
    r1, s1 = rg(g, "between_agedx_late"); r2, s2 = rg(g, "sibtrio_agedx_early")
    if np.isfinite(r1): a.errorbar(np.clip(r1, -1.5, 2), yy[k] + 0.15, xerr=1.96 * s1, fmt="o", color=ORANGE, ms=4, elinewidth=0.9)
    if np.isfinite(r2): a.errorbar(np.clip(r2, -1.5, 2), yy[k] - 0.15, xerr=1.96 * s2, fmt="s", color=AQUA, ms=3.5, elinewidth=0.9)
a.errorbar([], [], xerr=[], fmt="o", color=ORANGE, label="rg with late-diagnosis GWAS (>10 y)"); a.errorbar([], [], xerr=[], fmt="s", color=AQUA, label="rg with within-family autism signal (dx <6 y trios)")
a.axvline(1, color=MUTED, lw=0.5, ls="--"); a.axvline(0, color=MUTED, lw=0.5)
a.set_yticks(yy); a.set_yticklabels([l for _, l in GW], fontsize=6); a.set_xlim(-1.5, 2); a.set_xlabel("rg (95% CI)")
a.legend(frameon=False, fontsize=5.8, loc="lower center", bbox_to_anchor=(0.4, -0.45))
a.set_title("b  What the late-diagnosis GWAS tracks", loc="left", fontweight="bold", fontsize=7.5)
a.grid(axis="x", color=GRID, lw=0.4); a.set_axisbelow(True)
save(fig, "Fig2")

# ---------------- ED1: PC-AiR membership and projection shrinkage ----------------
ins = pd.read_csv(os.path.join(R, "eval_inset.tsv"), sep="\t"); shr = pd.read_csv(os.path.join(R, "eval_shrinkage.tsv"), sep="\t", index_col=0)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw=dict(wspace=0.45))
a = axs[0]; x = ins[ins.group.str.startswith("case")].copy()
x["st"] = x.group.str.split(", ").str[1]; x["np"] = x.group.str.split(", ").str[2]
for i, st in enumerate(["before 6", "6-10", "after 10"]):
    z = x[x.st == st]
    a.bar(np.arange(3) + (i - 1) * 0.26, z["mean"] * 100, width=0.25, color=[BLUE, YELLOW, ORANGE][i], label={"before 6": "dx <6 y", "6-10": "dx 6-10 y", "after 10": "dx >10 y"}[st])
a.set_xticks(range(3)); a.set_xticklabels(["0 parents", "1 parent", "2 parents"]); a.set_ylabel("% in PCA-defining (unrelated) set"); a.legend(frameon=False, fontsize=6)
a.set_title("a  Who defines the PCs (PC-AiR, Twist EUR)", loc="left", fontweight="bold", fontsize=7.5)
a = axs[1]; v = shr.iloc[0].dropna().astype(float); v = v[[c for c in v.index if c.startswith("PC")]]
a.bar(range(len(v)), v.values, color=BLUE); a.axhline(1, color=MUTED, lw=0.5, ls="--"); a.set_xticks(range(len(v))); a.set_xticklabels([c.replace("_AVG", "") for c in v.index], fontsize=6)
a.set_ylabel("SD of projected / SD of fitted PC scores"); a.set_ylim(0, 1.1)
a.set_title("b  Projected individuals are shrunk to the centre", loc="left", fontweight="bold", fontsize=7.5)
save(fig, "ED1")

# ---------------- ED2: external genetic correlations ----------------
E = pd.read_csv(os.path.join(R, "rg_external.tsv"), sep="\t")
order = ["between_agedx_early", "between_agedx_late", "between_adult_child", "between_adult_adult", "between_adi_low", "between_adi_high", "between_adhd_no", "between_adhd_yes",
         "between_adiE_low", "between_adiE_high", "between_adhdE_no", "between_adhdE_yes", "between_sexE_male", "between_sexE_female", "between_agedxC_early", "between_agedxC_late",
         "sibtrio_agedx_early", "sibtrio_adult_child", "sibtrio_adi_low", "sibtrio_adhd_no"]
TR = ["asd_ipsych_dbs", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018"]
M = E.pivot_table(index="p1", columns="trait", values="rg").reindex(index=[o for o in order if o in set(E.p1)], columns=TR)
Se = E.pivot_table(index="p1", columns="trait", values="se").reindex(index=M.index, columns=TR)
fig, a = plt.subplots(figsize=(5.2, 0.25 * len(M) + 1.2))
im = a.imshow(np.clip(M.values, -1, 1), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        r, s = M.values[i, j], Se.values[i, j]
        if np.isfinite(r): a.text(j, i, f"{r:.2f}\n({s:.2f})", ha="center", va="center", fontsize=4.8, color="white" if abs(r) > 0.6 else INK)
a.set_xticks(range(len(TR))); a.set_xticklabels(["Autism\n(iPSYCH)", "ADHD", "MDD", "Bipolar", "SCZ", "EA", "IQ"], fontsize=6)
a.set_yticks(range(len(M))); a.set_yticklabels([m.replace("between_", "between: ").replace("sibtrio_", "within-family: ").replace("_", " ") for m in M.index], fontsize=6)
plt.colorbar(im, ax=a, fraction=0.03, label="rg (clipped to ±1)")
a.set_title("Genetic correlation of each SPARK matched-design GWAS with external traits", fontsize=7.5, loc="left")
save(fig, "ED2")

# ---------------- ED3: rare-variant detail ----------------
sc = ["de novo", "inherited", "CNV"]; lab = {"de novo": "De novo SNV/indel", "inherited": "Inherited SNV/indel", "CNV": "Rare CNV (GSA v1 + exome)"}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.6), sharey=True, gridspec_kw=dict(wspace=0.08))
base = RD[RD.score == "rare SNV/indel (all)"].assign(c=lambda d: d.z.abs()).sort_values("c").split.tolist()
for a, s in zip(axs, sc):
    x = RD[RD.score == s].set_index("split").reindex(base)
    yy = np.arange(len(base))
    a.errorbar(x["diff"], yy, xerr=1.96 * x.se, fmt="o", ms=3, color=BLUE, ecolor=BLUE, elinewidth=0.7)
    k = base.index(AGE); a.errorbar(x["diff"].iloc[k], k, xerr=1.96 * x.se.iloc[k], fmt="o", ms=4, color=ORANGE, ecolor=ORANGE, elinewidth=1)
    a.axvline(0, color=MUTED, lw=0.5); a.set_title(lab[s], fontsize=7); a.set_xlabel("Difference in within-family excess\n(high minus low group, sibling SD)")
    a.grid(axis="x", color=GRID, lw=0.4); a.set_axisbelow(True)
axs[0].set_yticks(range(len(base))); axs[0].set_yticklabels([SHORT.get(s, s) for s in base], fontsize=6)
save(fig, "ED3")
print("figures written")
