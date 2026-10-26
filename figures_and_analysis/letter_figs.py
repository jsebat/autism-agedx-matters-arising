#!/usr/bin/env python3
"""Figures and sample-size tables for the Matters Arising on Zhang, Grove et al. Nature 2025.
Inputs (aggregate only): ../results/replication/* (their-pipeline replication, broad survey, random-split null, per-platform counts),
../results/survey_splits.tsv (split definitions, phi with late diagnosis). Outputs: ../figures/L_Fig1.*, L_Fig2.*, ../results/tables/*.tsv
DATA=all (default): all SPARK strata (GSA v1, GSA v2, Twist, iWGS v1.1, WGS 2026-08; results/full_data/all_*.tsv, all_ldsc.sbatch) -> L_Fig1, survey_table.tsv.
DATA=arrays: the arrays genotyped in Zhang et al.'s era (results/replication) -> L_Fig1_arrays, survey_table_arrays.tsv (supplementary sensitivity)."""
import os, glob
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results"); RP = os.path.join(R, "replication"); MS = os.path.join(H, "..", "..", "..", "Manuscripts", "Nature_LTE"); F = os.path.join(MS, "figures")
T = os.path.join(MS, "tables"); os.makedirs(T, exist_ok=True)
BLUE, ORANGE, AQUA, GREY, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#9a9a9a", "#222222", "#6b6b6b", "#dddddd"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, n): fig.savefig(os.path.join(F, n + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(F, n + ".pdf"), bbox_inches="tight"); plt.close(fig)

DATA = os.environ.get("DATA", "all"); SUF = "" if DATA == "all" else "_arrays"; FD = os.path.join(R, "full_data")
# ---------- counts ----------
if DATA == "all":
    ca = pd.read_csv(os.path.join(FD, "all_counts_all.tsv"), sep="\t")
    cs, cn, cr = ca[ca.gwas.str.match(r"s\d\d_")].copy(), ca[ca.gwas.str.match(r"(ageN|adiN)\d")].copy(), ca[ca.gwas.str.match(r"W_")].copy()
else:
    cnt = lambda pat: pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(RP, "counts", pat))])
    cs, cn, cr = cnt("counts_survey_*.tsv"), cnt("counts_null_*.tsv"), cnt("counts_spark_*.tsv")
sv = pd.read_csv(os.path.join(R, "survey_splits.tsv"), sep="\t").set_index("code")
SHORTNAME = {"s99": "design: no parent vs both parents genotyped"}
cs["code"] = cs.gwas.str[:3]; cs["half"] = np.where(cs.gwas.str.contains("_hi_"), "hi", "lo")
N = cs.groupby(["code", "half"]).agg(cases=("cases", "sum"), controls=("controls", "sum")).unstack()
N.columns = [f"{a}_{b}" for a, b in N.columns]
rg = pd.read_csv(os.path.join(FD, "all_survey_rg_all.tsv") if DATA == "all" else os.path.join(RP, "survey_rg.tsv"), sep="\t").set_index("code")
S = rg.join(N).join(sv[["split", "phi_late_dx", "n_high", "n_low"]].rename(columns={"n_high": "registered_hi", "n_low": "registered_lo"}))
S["split"] = S.split.fillna(S.index.to_series().map(SHORTNAME))
S["n_min"] = S[["cases_hi", "cases_lo"]].min(axis=1)
S["is_sparse"] = (S.h2z_hi < 3) | (S.h2z_lo < 3)
S.reset_index().to_csv(os.path.join(T, f"survey_table{SUF}.tsv"), sep="\t", index=False)
nul = pd.read_csv(os.path.join(FD, "all_null_rg_all.tsv"), sep="\t").rename(columns={"p1": "rep"}) if DATA == "all" else pd.read_csv(os.path.join(RP, "null_rg.tsv"), sep="\t")

# ---------- Figure 1 ----------
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw=dict(width_ratios=[0.8, 2.2], wspace=0.28))
a = axs[0]
if DATA == "all":
    pairs = pd.read_csv(os.path.join(FD, "all_rg_pairs_all.tsv"), sep="\t"); ag = pairs[pairs.pair.str.startswith("age split")].copy()
    ag["controls"] = ag.pair.str.replace("age split, ", "", regex=False); ag = ag.set_index("controls")
else:
    pairs = pd.read_csv(os.path.join(RP, "rg_pairs.tsv"), sep="\t"); ag = pairs[pairs.pair.str.startswith("W_under6")].set_index("controls")
rng = np.random.default_rng(3); z = nul[nul.set == "ageN"]
a.scatter(np.full(len(z), 0) + rng.uniform(-0.12, 0.12, len(z)), z.rg, s=10, color=GREY, alpha=0.8, linewidth=0, label="random splits, same sizes (n=16)")
for i, (lab, col) in enumerate([("shared controls", ORANGE), ("disjoint controls", ORANGE)]):
    r = ag.loc[lab]; a.errorbar(1 + i * 0.7, r.rg, yerr=1.96 * r.se, fmt="o" if i == 0 else "s", color=col, ms=4, elinewidth=0.9)
a.axhline(1, color=MUTED, lw=0.5, ls="--"); a.axhline(0, color=MUTED, lw=0.5)
a.set_xticks([0, 1, 1.7]); a.set_xticklabels(["random\nsplits", "age split\nshared\ncontrols", "age split\ndisjoint\ncontrols"], fontsize=6)
a.set_ylabel("Genetic correlation between the two halves (rg)"); a.set_ylim(-0.6, 2.0); a.set_xlim(-0.4, 2.2)
a.set_title("a  Replication of the age split", loc="left", fontweight="bold", fontsize=7.5)
a.text(0.02, 0.02, "Zhang et al.: rg 0.02 (0.13)", transform=a.transAxes, fontsize=6, color=MUTED)
a = axs[1]
ok = S.rg.notna() & S.se.notna() & (S.se < 1.5)
x = S[ok]
for _, r in x.iterrows():
    age = r.name == "s00"; strong = (abs(r.phi_late_dx) < 0.1) & (r.rg <= 0.35) & (not r["is_sparse"])
    col = ORANGE if age else (BLUE if strong else GREY)
    a.errorbar(r.n_min, r.rg, yerr=1.96 * r.se, fmt="o", color=col, mfc="white" if r["is_sparse"] else col, ms=4 if (age or strong) else 3, elinewidth=0.7, alpha=1 if (age or strong) else 0.7)
# null band (mean +- 2 SD of the random splits at their two sizes)
for st in ["ageN", "adiN"]:
    z = nul[nul.set == st]; c1n = cn.reset_index(drop=True)
    nn = min(c1n[c1n.gwas == f"{st}1_A_c2"].cases.sum(), c1n[c1n.gwas == f"{st}1_B_c1"].cases.sum())
    a.add_patch(plt.Rectangle((nn * 0.93, z.rg.mean() - 2 * z.rg.std()), nn * 0.14, 4 * z.rg.std(), color=GREY, alpha=0.25, lw=0))
    a.plot([nn * 0.93, nn * 1.07], [z.rg.mean()] * 2, color=GREY, lw=1)
a.axhline(1, color=MUTED, lw=0.5, ls="--"); a.axhline(0, color=MUTED, lw=0.5)
LAB = {"s00": "Age at diagnosis", "s11": "Maternal age", "s09": "Household income", "s10": "Mother's education", "s29": "ODD/conduct", "s27": "ADHD",
       "s28": "Mood/anxiety", "s08": "Deprivation", "s01": "Registered as adult", "s13": "Intellectual disability", "s32": "CBCL externalizing", "s31": "CBCL internalizing"}
OFF = {"s00": (-30, -14), "s11": (5, -8), "s09": (6, 2), "s10": (6, 4), "s29": (5, 2), "s27": (-14, 14), "s28": (5, 6), "s08": (-58, 6), "s01": (5, 2), "s13": (-80, 4),
       "s32": (5, 0), "s31": (5, -2)}
HA = {}
if DATA == "all":   # label placement for the all-strata positions (the 1-2 x 10^4 cluster)
    LAB.pop("s31")
    OFF.update({"s00": (6, 5), "s09": (-5, -3), "s11": (6, -10), "s10": (6, -4), "s08": (-5, 3), "s27": (-4, 9), "s28": (6, -1), "s01": (-5, 3), "s13": (-6, 4)})
    HA = {"s09": "right", "s08": "right", "s27": "right", "s01": "right", "s13": "right"}
for c, l in LAB.items():
    if c in x.index: a.annotate(l, (x.loc[c, "n_min"], x.loc[c, "rg"]), xytext=OFF[c], textcoords="offset points", fontsize=5.6, ha=HA.get(c, "left"), color=ORANGE if c == "s00" else INK)
a.set_xscale("log"); a.set_xlabel("Cases in the smaller half (summed over " + ("all SPARK genotyping and sequencing strata)" if DATA == "all" else "GSA v1, GSA v2, Twist)")); a.set_ylim(-1.2, 2.6)
a.set_title("b  Broad survey: 43 case characteristics, same pipeline", loc="left", fontweight="bold", fontsize=7.5)
a.scatter([], [], color=ORANGE, s=14, label="age at diagnosis"); a.scatter([], [], color=BLUE, s=14, label="|phi| < 0.1 with late diagnosis, rg <= 0.35, h2 z >= 3")
a.scatter([], [], facecolor="white", edgecolor=GREY, s=12, label="h2 z < 3 in a half (sparse)"); a.fill_between([], [], color=GREY, alpha=0.25, label="random-split null, mean +/- 2 SD")
a.legend(frameon=False, fontsize=5.6, loc="upper left")
a.grid(axis="y", color=GRID, lw=0.4); a.set_axisbelow(True)
save(fig, "L_Fig1" + SUF)

# ---------- sample-size tables ----------
cr = cr.assign(design="replication / robustness"); cr.to_csv(os.path.join(T, f"n_replication_by_platform{SUF}.tsv"), sep="\t", index=False)
cs.to_csv(os.path.join(T, f"n_survey_by_platform{SUF}.tsv"), sep="\t", index=False); cn.to_csv(os.path.join(T, f"n_null_by_platform{SUF}.tsv"), sep="\t", index=False)
print(S.sort_values("rg")[["split", "rg", "se", "h2z_hi", "h2z_lo", "cases_hi", "cases_lo", "controls_hi", "controls_lo", "phi_late_dx"]].round(2).head(16).to_string())
print("figures written")
