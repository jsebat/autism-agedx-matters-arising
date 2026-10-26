#!/usr/bin/env python3
"""Letter figures, part 2 (run after letter_figs.py, which builds the survey table S and Figure 1):
Figure 2 (what the halves capture: published-factor resemblance, external rg, our EFA) and Extended Data Figures 1-3."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results"); RP = os.path.join(R, "replication"); MS = os.path.join(H, "..", "..", "..", "Manuscripts", "Nature_LTE"); F = os.path.join(MS, "figures")
BLUE, ORANGE, AQUA, GREY, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#9a9a9a", "#222222", "#6b6b6b", "#dddddd"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, n): fig.savefig(os.path.join(F, n + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(F, n + ".pdf"), bbox_inches="tight"); plt.close(fig)
S = pd.read_csv(os.path.join(MS, "tables", "survey_table.tsv"), sep="\t").set_index("code")
LAB = {"s00": "Age at diagnosis", "s11": "Maternal age", "s09": "Household income", "s10": "Mother's education", "s27": "ADHD", "s28": "Mood/anxiety", "s08": "Deprivation"}

# ---------- Figure 2 ----------
HL = {"s00": ("dx > 10 y", "dx < 6 y"), "s09": ("low income", "high income"), "s10": ("mother no degree", "mother degree"), "s27": ("ADHD", "no ADHD"),
      "s28": ("mood/anxiety", "none"), "s08": ("high ADI", "low ADI"), "s11": ("older mothers", "younger mothers")}
Q = []
for c in HL:
    r = S.loc[c]
    for h in ["hi", "lo"]:
        Q.append(dict(code=c, half=h, d=r[f"{h}_rg_W_F2"] - r[f"{h}_rg_W_F1"], se=np.sqrt(r[f"{h}_se_W_F1"] ** 2 + r[f"{h}_se_W_F2"] ** 2)))
Q = pd.DataFrame(Q)
fig = plt.figure(figsize=(7.2, 6.6)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25], width_ratios=[1.05, 1], hspace=0.45, wspace=0.62)
a = fig.add_subplot(gs[0, 0])
for i, c in enumerate(HL):
    q = Q[Q.code == c].set_index("half"); y = len(HL) - i
    for h, dy, mk in [("hi", 0.15, "o"), ("lo", -0.15, "s")]:
        a.errorbar(q.loc[h, "d"], y + dy, xerr=1.96 * q.loc[h, "se"], fmt=mk, color=ORANGE if c == "s00" else BLUE, mfc="white" if h == "lo" else None, ms=3.5, elinewidth=0.7)
        a.text(1.22, y + dy, HL[c][0 if h == "hi" else 1], va="center", fontsize=5.3, color=MUTED)
a.axvline(0, color=MUTED, lw=0.5); a.set_yticks(range(1, len(HL) + 1)); a.set_yticklabels([LAB[c] for c in HL][::-1], fontsize=6)
a.set_xlabel("rg with Zhang F2 (later-diagnosed) minus\nrg with Zhang F1 (earlier-diagnosed)"); a.set_xlim(-1.3, 1.2)
a.set_title("a  Resemblance of each half to Zhang F1 / F2", loc="left", fontweight="bold", fontsize=7.2)
a.grid(axis="x", color=GRID, lw=0.4); a.set_axisbelow(True)
a = fig.add_subplot(gs[0, 1]); _X = pd.read_csv(os.path.join(R, "full_data", "all_rg_ext_all.tsv"), sep="\t").pivot_table(index="gwas", columns="trait", values="rg")
_S = pd.read_csv(os.path.join(R, "full_data", "all_rg_ext_all.tsv"), sep="\t").pivot_table(index="gwas", columns="trait", values="se")
_H = [("older mothers", "rep_s11_hi_c2"), ("younger mothers", "rep_s11_lo_c1"), ("low income", "rep_s09_hi_c2"), ("high income", "rep_s09_lo_c1"),
      ("diagnosed > 10 y", "rep_s00_hi_c2"), ("diagnosed < 6 y", "rep_s00_lo_c1")]
E = pd.DataFrame({f"{k}_{v}": [(_X if v == "rg" else _S).loc[g, t] for _, g in _H] for k, t in [("afb", "afb_2018"), ("edu", "edu_2022"), ("iq", "iq_2018"), ("adhd", "adhd_2022")] for v in ("rg", "se")},
                 index=[h for h, _ in _H])   # all SPARK strata (all_ldsc.sbatch); arrays-only version: results/afb_rg.tsv
TR = [("afb", "Age at first birth"), ("edu", "Educational attainment"), ("adhd", "ADHD")]; yl = []
for j, (k, lab) in enumerate(TR):
    base = (len(TR) - 1 - j) * (len(E) + 1.6)
    for i, h in enumerate(E.index):
        y = base + len(E) - i; col = ORANGE if "diagnosed" in h else BLUE
        a.errorbar(E.loc[h, f"{k}_rg"], y, xerr=1.96 * E.loc[h, f"{k}_se"], fmt="o" if i % 2 == 0 else "s", color=col, mfc="white" if i % 2 else None, ms=3.2, elinewidth=0.7)
        yl.append((y, h))
    a.text(-1.0, base + len(E) + 0.8, lab, fontsize=6.2, fontweight="bold")
a.set_yticks([y for y, _ in yl]); a.set_yticklabels([h for _, h in yl], fontsize=5.3); a.axvline(0, color=MUTED, lw=0.5)
a.set_xlim(-1.0, 1.0); a.set_xlabel("rg of the half-GWAS with the external GWAS")
a.set_title("b  The halves carry the stratifier's genetics", loc="left", fontweight="bold", fontsize=7.2)
a.grid(axis="x", color=GRID, lw=0.4); a.set_axisbelow(True)
a = fig.add_subplot(gs[1, :]); L3 = pd.read_csv(os.path.join(R, "gsem_all", "efa_k3_loadings.tsv"), sep="\t").set_index("Unnamed: 0")
L3["order"] = L3[["F1", "F2", "F3"]].abs().idxmax(axis=1)
L3 = pd.concat([L3[L3.order == f].sort_values(f, ascending=False) for f in ["F1", "F2", "F3"]])
M = L3[["F1", "F2", "F3"]].to_numpy()
im = a.imshow(M.T, cmap="RdBu_r", vmin=-1.1, vmax=1.1, aspect="auto")
for (i, j), v in np.ndenumerate(M.T): a.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=4.6, color="white" if abs(v) > 0.6 else INK)
NICE = {("s00", "hi"): "diagnosed > 10 y", ("s02", "hi"): "registered 2021+", ("s03", "hi"): "diagnosed 2019+", ("s04", "hi"): "multiplex family",
        ("s05", "hi"): "specialist diagnosis", ("s09", "hi"): "low income", ("s09", "lo"): "high income", ("s10", "hi"): "mother no degree", ("s10", "lo"): "mother degree",
        ("s11", "hi"): "older mothers", ("s11", "lo"): "younger mothers", ("s12", "hi"): "older fathers", ("s18", "hi"): "language regression", ("s21", "hi"): "motor delay",
        ("s23", "hi"): "high SCQ", ("s27", "hi"): "ADHD", ("s27", "lo"): "no ADHD", ("s28", "hi"): "mood/anxiety", ("s29", "hi"): "ODD/conduct",
        ("s32", "hi"): "high CBCL externalizing", ("s34", "hi"): "sleep disorder", ("s42", "hi"): "born before 2005", ("s99", "hi"): "no parent genotyped"}
def short(idx):
    p = idx.split("_"); return NICE.get((p[1], p[2]), idx)
a.set_xticks(range(len(L3))); a.set_xticklabels([short(i) for i in L3.index], rotation=60, ha="right", fontsize=5.2)
a.set_yticks(range(3)); a.set_yticklabels(["A: late diagnosis /\ndesign", "B: socioeconomic\nadvantage", "C: socioeconomic\ndisadvantage /\nexternalizing"], fontsize=6)
for k, idx in enumerate(L3.index):
    if idx == "rep_s00_hi_c2": a.get_xticklabels()[k].set_color(ORANGE); a.get_xticklabels()[k].set_fontweight("bold")
plt.colorbar(im, ax=a, fraction=0.02, pad=0.01, label="loading (promax)")
a.set_title("c  Exploratory genomic factor analysis of 23 split-half GWAS: at least three correlated factors", loc="left", fontweight="bold", fontsize=7.2)
save(fig, "L_Fig2")

# ---------- Extended Data Fig 1 ----------
comp = pd.read_csv(os.path.join(R, "composition.tsv"), sep="\t", index_col=0)
ins = pd.read_csv(os.path.join(R, "eval_inset.tsv"), sep="\t"); shr = pd.read_csv(os.path.join(R, "eval_shrinkage.tsv"), sep="\t", index_col=0)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=0.6, width_ratios=[1.2, 1, 1]))
a = axs[0]; rows = [("pct_0_parents", "No parent genotyped"), ("pct_adult_at_reg", "Registered as adult"), ("female_pct", "Female"), ("pct_dx_2019-26", "Diagnosed 2019-26"), ("pct_ID", "Intellectual disability")]
for k, (c, l) in enumerate(rows):
    y = len(rows) - k; e, lt = comp.loc[c, "before 6"], comp.loc[c, "after 10"]
    a.plot([e, lt], [y, y], color=GRID, lw=1.5); a.scatter([e], [y], color=BLUE, s=16, zorder=3); a.scatter([lt], [y], color=ORANGE, s=16, zorder=3)
a.set_yticks(range(1, len(rows) + 1)); a.set_yticklabels([l for _, l in rows][::-1], fontsize=6); a.set_xlim(0, 80); a.set_xlabel("% of autistic participants")
a.scatter([], [], color=BLUE, s=16, label=f"dx < 6 y (n={int(comp.loc['n','before 6']):,})"); a.scatter([], [], color=ORANGE, s=16, label=f"dx > 10 y (n={int(comp.loc['n','after 10']):,})")
a.legend(frameon=False, fontsize=5.5, loc="lower right"); a.set_title("a  Composition (SPARK EUR)", loc="left", fontweight="bold", fontsize=7)
a = axs[1]; x = ins[ins.group.str.startswith("case")].copy(); x["st"] = x.group.str.split(", ").str[1]
for i, st in enumerate(["before 6", "6-10", "after 10"]):
    z = x[x.st == st]; a.bar(np.arange(3) + (i - 1) * 0.26, z["mean"] * 100, width=0.25, color=[BLUE, GREY, ORANGE][i], label={"before 6": "dx < 6 y", "6-10": "dx 6-10 y", "after 10": "dx > 10 y"}[st])
a.set_xticks(range(3)); a.set_xticklabels(["0", "1", "2"]); a.set_xlabel("Parents genotyped"); a.set_ylabel("% in PC-defining unrelated set"); a.legend(frameon=False, fontsize=5.5)
a.set_title("b  Who defines the PCs", loc="left", fontweight="bold", fontsize=7)
a = axs[2]; v = shr.iloc[0].dropna().astype(float); v = v[[c for c in v.index if c.startswith("PC")]]
a.bar(range(len(v)), v.values, color=BLUE); a.axhline(1, color=MUTED, lw=0.5, ls="--"); a.set_xticks(range(len(v))); a.set_xticklabels([c.replace("_AVG", "").replace("PC", "") for c in v.index], fontsize=6)
a.set_xlabel("PC"); a.set_ylabel("SD projected / SD fitted"); a.set_ylim(0, 1.1); a.set_title("c  Projection shrinkage", loc="left", fontweight="bold", fontsize=7)
save(fig, "L_ED1")

# ---------- Extended Data Fig 2 ----------
gm = pd.read_csv(os.path.join(R, "full_data", "greml_all.tsv"), sep="\t"); gm = gm[gm["sample"] == "all SPARK strata"].rename(columns={"h2": "h2_ivw", "se": "se_ivw"})   # all SPARK strata (greml_combine_all.py)
gm["trait"] = gm.trait.map({"dx_age_u22": "age at diagnosis (registrants < 22 y; their sample restriction)", "late_dx": "diagnosed > 10 y vs < 6 y", "adi_pct": "area deprivation index percentile",
                            "adult_reg": "registered as adult", "adhd": "ADHD"}); gm = gm.dropna(subset=["trait"]).set_index("trait")
lit = [("Townsend deprivation, UK Biobank", 0.21, 0.005, "PMID 27818178"), ("Household income, UK Biobank", 0.11, 0.007, "PMID 27818178"),
       ("Research participation, UK Biobank", 0.125, np.nan, "PMID 37443256"), ("ADHD, PGC/iPSYCH", 0.14, 0.01, "PMID 36702997")]
order = [("age at diagnosis (registrants < 22 y; their sample restriction)", "Age at diagnosis"), ("diagnosed > 10 y vs < 6 y", "Diagnosed > 10 y vs < 6 y"),
         ("area deprivation index percentile", "Area deprivation (ADI)"), ("registered as adult", "Registered as adult"), ("ADHD", "ADHD")]
fig, a = plt.subplots(figsize=(4.8, 2.8)); ys = []
for i, (k, l) in enumerate(order):
    y = len(order) + len(lit) + 1 - i
    a.errorbar(gm.loc[k, "h2_ivw"], y, xerr=1.96 * gm.loc[k, "se_ivw"], fmt="o", color=ORANGE if "diagnos" in l.lower() else BLUE, ms=4, elinewidth=0.8)
    ys.append((y, f"{l} (SPARK cases, n={int(gm.loc[k, 'n']):,})"))
for i, (l, h, s, p) in enumerate(lit):
    y = len(lit) - i; a.errorbar(h, y, xerr=1.96 * s if np.isfinite(s) else None, fmt="D", color=GREY, ms=3.5, elinewidth=0.8); ys.append((y, f"{l} ({p})"))
a.set_yticks([y for y, _ in ys]); a.set_yticklabels([l for _, l in ys], fontsize=5.8); a.axvline(0, color=MUTED, lw=0.5)
a.set_xlabel("SNP heritability (GCTA-GREML in SPARK; published estimates in grey)"); a.grid(axis="x", color=GRID, lw=0.4); a.set_axisbelow(True)
save(fig, "L_ED2")

# ---------- Extended Data Fig 3 (full data, registered before 22) ----------
FD = os.path.join(R, "full_data")
P = pd.read_csv(os.path.join(FD, "ed3a_full.tsv"), sep="\t"); CND = pd.read_csv(os.path.join(FD, "ptdt_full.tsv"), sep="\t")
RC = pd.read_csv(os.path.join(FD, "rare_compare_full.tsv"), sep="\t")
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw=dict(wspace=1.05, width_ratios=[1.25, 0.9, 1]))
a = axs[0]; spl = ["maternal age (older vs younger third)", "household income (low vs high third)", "age at diagnosis (>10 vs <6 y)"]; pg = ["afb_2018", "edu_2022", "adhd_2022", "asd_2026"]
k = 0.0; yt = []
for s_ in spl:
    for p_ in pg:
        for m, mk, col in [("population", "o", BLUE), ("within-family pTDT", "s", AQUA)]:
            r = P[(P.split == s_) & (P.pgs == p_) & (P.measure == m)].iloc[0]
            a.errorbar(r["diff"], k + (0.15 if m == "population" else -0.15), xerr=1.96 * r.se, fmt=mk, color=col, ms=3, elinewidth=0.7)
        yt.append((k, f"{s_.split(' (')[0]}: {p_.split('_')[0].upper()}")); k -= 1
    k -= 0.6
a.set_yticks([y for y, _ in yt]); a.set_yticklabels([l for _, l in yt], fontsize=5.2); a.axvline(0, color=MUTED, lw=0.5)
a.errorbar([], [], fmt="o", color=BLUE, label="between families"); a.errorbar([], [], fmt="s", color=AQUA, label="within families (pTDT)"); a.legend(frameon=False, fontsize=5.3, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2)
a.set_xlabel("PGS difference between halves (SD)"); a.set_title("a  Stratifier PGS", loc="left", fontweight="bold", fontsize=7)
a = axs[1]; c2 = CND[(CND.term == "late") & (CND["sample"] == "registered before 22")]; mods = ["age only", "+ ADHD", "+ ADHD, adult registration, sex, era"]
PGL = [("adhd_2022", "ADHD PGS"), ("mdd_2025", "MDD PGS"), ("asd", "autism PGS")]
for i_, (p_, _) in enumerate(PGL):
    for j_, m in enumerate(mods):
        r = c2[(c2.pgs == p_) & (c2.model == m)].iloc[0]; a.errorbar(r.beta, -(i_ * 4 + j_), xerr=1.96 * r.se, fmt="o", color=[ORANGE, BLUE, GREY][j_], ms=3.5, elinewidth=0.8)
a.set_yticks([-(i_ * 4 + j_) for i_ in range(len(PGL)) for j_ in range(3)]); a.set_yticklabels([f"{l}: {m}" for _, l in PGL for m in ["age only", "+ ADHD", "+ ADHD, sex, era"]], fontsize=5.2); a.axvline(0, color=MUTED, lw=0.5)
a.set_xlabel("pTDT difference, dx > 10 vs < 6 y (SD)"); a.set_title("b  Transmission to late-diagnosed", loc="left", fontweight="bold", fontsize=7)
a = axs[2]; rr = RC.assign(c=lambda d: d.wald_chi2_4df).sort_values("c")
cl = lambda s_: ORANGE if s_.startswith("age at diagnosis") else BLUE if any(k_ in s_ for k_ in ["income", "deprivation", "degree", "maternal", "paternal"]) else AQUA if any(k_ in s_ for k_ in ["walking", "motor", "DCDQ", "onset", "intellectual", "Vineland", "first words", "language"]) else GREY
a.barh(range(len(rr)), rr.c, color=[cl(s_) for s_ in rr.split]); a.set_yticks(range(len(rr)))
a.set_yticklabels([{"age at diagnosis (top vs bottom third)": "age at diagnosis (thirds)", "age at diagnosis (>10 vs <6 y; paper)": "age at diagnosis (>10 vs <6 y)"}.get(s_, s_.split(" (")[0]) for s_ in rr.split], fontsize=4.8); a.axvline(9.49, color=MUTED, lw=0.5, ls=":")
a.set_xlabel("Wald chi-square, four rare-variant classes (4 df)"); a.set_title("c  Rare variants (13,282 trio probands)", loc="left", fontweight="bold", fontsize=7)
save(fig, "L_ED3")
print("part 2 figures written")
