#!/usr/bin/env python3
"""Main Figure 2. Order: (a) common variants, (b) rare variants, (c) rare-variant stratification ranked separately for inherited and de novo
classes, (d) common and rare variants on the same scale.
(a) Each split as an arrow in external-genetics space (x psychiatric: mean rg with ADHD and MDD; y EDU/SES: mean rg with age at first birth and EA);
    Zhang F2 - F1 for comparison. Input results/axis_vectors.tsv.
(b) Each split as an arrow in rare-variant space (probands, high/yes minus low/no): x de novo burden (De novo - NDD genes + De novo - other genes),
    y inherited burden (Inherited - NDD genes + Inherited - other genes). Input results/full_data/rare_compare_full.tsv (rare_full.py).
(c) Wald chi-square (2 df) for the two inherited classes and, separately, the two de novo classes. Same input.
(d) Out-of-sample canonical correlation of each axis linking 35 characteristics of the same 7,072 probands to polygenic scores alone, rare variants alone,
    and both; axes matched to the combined analysis by the similarity of their characteristic loadings; age-at-diagnosis loading (combined) under each.
    Inputs: results/full_data/cca_trio_{all,pgs,rare}_{cv,structure}.tsv (cca2.py, ANALYSIS=trio).
Output: Manuscripts/Nature_LTE/figures/L_Fig2main.(png|pdf)"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results"); MS = os.path.join(H, "..", "..", "..", "Manuscripts", "Nature_LTE")
BLUE, ORANGE, PURPLE, GREY, DGREY, INK, MUTED, GRID, GREEN = "#2a78d6", "#eb6834", "#6a3d9a", "#c4c4c4", "#7a7a7a", "#222222", "#6b6b6b", "#dddddd", "#1a9850"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(7.2, 7.4)); gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1], width_ratios=[1, 1], hspace=0.55, wspace=0.42)


# ---- directions of the Sebat factors (dashed lines through the origin) ----
TEAL = "#00798c"
def sebat_factor_directions():
    """Panel a: external-trait profile of the GWAS-based Sebat factors (unrotated principal-axis factors of 23 split-half GWAS), estimated by least squares
    of each half's rg with ADHD, MDD, age at first birth and EA on its factor loadings. Panel b: genetic structure of the rare-variant-only CCA axes."""
    import importlib.util, io, contextlib
    spec = importlib.util.spec_from_file_location("e", os.path.join(H, "efa_local.py")); e = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(e)
    L = e.paf(e.X, 3); keep = e.keep
    if L[:, 0].sum() < 0: L[:, 0] *= -1
    if L[keep.index("rep_s09_lo_c1"), 1] > 0: L[:, 1] *= -1      # orient EDU/SES factor to its disadvantaged pole
    if L[keep.index("rep_s00_hi_c2"), 2] < 0: L[:, 2] *= -1      # late-diagnosis factor to its late pole
    E = pd.read_csv(os.path.join(R, "full_data", "all_rg_ext_all.tsv"), sep="\t").pivot_table(index="gwas", columns="trait", values="rg")
    tr = ["adhd_2022", "mdd_2025", "afb_2018", "edu_2022"]; ok = [i for i, k in enumerate(keep) if k in E.index and E.loc[k, tr].notna().all()]
    prof = {x: np.linalg.lstsq(L[ok], E.loc[[keep[i] for i in ok], x].to_numpy(), rcond=None)[0] for x in tr}
    pa = {j: ((prof["adhd_2022"][j] + prof["mdd_2025"][j]) / 2, (prof["afb_2018"][j] + prof["edu_2022"][j]) / 2) for j in (1, 2)}
    s = pd.read_csv(os.path.join(R, "full_data", "cca_trio_rare_structure.tsv"), sep="\t"); s = s[s.block == "genetics"]
    pb = {}
    for k in (1, 2):
        x = s[s.component == k].set_index("variable").r; v = (x["de novo, NDD genes"] + x["de novo, other genes"], x["inherited, NDD genes"] + x["inherited, other genes"])
        if k == 1 and v[0] < 0: v = (-v[0], -v[1])
        if k == 2 and v[1] < 0: v = (-v[0], -v[1])
        pb[k] = v
    return pa, pb
PA, PB = sebat_factor_directions()
BOX = dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8)
def axisline(ax, v, length, label, ha="left", dy=0.0, dx=0.0, back=0.15, color=TEAL, lab_at=1.0):
    u = np.array(v) / np.hypot(*v)
    ax.plot([-u[0] * length * back, u[0] * length], [-u[1] * length * back, u[1] * length], ls=(0, (4, 3)), color=color, lw=0.9, alpha=0.55, zorder=1)
    ax.text(u[0] * length * lab_at + dx, u[1] * length * lab_at + dy, label, color=color, fontsize=6.4, fontweight="bold", ha=ha, va="center", zorder=6, bbox=BOX)
    ex, ey = u * length; (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim(); my = 0.06 * (y1 - y0)   # keep the line end and its label inside the axes
    ax.set_ylim(min(y0, ey - my), max(y1, ey + my)); ax.set_xlim(min(x0, ex - 0.02), max(x1, ex + 0.02))

# (a) common variants
V = pd.read_csv(os.path.join(R, "axis_vectors.tsv"), sep="\t")
LAB = {"s00": ("age at diagnosis split", ORANGE), "s09": ("low income", BLUE), "s10": ("mother no degree", BLUE),
       "s11": ("older mothers", BLUE), "s08": ("high ADI", BLUE), "s27": ("ADHD", DGREY)}
OFF = {"s00": (4, 3), "s27": (5, 2)}
LEADA = {"s11": ((0.85, -0.99), "right"), "s09": ((0.68, -1.09), "left"), "s10": ((0.50, -1.20), "left"), "s08": ((0.36, -0.97), "right")}
a = fig.add_subplot(gs[0, 0])
for r in V.itertuples():
    if r.code in LAB or r.code == "F2-F1": continue
    a.annotate("", xy=(r.d_psych, r.d_repro), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.6, mutation_scale=5))
for r in V.itertuples():
    if r.code not in LAB: continue
    lab, col = LAB[r.code]; a.annotate("", xy=(r.d_psych, r.d_repro), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2 if r.code == "s00" else 1.4, mutation_scale=8))
    if r.code in LEADA:   # crowded EDU/SES tips: label in free space (data coordinates), thin leader line to the arrow tip
        (lx, ly), ha_ = LEADA[r.code]
        a.annotate(lab, (r.d_psych, r.d_repro), xytext=(lx, ly), textcoords="data", ha=ha_, va="center", fontsize=6, color=col, zorder=6, bbox=BOX,
                   arrowprops=dict(arrowstyle="-", color=col, lw=0.4, shrinkA=0, shrinkB=2))
    else:
        a.annotate(lab, (r.d_psych, r.d_repro), xytext=OFF[r.code], textcoords="offset points", fontsize=6, color=col, fontweight="bold" if r.code == "s00" else None, zorder=6, bbox=BOX)
a.axhline(0, color=MUTED, lw=0.5); a.axvline(0, color=MUTED, lw=0.5); a.set_xlim(-0.12, 0.85); a.set_ylim(-1.25, 0.45)
zf = V[V.code == "F2-F1"].iloc[0]
LEN_AADX = float(np.hypot(zf.d_psych, zf.d_repro))
LEN_SES_A = float(np.mean([np.hypot(*V[V.code == c][["d_psych", "d_repro"]].iloc[0]) for c in ["s08", "s09", "s10", "s11", "s12"]]))
axisline(a, (zf.d_psych, zf.d_repro), LEN_AADX, "AADx", ha="left", dx=0.01, dy=0.0, color=INK, back=0.0)
axisline(a, PA[1], LEN_SES_A, "EDU/SES", ha="left", dx=0.01, dy=0.0, back=0.0)
# PSYCH: third axis of the combined CCA, drawn from its polygenic-score loadings on the same four traits (x: ADHD, MDD; y: age at first birth, EA);
# length = the ADHD split arrow, its closest split. LANG (axis 4) loads mainly on SCZ and height, outside this plane, and is not drawn.
gP = pd.read_csv(os.path.join(R, "full_data", "cca_trio_all_structure.tsv"), sep="\t"); gP = gP[(gP.block == "genetics") & (gP.component == 4)].set_index("variable").r
vP = ((gP["adhd_2022"] + gP["mdd_2025"]) / 2, (gP["afb_2018"] + gP["edu_2022"]) / 2)
if vP[0] < 0: vP = (-vP[0], -vP[1])
axisline(a, vP, float(np.hypot(*V[V.code == "s27"][["d_psych", "d_repro"]].iloc[0])), "PSYCH", ha="left", dx=0.01, dy=-0.06, back=0.0, color="#8c510a")
a.set_xlabel("psychiatric signal\n(difference in rg with ADHD and MDD)"); a.set_ylabel("EDU/SES signal (difference in rg with\nage at first birth and educational attainment)")
a.set_title("a  EDU/SES is a major polygenic axis;\n    AADx is secondary", loc="left", fontweight="bold", fontsize=7.2)
a.grid(color=GRID, lw=0.4); a.set_axisbelow(True)

# (b) rare variants
C = pd.read_csv(os.path.join(R, "full_data", "rare_compare_full.tsv"), sep="\t"); CP = C[C.outcome == "A probands, absolute"].set_index("split")
CP["dnm"] = CP["d_de novo, NDD genes"] + CP["d_de novo, other genes"]; CP["inh"] = CP["d_inherited, NDD genes"] + CP["d_inherited, other genes"]
DEV = ["late walking (top vs bottom third)", "motor delay", "motor coordination problems (DCDQ)", "early onset concern", "intellectual disability",
       "low Vineland ABC (bottom vs top third)", "late first words (top vs bottom third)"]
SES = ["household income (low vs high third)", "mother no degree vs degree", "deprivation (high vs low third)", "maternal age (older vs younger third)", "paternal age (older vs younger third)"]
AGE = ["age at diagnosis (top vs bottom third)"]; OTH = ["ADHD (yes vs no)", "SCQ (top vs bottom third)", "RBS-R (top vs bottom third)"]
SHORT = {"late walking (top vs bottom third)": "late walking", "motor delay": "motor delay", "motor coordination problems (DCDQ)": "DCDQ", "early onset concern": "early onset",
         "intellectual disability": "intellectual disability", "low Vineland ABC (bottom vs top third)": "low Vineland", "late first words (top vs bottom third)": "late first words",
         "household income (low vs high third)": "low income", "mother no degree vs degree": "mother no degree", "deprivation (high vs low third)": "high ADI",
         "maternal age (older vs younger third)": "older mothers", "paternal age (older vs younger third)": "older fathers", "age at diagnosis (top vs bottom third)": "late diagnosis",
         "age at diagnosis (>10 vs <6 y; paper)": "diagnosed >10 vs <6 y", "ADHD (yes vs no)": "ADHD", "SCQ (top vs bottom third)": "high SCQ", "RBS-R (top vs bottom third)": "high RBS-R"}
OFFB = {"motor delay": (-12, 6), "late walking (top vs bottom third)": (-28, -10), "high SCQ": (3, 4), "SCQ (top vs bottom third)": (3, 4), "RBS-R (top vs bottom third)": (-38, 4),
        "intellectual disability": (3, 3), "low Vineland ABC (bottom vs top third)": (3, -7), "age at diagnosis (top vs bottom third)": (-58, -14)}
LEAD = {"intellectual disability": ((0.40, 0.37), "center"), "motor coordination problems (DCDQ)": ((0.64, 0.30), "center"),
        "late walking (top vs bottom third)": ((0.92, 0.30), "center"), "motor delay": ((1.12, 0.24), "center"),
        "low Vineland ABC (bottom vs top third)": ((0.13, 0.30), "center"),
        "late first words (top vs bottom third)": ((0.52, -0.10), "left"), "early onset concern": ((0.72, -0.025), "left"),
        "age at diagnosis (top vs bottom third)": ((-0.13, 0.035), "right"), "ADHD (yes vs no)": ((-0.13, -0.04), "right"),
        "household income (low vs high third)": ((-0.52, 0.45), "left"), "deprivation (high vs low third)": ((-0.24, 0.40), "left"),
        "mother no degree vs degree": ((-0.52, 0.20), "left"), "RBS-R (top vs bottom third)": ((-0.38, 0.10), "center"),
        "maternal age (older vs younger third)": ((0.48, -0.25), "left"), "paternal age (older vs younger third)": ((0.62, -0.13), "left")}
b = fig.add_subplot(gs[0, 1])
for grp, col, lw in [(OTH, DGREY, 1.0), (DEV, GREEN, 1.5), (SES, BLUE, 1.5), (AGE, ORANGE, 2.4)]:
    for s in grp:
        x, y = CP.loc[s, "dnm"], CP.loc[s, "inh"]
        b.annotate("", xy=(x, y), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, mutation_scale=8))
        if s in LEAD:   # crowded developmental tips: label placed in free space, thin leader line to the arrow tip
            (lx, ly), ha_ = LEAD[s]
            b.plot([x, lx], [y, ly], color=col, lw=0.4, alpha=0.6, zorder=2)
            b.text(lx, ly, SHORT[s], fontsize=5.6, color=col, ha=ha_, va="center", zorder=6, bbox=BOX, fontweight="bold" if s in AGE else None)
        else:
            b.annotate(SHORT[s], (x, y), xytext=OFFB.get(s, (3, 2)), textcoords="offset points", fontsize=5.6, color=col, fontweight="bold" if s in AGE else None, zorder=6, bbox=BOX)
allx = CP.loc[DEV + SES + AGE + OTH, "dnm"]; ally = CP.loc[DEV + SES + AGE + OTH, "inh"]
b.set_xlim(min(allx.min(), 0) - 0.12, max(allx.max(), 0) + 0.25); b.set_ylim(min(ally.min(), 0) - 0.06, max(ally.max(), 0) + 0.06)
b.axhline(0, color=MUTED, lw=0.5); b.axvline(0, color=MUTED, lw=0.5)
LEN_DEV_B = float(np.mean([np.hypot(CP.loc[s, "dnm"], CP.loc[s, "inh"]) for s in DEV])); LEN_SES_B = float(np.mean([np.hypot(CP.loc[s, "dnm"], CP.loc[s, "inh"]) for s in SES[:3]]))
axisline(b, PB[1], LEN_DEV_B, "DEV", ha="left", dx=0.0, dy=0.0, back=0.0, lab_at=1.95)
axisline(b, PB[2], LEN_SES_B, "EDU/SES", ha="center", dx=0.0, dy=0.025, back=0.0)
b.set_xlabel("de novo burden difference (s.d.)\n(De novo - NDD genes + De novo - other genes)"); b.set_ylabel("inherited burden difference (s.d.)\n(Inherited - NDD + Inherited - other genes)")
b.set_title("b  Rare variants support 2 orthogonal\n    EDU/SES and DEV axes", loc="left", fontweight="bold", fontsize=7.2)
for col, lab in [(GREEN, "developmental / autistic traits"), (BLUE, "EDU/SES"), (ORANGE, "age at diagnosis"), (DGREY, "other")]:
    b.plot([], [], color=col, lw=2.5, label=lab)
b.legend(frameon=False, fontsize=5.6, loc="lower left", title="splits (panels a-c)", title_fontsize=5.8)
b.grid(color=GRID, lw=0.4); b.set_axisbelow(True)

# (c) rankings: inherited vs de novo
sub = gs[1, 0].subgridspec(1, 2, wspace=0.85)
ALL = DEV + SES + AGE + OTH + ["age at diagnosis (>10 vs <6 y; paper)"]
colof = lambda s: GREEN if s in DEV else BLUE if s in SES else ORANGE if s.startswith("age at") else DGREY
for k, (col_, ttl) in enumerate([("wald_inherited_2df", "Inherited"), ("wald_denovo_2df", "De novo")]):
    c = fig.add_subplot(sub[k]); order = sorted(ALL, key=lambda s: CP.loc[s, col_])
    c.barh(range(len(order)), CP.loc[order, col_], color=[colof(s) for s in order]); c.axvline(5.99, color=MUTED, lw=0.5, ls=":")
    c.set_yticks(range(len(order))); c.set_yticklabels([SHORT[s] for s in order], fontsize=5.4)
    for j, s in enumerate(order):
        if s.startswith("age at"): c.get_yticklabels()[j].set_color(ORANGE); c.get_yticklabels()[j].set_fontweight("bold")
    c.set_xlabel("Wald chi-square (2 df)", fontsize=6); c.set_title(ttl + " (NDD + other genes)", fontsize=6.6, fontweight="bold")
    c.grid(axis="x", color=GRID, lw=0.4); c.set_axisbelow(True)
    if k == 0: c.text(-0.75, 1.09, "c  EDU/SES and DEV axes are driven by inherited\n    and de novo rare variants, respectively", transform=c.transAxes, fontweight="bold", fontsize=7.2)

# (d) same scale: canonical correlation, PGS alone vs rare alone vs both (same 7,072 probands); (e) age-at-diagnosis loading on each axis
sub2 = gs[1, 1].subgridspec(2, 1, height_ratios=[1.3, 1.2], hspace=0.7)
d = fig.add_subplot(sub2[0]); e2 = fig.add_subplot(sub2[1])
def load(tag):
    cv = pd.read_csv(os.path.join(R, "full_data", f"cca_trio_{tag}_cv.tsv"), sep="\t"); st = pd.read_csv(os.path.join(R, "full_data", f"cca_trio_{tag}_structure.tsv"), sep="\t")
    st = st[(st.block == "characteristics") & ~st.variable.str.contains(r"\[missing\]")]
    return cv, st.pivot(index="variable", columns="component", values="r")
# expanded trio probands registered before 22 (Zhang et al. criterion); platform, sequencing cohort and family structure regressed out;
# characteristics exclude "genetic diagnosis" (reported from the same sequencing: circular for rare variants)
cvB, stB = load("all"); cvP, stP = load("pgs"); cvR, stR = load("rare")
AX = ["EDU/SES", "DEV", "PSYCH", "LANG"]   # display order
COMP = {"EDU/SES": 1, "DEV": 2, "LANG": 3, "PSYCH": 4}   # combined-analysis component of each axis
def matched(cv, st, nsig):
    out, lo = {}, {}
    for comp in range(1, min(nsig, st.shape[1]) + 1):
        cors = [abs(np.corrcoef(st[comp], stB.loc[st.index, COMP[a_]])[0, 1]) for a_ in AX]
        j = int(np.argmax(cors)); r = cv.cv_canonical_r.iloc[comp - 1]
        if AX[j] not in out or r > out[AX[j]]: out[AX[j]] = r; lo[AX[j]] = abs(st.loc["age at diagnosis", comp])
    return out, lo
nsig = lambda cv: int(((cv.p_perm <= 0.01) & (cv.cv_canonical_r > cv.null_95)).sum())
(mP, lP), (mR, lR) = matched(cvP, stP, nsig(cvP)), matched(cvR, stR, nsig(cvR))
mB = {a_: cvB.cv_canonical_r.iloc[COMP[a_] - 1] for a_ in AX}; lB = {a_: abs(stB.loc["age at diagnosis", COMP[a_]]) for a_ in AX}
xs = np.arange(4); w = 0.27
SER = [("#c51b7d", "rare variants alone"), ("#5e3c99", "polygenic scores alone"), ("#252525", "rare + polygenic")]
for k, (m, (col, lab)) in enumerate(zip([mR, mP, mB], SER)):
    vals = [m.get(a_, np.nan) for a_ in AX]; d.bar(xs + (k - 1) * w, np.nan_to_num(vals), width=w, color=col, label=lab)
    for x, v in zip(xs + (k - 1) * w, vals):
        if np.isfinite(v): d.text(x, v + 0.008, f"{v:.2f}", ha="center", fontsize=4.8, rotation=90, va="bottom")
d.set_xticks(xs); d.set_xticklabels(AX, fontsize=5.8); d.grid(axis="y", color=GRID, lw=0.4); d.set_axisbelow(True)
d.set_ylabel("how strongly genetics\npredicts the axis\n(out-of-sample r)", fontsize=6); d.set_ylim(0, 0.5)
d.legend(frameon=False, fontsize=5.4, loc="upper right")
d.set_title(f"d  Strength of each genetic axis\n    (same {int(cvB.n.iloc[0]):,} trio probands, {len(stB)} characteristics)", loc="left", fontweight="bold", fontsize=7.2)
# (e) what defines each axis (combined analysis): its top-loading characteristic vs age at diagnosis, |loading| and rank among all characteristics
w2 = 0.36; TOPS = []
SHORTX = {"deprivation (ADI)": "ADI", "CBCL externalizing": "externalizing", "minimal language": "language"}
for j in range(4):
    col_ = stB[COMP[AX[j]]].abs().sort_values(ascending=False); ad = col_["age at diagnosis"]; rk = list(col_.index).index("age at diagnosis") + 1
    oth = col_.drop("age at diagnosis"); top = oth.index[0]   # top-loading characteristic other than age at diagnosis
    e2.bar(j - w2 / 2, oth.iloc[0], width=w2, color=DGREY); e2.bar(j + w2 / 2, ad, width=w2, color=ORANGE)
    e2.text(j - w2 / 2, oth.iloc[0] / 2, SHORTX.get(top, top), ha="center", va="center", fontsize=4.8, rotation=90, color="white", fontweight="bold")
    e2.text(j + w2 / 2, ad + 0.02, f"{rk}/{len(col_)}", ha="center", va="bottom", fontsize=5.0, color=ORANGE, fontweight="bold")
    TOPS.append((AX[j], top, round(oth.iloc[0], 2), round(ad, 2), rk))
from matplotlib.patches import Patch
EH = [Patch(color=DGREY, label="top-loading other characteristic"), Patch(color=ORANGE, label="age at diagnosis loading on factor")]
e2.set_xticks(xs); e2.set_xticklabels(AX, fontsize=5.8); e2.grid(axis="y", color=GRID, lw=0.4); e2.set_axisbelow(True)
e2.set_ylabel("loading on the axis\n(|r| with the\ncharacteristics side)", fontsize=6); e2.set_ylim(0, 0.95)
e2.legend(handles=EH, frameon=False, fontsize=5.2, loc="upper right")
e2.set_title("e  What defines each axis (rare + polygenic)", loc="left", fontweight="bold", fontsize=7.2)
print("panel e:", TOPS)
fig.savefig(os.path.join(MS, "figures", "L_Fig2main.png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(MS, "figures", "L_Fig2main.pdf"), bbox_inches="tight")
print("written | PGS:", {k: round(v, 3) for k, v in mP.items()}, "| rare:", {k: round(v, 3) for k, v in mR.items()}, "| both:", {k: round(v, 3) for k, v in mB.items()})
print(CP.loc[ALL, ["wald_inherited_2df", "p_inherited", "wald_denovo_2df", "p_denovo"]].round(4).sort_values("wald_inherited_2df", ascending=False).to_string())
