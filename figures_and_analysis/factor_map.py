#!/usr/bin/env python3
"""Their factors vs our factors.
Our factors: unrotated principal-axis factors of the genetic correlation matrix of 23 SPARK split-half GWAS (efa_local.py):
  G   general (shared autism liability; every half loads 0.5-0.9)
  SES bipolar socioeconomic axis (+ advantaged halves: high income, mother degree, older parents; - disadvantaged/externalizing)
  LD  bipolar late-diagnosis/design axis (+ diagnosed > 10 y, recent registration/diagnosis, specialist, no parent genotyped;
      - early developmental: regression, motor delay, high SCQ)
Their factors: the published GenomicSEM factor GWAS (Zhang, Grove et al. 2025; F1 'earlier-diagnosed', F2 'later-diagnosed').
Estimate: for orthogonal factors, a half's rg with an external GWAS = sum_j loading_ij * rho_j, where rho_j is the genetic
correlation of that GWAS with our factor j. rho is estimated by least squares across the 23 halves (no intercept), SE by
bootstrap over halves (2,000). Implied rg(F1, F2) = sum_j rho_F1j * rho_F2j (published: 0.38).
Outputs: ../results/gsem_all (DATA=all, default) or ../results/gsem (DATA=arrays): their_vs_our_factors.tsv and the figure L_FactorMap (Manuscripts/Nature_LTE/figures)."""
import os, importlib.util, io, contextlib
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); MS = os.path.join(H, "..", "..", "..", "Manuscripts", "Nature_LTE")
spec = importlib.util.spec_from_file_location("e", os.path.join(H, "efa_local.py")); e = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(e)
X, keep = e.X, e.keep
L = e.paf(X, 3)
if L[:, 0].sum() < 0: L[:, 0] *= -1
if L[keep.index("rep_s09_lo_c1"), 1] < 0: L[:, 1] *= -1      # SES: high-income half positive
if L[keep.index("rep_s00_hi_c2"), 2] < 0: L[:, 2] *= -1      # LD: diagnosed > 10 y half positive
S = pd.read_csv(os.path.join(MS, "tables", "survey_table.tsv" if e.DATA == "all" else "survey_table_arrays.tsv"), sep="\t").set_index("code")
part = lambda k: (k.split("_")[1], k.split("_")[2])
Y = {f: np.array([S.loc[c, f"{h}_rg_W_{f}"] for c, h in map(part, keep)]) for f in ["F1", "F2"]}
FN = ["Sebat F1: general", "Sebat F2: EDU/SES (advantaged +)", "Sebat F3: late diagnosis / design (late +)"]
rng = np.random.default_rng(1); rows = []; RHO = {}
for f, y in Y.items():
    rho = np.linalg.lstsq(L, y, rcond=None)[0]; RHO[f] = rho
    bs = []
    for _ in range(2000):
        i = rng.integers(0, len(y), len(y)); bs.append(np.linalg.lstsq(L[i], y[i], rcond=None)[0])
    se = np.std(bs, axis=0)
    for j in range(3): rows.append(dict(their_factor=f, our_factor=FN[j], rho=rho[j], se=se[j]))
T = pd.DataFrame(rows); T.to_csv(os.path.join(e.G, "their_vs_our_factors.tsv"), sep="\t", index=False)
print(T.pivot(index="their_factor", columns="our_factor", values="rho").reindex(columns=FN).round(2).to_string())
print(T.pivot(index="their_factor", columns="our_factor", values="se").reindex(columns=FN).round(2).to_string())
print("implied rg(F1, F2) =", round(float(RHO["F1"] @ RHO["F2"]), 2), "(published 0.38)")
d = RHO["F2"] - RHO["F1"]; ang = np.degrees(np.arctan2(d[2], d[1]))
print("their F2-F1 contrast in our SES x LD plane: SES", round(d[1], 2), "LD", round(d[2], 2), "| share of contrast (excluding general) on SES:", round(d[1] ** 2 / (d[1] ** 2 + d[2] ** 2), 2))

# ---------------- figure ----------------
BLUE, ORANGE, GREEN, PURPLE, GREY, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#6a3d9a", "#9a9a9a", "#222222", "#6b6b6b"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False})
NICE = {("s00", "hi"): "diagnosed > 10 y", ("s02", "hi"): "registered 2021+", ("s03", "hi"): "diagnosed 2019+", ("s04", "hi"): "multiplex", ("s05", "hi"): "specialist dx",
        ("s09", "hi"): "low income", ("s09", "lo"): "high income", ("s10", "hi"): "mother no degree", ("s10", "lo"): "mother degree", ("s11", "hi"): "older mothers",
        ("s11", "lo"): "younger mothers", ("s12", "hi"): "older fathers", ("s18", "hi"): "regression", ("s21", "hi"): "motor delay", ("s23", "hi"): "high SCQ",
        ("s27", "hi"): "ADHD", ("s27", "lo"): "no ADHD", ("s28", "hi"): "mood/anxiety", ("s29", "hi"): "ODD/conduct", ("s32", "hi"): "high externalizing",
        ("s34", "hi"): "sleep disorder", ("s42", "hi"): "born < 2005", ("s99", "hi"): "no parent genotyped"}
# ---- extension: place halves that were not used to fit (e.g. diagnosed < 6 y, 61% of all cases) by their rg with the 23 indicators
Pinv = np.linalg.pinv(L); ALL = {k: L[i] for i, k in enumerate(keep)}
for t in e.R.index:
    if t not in ALL: ALL[t] = Pinv @ e.R.loc[t, keep].to_numpy()
d = RHO["F2"] - RHO["F1"]
cosf = lambda a, b: float(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
SPL = [("s00", "Age at diagnosis (>10 vs <6 y)", ORANGE), ("s27", "ADHD (yes vs no)", GREY), ("s09", "Household income (low vs high)", BLUE),
       ("s10", "Mother no degree vs degree", BLUE), ("s11", "Maternal age (younger vs older)", BLUE), ("s29", "ODD/conduct (yes vs no)", GREY), ("s21", "Motor delay (yes vs no)", GREY)]
rows = []
for c, nm, col in SPL:
    v = ALL[f"rep_{c}_hi_c2"] - ALL[f"rep_{c}_lo_c1"]
    if c == "s11": v = -v                                         # orient: younger minus older mothers
    cs = cosf(v[1:], d[1:]); rows.append(dict(split=nm, cos=cs, shared_variance=cs ** 2, angle=np.degrees(np.arccos(abs(cs))), color=col))
A = pd.DataFrame(rows); A.drop(columns="color").to_csv(os.path.join(e.G, "split_vs_their_axis.tsv"), sep="\t", index=False)
print(A.drop(columns="color").round(2).to_string(index=False))
fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1], width_ratios=[1.35, 1], hspace=0.38, wspace=0.42)
a = fig.add_subplot(gs[0, :])
for t, v in ALL.items():
    c, h = part(t)
    if (c, h) not in NICE and not (c == "s00" and h == "lo"): continue
    a.scatter(v[1], v[2], s=10, color="#c8c8c8", zorder=2)
for c, nm, col in SPL:
    hi, lo = ALL[f"rep_{c}_hi_c2"], ALL[f"rep_{c}_lo_c1"]; lw = 2.4 if c == "s00" else 1.3
    a.plot([lo[1], hi[1]], [lo[2], hi[2]], color=col, lw=lw, alpha=0.95, zorder=3, solid_capstyle="round")
    for p, lab in [(hi, NICE.get((c, "hi"), "")), (lo, NICE.get((c, "lo"), "diagnosed < 6 y" if c == "s00" else ""))]:
        a.scatter(p[1], p[2], s=22, color=col, zorder=4, edgecolor="white", linewidth=0.5)
        if c in ("s11",): lab = {"older mothers": "older mothers", "younger mothers": "younger mothers"}.get(lab, lab)
        if lab: a.annotate(lab, (p[1], p[2]), xytext={"low income": (-36, -9), "high income": (4, -6)}.get(lab, (4, 3)), textcoords="offset points", fontsize=5.8, color=ORANGE if c == "s00" else INK, fontweight="bold" if c == "s00" else None)
u = d[1:] / np.linalg.norm(d[1:]) * 0.85
a.annotate("", xy=(u[0], u[1]), xytext=(-u[0], -u[1]), arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=2.6, alpha=0.9))
a.text(-0.86, 0.71,"Zhang F1 (earlier-diagnosed)\n->  Zhang F2 (later-diagnosed)", color=PURPLE, fontsize=6.6, fontweight="bold", ha="left", va="top")
a.axhline(0, color=MUTED, lw=0.5); a.axvline(0, color=MUTED, lw=0.5)
a.set_xlabel("Sebat F2, the missing EDU/SES axis:   disadvantaged / externalizing   <-->   advantaged / older parents")
a.set_ylabel("Sebat F3, late-diagnosis axis:\nearly-developmental  <-->  late dx / recent / design")
a.set_xlim(-0.9, 0.95); a.set_ylim(-0.75, 0.75); a.set_aspect("equal", adjustable="box")
a.set_title("a  Each split is a line between its two halves; Zhang F1 -> F2 axis in purple (Sebat F1, general autism, removed)", loc="left", fontweight="bold", fontsize=7.2)
b = fig.add_subplot(gs[1, 0]); A2 = A.iloc[::-1]
b.barh(range(len(A2)), A2.shared_variance, color=A2.color)
for i, r in enumerate(A2.itertuples()): b.text(r.shared_variance + 0.02, i, f"{r.angle:.0f} deg", va="center", fontsize=5.8, color=MUTED)
b.set_yticks(range(len(A2))); b.set_yticklabels(A2.split, fontsize=6); b.set_xlim(0, 1.1); b.axvline(1, color=MUTED, lw=0.5, ls=":")
b.set_xlabel("overlap with the Zhang axis (cos^2 of the angle)")
b.set_title("b  How much of the Zhang F1 -> F2 axis each split recovers", loc="left", fontweight="bold", fontsize=7.2)
c3 = fig.add_subplot(gs[1, 1]); M = np.vstack([RHO["F1"], RHO["F2"]]); SE = T.pivot(index="their_factor", columns="our_factor", values="se").reindex(columns=FN).to_numpy()
im = c3.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
for (i, j), v in np.ndenumerate(M): c3.text(j, i, f"{v:.2f}\n({SE[i, j]:.2f})", ha="center", va="center", fontsize=6.3, color="white" if abs(v) > 0.6 else INK)
c3.set_xticks(range(3)); c3.set_xticklabels(["Sebat F1\ngeneral", "Sebat F2\nEDU/SES", "Sebat F3\nlate dx"], fontsize=6.2)
c3.set_yticks(range(2)); c3.set_yticklabels(["Zhang F1\n(earlier-dx)", "Zhang F2\n(later-dx)"], fontsize=6.3)
c3.set_title("c  rg, Zhang factors x Sebat factors", loc="left", fontweight="bold", fontsize=7.2)
fig.savefig(os.path.join(MS, "figures", "L_FactorMap.png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(MS, "figures", "L_FactorMap.pdf"), bbox_inches="tight")
print("figure written")
