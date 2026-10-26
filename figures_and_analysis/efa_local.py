#!/usr/bin/env python3
"""Our own genomic factor analysis of SPARK autism, from the GenomicSEM multivariable-LDSC genetic correlation matrix of the
broad-survey split-half GWAS (their fastGWA pipeline; results/gsem/rg.tsv, smoothed to the nearest PD matrix in R).
Indicators: halves with h2 z >= 3 (gsem_fa.R input) EXCLUDING near-total halves (> 60% of all genotyped cases, e.g. 'no seizures'),
which are essentially the all-case GWAS and would dominate the first factor.
EFA: principal-axis factoring (iterated communalities), 1-3 factors, promax rotation (kappa 4).
Interpretation: loadings vs the composition of each half's cases (late diagnosis, adult registration, female, high deprivation, ADHD,
ID, born before 2005, no parent genotyped) and vs each half's resemblance to the published factors (rg F2 - rg F1).
Outputs: results/gsem/efa_k{2,3}_loadings.tsv, efa_composition.tsv, efa_eigen.tsv
DATA=all (default): the same 23 indicators (chosen on the arrays, so the indicator set is fixed), with the genetic correlation matrix and published-factor rg from the
all-strata split GWAS (results/gsem_all, from gsem_fa_all.sbatch; results/full_data/all_survey_rg_all.tsv). Case composition columns remain the arrays' (descriptive only).
DATA=arrays: the arrays-only matrix (results/gsem)."""
import os
import numpy as np, pandas as pd
DATA = os.environ.get("DATA", "all")
H = os.path.dirname(os.path.abspath(__file__)); G = os.path.join(H, "..", "results", "gsem_all" if DATA == "all" else "gsem"); RP = os.path.join(H, "..", "results", "replication")
R = pd.read_csv(os.path.join(G, "rg.tsv"), sep="\t").set_index("trait"); R = R.loc[:, R.index]
comp = pd.read_csv(os.path.join(RP, "survey_composition.tsv"), sep="\t").set_index("gwas")
svr = pd.read_csv(os.path.join(H, "..", "results", "full_data", "all_survey_rg_all.tsv") if DATA == "all" else os.path.join(RP, "survey_rg.tsv"), sep="\t").set_index("code")
sv = pd.read_csv(os.path.join(H, "..", "results", "survey_splits.tsv"), sep="\t").set_index("code")
ncase_all = comp.loc["rep_s00_hi_c2", "n"] + comp.loc["rep_s00_lo_c1", "n"]
tot = 37000  # approx. genotyped EUR cases across GSA v1 / v2 / Twist in the replication (W_all: 7,054 + 16,298 + 14,124 = 37,476)
tot = 37476
keep = [t for t in R.index if comp.loc[t, "n"] <= 0.6 * tot]
print(f"indicators: {len(R)} -> {len(keep)} after dropping near-total halves:", sorted(set(R.index) - set(keep)))
X = R.loc[keep, keep].to_numpy()
ev = np.sort(np.linalg.eigvalsh(X))[::-1]; pd.DataFrame({"component": range(1, len(ev) + 1), "eigen": ev}).to_csv(os.path.join(G, "efa_eigen.tsv"), sep="\t", index=False)
print("eigenvalues:", np.round(ev[:6], 2), " Kaiser (>1):", int((ev > 1).sum()))

def paf(Rm, k, it=500):
    h = 1 - 1 / np.diag(np.linalg.pinv(Rm)); h = np.clip(h, 0.2, 0.99)
    for _ in range(it):
        Rr = Rm.copy(); np.fill_diagonal(Rr, h); w, v = np.linalg.eigh(Rr); idx = np.argsort(w)[::-1][:k]
        L = v[:, idx] * np.sqrt(np.clip(w[idx], 0, None)); hn = np.clip((L ** 2).sum(1), 0.05, 0.999)
        if np.max(np.abs(hn - h)) < 1e-6: break
        h = hn
    return L
def varimax(L, it=200):
    p, k = L.shape; Rt = np.eye(k); d = 0
    for _ in range(it):
        Lr = L @ Rt; u, s, vt = np.linalg.svd(L.T @ (Lr ** 3 - Lr @ np.diag((Lr ** 2).sum(0)) / p)); Rt = u @ vt
        if s.sum() < d * (1 + 1e-8): break
        d = s.sum()
    return L @ Rt
def promax(L, kappa=4):
    V = varimax(L); P = np.abs(V) ** (kappa - 1) * V
    U = np.linalg.lstsq(V, P, rcond=None)[0]; U = U @ np.diag(np.sqrt(np.diag(np.linalg.inv(U.T @ U))))
    Lp = V @ U; phi = np.linalg.inv(U.T @ U); return Lp, phi
code = lambda t: t.split("_")[1]; half = lambda t: t.split("_")[2]
lab = lambda t: f"{sv.loc[code(t), 'split'] if code(t) in sv.index else 'design: no parent vs both parents genotyped'} [{'high/yes' if half(t) == 'hi' else 'low/no'}]"
rows = []
for k in (2, 3):
    L = paf(X, k); Lp, phi = promax(L) if k > 1 else (L, np.eye(1))
    D = pd.DataFrame(Lp, index=keep, columns=[f"F{i+1}" for i in range(k)])
    D["label"] = [lab(t) for t in keep]
    D["F2minusF1_published"] = [svr.loc[code(t), f"{half(t)}_rg_W_F2"] - svr.loc[code(t), f"{half(t)}_rg_W_F1"] if code(t) in svr.index else np.nan for t in keep]
    D = D.join(comp[["n", "late_dx", "adult_reg", "female", "high_adi", "adhd", "id", "born_pre2005", "no_parent"]])
    D.to_csv(os.path.join(G, f"efa_k{k}_loadings.tsv"), sep="\t")
    print(f"\n=== EFA {k} factors (promax); factor correlations:\n", np.round(phi, 2))
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    print(D.sort_values("F1")[[f"F{i+1}" for i in range(k)] + ["label", "n", "late_dx", "F2minusF1_published"]].round(2).to_string())
    for f in [f"F{i+1}" for i in range(k)]:
        for v in ["late_dx", "adult_reg", "female", "high_adi", "adhd", "id", "born_pre2005", "no_parent", "F2minusF1_published"]:
            ok = D[[f, v]].dropna()
            rows.append(dict(k=k, factor=f, variable=v, n_ind=len(ok), r=np.corrcoef(ok[f], ok[v])[0, 1]))
C = pd.DataFrame(rows); C.to_csv(os.path.join(G, "efa_composition.tsv"), sep="\t", index=False)
print("\nloading vs composition (Pearson r across indicators)"); print(C.pivot_table(index=["k", "factor"], columns="variable", values="r").round(2).to_string())
