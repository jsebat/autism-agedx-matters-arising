#!/usr/bin/env python3
"""Simulation (review comment #4): can ONE autism liability plus a heritable stratifier of age at diagnosis reproduce Zhang et al.'s observations
(low rg between early- and late-diagnosed halves; two correlated factors preferred over one) without distinct developmental aetiologies?

Model (summary level, independent standardized SNPs):
  autism liability  L = gL + eL,  Var(gL) = h2L; stratifier S = gS + eS, Var(gS) = h2S, corr(SNP effects of gL, gS) = rho; Var(L) = Var(S) = 1.
  Cases: L > threshold (prevalence K). Age at diagnosis among cases D = a*S + sqrt(1-a^2)*noise (a = 0 -> D unrelated to S); no second liability.
  GWAS of a case group vs controls: expected per-SNP effect beta_j = (bL_j, bS_j) Sigma^-1 (mu_group - mu_controls), Sigma = Cov(L, S),
  mu = E[(L, S) | group] (Monte Carlo); observed z = beta * sqrt(Neff) + sampling noise, noise correlated between GWAS that share controls.
Design (approximate sizes; illustrative, not the original cohorts' exact counts):
  F1 indicators: PGC-2017-like (cases weighted toward early diagnosis), SPARK D < q_early, iPSYCH-like D < q_early
  F2 indicators: iPSYCH-like D > q_late, FinnGen-like (late-diagnosed cases), SPARK D > q_late
  Genetic covariances from z cross-products with the known noise covariance removed (LDSC with known intercept, no LD).
  Factor models fitted to the 6 x 6 genetic covariance by weighted least squares, weight = inverse sampling covariance of its 21 elements across
  replicates (GenomicSEM-style): 1 factor vs 2 correlated factors (Zhang assignment). Reported: rg(SPARK early, SPARK late), chi2 difference, AIC choice.
env SCEN (index into the grid, or 'all'), NREP (default 200). Output ~/rvh2/agedx/sim_twofactor[_<scen>].tsv"""
import os, itertools, json
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.stats import norm, chi2

OUT = os.path.expanduser("~/rvh2/agedx"); os.makedirs(OUT, exist_ok=True)
M, K, H2L = 20000, 0.02, 0.10
# cohort: (cases, controls, early quantile cut, late quantile cut); groups below
COH = {"SPARK": dict(nca=41695, nco=79123), "iPSYCH": dict(nca=19000, nco=40000),
       "PGC2017": dict(nca=18000, nco=28000), "FinnGen": dict(nca=2000, nco=300000)}
GWAS = ["PGC2017", "SPARK_early", "iPSYCH_early", "iPSYCH_late", "FinnGen", "SPARK_late"]   # F1 = first three, F2 = last three
GRID = [dict(h2S=h, rho=r, a=a) for h, r, a in itertools.product([0.0, 0.1, 0.3], [0.0, 0.3], [0.0, 0.3, 0.6])]


def group_means(h2S, rho, a, rng, n=2_000_000):
    """E[(L,S) | group] for every GWAS group and controls, by Monte Carlo; also case counts per group."""
    gl = rng.normal(size=n) * np.sqrt(H2L); gsr = rng.normal(size=n)
    # genetic components: corr(gL, gS) = rho (only meaningful if h2S > 0)
    gs = (rho * gl / np.sqrt(H2L) + np.sqrt(1 - rho ** 2) * gsr) * np.sqrt(h2S) if h2S > 0 else np.zeros(n)
    L = gl + rng.normal(size=n) * np.sqrt(1 - H2L); S = gs + rng.normal(size=n) * np.sqrt(1 - h2S)
    case = L > norm.ppf(1 - K); ctrl = ~case
    D = np.full(n, np.nan); D[case] = a * S[case] + np.sqrt(1 - a ** 2) * rng.normal(size=case.sum())
    mu = lambda m: np.array([L[m].mean(), S[m].mean()])
    out = {"controls": mu(ctrl)}
    for c, (qe, ql) in [("SPARK", (0.55, 0.75)), ("iPSYCH", (0.50, 0.70))]:   # early = bottom qe of diagnosis age, late = top (1 - ql)
        lo = np.nanquantile(D[case], qe); hi = np.nanquantile(D[case], ql)
        out[f"{c}_early"] = mu(case & (D <= lo)); out[f"{c}_late"] = mu(case & (D >= hi))
    # PGC-2017-like: cases weighted toward early diagnosis (sample 2:1 from below vs above the median)
    med = np.nanmedian(D[case]); w = np.where(D[case] <= med, 2.0, 1.0); out["PGC2017"] = np.array([np.average(L[case], weights=w), np.average(S[case], weights=w)])
    out["FinnGen"] = mu(case & (D >= np.nanquantile(D[case], 0.70)))
    Sig = np.cov(np.vstack([L, S]))
    return out, Sig


def sizes():
    sp, ip = COH["SPARK"], COH["iPSYCH"]
    N = {"SPARK_early": (28928, sp["nco"]), "SPARK_late": (12767, sp["nco"]), "iPSYCH_early": (11000, ip["nco"]), "iPSYCH_late": (5500, ip["nco"]),
         "PGC2017": (COH["PGC2017"]["nca"], COH["PGC2017"]["nco"]), "FinnGen": (COH["FinnGen"]["nca"], COH["FinnGen"]["nco"])}
    ctrlgrp = {"SPARK_early": "S", "SPARK_late": "S", "iPSYCH_early": "I", "iPSYCH_late": "I", "PGC2017": "P", "FinnGen": "F"}
    # sampling-noise correlation of two case-control mean differences sharing controls: (1/Nc) / sqrt((1/N1+1/Nc)(1/N2+1/Nc))
    R = np.eye(6)
    for i, j in itertools.combinations(range(6), 2):
        gi, gj = GWAS[i], GWAS[j]
        if ctrlgrp[gi] == ctrlgrp[gj]:
            (a1, c1), (a2, c2) = N[gi], N[gj]; R[i, j] = R[j, i] = (1 / c1) / np.sqrt((1 / a1 + 1 / c1) * (1 / a2 + 1 / c2))
    neff = np.array([1 / (1 / N[g][0] + 1 / N[g][1]) for g in GWAS])
    return neff, R


def simulate(sc, nrep, seed):
    rng = np.random.default_rng(seed); mus, Sig = group_means(sc["h2S"], sc["rho"], sc["a"], rng); neff, R = sizes(); Lc = np.linalg.cholesky(R)
    Sinv = np.linalg.inv(Sig); dmu = np.array([mus[g] - mus["controls"] for g in GWAS]) @ Sinv.T   # 6 x 2: loadings of each GWAS on (bL, bS)
    covs = []
    for r in range(nrep):
        bL = rng.normal(size=M) * np.sqrt(H2L / M)
        bS = ((sc["rho"] * bL / np.sqrt(H2L / M) + np.sqrt(1 - sc["rho"] ** 2) * rng.normal(size=M)) * np.sqrt(sc["h2S"] / M)) if sc["h2S"] > 0 else np.zeros(M)
        beta = dmu @ np.vstack([bL, bS])                                   # 6 x M
        z = beta * np.sqrt(neff)[:, None] + Lc @ rng.normal(size=(6, M))   # 6 x M
        C = (z @ z.T) / M - R                                              # genetic covariance on the z scale (known intercepts removed)
        C = C / np.sqrt(np.outer(neff, neff)) * M                          # per-SNP-h2 scale
        covs.append(C)
    return np.array(covs)


IU = np.triu_indices(6)
def fit_cfa(Sm, Vinv, nf):
    """WLS fit of a 1-factor or 2-correlated-factor model (factor variances fixed to 1, residual variances >= 0)."""
    s = Sm[IU]
    def implied(th):
        lam = th[:6]; psi = th[6:12] ** 2
        if nf == 1: Lm = lam[:, None]; Phi = np.eye(1)
        else:
            Lm = np.zeros((6, 2)); Lm[:3, 0] = lam[:3]; Lm[3:, 1] = lam[3:]; ph = np.tanh(th[12]); Phi = np.array([[1, ph], [ph, 1]])
        return (Lm @ Phi @ Lm.T + np.diag(psi))[IU]
    Lchol = np.linalg.cholesky(Vinv)
    res = lambda th: Lchol.T @ (s - implied(th))
    th0 = np.concatenate([np.sqrt(np.clip(np.diag(Sm), 1e-6, None)) * 0.9, np.sqrt(np.clip(np.diag(Sm), 1e-6, None)) * 0.3] + ([[0.5]] if nf == 2 else []))
    f = least_squares(res, th0, method="trf", max_nfev=5000)
    x2 = float(2 * f.cost); npar = 12 + (nf == 2); phi = float(np.tanh(f.x[12])) if nf == 2 else np.nan
    return x2, len(s) - npar, phi


def run(i, sc, nrep):
    covs = simulate(sc, nrep, seed=1000 + i)
    V = np.cov(np.array([c[IU] for c in covs]).T) + np.eye(len(IU[0])) * 1e-12; Vinv = np.linalg.inv(V)
    rows = []
    for c in covs:
        r_sp = c[1, 5] / np.sqrt(c[1, 1] * c[5, 5]) if c[1, 1] > 0 and c[5, 5] > 0 else np.nan
        x1, d1, _ = fit_cfa(c, Vinv, 1); x2, d2, phi = fit_cfa(c, Vinv, 2)
        rows.append(dict(rg_spark_early_late=r_sp, chi2_1f=x1, chi2_2f=x2, dchi2=x1 - x2, p_2f_better=float(chi2.sf(max(x1 - x2, 0), d1 - d2)), phi_2f=phi,
                         aic_prefers_2f=bool(x2 + 2 * 13 < x1 + 2 * 12)))
    R = pd.DataFrame(rows)
    return dict(scen=i, **sc, nrep=nrep, rg_spark_median=R.rg_spark_early_late.median(), rg_spark_q05=R.rg_spark_early_late.quantile(0.05),
                rg_spark_q95=R.rg_spark_early_late.quantile(0.95), frac_2f_p05=float((R.p_2f_better < 0.05).mean()), frac_aic_2f=float(R.aic_prefers_2f.mean()),
                phi_2f_median=R.phi_2f.median(), dchi2_median=R.dchi2.median())


if __name__ == "__main__":
    SC = os.environ.get("SCEN", "all"); NREP = int(os.environ.get("NREP", "200"))
    idx = range(len(GRID)) if SC == "all" else [int(SC)]
    T = pd.DataFrame([run(i, GRID[i], NREP) for i in idx])
    f = f"{OUT}/sim_twofactor{'' if SC == 'all' else '_' + SC}.tsv"; T.to_csv(f, sep="\t", index=False)
    pd.set_option("display.width", 200); print(T.round(3).to_string(index=False))
