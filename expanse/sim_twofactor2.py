#!/usr/bin/env python3
"""Extended simulation (review comment #4; follows sim_twofactor.py). ONE autism liability, no second aetiology, plus:
  * several heritable stratifiers S_1..S_K (heritabilities h2S_k, genetic correlation rho_k with autism, independent of each other);
  * age at diagnosis among cases D = sum_k a_k S_k + noise (a = 1 for one stratifier = splitting directly on the stratifier);
  * stratifier-dependent ascertainment: case (and optionally control) inclusion weight exp(gamma * S_asc) in chosen cohorts, optionally
    only for late-diagnosed cases (models self-registering adults entering a family study differently from children).
Same six-GWAS design, noise model and factor fits as sim_twofactor.py: F1 = PGC2017-like, SPARK early, iPSYCH early; F2 = iPSYCH late,
FinnGen-like late, SPARK late. Reported: rg(SPARK early, SPARK late) and how often two correlated factors beat one.
env SCEN (index or 'all'), NREP (default 200). Output ~/rvh2/agedx/sim_twofactor2[_<scen>].tsv"""
import os, itertools
import numpy as np, pandas as pd
from scipy.stats import norm, chi2
import sim_twofactor as S1   # sizes(), fit_cfa(), GWAS order, constants

OUT = os.path.expanduser("~/rvh2/agedx"); M, K_PREV, H2L = S1.M, S1.K, S1.H2L
GWAS = S1.GWAS


def scen(name, h2S, rho, a, asc=None):
    """h2S, rho, a: lists over stratifiers. asc: dict(var=index of stratifier, gamma=, cohorts=[...], late_only=bool, controls=bool)."""
    return dict(name=name, h2S=list(h2S), rho=list(rho), a=list(a), asc=asc)


GRID = [
    # A: split directly on one stratifier (a = 1)
    scen("A direct, h2S 0.3", [0.3], [0.0], [1.0]), scen("A direct, h2S 0.5", [0.5], [0.0], [1.0]),
    scen("A direct, h2S 0.3, rho 0.3", [0.3], [0.3], [1.0]), scen("A direct, h2S 0.5, rho 0.3", [0.5], [0.3], [1.0]),
    # B: several stratifiers jointly drive age at diagnosis
    scen("B three stratifiers, h2S 0.3 each", [0.3] * 3, [0.0] * 3, [0.55] * 3), scen("B three stratifiers, h2S 0.5 each", [0.5] * 3, [0.0] * 3, [0.55] * 3),
    # C: cases ascertained on the age-at-diagnosis stratifier in SPARK / PGC-like cohorts
    scen("C asc on S, cases, gamma 0.5", [0.3], [0.0], [0.6], dict(var=0, gamma=0.5, cohorts=["SPARK", "PGC2017"], late_only=False, controls=False)),
    scen("C asc on S, cases, gamma 1.0", [0.3], [0.0], [0.6], dict(var=0, gamma=1.0, cohorts=["SPARK", "PGC2017"], late_only=False, controls=False)),
    scen("C asc on S, cases + controls, gamma 1.0", [0.3], [0.0], [0.6], dict(var=0, gamma=1.0, cohorts=["SPARK", "PGC2017"], late_only=False, controls=True)),
    # D: only late-diagnosed SPARK cases ascertained, on a second heritable stratifier (adult self-registration)
    scen("D late-only asc on S2, gamma 0.5", [0.3, 0.3], [0.0, 0.0], [0.6, 0.0], dict(var=1, gamma=0.5, cohorts=["SPARK"], late_only=True, controls=False)),
    scen("D late-only asc on S2, gamma 1.0", [0.3, 0.3], [0.0, 0.0], [0.6, 0.0], dict(var=1, gamma=1.0, cohorts=["SPARK"], late_only=True, controls=False)),
    scen("D late-only asc on S2, gamma 1.0, a 0", [0.3, 0.3], [0.0, 0.0], [0.0, 0.0], dict(var=1, gamma=1.0, cohorts=["SPARK"], late_only=True, controls=False)),
    # combined: direct split + three stratifiers + late-only ascertainment
    scen("E combined: 3 stratifiers h2 0.5, late-only asc gamma 1.0", [0.5] * 3, [0.0] * 3, [0.55] * 3, dict(var=2, gamma=1.0, cohorts=["SPARK"], late_only=True, controls=False)),
]


def gcov(sc):
    k = len(sc["h2S"]); G = np.zeros((k + 1, k + 1)); G[0, 0] = H2L
    for i, (h, r) in enumerate(zip(sc["h2S"], sc["rho"]), start=1):
        G[i, i] = h; G[0, i] = G[i, 0] = r * np.sqrt(H2L * h)
    return G


def group_deltas(sc, rng, n=2_000_000):
    """Delta of E[(L, S_1..S_K)] between each GWAS case group and its controls (weighted for ascertainment), and Cov(L, S) in the population."""
    G = gcov(sc); k = G.shape[0]; g = rng.multivariate_normal(np.zeros(k), G, size=n)
    E = rng.normal(size=(n, k)) * np.sqrt(np.clip(1 - np.diag(G), 0, None)); P = g + E   # columns: L, S_1..S_K (unit variance)
    L = P[:, 0]; Sx = P[:, 1:]; case = L > norm.ppf(1 - K_PREV)
    a = np.array(sc["a"]); s2 = float((a ** 2).sum())
    D = Sx @ a + (np.sqrt(1 - s2) if s2 < 1 else 0.0) * rng.normal(size=n)
    if s2 == 0: D = rng.normal(size=n)
    Dc = D[case]; asc = sc["asc"]
    q = lambda p: np.quantile(Dc, p)
    groups = {"SPARK_early": ("SPARK", D <= q(0.55)), "SPARK_late": ("SPARK", D >= q(0.75)), "iPSYCH_early": ("iPSYCH", D <= q(0.50)),
              "iPSYCH_late": ("iPSYCH", D >= q(0.70)), "FinnGen": ("FinnGen", D >= q(0.70)), "PGC2017": ("PGC2017", np.ones(n, bool))}
    def w_for(cohort, sel_case, is_late):
        w = np.ones(n)
        if asc and cohort in asc["cohorts"]:
            ww = np.exp(asc["gamma"] * Sx[:, asc["var"]])
            if not asc["late_only"] or is_late: w = ww
        if cohort == "PGC2017": w = w * np.where(D <= np.median(Dc), 2.0, 1.0)   # PGC-like: cases weighted 2:1 toward early diagnosis
        return w
    out = {}
    for gname, (coh, sel) in groups.items():
        cm = case & sel; wc = w_for(coh, cm, gname.endswith("_late") or gname == "FinnGen")[cm]
        mu_c = (P[cm] * wc[:, None]).sum(0) / wc.sum()
        ctrl = ~case; wk = np.exp(asc["gamma"] * Sx[ctrl, asc["var"]]) if (asc and asc["controls"] and coh in asc["cohorts"]) else np.ones(ctrl.sum())
        mu_k = (P[ctrl] * wk[:, None]).sum(0) / wk.sum()
        out[gname] = mu_c - mu_k
    return np.array([out[gname] for gname in GWAS]), np.cov(P.T), G


def simulate(sc, nrep, seed):
    rng = np.random.default_rng(seed); dmu, Sig, G = group_deltas(sc, rng); neff, R = S1.sizes(); Lc = np.linalg.cholesky(R)
    W = dmu @ np.linalg.inv(Sig)                    # 6 x (K+1): each GWAS's loading on the SNP-effect vectors of (L, S_1..S_K)
    Gc = np.linalg.cholesky(G + np.eye(G.shape[0]) * 1e-12); covs = []
    for _ in range(nrep):
        B = (rng.normal(size=(M, G.shape[0])) @ Gc.T) / np.sqrt(M)   # per-SNP effects, total genetic covariance G
        z = (W @ B.T) * np.sqrt(neff)[:, None] + Lc @ rng.normal(size=(6, M))
        C = ((z @ z.T) / M - R) / np.sqrt(np.outer(neff, neff)) * M
        covs.append(C)
    return np.array(covs)


def run(i, sc, nrep):
    covs = simulate(sc, nrep, seed=2000 + i); IU = S1.IU
    V = np.cov(np.array([c[IU] for c in covs]).T) + np.eye(len(IU[0])) * 1e-12; Vinv = np.linalg.inv(V); rows = []
    for c in covs:
        r_sp = c[1, 5] / np.sqrt(c[1, 1] * c[5, 5]) if c[1, 1] > 0 and c[5, 5] > 0 else np.nan
        x1, d1, _ = S1.fit_cfa(c, Vinv, 1); x2, d2, phi = S1.fit_cfa(c, Vinv, 2)
        rows.append(dict(r=r_sp, p=float(chi2.sf(max(x1 - x2, 0), d1 - d2)), aic=bool(x2 + 26 < x1 + 24), phi=phi))
    Rr = pd.DataFrame(rows)
    return dict(scen=i, name=sc["name"], nrep=nrep, rg_spark_median=Rr.r.median(), rg_spark_q05=Rr.r.quantile(0.05), rg_spark_q95=Rr.r.quantile(0.95),
                frac_2f_p05=float((Rr.p < 0.05).mean()), frac_aic_2f=float(Rr.aic.mean()), phi_2f_median=Rr.phi.median())


if __name__ == "__main__":
    SC = os.environ.get("SCEN", "all"); NREP = int(os.environ.get("NREP", "200"))
    idx = range(len(GRID)) if SC == "all" else [int(SC)]
    T = pd.DataFrame([run(i, GRID[i], NREP) for i in idx])
    T.to_csv(f"{OUT}/sim_twofactor2{'' if SC == 'all' else '_' + SC}.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220); print(T.round(3).to_string(index=False), flush=True)
