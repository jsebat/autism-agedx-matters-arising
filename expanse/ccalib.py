#!/usr/bin/env python3
"""Ridge-regularized CCA with fully held-out evaluation (shared by cca2.py).

Changes from the first version (review, 2026-10-26):
  * every preprocessing step is learned on the training fold only: median / binary imputation, missing-value indicators (training missingness > 20%),
    residualization on the adjustment covariates Z, dropping constant columns, standardization;
  * the ridge penalty is chosen by inner family cross-validation inside each outer training fold (nested CV);
  * held-out canonical correlations are SIGNED; components are matched across folds to a reference genetic structure (Hungarian matching on the
    absolute correlation of genetic-block structure coefficients, then a joint sign flip of the pair, which leaves the held-out correlation's sign intact);
  * permutation null: Freedman-Lane (residuals of Y on Z permuted, fitted part kept), permuting whole families among families of the same size,
    with the penalty re-tuned inside every permutation; axes 2..K tested sequentially after deflating Y on the first j-1 genetic variates;
  * uncertainty: family bootstrap of the whole pipeline (paired across axes and single characteristics in the same resample).
"""
import numpy as np
from scipy.optimize import linear_sum_assignment

LAMS = (0.01, 0.05, 0.1, 0.3, 1.0)
NF = 5


def folds(fam, rng, nf=NF):
    uf, inv = np.unique(fam, return_inverse=True)
    return rng.integers(0, nf, len(uf))[inv]


def make(X, cont, Y, Z, tr, te, binfill=0.5, names=None):
    """Learn preprocessing on rows tr, apply to tr and te. X: n x p float with NaN; cont: bool p; Y: n x q (complete); Z: n x m (with intercept)."""
    Xt = X[tr]
    med = np.nanmedian(Xt, axis=0) if np.isnan(Xt).any() else np.median(Xt, axis=0)
    fill = np.where(cont, med, med if binfill is None else binfill); fill = np.where(np.isnan(fill), 0.0, fill)
    miss = np.isnan(X); ind = np.isnan(Xt).mean(0) > 0.2
    M = np.column_stack([np.where(miss, fill, X), miss[:, ind].astype(float)])
    bm = np.linalg.lstsq(Z[tr], M[tr], rcond=None)[0]; by = np.linalg.lstsq(Z[tr], Y[tr], rcond=None)[0]
    M = M - Z @ bm; Yr = Y - Z @ by
    keep = M[tr].std(0) > 1e-8; M = M[:, keep]
    mx, sx = M[tr].mean(0), M[tr].std(0); my, sy = Yr[tr].mean(0), Yr[tr].std(0) + 1e-12
    Mx, My = (M - mx) / sx, (Yr - my) / sy
    cols = None
    if names is not None:
        cols = [n for n, k in zip(list(names) + [f"{n} [missing]" for n, i in zip(names, ind) if i], keep) if k]
    return Mx[tr], Mx[te], My[tr], My[te], cols


def fit(Xm, Ym, lam, k):
    k = min(k, Ym.shape[1], Xm.shape[1]); n = len(Xm)
    Sxx = Xm.T @ Xm / n + lam * np.eye(Xm.shape[1]); Syy = Ym.T @ Ym / n + lam * np.eye(Ym.shape[1]); Sxy = Xm.T @ Ym / n
    ix = np.linalg.inv(np.linalg.cholesky(Sxx)); iy = np.linalg.inv(np.linalg.cholesky(Syy))
    U, d, Vt = np.linalg.svd(ix @ Sxy @ iy.T)
    return ix.T @ U[:, :k], iy.T @ Vt.T[:, :k], d[:k]


def gstruct(Ym, B):
    """Genetic-block structure coefficients: corr(Y column, genetic variate)."""
    V = Ym @ B; V = (V - V.mean(0)) / (V.std(0) + 1e-12)
    return (Ym.T @ V) / len(Ym)


def choose_lam(X, cont, Y, Z, fam, seed, binfill):
    """Inner family CV on the given (training) rows: penalty maximizing the held-out first canonical correlation."""
    f = folds(fam, np.random.default_rng(seed)); best = (-np.inf, LAMS[0])
    for L in LAMS:
        u = np.empty(len(Y)); v = np.empty(len(Y))
        for q in range(NF):
            a, b = f != q, f == q
            Xa, Xb, Ya, Yb, _ = make(X, cont, Y, Z, a, b, binfill)
            A, B, _ = fit(Xa, Ya, L, 1); u[b] = (Xb @ A)[:, 0]; v[b] = (Yb @ B)[:, 0]
        r = np.corrcoef(u, v)[0, 1]
        if r > best[0]: best = (r, L)
    return best[1]


def nested(X, cont, Y, Z, fam, K, ref=None, seed=0, binfill=0.5, lam=None):
    """Outer 5-fold family CV with inner penalty tuning. Returns signed held-out r (K), chosen penalties, held-out U, V (n x K).
    ref: q x K genetic structure that defines component identity and sign (default: the first fold's training fit)."""
    n = len(Y); f = folds(fam, np.random.default_rng(seed))
    U = np.full((n, K), np.nan); V = np.full((n, K), np.nan); lams = []
    for q in range(NF):
        tr, te = f != q, f == q
        L = lam if lam is not None else choose_lam(X[tr], cont, Y[tr], Z[tr], fam[tr], seed + 101 + q, binfill)
        lams.append(L)
        Xa, Xb, Ya, Yb, _ = make(X, cont, Y, Z, tr, te, binfill)
        A, B, _ = fit(Xa, Ya, L, K + 2); S = gstruct(Ya, B)
        if ref is None: ref = S[:, :K].copy()
        C = np.array([[np.corrcoef(S[:, i], ref[:, j])[0, 1] if S.shape[0] > 1 else np.sign(S[0, i] * ref[0, j]) for j in range(ref.shape[1])] for i in range(S.shape[1])])
        rows, cols = linear_sum_assignment(-np.abs(C))
        for i, j in zip(rows, cols):
            if j >= K: continue
            sg = np.sign(C[i, j]) or 1.0
            U[te, j] = sg * (Xb @ A[:, i]); V[te, j] = sg * (Yb @ B[:, i])
    r = np.array([np.corrcoef(U[:, j], V[:, j])[0, 1] if np.isfinite(U[:, j]).all() else np.nan for j in range(K)])
    return r, lams, U, V, ref


def full_fit(X, cont, Y, Z, fam, K, binfill=0.5, names=None, seed=0):
    """Fit on all rows (penalty by CV on all rows), for structure coefficients and the reference structure."""
    L = choose_lam(X, cont, Y, Z, fam, seed + 7, binfill)
    allr = np.ones(len(Y), bool)
    Xa, _, Ya, _, cols = make(X, cont, Y, Z, allr, allr, binfill, names)
    A, B, d = fit(Xa, Ya, L, K)
    U, V = Xa @ A, Ya @ B
    sx = np.array([[np.corrcoef(U[:, j], Xa[:, i])[0, 1] for j in range(A.shape[1])] for i in range(Xa.shape[1])])
    sy = np.array([[np.corrcoef(V[:, j], Ya[:, i])[0, 1] for j in range(B.shape[1])] for i in range(Ya.shape[1])])
    return dict(lam=L, sx=sx, sy=sy, cols=cols, V=V, U=U, d=d)


def family_perm(fam, rng):
    """Row index that moves each family's rows onto another family of the same size (order within family kept)."""
    fam = np.asarray(fam); order = np.argsort(fam, kind="stable")
    _, start, cnt = np.unique(fam[order], return_index=True, return_counts=True)
    new = np.empty(len(fam), int)
    for s in np.unique(cnt):
        g = np.where(cnt == s)[0]; pg = rng.permutation(g)
        src = np.concatenate([order[start[b]:start[b] + s] for b in pg]); dst = np.concatenate([order[start[a]:start[a] + s] for a in g])
        new[dst] = src
    return new


def fl_perm(Y, Z, fam, rng):
    """Freedman-Lane: keep Z-fitted part of Y, permute residuals among families of the same size."""
    bz = np.linalg.lstsq(Z, Y, rcond=None)[0]; Yh = Z @ bz
    return Yh + (Y - Yh)[family_perm(fam, rng)]


def boot_rows(fam, rng):
    """Family bootstrap: resample families with replacement; returns row index and new family labels."""
    uf, inv = np.unique(fam, return_inverse=True)
    rows_by = np.split(np.argsort(inv, kind="stable"), np.cumsum(np.bincount(inv))[:-1])
    draw = rng.integers(0, len(uf), len(uf))
    idx = np.concatenate([rows_by[d] for d in draw]); newf = np.concatenate([np.full(len(rows_by[d]), i) for i, d in enumerate(draw)])
    return idx, newf
