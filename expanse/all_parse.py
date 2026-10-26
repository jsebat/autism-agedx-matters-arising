#!/usr/bin/env python3
"""Parse the ALL-DATA LDSC logs (warrier_rep/combined_all/rg/*.log; all_ldsc.sbatch) into the same tables as the arrays-only analysis:
  rg_pairs_all.tsv   age-split pairs (shared / disjoint controls) and random splits rand1-4, with h2 z of both halves
  null_rg_all.tsv    random-split null (ageN / adiN replicates)
  survey_rg_all.tsv  per survey split: rg(high, low), h2 / h2 z of both halves, each half's rg with W_F1, W_F2, W_agedx, late/early GWAS
  rg_ext_all.tsv     external-trait rg profile of every half (AFB, EA, IQ, ADHD, MDD, SCZ) and of the published factors
Aggregate only. Output in combined_all/ and copies in ~/rvh2/agedx/all_*."""
import glob, io, os, re, shutil
import numpy as np, pandas as pd
C = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined_all"; OUT = os.path.expanduser("~/rvh2/agedx")
nm = lambda s: os.path.basename(str(s)).replace(".sumstats.gz", "")
def summ(f):
    t = open(f).read()
    if "Summary of Genetic Correlation Results" not in t: return None, None
    blk = t.split("Summary of Genetic Correlation Results")[1].strip().split("\n\n")[0]; x = pd.read_csv(io.StringIO(blk), sep=r"\s+")
    x["p1"] = x.p1.map(nm); x["p2"] = x.p2.map(nm)
    m = re.search(r"Heritability of phenotype 1\n-+\nTotal Observed scale h2: ([-\d.e]+) \(([-\d.e]+)\)", t)
    return x, ((float(m.group(1)), float(m.group(2))) if m else None)
pairs, nulls, surv, ext, h2 = [], [], [], [], {}
for f in sorted(glob.glob(f"{C}/rg/*.log")):
    kind, me = os.path.basename(f)[:-4].split("__", 1)
    x, h = summ(f)
    if x is None: continue
    if h: h2[me] = h
    if kind == "pair": pairs.append(x)
    elif kind == "null": nulls.append(x.iloc[[0]])
    elif kind == "survey": surv.append(x)
    elif kind == "ext": ext.append(x)
A = pd.concat(pairs + surv, ignore_index=True) if pairs or surv else pd.DataFrame()
H = pd.Series({k: v[0] / v[1] for k, v in h2.items()}, name="h2_z")
def get(a, b):
    x = A[(A.p1 == a) & (A.p2 == b)]
    return x.iloc[0] if len(x) else None
rows = []
for a, b, lab in [("rep_W_under6", "rep_W_over10", "age split, shared controls"), ("rep_W_under6_c1", "rep_W_over10_c2", "age split, disjoint controls")] + \
                 [(f"rep_rand{r}_A_c2", f"rep_rand{r}_B_c1", f"random split {r}") for r in range(1, 5)]:
    r = get(a, b)
    if r is not None: rows.append(dict(pair=lab, rg=r.rg, se=r.se, p=r.p, h2z_a=H.get(a, np.nan), h2z_b=H.get(b, np.nan)))
P = pd.DataFrame(rows); P.to_csv(f"{C}/rg_pairs_all.tsv", sep="\t", index=False)
N = pd.concat(nulls, ignore_index=True) if nulls else pd.DataFrame(columns=["p1", "rg", "se"])
N = N.assign(set=N.p1.str.extract(r"rep_(ageN|adiN)")[0])[["p1", "set", "rg", "se"]]; N.to_csv(f"{C}/null_rg_all.tsv", sep="\t", index=False)
out = []
for hi in sorted(k for k in h2 if k.endswith("_hi_c2")):
    lo = hi.replace("_hi_c2", "_lo_c1"); d = dict(code=hi.replace("rep_", "").replace("_hi_c2", ""))
    r = get(hi, lo)
    if r is not None: d.update(rg=r.rg, se=r.se, p=r.p)
    for k_, k in [("hi", hi), ("lo", lo)]:
        if k in h2: d[f"h2_{k_}"], d[f"h2se_{k_}"] = h2[k]; d[f"h2z_{k_}"] = h2[k][0] / h2[k][1]
        for ref in ["W_F1", "W_F2", "W_agedx", "rep_W_over10_c1", "rep_W_under6_c1"]:
            q = get(k, ref)
            if q is not None: d[f"{k_}_rg_{ref}"] = q.rg; d[f"{k_}_se_{ref}"] = q.se
    out.append(d)
S = pd.DataFrame(out); S.to_csv(f"{C}/survey_rg_all.tsv", sep="\t", index=False)
E = pd.concat(ext, ignore_index=True) if ext else pd.DataFrame()
if len(E): E = E.rename(columns={"p1": "gwas", "p2": "trait"}).assign(trait=lambda d: d.trait.str.replace("ext_", ""))[["gwas", "trait", "rg", "se"]]
E.to_csv(f"{C}/rg_ext_all.tsv", sep="\t", index=False)
for f in ["rg_pairs_all.tsv", "null_rg_all.tsv", "survey_rg_all.tsv", "rg_ext_all.tsv", "counts_all.tsv", "lambda.tsv"]:
    if os.path.exists(f"{C}/{f}"): shutil.copy(f"{C}/{f}", f"{OUT}/all_{f}")
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 200)
print(P.round(3).to_string(index=False)); print(N.groupby("set").rg.describe().round(3).to_string())
print(S[["code", "rg", "se", "h2z_hi", "h2z_lo"]].sort_values("rg").round(2).to_string(index=False))
