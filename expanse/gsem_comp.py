#!/usr/bin/env python3
"""What defines our genomic factors? For the 2- and 3-factor EFA of the SPARK split-half GWAS (gsem_fa.R), correlate each
indicator's loadings with the composition of its cases (survey_composition.tsv): fraction diagnosed > 10 y, registered as adult,
female, high deprivation, ADHD, ID, born before 2005, no parent genotyped. Weighted by nothing; Spearman and Pearson reported.
Also the loadings of the two age-split halves for reference. Output: gsem/loading_composition.tsv"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
C = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined"; G = f"{C}/gsem"
comp = pd.read_csv(f"{C}/survey_composition.tsv", sep="\t").set_index("gwas")
rows = []
for k in (2, 3):
    f = f"{G}/efa_loadings_k{k}.tsv"
    if not os.path.exists(f): continue
    L = pd.read_csv(f, sep="\t").set_index("trait"); X = L.join(comp, how="inner")
    for fac in L.columns:
        for v in [c for c in comp.columns if c != "n"]:
            ok = X[[fac, v]].dropna()
            if len(ok) < 6: continue
            rows.append(dict(k=k, factor=fac, composition=v, n_ind=len(ok), spearman=spearmanr(ok[fac], ok[v])[0], pearson=pearsonr(ok[fac], ok[v])[0], p_pearson=pearsonr(ok[fac], ok[v])[1]))
    print(f"\nEFA k={k}, loadings with composition"); print(X[list(L.columns) + ["late_dx", "adult_reg", "female", "high_adi", "adhd", "id"]].round(2).sort_values(L.columns[-1]).to_string())
T = pd.DataFrame(rows); T.to_csv(f"{G}/loading_composition.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(T.round(3).sort_values(["k", "factor", "p_pearson"]).to_string(index=False))
