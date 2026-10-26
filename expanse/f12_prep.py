#!/usr/bin/env python3
"""Factor PGS for Zhang, Grove et al. 2025: convert the published genomic-SEM factor GWAS (figshare 29566052 v2) to COJO .ma for
SBayesRC. Earlier_diagnosed_autism = F1, Later_diagnosed_autism = F2. A1 = effect allele; MAF column is the A1 frequency (figshare
README). N = N_hat (GenomicSEM implied sample size, the N on which Beta/SE are scaled). Duplicated SNP ids dropped.
Caveat recorded in outputs: both factors contain SPARK cases (SPARKbefore6 in F1, SPARKafter10 in F2), so scores in SPARK are
inflated for the cases that were in each discovery GWAS."""
import pandas as pd
W = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_sumstats"; O = "/expanse/projects/sebat1/jsebat/pgs_warrier/input"
import os; os.makedirs(O, exist_ok=True)
for src, nm in [("Earlier_diagnosed_autism.txt", "warrier_F1"), ("Later_diagnosed_autism.txt", "warrier_F2")]:
    x = pd.read_csv(f"{W}/{src}", sep=r"\s+", low_memory=False)
    x = x[x.SNP.astype(str).str.startswith("rs")].drop_duplicates("SNP", keep=False)
    m = pd.DataFrame({"SNP": x.SNP, "A1": x.A1, "A2": x.A2, "freq": x.MAF, "b": x.Beta, "se": x.SE, "p": x.P, "N": x.N_hat.round().astype(int)})
    m = m.dropna(); m = m[(m.se > 0) & (m.freq > 0) & (m.freq < 1)]
    m.to_csv(f"{O}/{nm}.ma", sep="\t", index=False); print(nm, len(m), "median N", m.N.median())
open("/expanse/projects/sebat1/jsebat/pgs_warrier/gwas_list.txt", "w").write("warrier_F1.ma\nwarrier_F2.ma\n")
