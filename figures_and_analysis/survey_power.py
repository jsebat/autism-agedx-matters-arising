#!/usr/bin/env python3
"""Power of each split GWAS to detect a drop in genetic correlation as large as the age split's (Supplementary Table 7, split-GWAS part).
Effect = mean rg of the random-split null (ageN replicates) minus the age split's rg; power = Phi(effect / s.e. - 1.96) for each split's own rg s.e.;
minimum detectable drop at 80% power = 2.8 x s.e. Inputs (aggregate): Manuscripts/Nature_LTE/tables/survey_table.tsv (letter_figs.py),
results/full_data/all_null_rg_all.tsv (all_parse.py). Output: results/full_data/survey_power_all.tsv"""
import os
import pandas as pd
from scipy.stats import norm
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..", "results", "full_data")
S = pd.read_csv(os.path.join(H, "..", "..", "..", "Manuscripts", "Nature_LTE", "tables", "survey_table.tsv"), sep="\t")
N = pd.read_csv(os.path.join(R, "all_null_rg_all.tsv"), sep="\t")
effect = N[N.set == "ageN"].rg.mean() - S.set_index("code").loc["s00", "rg"]
S["mde"] = 2.8 * S.se; S["power_age_effect"] = norm.cdf(effect / S.se - 1.96)
S.to_csv(os.path.join(R, "survey_power_all.tsv"), sep="\t", index=False)
print(f"effect {effect:.3f}"); print(S[["code", "split", "se", "power_age_effect"]].round(3).sort_values("power_age_effect").to_string(index=False))
