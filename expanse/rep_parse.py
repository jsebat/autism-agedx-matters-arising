#!/usr/bin/env python3
"""Replication step 7: parse LDSC rg logs (warrier_rep/combined/rg) into rg_all.tsv with each trait's observed-scale h2 and h2 z,
flag pairs where either trait has h2 z < 4 (LDSC guidance: rg unreliable), and print the key tables: within-pair rg (shared vs
disjoint controls), the random-split null, rg with the published factors F1/F2 and SPARK age at diagnosis."""
import glob, io, os
import numpy as np, pandas as pd
C = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined"; rows = []
for f in sorted(glob.glob(f"{C}/rg/*.log")):
    t = open(f).read()
    if "Summary of Genetic Correlation Results" not in t: print("no summary:", os.path.basename(f)); continue
    blk = t.split("Summary of Genetic Correlation Results")[1].strip().split("\n\n")[0]; rows.append(pd.read_csv(io.StringIO(blk), sep=r"\s+"))
A = pd.concat(rows, ignore_index=True); nm = lambda s: os.path.basename(str(s)).replace(".sumstats.gz", "")
A["p1"] = A.p1.map(nm); A["p2"] = A.p2.map(nm); A = A[A.p1 != A.p2]
H = A.drop_duplicates("p2").set_index("p2")[["h2_obs", "h2_obs_se"]]; H["h2_z"] = H.h2_obs / H.h2_obs_se
A["h2z_1"] = A.p1.map(H.h2_z); A["h2z_2"] = A.p2.map(H.h2_z); A["reliable"] = (A.h2z_1 >= 4) & (A.h2z_2 >= 4)
A.to_csv(f"{C}/rg_all.tsv", sep="\t", index=False); H.to_csv(f"{C}/h2.tsv", sep="\t")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print(H[H.index.str.startswith("rep_") | H.index.str.startswith("W_")].round(4).to_string())
def g(a, b):
    x = A[(A.p1 == a) & (A.p2 == b)]
    if not len(x): x = A[(A.p1 == b) & (A.p2 == a)]
    return x.iloc[0] if len(x) else None
PAIRS = [("W_under6", "W_over10"), ("adi_low", "adi_high"), ("adult_child", "adult_adult"), ("adhd_no", "adhd_yes"), ("adiE_low", "adiE_high"),
         ("adhdE_no", "adhdE_yes"), ("sexE_male", "sexE_female"), ("agedxC_early", "agedxC_late")]
out = []
for a, b in PAIRS:
    for lab, x, y in [("shared controls", f"rep_{a}", f"rep_{b}"), ("disjoint controls", f"rep_{a}_c1", f"rep_{b}_c2")]:
        r = g(x, y)
        if r is not None: out.append(dict(pair=f"{a} vs {b}", controls=lab, rg=r.rg, se=r.se, p=r.p, h2z_a=r.h2z_1 if r.p1 == x else r.h2z_2, h2z_b=r.h2z_2 if r.p1 == x else r.h2z_1, reliable=r.reliable))
for k in range(1, 5):
    r = g(f"rep_rand{k}_A_c2", f"rep_rand{k}_B_c1")
    if r is not None: out.append(dict(pair=f"random split {k} (sizes of over10 / under6)", controls="disjoint controls", rg=r.rg, se=r.se, p=r.p, h2z_a=r.h2z_1, h2z_b=r.h2z_2, reliable=r.reliable))
O = pd.DataFrame(out); O.to_csv(f"{C}/rg_pairs.tsv", sep="\t", index=False); print("\nwithin-pair rg"); print(O.round(3).to_string(index=False))
F = A[A.p2.isin(["W_F1", "W_F2", "W_agedx", "sibtrio_agedx_early"]) & A.p1.str.startswith("rep_")]
F.to_csv(f"{C}/rg_factors.tsv", sep="\t", index=False)
print("\nrg with published factors (rg / se)")
print(F.pivot_table(index="p1", columns="p2", values="rg").round(2).join(F.pivot_table(index="p1", columns="p2", values="se").round(2), rsuffix="_se").to_string())
cross = [("rep_W_over10_c1", f"rep_{x}") for x in ["adi_high_c2", "adult_adult_c2", "adhd_yes_c2", "adiE_high_c2", "adhdE_yes_c2", "sexE_female_c2", "agedxC_late_c2", "W_under6_c1"]]
print("\nlate-diagnosis GWAS (C1 controls) vs other GWAS (C2 controls)")
for a, b in cross:
    r = g(a, b)
    if r is not None: print(f"{b:24s} rg {r.rg:6.2f} ({r.se:.2f})  reliable={r.reliable}")
E = A[A.p2.str.startswith("ext_") & (A.p1.str.startswith("rep_") | A.p1.str.startswith("W_F"))]; E.to_csv(f"{C}/rg_external.tsv", sep="\t", index=False)
