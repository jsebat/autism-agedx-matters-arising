#!/usr/bin/env python3
"""Broad survey, post-processing (aggregate outputs only):
1) survey_rg.tsv - per split: rg(high half, low half) with SE, h2 and h2 z of both halves, and each half's rg with the published
   factors F1/F2, the SPARK age-at-diagnosis GWAS and the replicated late (>10 y) / early (<6 y) GWAS (from rg_survey/*.log).
2) gsem/fa_traits.txt - indicators for the factor analysis: halves with h2 z >= 3 (age-split halves always kept).
3) survey_composition.tsv - composition of every half's cases pooled over platforms (fraction diagnosed > 10 y among those with a
   diagnosis age, registered as adult, female, high-ADI third, ADHD, ID, born before 2005, no parent genotyped)."""
import glob, io, os, re
import numpy as np, pandas as pd
C = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined"; B = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep"
os.makedirs(f"{C}/gsem", exist_ok=True)
h2, rows = {}, []
for f in sorted(glob.glob(f"{C}/rg_survey/*.log")):
    t = open(f).read(); me = os.path.basename(f).replace(".log", "")
    m = re.search(r"Heritability of phenotype 1\n-+\nTotal Observed scale h2: ([-\d.e]+) \(([-\d.e]+)\)", t)
    if m: h2[me] = (float(m.group(1)), float(m.group(2)))
    if "Summary of Genetic Correlation Results" in t:
        blk = t.split("Summary of Genetic Correlation Results")[1].strip().split("\n\n")[0]; x = pd.read_csv(io.StringIO(blk), sep=r"\s+")
        x["p1"] = me; x["p2"] = x.p2.map(lambda s: os.path.basename(s).replace(".sumstats.gz", "")); rows.append(x)
A = pd.concat(rows, ignore_index=True)
get = lambda a, b: A[(A.p1 == a) & (A.p2 == b)]
out = []
for hi in sorted(k for k in h2 if k.endswith("_hi_c2")):
    lo = hi.replace("_hi_c2", "_lo_c1"); code = hi.replace("rep_", "").replace("_hi_c2", "")
    r = get(hi, f"munged/{lo}" if False else lo); d = dict(code=code)
    if len(r): d.update(rg=r.rg.iloc[0], se=r.se.iloc[0], p=r.p.iloc[0])
    for nm, k in [("hi", hi), ("lo", lo)]:
        if k in h2: d[f"h2_{nm}"], d[f"h2se_{nm}"] = h2[k]; d[f"h2z_{nm}"] = h2[k][0] / h2[k][1]
        for ref in ["W_F1", "W_F2", "W_agedx", "rep_W_over10_c1", "rep_W_under6_c1"]:
            q = get(k, ref)
            if len(q): d[f"{nm}_rg_{ref}"] = q.rg.iloc[0]; d[f"{nm}_se_{ref}"] = q.se.iloc[0]
    out.append(d)
T = pd.DataFrame(out); T["reliable"] = (T.h2z_hi >= 4) & (T.h2z_lo >= 4); T.to_csv(f"{C}/survey_rg.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(T[["code", "rg", "se", "h2z_hi", "h2z_lo", "hi_rg_W_F1", "hi_rg_W_F2", "lo_rg_W_F1", "lo_rg_W_F2"]].sort_values("rg").round(2).to_string(index=False))
keep = [k for k, (h, s) in h2.items() if (h / s >= 3) or k.startswith("rep_s00_")]
open(f"{C}/gsem/fa_traits.txt", "w").write("\n".join(sorted(keep)) + "\n"); print("factor-analysis indicators:", len(keep))
# composition
sv = pd.read_csv("/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
comp_vars = {"late_dx": "s00", "adult_reg": "s01", "female": "s07", "high_adi": "s08", "adhd": "s27", "id": "s13", "born_pre2005": "s42"}
recs = []
for src in ["spark_gsa_v1", "spark_gsa_v2", "spark_twist"]:
    W = f"{B}/{src}"; names = open(f"{W}/survey_names.txt").read().split()
    P = pd.read_csv(f"{W}/survey_pheno.txt", sep=" ", header=None, dtype=str); P.columns = ["FID", "IID"] + names
    d = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, usecols=["source", "IID", "sample_id"])
    sp = P.IID.map(d[d.source == src].drop_duplicates("IID").set_index("IID").sample_id)
    np99 = pd.to_numeric(P["s99_hi_c2"], errors="coerce").eq(1) if "s99_hi_c2" in P else None
    for k in names:
        m = pd.to_numeric(P[k], errors="coerce").eq(1); ids = sp[m]
        r = dict(source=src, gwas=f"rep_{k}", n=int(m.sum()))
        for nm, c in comp_vars.items():
            v = pd.to_numeric(ids.map(sv[c]), errors="coerce"); r[f"sum_{nm}"] = float(v.sum()); r[f"n_{nm}"] = int(v.notna().sum())
        if np99 is not None: r["sum_no_parent"] = float((m & np99).sum()); r["n_no_parent"] = int(m.sum())
        recs.append(r)
R = pd.DataFrame(recs).groupby("gwas").sum(numeric_only=True)
for nm in list(comp_vars) + ["no_parent"]:
    if f"sum_{nm}" in R: R[nm] = R[f"sum_{nm}"] / R[f"n_{nm}"].replace(0, np.nan)
R[["n"] + [c for c in list(comp_vars) + ["no_parent"] if c in R]].reset_index().to_csv(f"{C}/survey_composition.tsv", sep="\t", index=False)
print("composition written", len(R))
