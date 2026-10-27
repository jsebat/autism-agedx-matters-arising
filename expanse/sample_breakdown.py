#!/usr/bin/env python3
"""Sample breakdown by subcohort for the Matters Arising : where the CCA sample of 7,072 probands is lost.
Step counts per cohort and role (proband / unaffected sibling / parent where present) through the CCA filters:
  rare-variant table -> proband with >= 1 unaffected sibling -> polygenic scores available -> SPARK phenotype record (characteristics).
Also: age at diagnosis in the CCA sample vs all SPARK autistic registrants (range restriction check, C38), and its correlation with
the EDU/SES characteristics in both. Output (aggregate counts only): ~/rvh2/agedx/sample_breakdown.txt"""
import os, io, contextlib
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
# rare-variant tier definitions (SNV/indel), as in the lab's exploratory rare_survey.py
TI = ["lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4"]
cell = lambda tiers, modes, sets: [f"{t}__{m}__{g}" for t in tiers for m in modes for g in sets]
CAT = {"de novo, NDD genes": cell(TI, ["dnm"], ["ndd"]), "de novo, other genes": cell(TI, ["dnm"], ["nonndd"]),
       "de novo LoF": cell(["lof_t1", "lof_t2"], ["dnm"], ["ndd", "nonndd"]), "de novo missense": cell(TI[2:], ["dnm"], ["ndd", "nonndd"]),
       "inherited, NDD genes": cell(TI, ["inh"], ["ndd"]), "inherited, other genes": cell(TI, ["inh"], ["nonndd"]),
       "inherited LoF tier 1": cell(["lof_t1"], ["inh"], ["ndd", "nonndd"]), "inherited missense": cell(TI[2:], ["inh"], ["ndd", "nonndd"]),
       "all rare SNV/indel": cell(TI, ["dnm", "inh"], ["ndd", "nonndd"])}
def load(root, cats):
    a = pd.read_parquet(f"/expanse/projects/sebat1/jsebat/{root}/assembled_children.parquet"); a = a[a.cohort.str.startswith("spark")].copy()
    W = pd.read_csv(f"/expanse/projects/sebat1/jsebat/{root}/results/collider_weights.tsv", sep="\t", index_col=0).iloc[:, 0]
    for k, ts in cats.items():
        ts = [t for t in ts if t in a.columns and t in W.index]
        a[k] = (a[ts].fillna(0).to_numpy(float) * W[ts].to_numpy()).sum(1)
        ref = a[a.case.eq(0)].groupby("cohort")[k].agg(["mean", "std"]); a[k] = (a[k] - a.cohort.map(ref["mean"])) / a.cohort.map(ref["std"])
    return a
RC = ["de novo, NDD genes", "de novo, other genes", "inherited, NDD genes", "inherited, other genes", "inherited LoF tier 1"]
a = load("pgs_rare_h2_drift_full", {k: CAT[k] for k in RC})
buf = io.StringIO(); pr = lambda *x: print(*x, file=buf)
pr("columns:", list(a.columns)[:40])
pr("\n== rare-variant table: rows by cohort x case ==")
pr(a.groupby(["cohort", "case"]).size().unstack(fill_value=0).to_string())
if "role" in a.columns: pr(a.groupby(["cohort", "role"]).size().unstack(fill_value=0).to_string())
sibfam = a[a.case.eq(0)].groupby(["cohort", "family_key"]).size()
pro = a[a.case.eq(1)].copy()
pro["has_sib"] = pd.MultiIndex.from_frame(pro[["cohort", "family_key"]]).isin(sibfam.index)
cs = pd.read_parquet(f"{OUT}/cases.parquet").drop_duplicates("sample_id").set_index("sample_id")
pro["has_pgs"] = pro["sample"].map(cs["asd_2026"]).notna()
pro["has_ptdt"] = pro["sample"].map(cs["asd_2026_ptdt"]).notna()
reg = pd.read_csv(f"{P}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
pro["spark_pheno"] = pro["sample"].isin(reg.index)
pr("\n== probands: cumulative filters by cohort ==")
rows = []
for c, g in pro.groupby("cohort"):
    rows.append(dict(cohort=c, probands=len(g), with_sib=int(g.has_sib.sum()), with_pgs=int(g.has_pgs.sum()), sib_and_pgs=int((g.has_sib & g.has_pgs).sum()),
                     pgs_no_sib=int((g.has_pgs & ~g.has_sib).sum()), spark_pheno=int(g.spark_pheno.sum()), sib_pgs_pheno=int((g.has_sib & g.has_pgs & g.spark_pheno).sum())))
T = pd.DataFrame(rows); pr(T.to_string(index=False)); pr("total", T.drop(columns="cohort").sum().to_dict())
pr("\n== unaffected siblings in rare table, by cohort ==", a[a.case.eq(0)].groupby("cohort").size().to_dict())
pr("\n== cases.parquet (GWAS/PGS case table) ==", len(cs), "rows;", list(cs.columns)[:25])
# age at diagnosis: CCA-type sample vs all SPARK autistic registrants
dxa = pd.to_numeric(reg.diagnosis_age, errors="coerce").where(reg.asd.astype(str).str.upper().eq("TRUE")) / 12
env = pd.read_csv(f"{P}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
sv = pd.read_csv(f"{P}/spark_survey.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
cca_ids = pro.loc[pro.has_sib & pro.has_pgs & pro.spark_pheno, "sample"]
def desc(ids, nm):
    x = dxa.reindex(ids).dropna(); pr(f"{nm}: n={len(x)}, age at dx mean {x.mean():.2f} y, sd {x.std():.2f}, median {x.median():.1f}, >10 y {100*(x>10).mean():.1f}%, <6 y {100*(x<6).mean():.1f}%")
    for v, ser in [("ADI", env.adi_pct), ("maternal age", env.maternal_age), ("low income third", sv.s09), ("mother no degree", sv.s10)]:
        y = pd.to_numeric(ser.reindex(x.index), errors="coerce"); ok = y.notna()
        pr(f"   r(age at dx, {v}) = {np.corrcoef(x[ok], y[ok])[0, 1]:+.3f} (n={ok.sum()})")
pr("\n== age at diagnosis: range restriction check ==")
desc(dxa.dropna().index, "all SPARK autistic registrants"); desc(cca_ids, "CCA sample (sib + PGS + SPARK phenotype)")
desc(pro.loc[pro.has_pgs & pro.spark_pheno, "sample"], "probands with PGS + phenotype, sibling not required")
open(f"{OUT}/sample_breakdown.txt", "w").write(buf.getvalue()); print(buf.getvalue())
