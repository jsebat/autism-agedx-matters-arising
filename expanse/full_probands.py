#!/usr/bin/env python3
"""Full proband table for the age-at-diagnosis Matters Arising (use all data; Zhang-era SPARK as the key sensitivity
analysis; family structure 0/1/2 parents as a sensitivity analysis). See README_DATA_WARNING.md: pgs_results_v2 alone misses the new release.

Sources
  PGS, SPARK arrays + iWGS v1.1 (and SSC): s3/.../pgs_results_v2/all_cohorts_with_ped_full.tsv (+ ASD 2026 leave-SPARK-wave-out score,
      jsebat/pgs_asd2026_loo/asd2026_loo_scores.tsv). These participants were available to Zhang et al. -> zhang_era = True.
  PGS, SPARK-WGS-2026_08: j3guevar/.../plink/pgs/analyses/scoring-ancestry-20261012/05_scores/*.sscore (SDSM IIDs); ASD = Grove April 2026 FULL GWAS (jsebat/pgs_asd2026_loo/scores/spark_wgs_2026_08.asd2026_full), not asd_loSSC,
      ancestry 06_pca/ancestry/ancestry_probabilities.tsv, SDSM -> SP map rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv.
      The 212 already in pgs_results_v2 keep their older score (zhang_era = True).
  PGS standardisation: within source, to non-autistic relatives of the same source (as before).
  Pedigree: pheno/spark_pedigree_ids.tsv (SPARK); pheno/ssc_ids_sanders.csv (SSC, Sanders roles).
  Parents genotyped = parent has a PGS row (any source); parents sequenced = parent in iWES v3 / iWGS v1.1 ped or flagged sequenced in the new
      release metadata (SSC: all quads sequenced).
  Rare-variant classes: jsebat/pgs_rare_h2_drift/assembled_children.parquet (families with an unaffected sibling; trios).
Outputs (Expanse only for participant-level): ~/rvh2/agedx/full_probands.parquet; aggregate counts ~/rvh2/agedx/full_probands_counts.tsv and .txt"""
import os, io
import numpy as np, pandas as pd
S = "/expanse/projects/sebat1/s3/data/sebat"; J = "/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08"; P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift"
OUT = os.path.expanduser("~/rvh2/agedx"); buf = io.StringIO(); pr = lambda *x: print(*x, file=buf, flush=True)
PGS = ["asd", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
# ---- SPARK + SSC, older table
old = pd.read_csv(f"{S}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype={"sample_id": str, "FID": str, "IID": str}, low_memory=False)
loo = pd.read_csv("/expanse/projects/sebat1/jsebat/pgs_asd2026_loo/asd2026_loo_scores.tsv", sep="\t", dtype={"IID": str})[["source", "IID", "asd_2026"]]
old = old.merge(loo, on=["source", "IID"], how="left").rename(columns={"asd_2026": "asd"})
pr("older table rows by source:", old.source.value_counts().to_dict(), "| asd score present:", old.groupby("source").asd.apply(lambda x: int(x.notna().sum())).to_dict())
old = old[["source", "sample_id", "case_status", "predicted_ancestry"] + [p for p in PGS] + [f"PC{i}" for i in range(1, 11)]].copy()
old["zhang_era"] = old.source.str.startswith("spark")
# ---- new release
meta = pd.read_csv(f"{J}/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv", sep="\t", dtype=str)
anc = pd.read_csv(f"{J}/plink/pgs/analyses/scoring-ancestry-20261012/06_pca/ancestry/ancestry_probabilities.tsv", sep="\t", dtype={"#IID": str}).rename(columns={"#IID": "IID"})
D = f"{J}/plink/pgs/analyses/scoring-ancestry-20261012/05_scores"
new = meta[["sample_id", "ind_id", "asd"]].rename(columns={"sample_id": "IID", "ind_id": "sp", "asd": "asd_status"})
for p in PGS:
    if p == "asd":   # Grove April 2026 FULL GWAS (these participants are not in it); NOT the lab's asd_loSSC (older GWAS). EUR only (score_wgs2026_asd.sbatch)
        sc = pd.read_csv("/expanse/projects/sebat1/jsebat/pgs_asd2026_loo/scores/spark_wgs_2026_08.asd2026_full.sscore", sep="\t", dtype={"#IID": str}).rename(columns={"#IID": "IID"})
    else:
        sc = pd.read_csv(f"{D}/{p}.sscore", sep="\t", dtype={"#IID": str}).rename(columns={"#IID": "IID"})
    new = new.merge(sc[["IID", "SCORE1_AVG"]].rename(columns={"SCORE1_AVG": p}), on="IID", how="left")
new = new.merge(anc[["IID", "ANCESTRY"]], on="IID", how="left")
new = pd.DataFrame(dict(source="spark_wgs_2026_08", sample_id=new.sp, case_status=np.where(new.asd_status.str.lower().eq("true"), "case", "control"),
                        predicted_ancestry=new.ANCESTRY, **{p: new[p] for p in PGS}))
pr("new release scored:", len(new), "| ancestry:", new.predicted_ancestry.value_counts().to_dict())
new["zhang_era"] = False
in_old = new.sample_id.isin(old.loc[old.source.str.startswith("spark"), "sample_id"])
pr("new-release individuals also in older table (kept with older score):", int(in_old.sum()))
new = new[~in_old]
A = pd.concat([old, new], ignore_index=True)
# standardise within source to non-autistic relatives
for p in PGS:
    ref = A[~A.case_status.eq("case")].groupby("source")[p].agg(["mean", "std"])
    A[p] = (A[p] - A.source.map(ref["mean"])) / A.source.map(ref["std"])
# one row per SPARK individual: prefer WGS v1.1 > Twist > GSA v2 > GSA v1 among older sources
pri = {"spark_wgs": 0, "spark_twist": 1, "spark_gsa_v2": 2, "spark_gsa_v1": 3, "spark_wgs_2026_08": 4, "ssc": 0}
A["pri"] = A.source.map(pri); A = A.sort_values("pri").drop_duplicates(["sample_id"]).drop(columns="pri")
A["cohort"] = np.where(A.source.eq("ssc"), "SSC", "SPARK")
genotyped = set(A.sample_id)
# ---- pedigree, parents
ped = pd.read_csv(f"{P}/pheno/spark_pedigree_ids.tsv", sep="\t", dtype=str).drop_duplicates("subject_sp_id").set_index("subject_sp_id")
seq = set()
for f in [f"{S}/SPARK_iWES_v3/SPARK_iWES_v3.ped", f"{S}/SPARK_iWGS_v1.1/SPARK_iWGS_v1.1.ped"]:
    seq |= set(pd.read_csv(f, sep=r"\s+", header=None, dtype=str)[1])
seq |= set(meta.ind_id)
# SSC: Sanders roles are authoritative (CLAUDE.md); SS_RUID matches the PGS sample_id. relation: fa, mo, p1 (proband), s1.. (siblings)
ssc = pd.read_csv(f"{P}/pheno/ssc_ids_sanders.csv", dtype=str, encoding="utf-8-sig").drop_duplicates("SS_RUID")
sm = ssc.set_index("SS_RUID")
sp = A.cohort.eq("SPARK")
A.loc[sp, "mother"] = A.loc[sp, "sample_id"].map(ped.biomother_sp_id); A.loc[sp, "father"] = A.loc[sp, "sample_id"].map(ped.biofather_sp_id)
A.loc[sp, "family"] = A.loc[sp, "sample_id"].map(ped.family_sf_id)
A.loc[~sp, "family"] = A.loc[~sp, "sample_id"].map(sm.family); A.loc[~sp, "role"] = A.loc[~sp, "sample_id"].map(sm.relation)
fam_par = ssc[ssc.relation.isin(["fa", "mo"])].drop_duplicates(["family", "relation"]).pivot(index="family", columns="relation", values="SS_RUID")
A.loc[~sp, "mother"] = A.loc[~sp, "family"].map(fam_par["mo"]); A.loc[~sp, "father"] = A.loc[~sp, "family"].map(fam_par["fa"])
A.loc[~sp, "case_status"] = np.where(A.loc[~sp, "role"].astype(str).str.match(r"^p\d"), "case", np.where(A.loc[~sp, "role"].notna(), "control", A.loc[~sp, "case_status"]))
pr("SSC PGS rows matched to Sanders:", int(A.loc[~sp, "role"].notna().sum()), "of", int((~sp).sum()), "| roles:", A.loc[~sp, "role"].value_counts().head(8).to_dict())
for r in ["mother", "father"]:
    A[f"{r}_genotyped"] = A[r].isin(genotyped); A[f"{r}_sequenced"] = A[r].isin(seq) | (A.cohort.eq("SSC") & A[r].notna())
A["n_parents_genotyped"] = A.mother_genotyped.astype(int) + A.father_genotyped.astype(int)
A["n_parents_sequenced"] = A.mother_sequenced.astype(int) + A.father_sequenced.astype(int)
A["sequenced"] = A.sample_id.isin(seq) | A.cohort.eq("SSC")
# ---- rare classes, characteristics
ac = pd.read_parquet(f"{P}/assembled_children.parquet")
A["rare_classes"] = A.sample_id.isin(set(ac.ind_id.astype(str)))
reg = pd.read_csv(f"{P}/pheno/spark_registration.tsv", sep="\t", dtype={"subject_sp_id": str}).drop_duplicates("subject_sp_id").set_index("subject_sp_id")
A["dx_age_known"] = A.sample_id.map(pd.to_numeric(reg.diagnosis_age, errors="coerce")).notna()
A["case"] = A.case_status.eq("case") if "case_status" in A else False
A.to_parquet(f"{OUT}/full_probands.parquet", index=False)
# ---- counts
C = A[A.case].copy(); C["EUR"] = C.predicted_ancestry.eq("EUR")
g = C.groupby(["cohort", "source"])
T = pd.DataFrame(dict(autistic_with_PGS=g.size(), EUR=g.EUR.sum(), zhang_era=g.zhang_era.sum(),
                      parents_genotyped_0=g.n_parents_genotyped.apply(lambda x: int((x == 0).sum())), parents_genotyped_1=g.n_parents_genotyped.apply(lambda x: int((x == 1).sum())),
                      parents_genotyped_2=g.n_parents_genotyped.apply(lambda x: int((x == 2).sum())), sequenced=g.sequenced.sum(),
                      parents_sequenced_2=g.n_parents_sequenced.apply(lambda x: int((x == 2).sum())), rare_classes=g.rare_classes.sum(), dx_age_known=g.dx_age_known.sum())).reset_index()
T.loc[len(T)] = ["total", ""] + list(T.iloc[:, 2:].sum())
T.to_csv(f"{OUT}/full_probands_counts.tsv", sep="\t", index=False)
pr("\n== autistic participants with PGS (all ancestries) =="); pr(T.to_string(index=False))
E = C[C.EUR]; pr("\nEUR autistic by parents genotyped (SPARK):", E[E.cohort.eq("SPARK")].n_parents_genotyped.value_counts().sort_index().to_dict())
pr("non-autistic relatives with PGS by cohort:", A[~A.case].groupby("cohort").size().to_dict())
pr("unaffected siblings in rare table:", int((ac.case == 0).sum()), "| autistic in rare table:", int((ac.case == 1).sum()), "| matched to full table:", int(A.rare_classes.sum()))
open(f"{OUT}/full_probands_counts.txt", "w").write(buf.getvalue()); print(buf.getvalue())
