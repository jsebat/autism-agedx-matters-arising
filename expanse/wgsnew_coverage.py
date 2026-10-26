#!/usr/bin/env python3
"""New SPARK WGS release (SPARK-WGS-2026_08): how much of it reaches the age-at-diagnosis analyses, and why not.
SDSM -> SP map: rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv (45,178 rows).
PGS for the release: plink/pgs/analyses/scoring-ancestry-20261012/05_scores/*.sscore (SDSM IIDs); ancestry from 06_pca if present.
Old PGS table used by our pipeline: s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv (array + iWGS v1.1 + SSC).
Aggregate output only: ~/rvh2/agedx/wgsnew_coverage.txt"""
import os, io, glob
import pandas as pd
J = "/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08"
buf = io.StringIO(); pr = lambda *x: print(*x, file=buf, flush=True)
m = pd.read_csv(f"{J}/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv", sep="\t", dtype=str)
m["asd"] = m.asd.str.lower().eq("true")
pr("new release:", len(m), "samples;", m.ind_id.nunique(), "SP IDs;", m.family_id.nunique(), "families; autistic", int(m.asd.sum()))
pr("  with both parents sequenced:", int((m.father_sequenced.eq("True") & m.mother_sequenced.eq("True")).sum()), "| autistic with both parents sequenced:", int((m.asd & m.father_sequenced.eq("True") & m.mother_sequenced.eq("True")).sum()))
sc = pd.read_csv(f"{J}/plink/pgs/analyses/scoring-ancestry-20261012/05_scores/edu_2022.sscore", sep="\t", dtype={"#IID": str})
pr("Lab PGS (edu_2022) scored:", len(sc), "| in metadata:", int(sc["#IID"].isin(m.sample_id).sum()))
pca = glob.glob(f"{J}/plink/pgs/analyses/scoring-ancestry-20261012/06_pca/*")
pr("06_pca files:", [os.path.basename(f) for f in pca][:12])
for f in pca:
    if f.endswith((".tsv", ".txt", ".csv")):
        try:
            t = pd.read_csv(f, sep=None, engine="python", nrows=5); pr("  ", os.path.basename(f), list(t.columns)[:12])
        except Exception: pass
old = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, usecols=["source", "sample_id", "predicted_ancestry"])
old = old[old.source.str.startswith("spark")]
inold = m.ind_id.isin(old.sample_id)
pr("\nnew-release individuals already in the old PGS table (array or iWGS v1.1):", int(inold.sum()), f"({100 * inold.mean():.1f}%)")
pr("  autistic:", int((inold & m.asd).sum()), "of", int(m.asd.sum()))
pr("  NOT in old table:", int((~inold).sum()), f"({100 * (~inold).mean():.1f}%); autistic {int((~inold & m.asd).sum())}")
anc = old.drop_duplicates("sample_id").set_index("sample_id").predicted_ancestry
pr("  ancestry of those in old table:", m.loc[inold, "ind_id"].map(anc).value_counts().to_dict())
open(os.path.expanduser("~/rvh2/agedx/wgsnew_coverage.txt"), "w").write(buf.getvalue()); print(buf.getvalue())
