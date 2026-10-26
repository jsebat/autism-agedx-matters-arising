#!/usr/bin/env python3
"""Why do the spark_wgs_new probands lack polygenic scores ? Checks:
 1. ID format of spark_wgs_new samples in the rare-variant table vs the other cohorts.
 2. Sources and ID format in pgs_results_v2 (the PGS table).
 3. Any SDSM -> SP ID map now present on Expanse (search the lab's SPARK-WGS-2026_08 tree for manifest/metadata files).
 4. If a map exists: of the 45,178 new-WGS samples, how many have a PGS on any platform (and EUR).
Aggregate output only: ~/rvh2/agedx/wgsnew_pgs_check.txt"""
import os, io, glob, re
import pandas as pd
buf = io.StringIO(); pr = lambda *x: print(*x, file=buf, flush=True)
src = open(os.path.expanduser("~/rvh2/rare_survey.py")).read(); exec(src[src.index("TI = "):src.index("def excess")])
a = load("pgs_rare_h2_drift_full", {"x": CAT["de novo, NDD genes"]})
for c, g in a.groupby("cohort"): pr(f"rare table {c}: n={len(g)}; example sample IDs {list(g['sample'].astype(str).head(3))}; family_key {list(g.family_key.astype(str).head(2))}")
d = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, usecols=lambda c: c in ("source", "sample_id", "FID", "IID", "predicted_ancestry", "case_status", "edu_2022"))
pr("\nPGS table sources:", d.source.value_counts().to_dict())
for s, g in d.groupby("source"): pr(f"  {s}: example sample_id {list(g.sample_id.head(2))} IID {list(g.IID.head(2))}")
# search for an SDSM -> SP map
root = "/expanse/projects/sebat1"
cands = []
for pat in ["*/SPARK-WGS-2026_08/*", "*/SPARK-WGS-2026_08/*/*", "j3guevar/*/SPARK-WGS-2026_08*", "s3/data/sebat/SPARK*2026*/*", "s3/data/sebat/SPARK*2026*/*/*"]:
    for f in glob.glob(f"{root}/{pat}"):
        if re.search(r"manifest|metadata|sample.?map|id.?map|crosswalk|key|ped|subject|barcode", os.path.basename(f), re.I) and os.path.isfile(f): cands.append(f)
pr("\ncandidate ID-map files:", len(cands)); [pr("  ", f) for f in sorted(set(cands))[:40]]
newids = set(a.loc[a.cohort.eq("spark_wgs_new"), "sample"].astype(str))
for f in sorted(set(cands))[:40]:
    try:
        t = pd.read_csv(f, sep=None, engine="python", dtype=str, nrows=200000)
    except Exception as e:
        continue
    cols = list(t.columns); hit = [c for c in cols if t[c].astype(str).isin(newids).any()]
    spc = [c for c in cols if t[c].astype(str).str.match(r"^SP\d{7}$").mean() > 0.5]
    if hit or spc: pr(f"  {f}: rows {len(t)}; cols matching rare-table new IDs {hit}; SP-ID cols {spc}")
open(os.path.expanduser("~/rvh2/agedx/wgsnew_pgs_check.txt"), "w").write(buf.getvalue()); print(buf.getvalue())
