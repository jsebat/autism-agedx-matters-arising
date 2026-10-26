#!/usr/bin/env python3
"""What the diagnostic-drift project assembled (pgs_rare_h2_drift/assembled_children.parquet): rows by cohort, role/case,
coverage of PGS, rare-variant classes, age at diagnosis and SES. Aggregate output: ~/rvh2/agedx/drift_children.txt"""
import os, io
import pandas as pd
buf = io.StringIO(); pr = lambda *x: print(*x, file=buf, flush=True)
d = pd.read_parquet("/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/assembled_children.parquet")
pr("rows", len(d)); pr("columns:", list(d.columns))
key = [c for c in ["cohort", "source", "case", "asd", "role", "case_status"] if c in d.columns]
pr("key cols:", key)
for k in key[:3]: pr(f"\n{k}:", d[k].value_counts(dropna=False).to_dict())
if len(key) >= 2: pr(d.groupby(key[:2]).size().unstack(fill_value=0).to_string())
nn = d.notna().mean().round(3); pr("\nnon-missing fraction:"); pr(nn.to_string())
gcol = key[0] if key else None
if gcol:
    pr("\nnon-missing by", gcol); pr(d.groupby(gcol).apply(lambda g: g.notna().sum()).T.to_string())
for f in sorted(os.listdir("/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno")): pr("pheno/", f)
for f in sorted(os.listdir("/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/results"))[:60]: pr("results/", f)
open(os.path.expanduser("~/rvh2/agedx/drift_children.txt"), "w").write(buf.getvalue()); print(buf.getvalue())
