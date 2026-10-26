import glob, io, os, re, pandas as pd
C = "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined"; rows = []
for f in glob.glob(f"{C}/rg_ext/*.log"):
    t = open(f).read(); me = os.path.basename(f)[:-4]
    if "Summary of Genetic Correlation Results" not in t: continue
    blk = t.split("Summary of Genetic Correlation Results")[1].strip().split("\n\n")[0]; x = pd.read_csv(io.StringIO(blk), sep=r"\s+")
    for r in x.itertuples(): rows.append(dict(gwas=me, trait=os.path.basename(r.p2).replace(".sumstats.gz", "").replace("ext_", ""), rg=r.rg, se=r.se))
D = pd.DataFrame(rows); D.to_csv(f"{C}/rg_ext_all.tsv", sep="\t", index=False); print(len(D), D.gwas.nunique())
