import glob, os, pandas as pd, numpy as np
C = "/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/combined"
for f in sorted(glob.glob(f"{C}/*sumstats.txt")):
    x = pd.read_csv(f, sep="\t"); z = x.Z
    print(f"{os.path.basename(f):35s} n={len(x):8d} medianZ={z.median():+.3f} meanZ={z.mean():+.3f} lambda={np.median(z**2)/0.4549:.3f} medN={x.N.median():.0f}")
for s in ["spark_gsa_v1", "spark_gsa_v2", "spark_twist"]:
    f = f"/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/{s}/tdt_agedx_early.tdt"
    if os.path.exists(f):
        x = pd.read_csv(f, sep=r"\s+"); x = x[(x["T"] + x["U"]) > 0]; z = (x["T"] - x["U"]) / np.sqrt(x["T"] + x["U"])
        print(s, "TDT agedx early: median Z", round(z.median(), 3), "mean T-U", round((x["T"] - x["U"]).mean(), 2), "A1 freq-like: share T<U", round((x["T"] < x["U"]).mean(), 3))
