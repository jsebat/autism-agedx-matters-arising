import pandas as pd
d = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "FID", "IID", "sample_id", "predicted_ancestry"])
e = d[(d.source == "spark_twist") & (d.predicted_ancestry == "EUR")].drop_duplicates("IID")
psam = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_wf_spark_twist/03_missingness/cohort.psam", sep="\t", dtype=str)
print(psam.columns.tolist(), psam.head(2).to_string())
idc = [c for c in psam.columns if "IID" in c][0]; fidc = [c for c in psam.columns if "FID" in c]
keep = psam[psam[idc].isin(set(e.IID))]
cols = (fidc + [idc]) if fidc else [idc]
keep[cols].to_csv("/expanse/projects/sebat1/jsebat/cnv_drift/pcproj/eur_keep.txt", sep="\t", index=False, header=False)
print("EUR twist samples kept:", len(keep))
