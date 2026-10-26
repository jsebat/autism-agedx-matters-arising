#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --array=1-3
#SBATCH --job-name=mgwas_sib
#SBATCH -o /home/jsebat/rvh2/logs/mgwas_sib_%A_%a.out
# Sibling-null TDT: transmissions from the same trio parents to NON-autistic children (coded affected for plink --tdt), per
# platform; gives the per-SNP null transmission rate that absorbs genotyping/imputation transmission bias.
set -e
SRCS=(spark_gsa_v1 spark_gsa_v2 spark_twist); SRC=${SRCS[$((SLURM_ARRAY_TASK_ID-1))]}
W=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/$SRC; cd $W
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python - <<PY
import pandas as pd, numpy as np
d = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "case_status"])
d = d[d.source == "$SRC"].drop_duplicates("IID"); nonaut = set(d[d.case_status != "case"].IID)
fam = pd.read_csv("trio.fam", sep=" ", header=None, dtype=str); fam.columns = ["FID", "IID", "PAT", "MAT", "SEX", "PHE"]
kid = (fam.PAT != "0") & (fam.MAT != "0") & fam.IID.isin(nonaut)
tp = fam[["FID", "IID"]].copy(); tp["Y"] = np.where(kid, 2, 1); tp.to_csv("tdt_sibnull.pheno", sep=" ", index=False, header=False)
print("non-autistic trio offspring:", int(kid.sum()))
PY
/expanse/projects/sebat1/tools/plink1.90b6.15 --bed hm3_bed.bed --bim hm3_bed.bim --fam trio.fam --extract /expanse/projects/sebat1/jsebat/tools/ldref/w_hm3.snplist --pheno tdt_sibnull.pheno --tdt --threads 16 --memory 44000 --out tdt_sibnull --allow-no-sex
echo DONE
