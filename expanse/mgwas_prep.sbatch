#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 05:00:00
#SBATCH --array=1-3
#SBATCH --job-name=mgwas_prep
#SBATCH -o /home/jsebat/rvh2/logs/mgwas_prep_%A_%a.out
# Matched-design GWAS, step 1 per platform: EUR samples, HapMap3 SNPs (SNPs of the PGS-pipeline .ma files) with MAF >= 1% and
# missingness <= 5%; KING unrelated set (3rd degree); PCA (10) computed IN the unrelated set (no projection). Also a bed copy
# of all EUR samples (families intact) for the trio TDT.
set -e
PL=(spark_gsa_v1:pgs_wf_spark_gsa_v1/03_missingness/cohort spark_gsa_v2:pgs_wf_spark_gsa_v2/03_missingness/cohort spark_twist:pgs_wf_spark_twist/03_missingness/cohort)
IFS=: read -r SRC PF <<< "${PL[$((SLURM_ARRAY_TASK_ID-1))]}"
B=/expanse/projects/sebat1/s3/data/sebat; P2=/home/j3guevar/bin/plink2; W=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/$SRC; mkdir -p $W; cd $W
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python - <<EOF
import pandas as pd
d = pd.read_csv("$B/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "predicted_ancestry"])
e = set(d[(d.source == "$SRC") & (d.predicted_ancestry == "EUR")].IID)
ps = pd.read_csv("$B/$PF.psam", sep="\t", dtype=str); ic = [c for c in ps.columns if "IID" in c][0]
ps[ps[ic].isin(e)][[ic]].to_csv("eur_keep.txt", sep="\t", index=False, header=False)
pd.read_csv("$B/pgs_pipeline/input/asd_ipsych_dbs.ma", sep="\t", usecols=["SNP"]).SNP.to_csv("hm3.txt", index=False, header=False)
print("EUR", len(e))
EOF
printf "6 25000000 34000000 MHC\n8 8000000 12000000 inv8p23\n" > longld.txt
$P2 --pfile $B/$PF --keep eur_keep.txt --extract hm3.txt --autosome --snps-only just-acgt --maf 0.01 --geno 0.05 --make-pgen --threads 16 --memory 88000 --out hm3
$P2 --pfile hm3 --maf 0.05 --exclude bed1 longld.txt --indep-pairwise 1000kb 1 0.1 --threads 16 --memory 88000 --out prune
$P2 --pfile hm3 --extract prune.prune.in --king-cutoff 0.0442 --threads 16 --memory 88000 --out king
$P2 --pfile hm3 --keep king.king.cutoff.in.id --extract prune.prune.in --pca 10 approx --threads 16 --memory 88000 --out pca_unrel
$P2 --pfile hm3 --make-bed --threads 16 --memory 88000 --out hm3_bed
echo DONE
