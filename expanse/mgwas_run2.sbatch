#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --array=1-3
#SBATCH --job-name=mgwas_run2
#SBATCH -o /home/jsebat/rvh2/logs/mgwas_run2_%A_%a.out
# Matched-design GWAS round 2 per platform: between-family logistic GWAS for the orthogonal splits (mgwas_pheno2.py).
set -e
SRCS=(spark_gsa_v1 spark_gsa_v2 spark_twist); export SRC=${SRCS[$((SLURM_ARRAY_TASK_ID-1))]}
W=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/$SRC; cd $W
P2=/home/j3guevar/bin/plink2
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/mgwas_pheno2.py
for f in pheno2_*_*.txt; do t=${f#pheno2_}; t=${t%.txt}
  $P2 --pfile hm3 --extract /expanse/projects/sebat1/jsebat/tools/ldref/w_hm3.snplist --pheno $f --pheno-name Y --covar covar.txt --covar-variance-standardize --glm firth-fallback hide-covar cols=+a1freq,+nobs \
      --threads 16 --memory 56000 --out bf_$t
done
echo DONE
