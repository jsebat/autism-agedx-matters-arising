#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --array=1-3
#SBATCH --job-name=mgwas_run
#SBATCH -o /home/jsebat/rvh2/logs/mgwas_run_%A_%a.out
# Matched-design GWAS, step 3 per platform: between-family logistic GWAS (plink2, Firth fallback, PC1-10 + sex) and trio TDT
# (plink 1.9 --tdt) for 4 splits x 2 halves.
set -e
SRCS=(spark_gsa_v1 spark_gsa_v2 spark_twist); export SRC=${SRCS[$((SLURM_ARRAY_TASK_ID-1))]}
W=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/$SRC; cd $W
P2=/home/j3guevar/bin/plink2; P1=/expanse/projects/sebat1/tools/plink1.90b6.15
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/mgwas_pheno.py
for f in pheno_*_*.txt; do t=${f#pheno_}; t=${t%.txt}
  $P2 --pfile hm3 --extract /expanse/projects/sebat1/jsebat/tools/ldref/w_hm3.snplist --pheno $f --pheno-name Y --covar covar.txt --covar-variance-standardize --glm firth-fallback hide-covar cols=+a1freq,+nobs \
      --threads 16 --memory 56000 --out bf_$t
done
for f in tdt_*_*.pheno; do t=${f#tdt_}; t=${t%.pheno}
  $P1 --bed hm3_bed.bed --bim hm3_bed.bim --fam trio.fam --extract /expanse/projects/sebat1/jsebat/tools/ldref/w_hm3.snplist --pheno $f --tdt --threads 16 --memory 56000 --out tdt_$t --allow-no-sex
done
echo DONE
