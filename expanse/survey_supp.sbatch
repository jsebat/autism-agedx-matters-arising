#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --array=1-3
#SBATCH --job-name=survey_supp
#SBATCH -o /home/jsebat/rvh2/logs/survey_supp_%A_%a.out
# Random-split null GWAS (their pipeline) per platform: survey_pheno.py then fastGWA-mlm-binary for each random set.
set -e
SRCS=(spark_gsa_v1 spark_gsa_v2 spark_twist); export SRC=${SRCS[$((SLURM_ARRAY_TASK_ID-1))]}
W=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/$SRC; cd $W; G=/expanse/projects/sebat1/jsebat/tools/gcta/gcta64
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/survey_pheno.py
while read nm; do k=$(grep -n -x "$nm" survey_names.txt | cut -d: -f1)
  $G --bfile hm3w --grm-sparse grm_sp --fastGWA-mlm-binary --pheno survey_pheno.txt --mpheno $k --qcovar qcovar.txt --covar covar.txt --maf 0.01 --thread-num 16 --out fg_$nm > fg_$nm.stdout 2>&1 || echo "failed $nm"
done < <(grep -E "^s4[0-2]_|^s99_" survey_names.txt)
echo DONE
