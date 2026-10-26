#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 10:00:00
#SBATCH --array=0-17
#SBATCH --job-name=allstrata_gwas
#SBATCH -o /home/jsebat/rvh2/logs/allstrata_gwas_%A_%a.out
# fastGWA-mlm-binary (their pipeline) for every split GWAS of the new strata; 3 strata x 6 chunks of allpheno_names.txt.
set -e
SRCS=(spark_wgs spark_wgs_2026_08 ssc); SRC=${SRCS[$((SLURM_ARRAY_TASK_ID / 6))]}; CH=$((SLURM_ARRAY_TASK_ID % 6))
W=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/$SRC; cd $W; G=/expanse/projects/sebat1/jsebat/tools/gcta/gcta64
k=0; while read nm; do k=$((k+1))
  [ $(( (k-1) % 6 )) -eq $CH ] || continue
  [ -s fg_$nm.fastGWA ] && continue
  $G --bfile hm3w --grm-sparse grm_sp --fastGWA-mlm-binary --pheno allpheno.txt --mpheno $k --qcovar qcovar.txt --covar covar.txt --maf 0.01 --thread-num 16 --out fg_$nm > fg_$nm.stdout 2>&1 || echo "failed $nm"
done < allpheno_names.txt
echo DONE $SRC chunk $CH
