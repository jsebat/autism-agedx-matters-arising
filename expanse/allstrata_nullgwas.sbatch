#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 10:00:00
#SBATCH --array=0-7
#SBATCH --job-name=allstrata_null
#SBATCH -o /home/jsebat/rvh2/logs/allstrata_null_%A_%a.out
# Random-split null GWAS for the new SPARK strata (allstrata_null.py): 2 strata x 4 chunks.
set -e
SRCS=(spark_wgs spark_wgs_2026_08); SRC=${SRCS[$((SLURM_ARRAY_TASK_ID / 4))]}; CH=$((SLURM_ARRAY_TASK_ID % 4))
W=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/$SRC; cd $W; G=/expanse/projects/sebat1/jsebat/tools/gcta/gcta64
[ $CH -eq 0 ] && SRC=$SRC /home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/allstrata_null.py
for i in $(seq 1 60); do [ -s allnull_names.txt ] && break; sleep 10; done
k=0; while read nm; do k=$((k+1))
  [ $(( (k-1) % 4 )) -eq $CH ] || continue
  [ -s fg_$nm.fastGWA ] && continue
  $G --bfile hm3w --grm-sparse grm_sp --fastGWA-mlm-binary --pheno allnull_pheno.txt --mpheno $k --qcovar qcovar.txt --covar covar.txt --maf 0.01 --thread-num 16 --out fg_$nm > fg_$nm.stdout 2>&1 || echo "failed $nm"
done < allnull_names.txt
echo DONE $SRC chunk $CH
