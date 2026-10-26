#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --array=1-3
#SBATCH -o /home/jsebat/rvh2/logs/ccassc_%A_%a.out
# SSC replication CCA (pgs; pgs + rare) and SPARK on the same harmonised characteristics
export OPENBLAS_NUM_THREADS=2
case $SLURM_ARRAY_TASK_ID in
  1) export COHORT=ssc YSET=pgs ;;
  2) export COHORT=ssc YSET=all ;;
  3) export COHORT=spark YSET=pgs ;;
esac
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_ssc.py
