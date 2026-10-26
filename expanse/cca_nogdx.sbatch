#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --array=1-3
#SBATCH -o /home/jsebat/rvh2/logs/cca_nogdx_%A_%a.out
# CCA re-run without "genetic diagnosis" among the characteristics (reported from the same sequencing: circular for rare variants)
export DROP="genetic diagnosis" OPENBLAS_NUM_THREADS=1
case $SLURM_ARRAY_TASK_ID in
  1) export YSET=all ;;
  2) export YSET=pgs ;;
  3) export YSET=rare COMMON=1 ;;
esac
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_axes.py
