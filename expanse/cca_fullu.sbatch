#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --array=1-3
#SBATCH -o /home/jsebat/rvh2/logs/cca_fullu_%A_%a.out
# Polygenic-score CCA on the full SPARK sample and sensitivity subsets (Zhang-era, new release only, 0/1/2 parents genotyped)
export OPENBLAS_NUM_THREADS=2 PARTIAL=source
SUBS=(x u22 p0u22 zhangu22)
export SUBSET=${SUBS[$SLURM_ARRAY_TASK_ID]}
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_full.py
