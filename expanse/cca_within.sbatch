#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --array=1-2
#SBATCH -o /home/jsebat/rvh2/logs/ccawithin_%A_%a.out
# Within-family (pTDT) polygenic CCA on the full SPARK data: 1 = registered before 22 (main), 2 = all ages
export OPENBLAS_NUM_THREADS=2 WITHIN=1 PARTIAL=source
[ "$SLURM_ARRAY_TASK_ID" = "1" ] && export SUBSET=u22 || export SUBSET=all
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_full.py
