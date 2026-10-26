#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --array=1-2
#SBATCH -o /home/jsebat/rvh2/logs/rarefull_%A_%a.out
# Expanded-trio rare-variant comparison: 1 = registered before 22 (main, Zhang criterion), 2 = all ages (supplement)
export OPENBLAS_NUM_THREADS=2
[ "$SLURM_ARRAY_TASK_ID" = "1" ] && export U22=1 || export U22=0
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/rare_full.py
