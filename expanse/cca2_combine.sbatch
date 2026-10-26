#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --array=1-26
#SBATCH --job-name=cca2_comb
#SBATCH -o /home/jsebat/rvh2/logs/cca2_comb_%A_%a.out
export OPENBLAS_NUM_THREADS=1
export $(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/rvh2/cca2_configs.cfg) TASK=combine
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca2.py
