#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 04:00:00
#SBATCH --array=1-26%10
#SBATCH --job-name=cca2_main
#SBATCH -o /home/jsebat/rvh2/logs/cca2_main_%A_%a.out
# Held-out CCA (cca2.py / ccalib.py): observed fits for every configuration in cca2_configs.cfg (line = array index)
export OPENBLAS_NUM_THREADS=2
export $(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/rvh2/cca2_configs.cfg) TASK=main
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca2.py
