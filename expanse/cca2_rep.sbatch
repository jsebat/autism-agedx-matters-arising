#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 06:00:00
#SBATCH --job-name=cca2_rep
#SBATCH -o /home/jsebat/rvh2/logs/cca2_rep_%A_%a.out
# Permutation (TASK=null) or family-bootstrap (TASK=boot) chunk for config CFG (line of cca2_configs.cfg); CHUNK = array index, NCHUNK set at submission
export OPENBLAS_NUM_THREADS=2
export $(sed -n "${CFG}p" ~/rvh2/cca2_configs.cfg) CHUNK=$SLURM_ARRAY_TASK_ID
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca2.py
