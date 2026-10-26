#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --job-name=ccatriosgl
#SBATCH -o /home/jsebat/rvh2/logs/ccatriosingle_%j.out
# Strength of age at diagnosis alone in the combined rare + polygenic trio sample (registered before 22; platform, sequencing cohort, family structure adjusted)
export OPENBLAS_NUM_THREADS=2 YSET=all U22=1 SINGLE=1
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_trio.py
