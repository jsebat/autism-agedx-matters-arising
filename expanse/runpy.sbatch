#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -o /home/jsebat/rvh2/logs/runpy_%x_%j.out
# Run one analysis script with the caller's environment (sbatch --export=ALL,VAR=... runpy.sbatch script.py)
export OPENBLAS_NUM_THREADS=2
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/$1
