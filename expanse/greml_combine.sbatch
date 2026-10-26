#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH -o /home/jsebat/rvh2/logs/greml_combine_%j.out
export OPENBLAS_NUM_THREADS=1
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/greml_combine_all.py
