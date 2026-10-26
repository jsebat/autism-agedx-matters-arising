#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /home/jsebat/rvh2/logs/ccapow_%j.out
# Power analysis for the main polygenic CCA (registered before 22, platform + family structure adjusted)
export OPENBLAS_NUM_THREADS=2 SUBSET=u22 PARTIAL=source_family POWER=1
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_full.py
