#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --job-name=ccasingle
#SBATCH -o /home/jsebat/rvh2/logs/ccasingle_all_%j.out
# Strength of single characteristics (age at diagnosis and comparators) in the main polygenic CCA sample (ALL ages incl. adult registrants, platform + family structure adjusted)
export OPENBLAS_NUM_THREADS=2 SUBSET=all PARTIAL=source_family SINGLE=1
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_full.py
