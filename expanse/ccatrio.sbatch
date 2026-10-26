#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --array=1-4
#SBATCH -o /home/jsebat/rvh2/logs/ccatrio_%A_%a.out
# Combined CCA on expanded trio probands: 1-3 = all/pgs/rare registered before 22 (main); 4 = all blocks, all ages (supplement)
export OPENBLAS_NUM_THREADS=2
case $SLURM_ARRAY_TASK_ID in
  1) export YSET=all U22=1 ;;
  2) export YSET=pgs U22=1 ;;
  3) export YSET=rare U22=1 ;;
  4) export YSET=all U22=0 ;;
esac
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_trio.py
