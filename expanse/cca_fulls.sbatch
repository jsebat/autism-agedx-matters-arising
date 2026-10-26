#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --array=1-10
#SBATCH -o /home/jsebat/rvh2/logs/cca_fulls_%A_%a.out
# Sensitivity set for the supplementary section: unadjusted, platform + family structure adjusted, and per-platform strata
export OPENBLAS_NUM_THREADS=2
case $SLURM_ARRAY_TASK_ID in
  1) export SUBSET=all PARTIAL=source_family ;;
  2) export SUBSET=u22 PARTIAL=source_family ;;
  3) export SUBSET=zhangu22 PARTIAL=source_family ;;
  4) export SUBSET=u22 PARTIAL= ;;
  5) export SUBSET=src_spark_gsa_v1 PARTIAL= ;;
  6) export SUBSET=src_spark_gsa_v2 PARTIAL= ;;
  7) export SUBSET=src_spark_twist PARTIAL= ;;
  8) export SUBSET=src_spark_wgs PARTIAL= ;;
  9) export SUBSET=src_spark_wgs_2026_08 PARTIAL= ;;
  10) export SUBSET=zhang PARTIAL=source_family ;;
esac
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/cca_full.py
