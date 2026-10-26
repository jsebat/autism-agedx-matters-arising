#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH -o /home/jsebat/rvh2/logs/score_wgs2026_asd_%j.out
# Autism PGS for SPARK-WGS-2026_08 with the Grove April 2026 FULL GWAS weights (these participants are not in the GWAS).
# Replaces the lab's asd_loSSC (older GWAS) for this cohort so that every SPARK subcohort uses the 2026 Grove summary statistics.
# Genotypes: EUR, rsID-renamed pfile from allgeno_prep.sbatch (cnv_drift/mgwas/spark_wgs_2026_08/hm3). Same plink2 flags as score_asd2026_loo.sh.
set -euo pipefail
WD=/expanse/projects/sebat1/jsebat/pgs_asd2026_loo; M=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/spark_wgs_2026_08
/home/j3guevar/bin/plink2 --pfile $M/hm3 --score $WD/output/weights/asd2026_full.snpRes 2 5 8 header center list-variants no-mean-imputation \
  --threads 4 --memory 28000 --out $WD/scores/spark_wgs_2026_08.asd2026_full
grep -c . $WD/scores/spark_wgs_2026_08.asd2026_full.sscore.vars
echo DONE
