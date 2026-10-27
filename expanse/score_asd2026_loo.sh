#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:30:00
#SBATCH --array=1-7
#SBATCH -o /expanse/projects/sebat1/jsebat/pgs_asd2026_loo/logs/score_%A_%a.out
#SBATCH -e /expanse/projects/sebat1/jsebat/pgs_asd2026_loo/logs/score_%A_%a.err
#SBATCH --job-name=score_asdloo
# Score the April 2026 ASD SBayesRC weights on each genotype set with the matching leave-one-out weights
# (README: /expanse/projects/sebat1/jsebat/sumstats/README.md). iWGS is scored with full, W1-4 and W5 weights so the right
# score can be chosen per person (iWGS-only -> full; also on GSA -> W1-4; also on Twist -> W5).
# Same pfiles and plink2 flags as the lab's asd_ipsych_dbs / asd_loSSC scoring.
set -euo pipefail
B=/expanse/projects/sebat1/s3/data/sebat; WD=/expanse/projects/sebat1/jsebat/pgs_asd2026_loo
JOBS=(spark_gsa_v1:pgs_wf_spark_gsa_v1/03_missingness/cohort:asd2026_loSPARKW1to4
      spark_gsa_v2:pgs_wf_spark_gsa_v2/03_missingness/cohort:asd2026_loSPARKW1to4
      spark_twist:pgs_wf_spark_twist/03_missingness/cohort:asd2026_loSPARKW5
      ssc:pgs_wf_ssc/05_summary/ssc:asd2026_loSSC
      spark_wgs:pgs_wf_spark_wgs/05_summary/spark_wgs:asd2026_full
      spark_wgs:pgs_wf_spark_wgs/05_summary/spark_wgs:asd2026_loSPARKW1to4
      spark_wgs:pgs_wf_spark_wgs/05_summary/spark_wgs:asd2026_loSPARKW5)
IFS=: read -r SET PF W <<< "${JOBS[$((SLURM_ARRAY_TASK_ID-1))]}"
echo "[$(date)] scoring $SET ($B/$PF) with $W"
/home/j3guevar/bin/plink2 --pfile $B/$PF --score $WD/output/weights/$W.snpRes 2 5 8 header center list-variants no-mean-imputation \
  --threads 4 --memory 28000 --out $WD/scores/${SET}.${W}
echo "[$(date)] DONE"
