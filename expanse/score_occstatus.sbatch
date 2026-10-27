#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --array=1-6
#SBATCH --job-name=score_occ
#SBATCH -o /home/jsebat/rvh2/logs/score_occ_%A_%a.out
# SES polygenic score: occupational status (CAMSIS; Akimova et al. 2025, Nat Hum Behav, doi:10.1038/s41562-024-02076-3; UK Biobank, N = 271,197),
# SBayesRC weights from the lab collection (a1sriniv/pgs_sumstats_organized/processed/txt_weights/occupational_status_cojo.txt: SNP A1 BETA ...).
# Same genotype sets and plink2 flags as the autism scores (score_asd2026_loo.sh, score_wgs2026_asd.sbatch); SPARK and SSC are not in UK Biobank.
set -euo pipefail
B=/expanse/projects/sebat1/s3/data/sebat; W=/expanse/projects/sebat1/a1sriniv/pgs_sumstats_organized/processed/txt_weights/occupational_status_cojo.txt
OUT=/expanse/projects/sebat1/jsebat/pgs_occstatus; mkdir -p $OUT
JOBS=(spark_gsa_v1:$B/pgs_wf_spark_gsa_v1/03_missingness/cohort spark_gsa_v2:$B/pgs_wf_spark_gsa_v2/03_missingness/cohort
      spark_twist:$B/pgs_wf_spark_twist/03_missingness/cohort ssc:$B/pgs_wf_ssc/05_summary/ssc spark_wgs:$B/pgs_wf_spark_wgs/05_summary/spark_wgs
      spark_wgs_2026_08:/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/spark_wgs_2026_08/hm3)
IFS=: read -r SET PF <<< "${JOBS[$((SLURM_ARRAY_TASK_ID-1))]}"
echo "[$(date)] scoring $SET ($PF)"
/home/j3guevar/bin/plink2 --pfile $PF --score $W 1 2 3 header center list-variants no-mean-imputation --threads 4 --memory 28000 --out $OUT/${SET}.occstatus
echo "variants used: $(grep -c . $OUT/${SET}.occstatus.sscore.vars)"
echo "[$(date)] DONE"
