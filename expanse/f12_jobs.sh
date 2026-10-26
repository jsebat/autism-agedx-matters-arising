#!/bin/bash
# Factor PGS chain: convert -> SBayesRC (lab pipeline, weights) -> score GSA v1/v2/Twist -> driver analysis
export PATH=/cm/shared/apps/slurm/current/bin:$PATH
WD=/expanse/projects/sebat1/jsebat/pgs_warrier; mkdir -p $WD/logs $WD/scores $WD/output
sed -e "s#WORK_DIR=\"/expanse/projects/sebat1/jsebat/pgs_asd2026_loo\"#WORK_DIR=\"$WD\"#" -e "s#--array=1-4#--array=1-2#" -e "s#pgs_asd2026_loo/logs#pgs_warrier/logs#g" -e "s#sbrc_asdloo#sbrc_warrier#" /expanse/projects/sebat1/jsebat/pgs_asd2026_loo/run_sbayesrc.sh > $WD/run_sbayesrc.sh
cat > $WD/prep.sbatch <<EOF
#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH -o $WD/logs/prep_%j.out
export OPENBLAS_NUM_THREADS=1
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/f12_prep.py
EOF
cat > $WD/score.sbatch <<'EOF'
#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH --array=1-6
#SBATCH -o /expanse/projects/sebat1/jsebat/pgs_warrier/logs/score_%A_%a.out
# Score the factor SBayesRC weights on each SPARK genotype set (same pfiles and plink2 flags as the ASD 2026 LOO scoring).
set -euo pipefail
B=/expanse/projects/sebat1/s3/data/sebat; WD=/expanse/projects/sebat1/jsebat/pgs_warrier
JOBS=(spark_gsa_v1:pgs_wf_spark_gsa_v1/03_missingness/cohort:warrier_F1 spark_gsa_v1:pgs_wf_spark_gsa_v1/03_missingness/cohort:warrier_F2
      spark_gsa_v2:pgs_wf_spark_gsa_v2/03_missingness/cohort:warrier_F1 spark_gsa_v2:pgs_wf_spark_gsa_v2/03_missingness/cohort:warrier_F2
      spark_twist:pgs_wf_spark_twist/03_missingness/cohort:warrier_F1 spark_twist:pgs_wf_spark_twist/03_missingness/cohort:warrier_F2)
IFS=: read -r SET PF W <<< "${JOBS[$((SLURM_ARRAY_TASK_ID-1))]}"
/home/j3guevar/bin/plink2 --pfile $B/$PF --score $WD/output/weights/$W.snpRes 2 5 8 header center list-variants no-mean-imputation --threads 4 --memory 28000 --out $WD/scores/${SET}.${W}
EOF
cat > $WD/driver.sbatch <<EOF
#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -o $WD/logs/driver_%j.out
export OPENBLAS_NUM_THREADS=1
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/f12_driver.py
EOF
J1=$(sbatch --parsable $WD/prep.sbatch); J2=$(sbatch --parsable --dependency=afterok:$J1 $WD/run_sbayesrc.sh); J3=$(sbatch --parsable --dependency=afterok:$J2 $WD/score.sbatch); J4=$(sbatch --parsable --dependency=afterok:$J3 $WD/driver.sbatch)
echo prep=$J1 sbrc=$J2 score=$J3 driver=$J4; grep -n "WORK_DIR=\|array\|logs" $WD/run_sbayesrc.sh | head -5
