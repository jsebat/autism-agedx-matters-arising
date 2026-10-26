#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 03:00:00
#SBATCH --array=1-8
#SBATCH --job-name=ext_rg
#SBATCH -o /home/jsebat/rvh2/logs/ext_rg_%A_%a.out
# External-trait profile of every survey half-GWAS (and the replicated age-split GWAS): LDSC rg with age at first birth, educational
# attainment, IQ, ADHD, MDD, SCZ. Array task k processes every 8th half. Output: warrier_rep/combined/rg_ext/<half>.log
set -e
export OPENBLAS_NUM_THREADS=1
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/
C=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined; cd $C; mkdir -p rg_ext
EXT=munged/ext_afb_2018.sumstats.gz,munged/ext_edu_2022.sumstats.gz,munged/ext_iq_2018.sumstats.gz,munged/ext_adhd_2022.sumstats.gz,munged/ext_mdd_2025.sumstats.gz,munged/ext_scz_2022.sumstats.gz
ls munged/rep_s*_hi_c2.sumstats.gz munged/rep_s*_lo_c1.sumstats.gz munged/W_F1.sumstats.gz munged/W_F2.sumstats.gz | sort > rg_ext_list.txt
i=0; while read f; do i=$((i+1)); [ $(( (i-1) % 8 + 1 )) -eq $SLURM_ARRAY_TASK_ID ] || continue
  b=$(basename $f .sumstats.gz); [ -s rg_ext/$b.log ] && grep -q "Summary of Genetic" rg_ext/$b.log && continue
  $PY $L/ldsc.py --rg $f,$EXT --ref-ld-chr $REF --w-ld-chr $REF --out rg_ext/$b > /dev/null || echo "failed $b"
done < rg_ext_list.txt
echo DONE
