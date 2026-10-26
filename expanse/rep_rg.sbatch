#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH --array=1-48
#SBATCH --job-name=rep_rg
#SBATCH -o /home/jsebat/rvh2/logs/rep_rg_%A_%a.out
# Replication step 6 (array over rg_targets.txt): LDSC rg of one GWAS vs every replication GWAS, the published factors (F1, F2,
# SPARK age at diagnosis), the matched-design within-family signal and external traits.
set -e
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/
C=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined; cd $C
f=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rg_targets.txt); [ -z "$f" ] && { echo "no target"; exit 0; }
b=$(basename $f .sumstats.gz)
ALL=$(ls munged/rep_*.sumstats.gz | tr "\n" ","); OTH="munged/W_F1.sumstats.gz,munged/W_F2.sumstats.gz,munged/W_agedx.sumstats.gz,munged/sibtrio_agedx_early.sumstats.gz,"
EXT=$(ls munged/ext_*.sumstats.gz | tr "\n" ",")
$PY $L/ldsc.py --rg ${f},${ALL}${OTH}${EXT%,} --ref-ld-chr $REF --w-ld-chr $REF --out rg/$b > /dev/null
echo DONE $b
