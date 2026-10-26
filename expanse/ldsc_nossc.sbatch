#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --job-name=ldsc_nossc
#SBATCH -o /home/jsebat/rvh2/logs/ldsc_nossc_%j.out
# Review #5: SES splits (maternal age s11, income s09, mother's education s10) on SPARK strata only (SSC excluded, as in the age split), same pipeline
set -e
export OPENBLAS_NUM_THREADS=1
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/; HM=$T/ldref/w_hm3.snplist
EXCLUDE=ssc ONLY='^s(09|10|11)_' /home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/combine_all.py
C=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined_all_nossc; cd $C; mkdir -p munged rg
for f in rep_s0[9]_*_c[12].sumstats.txt rep_s1[01]_*_c[12].sumstats.txt; do b=${f%.sumstats.txt}; $PY $L/munge_sumstats.py --sumstats $f --merge-alleles $HM --signed-sumstats Z,0 --out munged/$b > /dev/null; done
for s in s09 s10 s11; do $PY $L/ldsc.py --rg munged/rep_${s}_hi_c2.sumstats.gz,munged/rep_${s}_lo_c1.sumstats.gz --ref-ld-chr $REF --w-ld-chr $REF --out rg/pair_$s > /dev/null; done
grep -h -A2 "Summary of Genetic" rg/pair_*.log
