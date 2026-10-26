#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --job-name=rep_ldsc
#SBATCH -o /home/jsebat/rvh2/logs/rep_ldsc_%j.out
# Replication step 5: combine platforms, munge, LDSC rg of every replication GWAS vs every other, vs the published factors of
# Zhang, Grove et al. 2025 (figshare 29566052 v2: Earlier_diagnosed_autism = F1, Later_diagnosed_autism = F2, SPARK age at diagnosis)
# and vs external traits; matched-design within-family signal (sibtrio_agedx_early) included for reference.
# This job: combine + munge and the target list; rg runs in rep_rg.sbatch (array), parsing in rep_parse.py.
set -e
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/; HM=$T/ldref/w_hm3.snplist
C=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined; WS=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_sumstats
MG=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/combined/munged
mkdir -p $C/munged $C/rg; cd $C
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/rep_combine.py
for f in rep_*.sumstats.txt; do b=${f%.sumstats.txt}; $PY $L/munge_sumstats.py --sumstats $f --merge-alleles $HM --signed-sumstats Z,0 --out munged/$b > /dev/null; done
[ -s munged/W_F1.sumstats.gz ] || $PY $L/munge_sumstats.py --sumstats $WS/Earlier_diagnosed_autism.txt --merge-alleles $HM --signed-sumstats Beta,0 --N-col N_effective --a1 A1 --a2 A2 --out munged/W_F1 > /dev/null
[ -s munged/W_F2.sumstats.gz ] || $PY $L/munge_sumstats.py --sumstats $WS/Later_diagnosed_autism.txt --merge-alleles $HM --signed-sumstats Beta,0 --N-col N_effective --a1 A1 --a2 A2 --out munged/W_F2 > /dev/null
[ -s munged/W_agedx.sumstats.gz ] || $PY $L/munge_sumstats.py --sumstats $WS/SPARK_agediagnosis_meta_cleanedsumstats.txt --merge-alleles $HM --signed-sumstats BETA,0 --N-col N --a1 A1 --a2 A2 --out munged/W_agedx > /dev/null
cp -n $MG/ext_*.sumstats.gz munged/ ; cp -n $MG/sibtrio_agedx_early.sumstats.gz munged/
ls munged/rep_*.sumstats.gz munged/W_F1.sumstats.gz munged/W_F2.sumstats.gz > rg_targets.txt; wc -l rg_targets.txt
echo DONE
