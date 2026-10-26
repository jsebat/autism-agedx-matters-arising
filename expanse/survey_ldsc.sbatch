#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH --job-name=survey_ldsc
#SBATCH -o /home/jsebat/rvh2/logs/survey_ldsc_%j.out
# Broad survey: combine platforms, munge, and for every split: rg(high half, low half) and the h2 of each half; each half's rg with
# the published factors (F1, F2), the SPARK age-at-diagnosis GWAS and the replicated late / early GWAS. Output rg_survey/*.log.
set -e
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/; HM=$T/ldref/w_hm3.snplist
C=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined; cd $C; mkdir -p rg_survey
COUNTS=counts_survey.tsv /home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/rep_combine.py > /dev/null
REFS="munged/W_F1.sumstats.gz,munged/W_F2.sumstats.gz,munged/W_agedx.sumstats.gz,munged/rep_W_over10_c1.sumstats.gz,munged/rep_W_under6_c1.sumstats.gz"
for f in rep_s*_hi_c2.sumstats.txt; do a=${f%.sumstats.txt}; b=${a/_hi_c2/_lo_c1}
  for x in $a $b; do [ -s $x.sumstats.txt ] && { [ -s munged/$x.sumstats.gz ] || $PY $L/munge_sumstats.py --sumstats $x.sumstats.txt --merge-alleles $HM --signed-sumstats Z,0 --out munged/$x > /dev/null; }; done
  [ -s munged/$b.sumstats.gz ] && $PY $L/ldsc.py --rg munged/$a.sumstats.gz,munged/$b.sumstats.gz,$REFS --ref-ld-chr $REF --w-ld-chr $REF --out rg_survey/$a > /dev/null || echo "rg failed $a"
  [ -s munged/$b.sumstats.gz ] && $PY $L/ldsc.py --rg munged/$b.sumstats.gz,$REFS --ref-ld-chr $REF --w-ld-chr $REF --out rg_survey/$b > /dev/null || echo "rg failed $b"
done
echo DONE
