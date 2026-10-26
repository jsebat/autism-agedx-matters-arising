#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --job-name=mgwas_ldsc
#SBATCH -o /home/jsebat/rvh2/logs/mgwas_ldsc_%j.out
# Matched-design GWAS, step 5: combine platforms, munge to HapMap3, LDSC h2 and rg (CBIIT python3 ldsc; eur_w_ld_chr).
#  rg: each split-half GWAS vs every other split-half GWAS (both designs) and vs external sumstats (pgs_pipeline inputs).
set -e
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/; HM=$T/ldref/w_hm3.snplist
C=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/combined; mkdir -p $C/munged $C/rg; cd $C
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/mgwas_combine.py
for f in between_*.sumstats.txt sibtrio_*.sumstats.txt; do b=${f%.sumstats.txt}
  $PY $L/munge_sumstats.py --sumstats $f --merge-alleles $HM --signed-sumstats Z,0 --out munged/$b > /dev/null
done
X=/expanse/projects/sebat1/s3/data/sebat/pgs_pipeline/input
for t in asd_ipsych_dbs adhd_2022 mdd_2025 bip_2024 scz_2022 edu_2022 iq_2018; do
  [ -s munged/ext_$t.sumstats.gz ] || $PY $L/munge_sumstats.py --sumstats $X/$t.ma --merge-alleles $HM --signed-sumstats b,0 --N-col N --out munged/ext_$t > /dev/null
done
ALL=$(ls munged/between_*.sumstats.gz munged/sibtrio_*.sumstats.gz | tr "\n" "," )
EXT=$(ls munged/ext_*.sumstats.gz | tr "\n" ",")
for f in munged/between_*.sumstats.gz munged/sibtrio_*.sumstats.gz; do b=$(basename $f .sumstats.gz)
  $PY $L/ldsc.py --rg ${f},${ALL}${EXT%,} --ref-ld-chr $REF --w-ld-chr $REF --out rg/$b > /dev/null || echo "rg failed $b"
done
echo DONE
