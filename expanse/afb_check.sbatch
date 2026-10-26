#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /home/jsebat/rvh2/logs/afb_check_%j.out
# PGS check + LDSC rg of the maternal-age / income halves with the AFB 2018 and EA 2022 GWAS.
export OPENBLAS_NUM_THREADS=1
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/afb_check.py
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/; HM=$T/ldref/w_hm3.snplist
C=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined; cd $C; mkdir -p rg_afb
[ -s munged/ext_afb_2018.sumstats.gz ] || $PY $L/munge_sumstats.py --sumstats /expanse/projects/sebat1/s3/data/sebat/pgs_pipeline/input/afb_2018.ma --merge-alleles $HM --signed-sumstats b,0 --N-col N --out munged/ext_afb_2018 > /dev/null
for x in rep_s11_hi_c2 rep_s11_lo_c1 rep_s09_hi_c2 rep_s09_lo_c1 rep_s00_hi_c2 rep_s00_lo_c1; do
  $PY $L/ldsc.py --rg munged/$x.sumstats.gz,munged/ext_afb_2018.sumstats.gz,munged/ext_edu_2022.sumstats.gz,munged/ext_iq_2018.sumstats.gz,munged/ext_adhd_2022.sumstats.gz --ref-ld-chr $REF --w-ld-chr $REF --out rg_afb/$x > /dev/null
  echo "== $x"; sed -n '/Summary of Genetic Correlation Results/,$p' rg_afb/$x.log | awk 'NR>1 && NF>4 {n=split($2,a,"/"); sub(".sumstats.gz","",a[n]); printf "%-18s rg %6.2f (%.2f)\n", a[n], $3, $4}'
done
