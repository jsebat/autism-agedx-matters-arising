#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --job-name=mgwas_ldsc2
#SBATCH -o /home/jsebat/rvh2/logs/mgwas_ldsc2_%j.out
# Round 2: combine platforms (all splits), munge the new between-family GWAS, rg of each new GWAS vs all between-family GWAS,
# the within-family early-diagnosis signal, and external traits; then parse everything.
set -e
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/; HM=$T/ldref/w_hm3.snplist
C=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/combined; cd $C
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/mgwas_combine.py
NEW="adiE_high adiE_low adhdE_yes adhdE_no sexE_female sexE_male agedxC_late agedxC_early"
for s in $NEW; do $PY $L/munge_sumstats.py --sumstats between_$s.sumstats.txt --merge-alleles $HM --signed-sumstats Z,0 --out munged/between_$s > /dev/null; done
ALL=$(ls munged/between_*.sumstats.gz | tr "\n" ","); EXT=$(ls munged/ext_*.sumstats.gz | tr "\n" ",")
for s in $NEW; do
  $PY $L/ldsc.py --rg munged/between_$s.sumstats.gz,${ALL}munged/sibtrio_agedx_early.sumstats.gz,${EXT%,} --ref-ld-chr $REF --w-ld-chr $REF --out rg/between_$s > /dev/null || echo "rg failed $s"
done
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/mgwas_parse.py > /dev/null
echo DONE
