#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --array=1-3
#SBATCH --job-name=greml
#SBATCH -o /home/jsebat/rvh2/logs/greml_%A_%a.out
# SNP heritability (GCTA-GREML, as the paper did for age at diagnosis) of each case characteristic, among unrelated autistic
# participants per platform; dense GRM on HapMap3 SNPs (hm3w from rep_prep), --grm-cutoff 0.05; binary traits on the liability
# scale with K = sample prevalence (--prevalence).
set -e
SRCS=(spark_gsa_v1 spark_gsa_v2 spark_twist); export SRC=${SRCS[$((SLURM_ARRAY_TASK_ID-1))]}
W=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/$SRC; cd $W; G=/expanse/projects/sebat1/jsebat/tools/gcta/gcta64
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/greml_pheno.py
[ -s grm_cases.grm.bin ] || $G --bfile hm3w --keep greml_keep.txt --maf 0.01 --make-grm --thread-num 16 --out grm_cases > /dev/null
k=0; while read nm; do k=$((k+1))
  PV=$(awk -v n=$nm '$1==n{print "--prevalence "$2}' greml_prev.txt)
  $G --grm grm_cases --grm-cutoff 0.05 --reml --pheno greml_pheno.txt --mpheno $k --qcovar greml_qcovar.txt --covar greml_covar.txt $PV --reml-no-lrt --thread-num 16 --out greml_$nm > /dev/null 2>&1 || echo "failed $nm"
  [ -s greml_$nm.hsq ] && { echo "== $nm"; grep -E "V\(G\)/Vp|n\b" greml_$nm.hsq; }
done < greml_names.txt
echo DONE
