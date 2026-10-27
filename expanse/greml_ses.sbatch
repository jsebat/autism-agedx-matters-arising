#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --array=1-5
#SBATCH --job-name=greml_ses
#SBATCH -o /home/jsebat/rvh2/logs/greml_ses_%A_%a.out
# SNP heritability among unrelated autistic participants of the SES stratifiers named in the text (maternal and paternal age, low household
# income, mother without degree), every SPARK stratum; same design as greml_all.sbatch, existing case GRMs reused; only the new traits are fitted.
set -e
SRCS=(spark_gsa_v1 spark_gsa_v2 spark_twist spark_wgs spark_wgs_2026_08); export SRC=${SRCS[$((SLURM_ARRAY_TASK_ID-1))]}
W=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/$SRC; cd $W; G=/expanse/projects/sebat1/jsebat/tools/gcta/gcta64
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/greml_pheno_all.py
[ -s grm_cases.grm.bin ] || $G --bfile hm3w --keep greml_keep.txt --maf 0.01 --make-grm --thread-num 8 --out grm_cases > /dev/null
k=0; while read nm; do k=$((k+1))
  case $nm in mat_age|pat_age|low_income|mother_nodeg) ;; *) continue ;; esac
  PV=$(awk -v n=$nm '$1==n{print "--prevalence "$2}' greml_prev.txt)
  $G --grm grm_cases --grm-cutoff 0.05 --reml --pheno greml_pheno.txt --mpheno $k --qcovar greml_qcovar.txt --covar greml_covar.txt $PV --reml-no-lrt --thread-num 8 --out greml_$nm > /dev/null 2>&1 || echo "failed $nm"
  [ -s greml_$nm.hsq ] && { echo "== $nm"; grep -E "V\(G\)/Vp|n\b" greml_$nm.hsq; }
done < greml_names.txt
echo DONE
