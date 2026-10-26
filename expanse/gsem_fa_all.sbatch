#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --job-name=gsem_all
#SBATCH -o /home/jsebat/rvh2/logs/gsem_all_%j.out
# GenomicSEM multivariable LDSC + EFA on the ALL-STRATA split-half GWAS, same 23 indicators as the arrays-only run (for Fig 2a).
set -e
export OPENBLAS_NUM_THREADS=1
B=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep; C=$B/combined_all
mkdir -p $C/gsem; cp $B/combined/gsem/fa_traits.txt $C/gsem/
for t in $(cat $C/gsem/fa_traits.txt); do [ -s $C/munged/$t.sumstats.gz ] || { echo "missing $t"; exit 1; }; done
GSEM_C=$C /expanse/projects/sebat1/jsebat/tools/envs/gsem/bin/Rscript ~/rvh2/gsem_fa_all.R
