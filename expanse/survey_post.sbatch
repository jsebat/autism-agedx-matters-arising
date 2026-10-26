#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH --job-name=survey_post
#SBATCH -o /home/jsebat/rvh2/logs/survey_post_%j.out
# Broad survey post-processing + our own genomic factor analysis (GenomicSEM multivariable LDSC, parallel analysis, EFA) and the
# composition analysis of the factor loadings.
set -e
export OPENBLAS_NUM_THREADS=1
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/survey_post.py
/expanse/projects/sebat1/jsebat/tools/envs/gsem/bin/Rscript ~/rvh2/gsem_fa.R
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/gsem_comp.py
echo DONE
