#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH -o /home/jsebat/rvh2/logs/pcproj2_%j.out
set -e
cd /expanse/projects/sebat1/jsebat/cnv_drift/pcproj
/home/j3guevar/bin/plink2 --pfile base --read-freq pca_unrel.afreq --score pca_unrel.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --threads 8 --memory 28000 --out proj_all
