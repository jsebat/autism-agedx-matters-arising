#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 05:00:00
#SBATCH --job-name=pcproj
#SBATCH -o /home/jsebat/rvh2/logs/pcproj_%j.out
# PC-AiR-style PCs on SPARK Twist (EUR): PCA in a KING-unrelated set, then project every sample with the unrelated-set
# allele weights and frequencies (plink2). Output: projected PCs for all samples + in-set membership (Expanse only).
set -e
P=/home/j3guevar/bin/plink2; PF=/expanse/projects/sebat1/s3/data/sebat/pgs_wf_spark_twist/03_missingness/cohort
W=/expanse/projects/sebat1/jsebat/cnv_drift/pcproj; mkdir -p $W; cd $W
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/pcproj_ids.py
printf "6 25000000 34000000 MHC\n8 8000000 12000000 inv8p23\n" > longld.txt
$P --pfile $PF --keep eur_keep.txt --autosome --maf 0.05 --geno 0.02 --snps-only just-acgt --exclude bed1 longld.txt --indep-pairwise 1000kb 1 0.1 --threads 16 --memory 88000 --out prune
$P --pfile $PF --keep eur_keep.txt --extract prune.prune.in --make-pgen --threads 16 --memory 88000 --out base
$P --pfile base --king-cutoff 0.0442 --threads 16 --memory 88000 --out king
$P --pfile base --keep king.king.cutoff.in.id --freq --pca 10 approx allele-wts --threads 16 --memory 88000 --out pca_unrel
$P --pfile base --read-freq pca_unrel.afreq --score pca_unrel.eigenvec.allele 2 6 header-read no-mean-imputation variance-standardize --score-col-nums 7-16 --threads 16 --memory 88000 --out proj_all
echo DONE
