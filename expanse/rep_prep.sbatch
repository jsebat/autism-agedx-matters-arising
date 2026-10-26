#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 06:00:00
#SBATCH --array=1-3
#SBATCH --job-name=rep_prep
#SBATCH -o /home/jsebat/rvh2/logs/rep_prep_%A_%a.out
# Replication of Zhang, Grove et al. 2025 SPARK stratified GWAS (vwarrier/autism_agediagnosis SPARK_replicationGWAS.R), step 1 per
# platform, on the matched-design EUR HapMap3 files (cnv_drift/mgwas/<src>/hm3*):
#   PC-AiR-style PCs: PCA (10, allele weights) in the KING-unrelated set, projected onto every sample (their PCs: GENESIS PC-AiR)
#   sparse GRM as in their script: GCTA --make-grm (MAF 0.01; here on the LD-pruned SNP set, enough to find relatives for the 0.05 cut) -> --make-bK-sparse 0.05
# Output: cnv_drift/warrier_rep/<src>/{pcair.sscore, grm_sp.*}
set -e
SRCS=(spark_gsa_v1 spark_gsa_v2 spark_twist); SRC=${SRCS[$((SLURM_ARRAY_TASK_ID-1))]}
M=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/$SRC; W=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/$SRC; mkdir -p $W; cd $W
P2=/home/j3guevar/bin/plink2; G=/expanse/projects/sebat1/jsebat/tools/gcta/gcta64; HM=/expanse/projects/sebat1/jsebat/tools/ldref/w_hm3.snplist
$P2 --pfile $M/hm3 --keep $M/king.king.cutoff.in.id --extract $M/prune.prune.in --freq --pca 10 approx allele-wts --threads 16 --memory 88000 --out pca_unrel
$P2 --pfile $M/hm3 --extract $M/prune.prune.in --read-freq pca_unrel.afreq --score pca_unrel.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --threads 16 --memory 88000 --out pcair
$P2 --bfile $M/hm3_bed --extract $HM --maf 0.01 --make-bed --threads 16 --memory 88000 --out hm3w
for i in $(seq 1 8); do $G --bfile hm3w --extract $M/prune.prune.in --make-grm-part 8 $i --maf 0.01 --thread-num 16 --out grm > /dev/null; done
cat grm.part_8_*.grm.id > grm_full.grm.id; cat grm.part_8_*.grm.bin > grm_full.grm.bin; cat grm.part_8_*.grm.N.bin > grm_full.grm.N.bin
$G --grm grm_full --make-bK-sparse 0.05 --thread-num 16 --out grm_sp
rm -f grm.part_8_* grm_full.grm.bin grm_full.grm.N.bin
wc -l pcair.sscore grm_sp.grm.id
echo DONE
