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
#SBATCH --job-name=allgeno_prep
#SBATCH -o /home/jsebat/rvh2/logs/allgeno_prep_%A_%a.out
# Split-GWAS strata beyond the three array platforms (use all data): SPARK iWGS v1.1, SPARK-WGS-2026_08, SSC.
# Same steps as mgwas_prep.sbatch: EUR samples, HapMap3 SNPs (.ma SNP list) MAF >= 1%, missingness <= 5%; KING unrelated set; 10 PCs in it; bed copy.
# SPARK-WGS-2026_08 is GRCh38 with chr:pos:ref:alt IDs: renamed to rsIDs with pgs_wf_spark_wgs/rsid_maps (GRCh38) before extraction.
set -e
SRCS=(spark_wgs spark_wgs_2026_08 ssc); SRC=${SRCS[$((SLURM_ARRAY_TASK_ID-1))]}
B=/expanse/projects/sebat1/s3/data/sebat; J=/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08; P2=/home/j3guevar/bin/plink2
W=/expanse/projects/sebat1/jsebat/cnv_drift/mgwas/$SRC; mkdir -p $W; cd $W
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
$PY - <<EOF
import pandas as pd, glob
src = "$SRC"; B = "$B"; J = "$J"
hm3 = set(pd.read_csv(f"{B}/pgs_pipeline/input/asd_ipsych_dbs.ma", sep="\t", usecols=["SNP"]).SNP)
if src == "spark_wgs_2026_08":
    a = pd.read_csv(f"{J}/plink/pgs/analyses/scoring-ancestry-20261012/06_pca/ancestry/ancestry_probabilities.tsv", sep="\t", dtype=str)
    a[a.ANCESTRY.eq("EUR")][["#IID"]].to_csv("eur_keep.txt", sep="\t", index=False, header=False)
    rows = []
    for f in sorted(glob.glob(f"{B}/pgs_wf_spark_wgs/rsid_maps/chr*.map")):
        m = pd.read_csv(f, sep=r"\s+", header=None, names=["key", "rs"], dtype=str); m = m[m.rs.isin(hm3)]
        k = m.key.str.split(":", expand=True); swap = k[0] + ":" + k[1] + ":" + k[3] + ":" + k[2]   # allow REF/ALT swapped keys
        rows += [m[["key", "rs"]], pd.DataFrame({"key": swap, "rs": m.rs})]
    r = pd.concat(rows).drop_duplicates("key")
    r.to_csv("rename.txt", sep="\t", index=False, header=False); r.key.to_csv("hm3keys.txt", index=False, header=False)
    print("EUR", (a.ANCESTRY == "EUR").sum(), "| HM3 keys", len(r))
else:
    d = pd.read_csv(f"{B}/pgs_results_v2/all_cohorts_with_ped_full.tsv", sep="\t", dtype=str, low_memory=False, usecols=["source", "IID", "predicted_ancestry"])
    e = set(d[(d.source == src) & (d.predicted_ancestry == "EUR")].IID)
    ps = pd.read_csv(f"{B}/pgs_wf_{src}/05_summary/{src}.psam", sep="\t", dtype=str); ic = [c for c in ps.columns if "IID" in c][0]
    ps[ps[ic].isin(e)][[ic]].to_csv("eur_keep.txt", sep="\t", index=False, header=False)
    print("EUR", len(e), "| in psam", int(ps[ic].isin(e).sum()))
pd.Series(sorted(hm3)).to_csv("hm3.txt", index=False, header=False)
EOF
printf "6 25000000 34000000 MHC\n8 8000000 12000000 inv8p23\n" > longld.txt
if [ "$SRC" = "spark_wgs_2026_08" ]; then
  $P2 --pfile $J/plink/pgs/results/04_score_input/score_input vzs --keep eur_keep.txt --extract hm3keys.txt --autosome --snps-only just-acgt --make-pgen --threads 16 --memory 88000 --out keyed
  $P2 --pfile keyed --update-name rename.txt --make-pgen --threads 16 --memory 88000 --out renamed
  IN=renamed
else
  IN=$B/pgs_wf_$SRC/05_summary/$SRC
fi
$P2 --pfile $IN --keep eur_keep.txt --extract hm3.txt --autosome --snps-only just-acgt --maf 0.01 --geno 0.05 --rm-dup exclude-all --make-pgen --threads 16 --memory 88000 --out hm3
$P2 --pfile hm3 --maf 0.05 --exclude bed1 longld.txt --indep-pairwise 1000kb 1 0.1 --threads 16 --memory 88000 --out prune
$P2 --pfile hm3 --extract prune.prune.in --king-cutoff 0.0442 --threads 16 --memory 88000 --out king
$P2 --pfile hm3 --keep king.king.cutoff.in.id --extract prune.prune.in --pca 10 approx --threads 16 --memory 88000 --out pca_unrel
$P2 --pfile hm3 --make-bed --threads 16 --memory 88000 --out hm3_bed
rm -f keyed.pgen renamed.pgen
echo DONE $(wc -l < hm3_bed.fam) samples $(wc -l < hm3_bed.bim) SNPs
