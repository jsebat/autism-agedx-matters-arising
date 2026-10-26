#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --job-name=null_ldsc
#SBATCH -o /home/jsebat/rvh2/logs/null_ldsc_%j.out
# Random-split null: combine platforms, munge, pairwise rg A vs B for each replicate; summary table null_rg.tsv.
set -e
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/; HM=$T/ldref/w_hm3.snplist
C=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined; cd $C; mkdir -p rg_null
COUNTS=counts_null.tsv /home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/rep_combine.py > /dev/null
for f in rep_ageN*_A_c2.sumstats.txt rep_adiN*_A_c2.sumstats.txt; do a=${f%.sumstats.txt}; b=${a/_A_c2/_B_c1}
  for x in $a $b; do [ -s munged/$x.sumstats.gz ] || $PY $L/munge_sumstats.py --sumstats $x.sumstats.txt --merge-alleles $HM --signed-sumstats Z,0 --out munged/$x > /dev/null; done
  $PY $L/ldsc.py --rg munged/$a.sumstats.gz,munged/$b.sumstats.gz --ref-ld-chr $REF --w-ld-chr $REF --out rg_null/$a > /dev/null || echo "rg failed $a"
done
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python - <<'EOF'
import glob, io, os, pandas as pd
rows = []
for f in sorted(glob.glob("rg_null/*.log")):
    t = open(f).read()
    if "Summary of Genetic Correlation Results" not in t: continue
    blk = t.split("Summary of Genetic Correlation Results")[1].strip().split("\n\n")[0]; x = pd.read_csv(io.StringIO(blk), sep=r"\s+").iloc[0]
    rows.append(dict(rep=os.path.basename(f).replace(".log", ""), rg=x.rg, se=x.se, h2_A=x.h2_obs, h2_A_se=x.h2_obs_se))
N = pd.DataFrame(rows); N["set"] = N.rep.str.extract(r"rep_(ageN|adiN)")[0]; N.to_csv("null_rg.tsv", sep="\t", index=False)
print(N.round(3).to_string(index=False))
print(N.groupby("set").rg.describe().round(3).to_string())
EOF
echo DONE
