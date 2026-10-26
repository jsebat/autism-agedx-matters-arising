#!/bin/bash
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --array=0-8
#SBATCH --job-name=all_ldsc
#SBATCH -o /home/jsebat/rvh2/logs/all_ldsc_%A_%a.out
# ALL-DATA split GWAS, LDSC. Task 0: combine six strata (combine_all.py), munge, write the rg task list. Tasks 1-8 (wait for task 0): every
# rg in the list (age-split pairs, random-split null, survey halves vs each other and vs the published factors, external-trait profiles).
set -e
export OPENBLAS_NUM_THREADS=1
T=/expanse/projects/sebat1/jsebat/tools; PY=$T/envs/ldsc3/bin/python; L=$T/ldsc_py3; REF=$T/ldref/eur_w_ld_chr/; HM=$T/ldref/w_hm3.snplist
C=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined_all; OLD=/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined
mkdir -p $C/munged $C/rg; cd $C
if [ "$SLURM_ARRAY_TASK_ID" = "0" ]; then
  rm -f tasks.ready
  /home/jsebat/micromamba/envs/dnmt3a-py/bin/python ~/rvh2/combine_all.py
  for f in rep_*.sumstats.txt; do b=${f%.sumstats.txt}; $PY $L/munge_sumstats.py --sumstats $f --merge-alleles $HM --signed-sumstats Z,0 --out munged/$b > /dev/null || echo "munge failed $b"; done
  cp -n $OLD/munged/W_F1.sumstats.gz $OLD/munged/W_F2.sumstats.gz $OLD/munged/W_agedx.sumstats.gz $OLD/munged/ext_*.sumstats.gz munged/
  /home/jsebat/micromamba/envs/dnmt3a-py/bin/python - <<'EOF'
import os, glob
M = "munged"; ex = lambda n: os.path.exists(f"{M}/{n}.sumstats.gz"); p = lambda n: f"{M}/{n}.sumstats.gz"
REFS = [x for x in ["W_F1", "W_F2", "W_agedx", "rep_W_over10_c1", "rep_W_under6_c1"] if ex(x)]
EXT = [os.path.basename(x)[:-12] for x in sorted(glob.glob(f"{M}/ext_*.sumstats.gz"))]
T = []
for a, b in [("rep_W_under6", "rep_W_over10"), ("rep_W_under6_c1", "rep_W_over10_c2")] + [(f"rep_rand{r}_A_c2", f"rep_rand{r}_B_c1") for r in range(1, 5)]:
    if ex(a) and ex(b): T.append((f"pair__{a}", [a, b] + REFS))
for a in sorted(glob.glob(f"{M}/rep_*N*_A_c2.sumstats.gz")):
    a = os.path.basename(a)[:-12]; b = a.replace("_A_c2", "_B_c1")
    if ex(b): T.append((f"null__{a}", [a, b]))
for a in sorted(glob.glob(f"{M}/rep_s*_hi_c2.sumstats.gz")):
    a = os.path.basename(a)[:-12]; b = a.replace("_hi_c2", "_lo_c1")
    if ex(b): T.append((f"survey__{a}", [a, b] + REFS)); T.append((f"survey__{b}", [b] + REFS))
for a in sorted(glob.glob(f"{M}/rep_s*_hi_c2.sumstats.gz")) + sorted(glob.glob(f"{M}/rep_s*_lo_c1.sumstats.gz")) + [p("rep_W_over10_c1"), p("rep_W_under6_c1"), p("W_F1"), p("W_F2")]:
    n = os.path.basename(a)[:-12]
    if ex(n): T.append((f"ext__{n}", [n] + EXT))
open("rg_tasks.txt", "w").write("\n".join(o + "|" + ",".join(p(x) for x in xs) for o, xs in T) + "\n"); print(len(T), "rg tasks")
EOF
  touch tasks.ready; echo DONE task0
else
  for i in $(seq 1 720); do [ -f tasks.ready ] && break; sleep 30; done
  k=0; while IFS='|' read -r out files; do k=$((k+1)); [ $(( (k-1) % 8 + 1 )) -eq $SLURM_ARRAY_TASK_ID ] || continue
    [ -s rg/$out.log ] && grep -q "Summary of Genetic" rg/$out.log && continue
    $PY $L/ldsc.py --rg $files --ref-ld-chr $REF --w-ld-chr $REF --out rg/$out > /dev/null || echo "rg failed $out"
  done < rg_tasks.txt
  echo DONE task $SLURM_ARRAY_TASK_ID
fi
