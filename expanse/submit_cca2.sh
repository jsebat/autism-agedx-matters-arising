#!/bin/bash
# Held-out CCA rerun (review, 2026-10-26): main fits -> permutation / bootstrap chunks -> combine
S=/cm/shared/apps/slurm/current/bin/sbatch; cd ~/rvh2
M=$($S --parsable cca2_main.sbatch); echo "main $M"
DEPS=""
sub() {  # cfg task total nchunk
  j=$($S --parsable --dependency=afterok:${M}_$1 --array=0-$(($4-1))%8 --export=ALL,CFG=$1,TASK=$2,NPERM=$3,NBOOT=$3,NCHUNK=$4 -J cca2_$2_$1 cca2_rep.sbatch); echo "cfg $1 $2 $j"; DEPS="$DEPS:$j"; }
for c in 1 2 3; do sub $c null 200 20; sub $c boot 200 20; done
for c in 4 5 6 7 8 9; do sub $c null 100 10; done
C=$($S --parsable --dependency=afterany:${M}${DEPS} cca2_combine.sbatch); echo "combine $C"
