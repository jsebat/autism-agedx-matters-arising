#!/bin/bash
# After the Grove-2026-full rescoring of SPARK-WGS-2026_08 (job $1): rebuild full_probands, then every analysis that uses the autism PGS.
S=/cm/shared/apps/slurm/current/bin; D=$1
A=$($S/sbatch --parsable --dependency=afterok:$D ~/rvh2/fullpro.sbatch); echo fullpro $A
for f in cca_fulls.sbatch cca_fullp.sbatch cca_fullu.sbatch cca_within.sbatch ccapow.sbatch ccamiss.sbatch ptdtfull.sbatch supptab.sbatch ccassc.sbatch rarefull.sbatch; do
  K=$($S/sbatch --parsable --dependency=afterok:$A ~/rvh2/$f); echo $f $K
  [ "$f" = "rarefull.sbatch" ] && R=$K
done
T=$($S/sbatch --parsable --dependency=afterok:$R ~/rvh2/ccatrio.sbatch); echo ccatrio $T
