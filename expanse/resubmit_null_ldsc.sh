#!/bin/bash
S=/cm/shared/apps/slurm/current/bin
$S/scancel 55279927
N=$($S/sbatch --parsable ~/rvh2/allstrata_nullgwas.sbatch); echo null $N
L=$($S/sbatch --parsable --dependency=afterok:$N ~/rvh2/all_ldsc.sbatch); echo ldsc $L
G=$($S/sbatch --parsable ~/rvh2/greml_combine.sbatch); echo greml_combine $G
