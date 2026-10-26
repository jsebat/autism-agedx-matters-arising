#!/bin/bash
# Review fixes (2026-10-26): allele audit, pTDT (all trios; same-platform trios), rare-variant comparison (pooled; within cohort; + parental age)
S=/cm/shared/apps/slurm/current/bin/sbatch; cd ~/rvh2
$S -J allele_audit --export=ALL,AUDIT=1 runpy.sbatch combine_all.py
$S -J ptdt --export=ALL runpy.sbatch ptdt_full.py
$S -J ptdt_same --export=ALL,SAMESRC=1 runpy.sbatch ptdt_full.py
$S -J rare_pooled --export=ALL,ADJ=pooled runpy.sbatch rare_full.py
$S -J rare_coh --export=ALL,ADJ=cohort runpy.sbatch rare_full.py
$S -J rare_par --export=ALL,ADJ=cohort_parage runpy.sbatch rare_full.py
