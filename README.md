# Axes of genetic heterogeneity in autism: code for a Matters Arising

Analysis code for a Matters Arising on Zhang, Grove et al., "Polygenic and developmental profiles of autism differ by age at diagnosis", *Nature* 646, 1146–1155 (2025).

The code does four things:
1. Replicates the paper's SPARK case-stratified GWAS.
2. Applies the same design to many other splits of autistic participants.
3. Ranks the axes that link characteristics of autistic people to their polygenic scores and rare variants.
4. Runs sensitivity analyses for genotyping platform, family structure and registration age.

No participant-level data, study identifiers or results are included.

## Layout

| Folder | Contents |
|---|---|
| `expanse/` | Cluster code (SLURM, SDSC Expanse): genotype preparation, split GWAS (GCTA fastGWA-mlm-binary), meta-analysis, LDSC, GREML, polygenic-score and rare-variant analyses, canonical correlation analyses |
| `figures_and_analysis/` | Local code: genomic factor analysis, survey assembly, figures |
| `scripts/id_scan.sh` | Pre-push check that the repository holds no identifiers or data |
| `README_DATA_WARNING.md` | Which polygenic-score and genotype tables to use (and which not to) |

## Pipeline (in order)

1. **Genotypes per stratum.** `mgwas_prep.sbatch` (SPARK GSA v1, GSA v2, Twist) and `allgeno_prep.sbatch` (SPARK iWGS v1.1, SPARK whole-genome release 2026-08, SSC). Each keeps European ancestry and HapMap3 SNPs, and builds a KING-unrelated set and principal components. One person per stratum.
2. **PCs, sparse GRM, phenotypes.** `rep_prep.sbatch` / `allstrata_prep.sbatch`, then `rep_pheno.py`, `survey_pheno.py`, `null_pheno.py`, `allstrata_pheno.py`, `allstrata_null.py`. Cases are compared with all non-autistic relatives, using disjoint control halves split by family.
3. **Split GWAS.** `rep_gwas.sbatch`, `survey_run.sbatch`, `null_run.sbatch`, `allstrata_gwas.sbatch`, `allstrata_nullgwas.sbatch`.
4. **Meta-analysis and LDSC.**
   - Arrays only (replication of the paper's data): `rep_combine.py`, `rep_ldsc.sbatch`, `rep_rg.sbatch`, `survey_ldsc.sbatch`, `null_ldsc.sbatch`, `ext_rg.sbatch`, parsers.
   - All strata: `combine_all.py`, `all_ldsc.sbatch`, `all_parse.py`.
   - Genomic factor analysis on the all-strata split GWAS: `gsem_fa_all.sbatch` (`gsem_fa_all.R`).
5. **SNP heritability.** `greml.sbatch`, `greml_all.sbatch`, `greml_combine_all.py`.
6. **Polygenic scores.** `full_probands.py` assembles scores for all subcohorts and records family structure. Autism scores follow the leave-one-out rule: `score_wgs2026_asd.sbatch`.
7. **Canonical correlation analyses.** Current version: `cca2.py` with the shared module `ccalib.py` (configurations in `cca2_configs.cfg`; `submit_cca2.sh` runs main fits, permutation and bootstrap chunks, then `cca2_combine.sbatch`). Evaluation is fully held out: nested family cross-validation with imputation, covariate adjustment and penalty choice inside the training folds; signed held-out correlations with axes matched across folds; Freedman-Lane family permutations with sequential tests for later axes; family bootstrap of the whole procedure, including paired fold differences between axes and single characteristics. The scripts below are the first version, kept for reference:
   - `cca_full.py`: polygenic scores; main sample, sensitivity subsets, within-family, power, missingness by family structure.
   - `cca_trio.py`: polygenic + rare variants.
   - `cca_ssc.py`: SSC replication.
   - `SINGLE=1` (`ccasingle.sbatch`, `ccasingle_all.sbatch`, `ccatriosingle.sbatch`): strength of single characteristics, including age at diagnosis, on the same scale as the axes.
8. **Rare variants.** `build_trio_rare.py` (de novo / inherited × NDD / other genes, SNV/indel), then `rare_full.py` (split comparison) and `rare_between_full.py` (probands vs siblings, parental age).
9. **Within-family transmission.** `ptdt_full.py`.
10. **Factor analysis and figures.** `figures_and_analysis/`.
11. **Simulation.** `sim_twofactor.py`: one autism liability plus a heritable stratifier of age at diagnosis, passed through a six-GWAS design like the original study's, with one- and two-factor models fitted by weighted least squares.

Paths in the scripts point to the authors' cluster storage. Adapt them to your own copies of the inputs.

## Data needed (not distributed here)

- **SPARK and the Simons Simplex Collection** genotypes, sequence data and phenotypes: from SFARI Base, by application.
- **Published factor GWAS** of Zhang, Grove et al.: figshare, doi:10.6084/m9.figshare.29566052.
- **Public GWAS summary statistics** for the polygenic scores: listed with references in the manuscript's Supplementary Table 5.
- **Autism polygenic scores** use pre-publication leave-one-out summary statistics from the PGC Autism Working Group. These are available from the working group, not from this repository.

## Software

- GCTA 1.94.1 and PLINK 2.
- LDSC (Python 3 port) and GenomicSEM.
- SBayesRC.
- Python 3.11: pandas, numpy, scipy, statsmodels, matplotlib.

## Licence

MIT
