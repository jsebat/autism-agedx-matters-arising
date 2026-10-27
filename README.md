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

## Figures and tables → code

Numbering follows the manuscript (items numbered by first mention in the main text). Cluster scripts are in `expanse/`, figure scripts in `figures_and_analysis/`.

### Main figures

| Item | Content | Figure script | Analysis scripts |
|---|---|---|---|
| Fig. 1 | Age split and 43 case splits: genetic correlation between halves | `letter_figs.py` | `allgeno_prep.sbatch`, `allstrata_prep.sbatch`, `allstrata_pheno.py`, `allstrata_null.py`, `allstrata_gwas.sbatch`, `allstrata_nullgwas.sbatch`, `combine_all.py`, `all_ldsc.sbatch`, `all_parse.py` |
| Fig. 2 | a, split arrows in external-trait space; b–c, rare-variant stratification; d–e, axis strength and loadings | `fig2_main.py`, `axis_figure.py` (panel a vectors), `efa_local.py` (EDU/SES direction) | `all_ldsc.sbatch` + `all_parse.py` (a); `build_trio_rare.py`, `rare_full.py` (b, c); `cca2.py` + `ccalib.py`, trio configurations in `cca2_configs.cfg` (d, e) |

### Extended Data

| Item | Content | Figure script | Analysis scripts |
|---|---|---|---|
| Extended Data Fig. 1 | Strength and structure of the canonical axes (between families, within families, trio probands) | `cca_figure.py` | `cca2.py`, `ccalib.py`, `submit_cca2.sh` |
| Extended Data Fig. 2 | Zhang factors in the space of our common-variant genomic factors | `factor_map.py` (uses `efa_local.py`) | `gsem_fa_all.sbatch`, `gsem_fa_all.R`; external-trait correlations from `all_ldsc.sbatch` |
| Extended Data Fig. 3 | Within-family analyses | `letter_figs2.py` | `ed3_full.py`, `ptdt_full.py`, `rare_full.py` |

### Supplementary Information

| Item | Content | Figure / table script | Analysis scripts |
|---|---|---|---|
| Supplementary Section 1 | Platform, family structure, missingness, adult registrants | (text) | `cca2.py` sensitivity configurations in `cca2_configs.cfg`; `rare_full.py` (`U22=0`) |
| Supplementary Table 1 | Data used in every analysis, by subcohort and role | manuscript builder (not in this repository) | `full_probands.py`, `supp_sample_table.py`, `combine_all.py` (split-GWAS counts), `greml_combine_all.py` |
| Supplementary Table 2 | Broad survey: within-pair genetic correlation for every split | `letter_figs.py` (survey table) | as Fig. 1 |
| Supplementary Table 3 | Summary statistics used for polygenic scores | manuscript builder | `score_asd2026_loo.sh`, `score_wgs2026_asd.sbatch` (autism score); other weights from the lab collection (see `README_DATA_WARNING.md`) |
| Supplementary Table 4 | Polygenic axes by sample, platform, family structure, genetic PCs; SSC replication | manuscript builder | `cca2.py` (configurations in `cca2_configs.cfg`, including `PCS=1` and SSC) |
| Supplementary Table 5 | Strength of each axis relative to age at diagnosis | manuscript builder | `cca2.py` (single-characteristic models, paired bootstrap; `*_single.tsv`, `*_ratio.tsv`) |
| Supplementary Table 6 | Within-family associations of the genetic predictors with autism | manuscript builder | `within_family_table.py` |
| Supplementary Table 7 | Power | manuscript builder; `survey_power.py` (split-GWAS power) | `cca2.py` (`POWER=1`) |
| Supplementary Table 8 | Recording by family structure and its effect on the EDU/SES associations | manuscript builder | `cca2.py` (`MISSBIAS=1`) |
| Supplementary Fig. 1 | Two factors from one autism liability (simulation) | `sim_figure.py` | `sim_twofactor.py`, `sim_twofactor2.py` |
| Supplementary Fig. 2 | Fig. 1 restricted to the arrays (data available to Zhang et al.) | `letter_figs.py` (`DATA=arrays`) | `mgwas_prep.sbatch`, `rep_prep.sbatch`, `rep_pheno.py`, `survey_pheno.py`, `null_pheno.py`, `rep_gwas.sbatch`, `survey_run.sbatch`, `null_run.sbatch`, `rep_combine.py`, `rep_ldsc.sbatch`, `survey_ldsc.sbatch`, `null_ldsc.sbatch` |
| Supplementary Fig. 3 | SNP heritability of the stratifiers among autistic participants | `letter_figs2.py` | `greml_all.sbatch` (runs `greml_pheno_all.py`), `greml_combine.sbatch` (`greml_combine_all.py`) |
| Supplementary Fig. 4 | What the split halves capture | `letter_figs2.py` | `all_ldsc.sbatch` (external-trait correlations), `gsem_fa_all.sbatch` (rotated factor solution) |
| Supplementary Fig. 5 | Design differences between the age strata | `letter_figs2.py` | `agedx_strata.py`, `agedx_decomp.py`, `pcproj_eval.py` |

The manuscript builders (text, table assembly, docx/PDF) are not included; every number they print comes from the aggregate outputs of the scripts above. Exploratory analyses that do not feed any manuscript item, including an earlier version of the canonical correlation analysis, are not included.

No copy-number variants are used anywhere in this analysis: every rare-variant result is from coding SNVs and indels (de novo and inherited).


## Pipeline (in order)

1. **Genotypes per stratum.** `mgwas_prep.sbatch` (SPARK GSA v1, GSA v2, Twist) and `allgeno_prep.sbatch` (SPARK iWGS v1.1, SPARK whole-genome release 2026-08, SSC). Each keeps European ancestry and HapMap3 SNPs, and builds a KING-unrelated set and principal components. One person per stratum.
2. **PCs, sparse GRM, phenotypes.** `rep_prep.sbatch` / `allstrata_prep.sbatch`, then `rep_pheno.py`, `survey_pheno.py`, `null_pheno.py`, `allstrata_pheno.py`, `allstrata_null.py`. Cases are compared with all non-autistic relatives, using disjoint control halves split by family.
3. **Split GWAS.** `rep_gwas.sbatch`, `survey_run.sbatch`, `null_run.sbatch`, `allstrata_gwas.sbatch`, `allstrata_nullgwas.sbatch`.
4. **Meta-analysis and LDSC.**
   - Arrays only (replication of the paper's data): `rep_combine.py`, `rep_ldsc.sbatch`, `rep_rg.sbatch`, `survey_ldsc.sbatch`, `null_ldsc.sbatch`, `ext_rg.sbatch`, parsers.
   - All strata: `combine_all.py`, `all_ldsc.sbatch`, `all_parse.py`.
   - Genomic factor analysis on the all-strata split GWAS: `gsem_fa_all.sbatch` (`gsem_fa_all.R`).
5. **SNP heritability.** `greml_all.sbatch` (every SPARK stratum; traits from `greml_pheno_all.py`), then `greml_combine.sbatch` (`greml_combine_all.py`).
6. **Polygenic scores.** `full_probands.py` assembles scores for all subcohorts and records family structure. Autism scores follow the leave-one-out rule: `score_wgs2026_asd.sbatch`.
7. **Canonical correlation analyses.** Current version: `cca2.py` with the shared module `ccalib.py` (configurations in `cca2_configs.cfg`; `submit_cca2.sh` runs main fits, permutation and bootstrap chunks, then `cca2_combine.sbatch`). Evaluation is fully held out: nested family cross-validation with imputation, covariate adjustment and penalty choice inside the training folds; signed held-out correlations with axes matched across folds; Freedman-Lane family permutations with sequential tests for later axes; family bootstrap of the whole procedure, including paired fold differences between axes and single characteristics. `ANALYSIS=full` fits polygenic scores (main sample, sensitivity subsets, within-family, power, missingness by family structure), `ANALYSIS=trio` polygenic scores plus rare variants, `ANALYSIS=ssc` the SSC replication; unless `NOSINGLE=1`, a run also measures the strength of single characteristics, including age at diagnosis, on the same scale as the axes.
8. **Rare variants.** `build_trio_rare.py` (de novo / inherited × NDD / other genes, SNV/indel), then `rare_full.py` (split comparison) and `rare_between_full.py` (probands vs siblings, parental age).
9. **Within-family transmission.** `ptdt_full.py`.
10. **Factor analysis and figures.** `figures_and_analysis/`.
11. **Simulation.** `sim_twofactor.py` and `sim_twofactor2.py` (direct splits, several stratifiers, ascertainment); figure `figures_and_analysis/sim_figure.py`. `sim_twofactor.py`: one autism liability plus a heritable stratifier of age at diagnosis, passed through a six-GWAS design like the original study's, with one- and two-factor models fitted by weighted least squares.

Paths in the scripts point to the authors' cluster storage. Adapt them to your own copies of the inputs.

## Data needed (not distributed here)

- **SPARK and the Simons Simplex Collection** genotypes, sequence data and phenotypes: from SFARI Base, by application.
- **Published factor GWAS** of Zhang, Grove et al.: figshare, doi:10.6084/m9.figshare.29566052.
- **Public GWAS summary statistics** for the polygenic scores: listed with references in the manuscript's Supplementary Table 5.
- **Autism polygenic scores** use pre-publication leave-one-out summary statistics from the PGC Autism Working Group. These are available from the working group, not from this repository.

## Software

- GCTA 1.94.1 and PLINK 2.
- LDSC (Python 3 port) and GenomicSEM 0.0.5 (R 4.5.2, lavaan 0.7.2; environment built by `expanse/inst_gsem.sh`).
- SBayesRC.
- Python 3.11: pandas, numpy, scipy, statsmodels, matplotlib.

## Licence

MIT
