# WARNING: `pgs_results_v2` does not include the new SPARK WGS release

`/expanse/projects/sebat1/s3/data/sebat/pgs_results_v2/all_cohorts_with_ped_full.tsv` was assembled in March 2026. It covers SPARK GSA v1, GSA v2, Twist, iWGS v1.1 and SSC.

It **does not include SPARK-WGS-2026_08**. Only 212 of that release's 45,178 genomes (0.5%) are in it. An analysis built on this table alone silently loses about 20,800 autistic SPARK participants.

## Use instead, or in addition

| What | Where (`/expanse/projects/sebat1/...`) |
|---|---|
| PGS, all 45,178 new genomes | `j3guevar/SPARK-WGS-2026_08/plink/pgs/analyses/scoring-ancestry-20261012/05_scores/*.sscore` |
| Ancestry, new genomes | `j3guevar/SPARK-WGS-2026_08/plink/pgs/analyses/scoring-ancestry-20261012/06_pca/` |
| SDSM barcode → SP ID, family, parents, sequenced parents, ASD, age | `j3guevar/SPARK-WGS-2026_08/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv` |
| SPARK + SSC children with PGS and rare variants (21,291) | `jsebat/pgs_rare_h2_drift/assembled_children.parquet` |

The new release is keyed by **SDSM barcodes** (`SDSM-xxxx`), not SP IDs. A join on SP IDs alone drops every new genome without an error.

## Check before every analysis

Tabulate the sample size per subcohort (arrays, WES IDT, WES Twist, iWGS v1.1, WGS-2026_08, SSC) and per role (proband, sibling, parent). A subcohort with zero or near-zero rows means a broken ID join, not missing data.

Use all available data by default. Restrictions (trios only, siblings required, overlap with a published study) belong in sensitivity analyses.

Written 2026-10-25 after the age-at-diagnosis Matters Arising used this table alone (J. Sebat lab).
