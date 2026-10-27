# GenomicSEM environment used for the genomic factor analyses (gsem_fa_all.R, gsem_fa.R). Compiled R dependencies come from conda-forge;
# simsalapar from CRAN; GenomicSEM from GitHub at the commit used (0.0.5, da95d43). Environment as used: R 4.5.2, lavaan 0.7.2, mgsub 2.0.0,
# simsalapar 1.0.13.
export OPENBLAS_NUM_THREADS=1
MM=/home/jsebat/bin/micromamba; P=/expanse/projects/sebat1/jsebat/tools/envs/gsem
$MM create -y -p $P -c conda-forge r-base=4.5 r-lavaan r-gdata r-corpcor r-data.table r-matrix r-remotes r-doparallel r-foreach r-stringr \
  r-splitstackshape r-r.utils r-mnormt r-gtools r-rcpp r-psych r-plyr r-mgcv r-e1071 r-readr r-dplyr r-colorspace r-sfsmisc r-mgsub
$P/bin/Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("simsalapar")
  remotes::install_github("GenomicSEM/GenomicSEM@da95d431a4709693d07eae55e172ade55526948b", upgrade = "never", dependencies = FALSE)
  library(GenomicSEM); cat("GenomicSEM", as.character(packageVersion("GenomicSEM")), "OK\n")'
