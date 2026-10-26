export OPENBLAS_NUM_THREADS=1; P=/expanse/projects/sebat1/jsebat/tools/envs/gsem
/home/jsebat/bin/micromamba install -y -p $P -c conda-forge r-e1071 r-readr r-dplyr r-colorspace r-sfsmisc >> ~/rvh2/logs/gsem_env.log 2>&1
$P/bin/Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages(c("mgsub","simsalapar")); remotes::install_github("GenomicSEM/GenomicSEM", upgrade = "never", dependencies = FALSE); library(GenomicSEM); cat("GenomicSEM", as.character(packageVersion("GenomicSEM")), "OK\n")' >> ~/rvh2/logs/gsem_env.log 2>&1
tail -3 ~/rvh2/logs/gsem_env.log; grep -i "ERROR" ~/rvh2/logs/gsem_env.log | tail -3
