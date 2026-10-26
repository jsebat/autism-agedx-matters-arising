export OPENBLAS_NUM_THREADS=1; P=/expanse/projects/sebat1/jsebat/tools/envs/gsem
/home/jsebat/bin/micromamba install -y -p $P -c conda-forge r-mgsub >> ~/rvh2/logs/gsem_env.log 2>&1 || /home/jsebat/bin/micromamba install -y -p $P -c conda-forge compilers make >> ~/rvh2/logs/gsem_env.log 2>&1
$P/bin/Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); if (!requireNamespace("mgsub", quietly=TRUE)) install.packages("mgsub"); if (!requireNamespace("simsalapar", quietly=TRUE)) install.packages("simsalapar"); remotes::install_github("GenomicSEM/GenomicSEM", upgrade = "never", dependencies = FALSE); library(GenomicSEM); cat("GenomicSEM", as.character(packageVersion("GenomicSEM")), "OK\n")' >> ~/rvh2/logs/gsem_env.log 2>&1
tail -2 ~/rvh2/logs/gsem_env.log
