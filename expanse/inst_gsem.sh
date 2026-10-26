# GenomicSEM env: conda-forge binaries for all compiled dependencies (the dnmt3a-r env has no compiler), GenomicSEM itself from GitHub
export OPENBLAS_NUM_THREADS=1
MM=/home/jsebat/micromamba/bin/micromamba; [ -x $MM ] || MM=$(ls /home/jsebat/bin/micromamba 2>/dev/null)
P=/expanse/projects/sebat1/jsebat/tools/envs/gsem
$MM create -y -p $P -c conda-forge r-base=4.3 r-lavaan r-gdata r-corpcor r-data.table r-matrix r-remotes r-doparallel r-foreach r-stringr r-splitstackshape r-r.utils r-mnormt r-gtools r-rcpp r-psych r-plyr r-mgcv > ~/rvh2/logs/gsem_env.log 2>&1
$P/bin/Rscript -e 'remotes::install_github("GenomicSEM/GenomicSEM", upgrade = "never", dependencies = FALSE); library(GenomicSEM); cat("GenomicSEM", as.character(packageVersion("GenomicSEM")), "OK\n")' >> ~/rvh2/logs/gsem_env.log 2>&1
tail -3 ~/rvh2/logs/gsem_env.log
