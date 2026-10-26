# Our own genomic factor analysis of SPARK autism (Matters Arising on Zhang, Grove et al. 2025).
# Indicators: the split-half GWAS of the broad survey (their fastGWA pipeline, disjoint controls) that pass an h2 filter, chosen by
# gsem_select.py (fa_traits.txt). GenomicSEM multivariable LDSC (Grotzinger 2019) -> genetic covariance S and sampling covariance V.
# 1) eigen-decomposition of the genetic correlation matrix and a parallel analysis against S + V-based sampling noise;
# 2) exploratory factor analysis (1-3 factors, promax) of the smoothed genetic correlation matrix;
# 3) (gsem_comp.py) what predicts each indicator's factor loadings: the composition of its cases (late diagnosis, adult
#    registration, female, deprivation, ADHD, ID, birth cohort, family structure).
# Output: gsem/S.tsv, rg.tsv, eigen.tsv, efa_loadings_k*.tsv
suppressMessages({library(GenomicSEM); library(data.table)})
C <- "/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/combined"; O <- file.path(C, "gsem"); dir.create(O, showWarnings = FALSE)
setwd(C); REF <- "/expanse/projects/sebat1/jsebat/tools/ldref/eur_w_ld_chr/"
tr <- readLines(file.path(O, "fa_traits.txt"))
files <- file.path("munged", paste0(tr, ".sumstats.gz"))
if (!file.exists(file.path(O, "ldsc.rds"))) {
  L <- ldsc(traits = files, sample.prev = rep(NA, length(tr)), population.prev = rep(NA, length(tr)), ld = REF, wld = REF,
            trait.names = tr, ldsc.log = file.path(O, "gsem_ldsc"))
  saveRDS(L, file.path(O, "ldsc.rds"))
} else L <- readRDS(file.path(O, "ldsc.rds"))
S <- L$S; dimnames(S) <- list(tr, tr); fwrite(data.table(trait = tr, S), file.path(O, "S.tsv"), sep = "\t")
R <- cov2cor(as.matrix(Matrix::nearPD(S, corr = FALSE)$mat)); dimnames(R) <- list(tr, tr)
fwrite(data.table(trait = tr, R), file.path(O, "rg.tsv"), sep = "\t")
# 1) eigenvalues + parallel analysis: draw S* ~ N(S0, V) under a one-factor null (S0 = rank-1 approximation of S)
ev <- eigen(R, symmetric = TRUE)$values
e1 <- eigen(S, symmetric = TRUE); S0 <- e1$values[1] * tcrossprod(e1$vectors[, 1])
k <- length(tr); idx <- which(lower.tri(S, diag = TRUE)); V <- as.matrix(Matrix::nearPD(L$V)$mat)
set.seed(1); nul <- replicate(500, { s <- S0[idx] + as.vector(t(chol(V)) %*% rnorm(length(idx))); M <- matrix(0, k, k); M[idx] <- s
  M <- M + t(M) - diag(diag(M)); M <- as.matrix(Matrix::nearPD(M)$mat); eigen(cov2cor(M), symmetric = TRUE, only.values = TRUE)$values })
pa <- data.table(component = seq_len(k), eigen_obs = ev, null_mean = rowMeans(nul), null_95 = apply(nul, 1, quantile, 0.95))
fwrite(pa, file.path(O, "eigen.tsv"), sep = "\t"); print(pa[1:min(6, k)])
# 2) EFA
for (nf in 1:3) {
  if (nf >= k / 2) next
  fa <- tryCatch(factanal(covmat = R, factors = nf, rotation = if (nf > 1) "promax" else "none", n.obs = 5000), error = function(e) NULL)
  if (is.null(fa)) next
  ld <- unclass(fa$loadings); fwrite(data.table(trait = tr, ld), file.path(O, sprintf("efa_loadings_k%d.tsv", nf)), sep = "\t")
  cat("\nEFA", nf, "factors\n"); print(round(ld, 2))
  if (!is.null(fa$rotmat)) { phi <- solve(t(fa$rotmat) %*% fa$rotmat); cat("factor correlations\n"); print(round(phi, 2)) }
}
cat("DONE\n")
