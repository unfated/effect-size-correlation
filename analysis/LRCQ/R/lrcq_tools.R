# Helpers for the LRCQ paper analyses, on top of the lrcpq package.
# Estimation itself uses lrcpq::lrcq_window(); this file adds: reading the
# UKB HapMap3 LD windows, LD pruning, clump-level LRCQ (rectangular design),
# and the comparator per-SNP scores.

suppressPackageStartupMessages({library(lrcpq); library(Matrix)})

# Read a window written by data/10_ukb_ld_subset.py.
read_window <- function(prefix, fix_psd = TRUE, eig_floor = 1e-4) {
  snps <- read.delim(paste0(prefix, ".snps.tsv"))
  m <- nrow(snps)
  R <- matrix(readBin(paste0(prefix, ".R.f32"), "numeric", m * m, size = 4), m, m)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  if (fix_psd) {
    e <- eigen(R, symmetric = TRUE)
    if (min(e$values) < eig_floor) {
      v <- pmax(e$values, eig_floor)
      R <- e$vectors %*% (v * t(e$vectors))
      d <- sqrt(diag(R)); R <- R / outer(d, d)
    }
  }
  list(snps = snps, R = R)
}

# Greedy LD pruning: keep SNPs in decreasing order of `priority` so that no
# two kept SNPs have r^2 >= thr. Returns kept indices (sorted) and, for every
# SNP, the kept SNP it is assigned to (highest r^2).
prune_tags <- function(R, thr = 0.5, priority = rowSums(R^2)) {
  m <- nrow(R)
  keep <- integer(0)
  blocked <- rep(FALSE, m)
  for (i in order(-priority)) {
    if (blocked[i]) next
    keep <- c(keep, i)
    blocked <- blocked | (R[i, ]^2 >= thr)
  }
  keep <- sort(keep)
  assign <- keep[max.col(R[, keep, drop = FALSE]^2, ties.method = "first")]
  list(keep = keep, assign = assign)
}

# Clump-level LRCQ in one window: responses z^2 - C for all m SNPs, parameters
# only for tag SNPs P: E[z_ka^2 - C_aa] = s_a (D[, P] w_P)_k. OLS or one-step
# WLS with weights 1 / (2 mu^2) (Isserlis, independent traits approximation
# for the weights only). Returns w_P.
lrcq_clump <- function(Z, R, P, n, h2, M, Cd, method = c("wls", "ols"),
                       ridge = 0, N_ref = Inf) {
  method <- match.arg(method)
  m <- nrow(Z); q <- ncol(Z)
  D <- ld_r2(R, N_ref)
  X <- D[, P, drop = FALSE]
  s <- n * h2 / M
  Y <- Z^2 - matrix(Cd, m, q, byrow = TRUE)
  solve1 <- function(V) {
    Sm <- matrix(s, m, q, byrow = TRUE)
    omega <- rowSums(Sm^2 / V)
    u <- rowSums(Y * Sm / V)
    A <- crossprod(X, X * omega)
    if (ridge > 0) A <- A + diag(ridge * mean(diag(A)), ncol(X))
    as.vector(solve(A, crossprod(X, u)))
  }
  w <- solve1(matrix(1, m, q))
  if (method == "wls") {
    mu <- outer(pmax(as.vector(X %*% w), 0), s) + matrix(Cd, m, q, byrow = TRUE)
    w <- solve1(2 * mu^2)
  }
  w
}

# Comparator per-SNP scores across traits.
#  mean_chi2   : mean over traits of z^2 - C_aa
#  ldnorm_chi2 : mean_chi2 / window LD score (LDSC-style normalisation)
#  n_sig       : number of traits with p < alpha
#  omnibus     : z' C^{-1} z / q, multi-trait Wald statistic using the
#                intercept (null correlation) matrix
comparators <- function(Z, Cmat, ldscore, alpha = 5e-8) {
  q <- ncol(Z)
  Cd <- diag(Cmat)
  mean_chi2 <- rowMeans(Z^2 - matrix(Cd, nrow(Z), q, byrow = TRUE))
  thr <- qnorm(alpha / 2, lower.tail = FALSE)
  Ci <- solve(Cmat)
  data.frame(mean_chi2 = mean_chi2,
             ldnorm_chi2 = mean_chi2 / ldscore,
             n_sig = rowSums(abs(Z) > thr),
             omnibus = rowSums((Z %*% Ci) * Z) / q)
}

# Rank-based AUC of score for a binary truth.
auc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth) || all(truth)) return(NA_real_)
  r <- rank(score)
  n1 <- sum(truth); n0 <- sum(!truth)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
