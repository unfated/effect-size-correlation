# S7a: non-separable effects under the null. Background as S1 (clustered traits, overlap),
# Rb = I, plus one "locus-specific" SNP per region whose effects fall only on the traits of
# one cluster (independent between the two regions, so true rho = 0). Compares calibration of
# the lrcp_distal GLS null-SE test and the lrcp_gene conditional test for the locus pair.
suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
args <- commandArgs(TRUE); out <- args[1]; reps <- as.integer(args[2]); set.seed(as.integer(args[3]))
one <- function(q, zloc, same_cluster) {
  blocks <- make_ld_ar1(rep(100, 4), rho = 0.6); m <- 400; reg1 <- 1:200; reg2 <- 201:400
  w <- make_w(m, prop_nonzero = 0.1, dist = "lnorm")
  gc <- make_gcov(q, type = "cluster", rg = 0.6, n_clusters = 10)
  C <- make_intercept(q, overlap = 1, rp = 0.3 * (gc$Rg > 0))
  sim <- simulate_lrcpq(blocks, n = rep(3e5, q), w = w, Rb = diag(m), gcov = gc$gcov, intercept = C, M = 2e5)
  R1 <- as.matrix(bdiag(blocks[1:2])); R2 <- as.matrix(bdiag(blocks[3:4]))
  cl <- rep(1:10, length.out = q)  # cluster labels used by make_gcov are block-contiguous; derive from Rg
  cl <- cutree(hclust(as.dist(1 - gc$Rg)), h = 0.5)
  k1 <- sample(which(w[reg1] > 0), 1); k2 <- sample(which(w[reg2] > 0), 1)
  tr1 <- cl == 1; tr2 <- if (same_cluster) cl == 1 else cl == 2
  # locus effects: shared within the cluster (r = 0.8 among its traits), independent across regions
  mk <- function(tr) { b <- numeric(q); nt <- sum(tr); S <- matrix(0.8, nt, nt); diag(S) <- 1
    b[tr] <- zloc * as.vector(t(chol(S)) %*% rnorm(nt)); b }
  Z1 <- sim$Z[reg1, ] + outer(R1[, k1], mk(tr1)); Z2 <- sim$Z[reg2, ] + outer(R2[, k2], mk(tr2))
  w1 <- w[reg1]; w2 <- w[reg2]; S1 <- which(w1 > 0); S2 <- which(w2 > 0)
  f <- lrcp_distal(Z1, Z2, R1, R2, w1, w2, sim$n, gc$gcov, 2e5, C, S1 = S1, S2 = S2, method = "gls", denominators = "R")
  g <- lrcp_gene(Z1[S1, ], Z2[S2, ], R1[S1, S1], R2[S2, S2], n = sim$n, gcov = gc$gcov, M = 2e5, intercept = C, n_sim = 0)
  i <- match(k1, S1); j <- match(k2, S2)
  data.frame(q = q, zloc = zloc, same = same_cluster, z_gls = (f$rho / f$se)[i, j], z_cond = g$z[i, j],
             rej_gls_other = mean(abs(f$rho / f$se)[-i, -j] > 1.96, na.rm = TRUE), rej_cond_other = mean(abs(g$z[-i, -j]) > 1.96))
}
grid <- expand.grid(q = c(100, 300), zloc = c(3, 6), same = c(TRUE, FALSE))
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid) * reps), function(k) {
  g <- grid[(k - 1) %% nrow(grid) + 1, ]
  tryCatch(suppressWarnings(one(g$q, g$zloc, g$same)), error = function(e) NULL) }, mc.cores = 3))
saveRDS(res, out)
a <- aggregate(cbind(gls = abs(z_gls) > 1.96, cond = abs(z_cond) > 1.96, rej_gls_other, rej_cond_other) ~ q + zloc + same, res, mean)
a$sd_z_gls <- aggregate(z_gls ~ q + zloc + same, res, sd)$z_gls; a$sd_z_cond <- aggregate(z_cond ~ q + zloc + same, res, sd)$z_cond
print(a, digits = 3)
