# S7b-c: misspecified locus effects. Background as S7a. One added locus per region:
#  S7b "sparse": effects on a random 10% of traits, iid, independent between regions (true rho = 0);
#  S7c "mediated": effects proportional to the genetic-correlation row of one mediator trait,
#      b = gamma * Rg[, med] (gamma ~ +-zloc); same mediator in both regions (|rho| = 1, sign of
#      gamma1 * gamma2) or different mediators (rho = Rg-implied alignment, mostly small).
# Reports the GLS rho-hat for the locus pair and rejection rates of the GLS and conditional tests.
# usage: s7bc_misspec.R <out.rds> <reps> <seed>
if (nzchar(Sys.getenv("RLIB"))) .libPaths(c(Sys.getenv("RLIB"), .libPaths()))
suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
args <- commandArgs(TRUE); out <- args[1]; reps <- as.integer(args[2]); set.seed(as.integer(args[3]))
one <- function(q, zloc, kind) {
  blocks <- make_ld_ar1(rep(100, 4), rho = 0.6); m <- 400; reg1 <- 1:200; reg2 <- 201:400
  w <- make_w(m, prop_nonzero = 0.1, dist = "lnorm")
  gc <- make_gcov(q, type = "cluster", rg = 0.6, n_clusters = 10)
  C <- make_intercept(q, overlap = 1, rp = 0.3 * (gc$Rg > 0))
  sim <- simulate_lrcpq(blocks, n = rep(3e5, q), w = w, Rb = diag(m), gcov = gc$gcov, intercept = C, M = 2e5)
  R1 <- as.matrix(bdiag(blocks[1:2])); R2 <- as.matrix(bdiag(blocks[3:4]))
  k1 <- sample(which(w[reg1] > 0), 1); k2 <- sample(which(w[reg2] > 0), 1)
  sparse <- function() { b <- numeric(q); s <- sample(q, ceiling(0.1 * q)); b[s] <- zloc * rnorm(length(s)); b }
  medv <- function(t) { u <- gc$Rg[, t]; u / sqrt(mean(u^2)) }
  g1 <- sample(c(-1, 1), 1); g2 <- sample(c(-1, 1), 1); truth <- 0
  if (kind == "sparse") { b1 <- sparse(); b2 <- sparse() }
  if (kind == "med_same") { t <- sample(q, 1); b1 <- g1 * zloc * medv(t); b2 <- g2 * zloc * medv(t); truth <- g1 * g2 }
  if (kind == "med_diff") { t <- sample(q, 2); b1 <- g1 * zloc * medv(t[1]); b2 <- g2 * zloc * medv(t[2])
    truth <- g1 * g2 * cor(medv(t[1]), medv(t[2])) }
  Z1 <- sim$Z[reg1, ] + outer(R1[, k1], b1); Z2 <- sim$Z[reg2, ] + outer(R2[, k2], b2)
  w1 <- w[reg1]; w2 <- w[reg2]; S1 <- which(w1 > 0); S2 <- which(w2 > 0)
  f <- lrcp_distal(Z1, Z2, R1, R2, w1, w2, sim$n, gc$gcov, 2e5, C, S1 = S1, S2 = S2, method = "gls", denominators = "R")
  g <- lrcp_gene(Z1[S1, ], Z2[S2, ], R1[S1, S1], R2[S2, S2], n = sim$n, gcov = gc$gcov, M = 2e5, intercept = C, n_boot = 0, n_sim = 0)
  i <- match(k1, S1); j <- match(k2, S2)
  data.frame(q = q, zloc = zloc, kind = kind, truth = truth, rho = f$rho[i, j], aligned = sign(f$rho[i, j]) == sign(truth),
             z_gls = (f$rho / f$se)[i, j], z_cond = g$z[i, j])
}
grid <- expand.grid(q = c(100, 300), zloc = c(3, 6), kind = c("sparse", "med_same", "med_diff"), stringsAsFactors = FALSE)
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid) * reps), function(k) {
  g <- grid[(k - 1) %% nrow(grid) + 1, ]
  tryCatch(suppressWarnings(one(g$q, g$zloc, g$kind)), error = function(e) NULL) }, mc.cores = 3))
saveRDS(res, out)
res$abs_truth <- abs(res$truth); res$rho_signed <- res$rho * sign(ifelse(res$truth == 0, 1, res$truth))
a <- aggregate(cbind(rej_gls = abs(z_gls) > 1.96, rej_cond = abs(z_cond) > 1.96, rho_signed, abs_truth) ~ q + zloc + kind, res, mean)
a$n <- aggregate(rho ~ q + zloc + kind, res, length)$rho
print(a, digits = 3)
