# S2: local pairs. One AR(1) LD block of 200 SNPs (neighbour r = 0.9); 10% causal SNPs plus a
# causal pair (i, j) whose LD r_ij = 0.9^k is set by its distance k; true rho in {0, 0.5}.
# LRCP-local (lrcp_local, oracle screen) vs naive Z correlation across traits. No SEs (empirical SD).
.libPaths(c("/home/user/rlib2", .libPaths()))
suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
args <- commandArgs(TRUE); out <- args[1]; reps <- as.integer(args[2]); set.seed(as.integer(args[3]))
one <- function(q, rho, r_ld) {
  m <- 200; blocks <- make_ld_ar1(m, rho = 0.9); R <- as.matrix(blocks[[1]])
  k <- if (r_ld == 0) 80 else round(log(r_ld) / log(0.9))
  i <- sample(20:(m - k - 20), 1); j <- i + k
  w <- make_w(m, prop_nonzero = 0.1, dist = "lnorm"); w[c(i, j)] <- mean(w[w > 0]); w <- w / mean(w)
  Rb <- make_Rb(m, "pairs", pairs = cbind(i, j), pair_rho = rho)
  gc <- make_gcov(q, type = "cluster", rg = 0.6, n_clusters = 10); C <- make_intercept(q, overlap = 1, rp = 0.3 * (gc$Rg > 0))
  sim <- simulate_lrcpq(blocks, n = rep(3e5, q), w = w, Rb = Rb, gcov = gc$gcov, intercept = C, M = 2e5)
  S <- which(w > 0)
  f <- suppressWarnings(lrcp_local(sim$Z, R, w, sim$n, gc$gcov, 2e5, C, S = S))
  est <- f$rho[(f$i == i & f$j == j) | (f$i == j & f$j == i)]
  if (!length(est)) { ii <- match(i, S); jj <- match(j, S); est <- f$rho[(f$i == ii & f$j == jj) | (f$i == jj & f$j == ii)] }
  data.frame(q = q, rho = rho, r_ld = round(R[i, j], 2), lrcp = est[1], naive = cor(sim$Z[i, ], sim$Z[j, ]))
}
grid <- expand.grid(q = 300, rho = c(0, 0.5), r_ld = c(0, 0.2, 0.5, 0.8))
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid) * reps), function(k) { g <- grid[(k - 1) %% nrow(grid) + 1, ]
  tryCatch(one(g$q, g$rho, g$r_ld), error = function(e) { message(conditionMessage(e)); NULL }) }, mc.cores = 3))
saveRDS(res, out)
a <- aggregate(cbind(lrcp, naive) ~ q + rho + r_ld, res, mean)
a$sd_lrcp <- aggregate(lrcp ~ q + rho + r_ld, res, sd)$lrcp; a$sd_naive <- aggregate(naive ~ q + rho + r_ld, res, sd)$naive
a$n <- aggregate(lrcp ~ q + rho + r_ld, res, length)$lrcp
print(a, digits = 3)
