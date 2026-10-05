# S6: sensitivity to stage 1-2 plug-in errors. S1 generator (q traits, 10 clusters, rg 0.6,
# full overlap); estimation with (a) true plug-ins, (b) h2 scaled by exp(N(0, 0.2)) per trait,
# (c) intercept off-diagonals perturbed by U(-0.05, 0.05), (d) both, (e) genetic correlations
# ignored (diagonal gcov). GLS estimate / null test and conditional test.
.libPaths(c("/home/user/rlib2", .libPaths()))
suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
args <- commandArgs(TRUE); out <- args[1]; reps <- as.integer(args[2]); set.seed(as.integer(args[3]))
psd <- function(A) { e <- eigen((A + t(A)) / 2, symmetric = TRUE); e$vectors %*% (pmax(e$values, 1e-6) * t(e$vectors)) }
one <- function(q, rho) {
  blocks <- make_ld_ar1(rep(100, 4), rho = 0.6); m <- 400; reg1 <- 1:200; reg2 <- 201:400
  w <- make_w(m, prop_nonzero = 0.1, dist = "lnorm")
  nz1 <- intersect(which(w > 0), reg1); nz2 <- intersect(which(w > 0), reg2); pr <- cbind(sample(nz1, 3), sample(nz2, 3))
  Rb <- make_Rb(m, "pairs", pairs = pr, pair_rho = c(rho, rho, -rho))
  gc <- make_gcov(q, type = "cluster", rg = 0.6, n_clusters = 10); C <- make_intercept(q, overlap = 1, rp = 0.3 * (gc$Rg > 0))
  sim <- simulate_lrcpq(blocks, n = rep(3e5, q), w = w, Rb = Rb, gcov = gc$gcov, intercept = C, M = 2e5)
  R1 <- as.matrix(bdiag(blocks[1:2])); R2 <- as.matrix(bdiag(blocks[3:4])); Z1 <- sim$Z[reg1, ]; Z2 <- sim$Z[reg2, ]
  S1 <- which(w[reg1] > 0); S2 <- which(w[reg2] > 0); ti <- match(pr[, 1], S1); tj <- match(pr[, 2] - 200, S2)
  sc <- exp(rnorm(q, 0, 0.2)); E <- matrix(runif(q * q, -0.05, 0.05), q); E <- (E + t(E)) / 2; diag(E) <- 0
  G_h2 <- gc$gcov * sqrt(outer(sc, sc)); C_ic <- psd(C + E)
  variants <- list(true = list(gc$gcov, C), h2 = list(G_h2, C), intercept = list(gc$gcov, C_ic),
                   both = list(G_h2, C_ic), diag_gcov = list(diag(diag(gc$gcov)), C))
  do.call(rbind, lapply(names(variants), function(v) {
    G <- variants[[v]][[1]]; Ci <- variants[[v]][[2]]
    f <- suppressWarnings(lrcp_distal(Z1, Z2, R1, R2, w[reg1], w[reg2], sim$n, G, 2e5, Ci, S1 = S1, S2 = S2, method = "gls"))
    g <- suppressWarnings(lrcp_gene(Z1[S1, ], Z2[S2, ], R1[S1, S1], R2[S2, S2], n = sim$n, gcov = G, M = 2e5, intercept = Ci, n_sim = 0))
    nul <- matrix(TRUE, length(S1), length(S2)); nul[cbind(ti, tj)] <- FALSE
    data.frame(variant = v, q = q, rho = rho, est = mean(f$rho[cbind(ti, tj)] * c(1, 1, -1)),
               rej_gls = mean(abs(f$rho / f$se)[nul] > 1.96, na.rm = TRUE), rej_cond = mean(abs(g$z)[nul] > 1.96, na.rm = TRUE))
  }))
}
grid <- expand.grid(q = c(100, 300), rho = c(0.5))
R <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid) * reps), function(k) { g <- grid[(k - 1) %% nrow(grid) + 1, ]
  tryCatch(one(g$q, g$rho), error = function(e) { message(conditionMessage(e)); NULL }) }, mc.cores = 3))
saveRDS(R, out)
a <- aggregate(cbind(est, rej_gls, rej_cond) ~ variant + q, R, mean); a$sd_est <- aggregate(est ~ variant + q, R, sd)$est
print(a, digits = 3)
