# S1: distal pairs, main simulation (notes/simulation-design.md).
# Two LD regions on "different chromosomes" (2 AR(1) blocks of 100 SNPs each),
# 3 cross-region causal pairs with correlation rho, oracle prescreen (w > 0)
# or a top-k screen; estimators: lrcp_distal GLS / OLS, lrcp_mle (small
# screened sets), naive Z correlation (Fisher test on q-3 df), oracle.
suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
args <- commandArgs(TRUE)
out <- args[1]; reps <- as.integer(args[2]); seed <- as.integer(args[3])
cells <- if (length(args) > 3) args[4] else "main"
set.seed(seed)

one_rep <- function(q, rho, traits, ld_r, M, do_mle = FALSE) {
  blocks <- make_ld_ar1(c(100, 100, 100, 100), rho = ld_r)
  m <- 400; reg1 <- 1:200; reg2 <- 201:400
  w <- make_w(m, prop_nonzero = 0.1, dist = "lnorm")
  nz1 <- intersect(which(w > 0), reg1); nz2 <- intersect(which(w > 0), reg2)
  pr <- cbind(sample(nz1, 3), sample(nz2, 3))
  Rb <- make_Rb(m, "pairs", pairs = pr, pair_rho = c(rho, rho, -rho))
  if (traits == "indep") {
    gc <- make_gcov(q); C <- diag(q)
  } else {
    gc <- make_gcov(q, type = "cluster", rg = 0.6, n_clusters = 10)
    C <- make_intercept(q, overlap = 1, rp = 0.3 * (gc$Rg > 0))
  }
  sim <- simulate_lrcpq(blocks, n = rep(3e5, q), w = w, Rb = Rb, gcov = gc$gcov,
                        intercept = C, M = M)
  R1 <- as.matrix(bdiag(blocks[1:2])); R2 <- as.matrix(bdiag(blocks[3:4]))
  Z1 <- sim$Z[reg1, ]; Z2 <- sim$Z[reg2, ]
  w1 <- w[reg1]; w2 <- w[reg2]
  S1 <- which(w1 > 0); S2 <- which(w2 > 0)
  args <- list(Z1 = Z1, Z2 = Z2, R1 = R1, R2 = R2, w1 = w1, w2 = w2, n = sim$n,
               gcov = gc$gcov, M = M, intercept = C, S1 = S1, S2 = S2)
  fg <- do.call(lrcp_distal, c(args, method = "gls"))
  fo <- do.call(lrcp_distal, c(args, method = "ols"))
  ti <- match(pr[, 1], S1); tj <- match(pr[, 2] - 200, S2)
  truth <- c(rho, rho, -rho)
  isnull <- matrix(TRUE, length(S1), length(S2)); isnull[cbind(ti, tj)] <- FALSE
  naive <- cor(t(Z1[S1, , drop = FALSE]), t(Z2[S2, , drop = FALSE]))
  nz <- atanh(naive) * sqrt(q - 3)
  B <- as.matrix(sim$B)
  orc <- sapply(1:3, function(k) cor(B[pr[k, 1], ], B[pr[k, 2], ]))
  res <- data.frame(q = q, rho = rho, traits = traits, ld_r = ld_r, M = M, pair = 1:3, truth = truth,
                    gls = fg$rho[cbind(ti, tj)], gls_se = fg$se[cbind(ti, tj)],
                    ols = fo$rho[cbind(ti, tj)], ols_se = fo$se[cbind(ti, tj)],
                    naive = naive[cbind(ti, tj)], oracle = orc)
  if (do_mle) {
    # MLE on a small screened set: the 3 pair SNPs + 2 random other candidates per region
    s1 <- sort(unique(c(pr[, 1], sample(setdiff(nz1, pr[, 1]), 2))))
    s2 <- sort(unique(c(pr[, 2] - 200, sample(setdiff(nz2 - 200, pr[, 2] - 200), 2))))
    a2 <- args; a2$S1 <- s1; a2$S2 <- s2
    fm <- do.call(lrcp_mle, a2)
    res$mle <- fm$rho[cbind(match(pr[, 1], s1), match(pr[, 2] - 200, s2))]
    res$mle_se <- fm$se[cbind(match(pr[, 1], s1), match(pr[, 2] - 200, s2))]
  } else { res$mle <- NA; res$mle_se <- NA }
  nulls <- data.frame(q = q, rho = rho, traits = traits, ld_r = ld_r, M = M,
                      z_gls = fg$z[isnull], z_ols = fo$z[isnull], z_naive = nz[isnull])
  list(res = res, nulls = nulls)
}

grid <- if (cells == "main") {
  expand.grid(q = c(30, 100, 300), rho = c(0, 0.5, 0.8), traits = c("indep", "clustered"),
              ld_r = c(0.3, 0.6, 0.9), M = 2e5, stringsAsFactors = FALSE)
} else {
  expand.grid(q = c(30, 100, 300), rho = c(0.5), traits = c("clustered"), ld_r = 0.6, M = c(2e5, 1e6),
              stringsAsFactors = FALSE)
}
run_cell <- function(g) {
  p <- grid[g, ]
  set.seed(seed * 1000 + g)
  R <- list(); N <- list(); t0 <- Sys.time()
  for (r in seq_len(reps)) {
    o <- tryCatch(one_rep(p$q, p$rho, p$traits, p$ld_r, p$M, do_mle = (Sys.getenv("MLE") == "1" && r <= 30)),
                  error = function(e) { message("rep error: ", conditionMessage(e)); NULL })
    if (is.null(o)) next
    o$res$rep <- r; R[[length(R) + 1]] <- o$res
    if (r <= 100) N[[length(N) + 1]] <- o$nulls
  }
  message(sprintf("cell %d/%d q=%d rho=%.1f %s ld=%.1f M=%g done in %.0fs", g, nrow(grid), p$q, p$rho,
                  p$traits, p$ld_r, p$M, as.numeric(Sys.time() - t0, units = "secs")))
  list(res = do.call(rbind, R), nulls = do.call(rbind, N))
}
ncores <- as.integer(Sys.getenv("NCORES", "3"))
out_list <- parallel::mclapply(seq_len(nrow(grid)), run_cell, mc.cores = ncores, mc.preschedule = FALSE)
saveRDS(list(res = do.call(rbind, lapply(out_list, `[[`, "res")),
             nulls = do.call(rbind, lapply(out_list, `[[`, "nulls")), grid = grid), out)
