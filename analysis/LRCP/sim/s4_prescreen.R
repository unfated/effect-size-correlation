# S4: prescreening on estimated enrichment (stage 3 LRCQ w-hat) before stage 4.
# Same generator as S1 (clustered traits, full overlap); w estimated per region
# with lrcq(); screen at several thresholds t; report screening sensitivity,
# bias of rho-hat for true pairs, null calibration and BH-FDR among screened pairs.
suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
args <- commandArgs(TRUE)
out <- args[1]; reps <- as.integer(args[2]); seed <- as.integer(args[3])
set.seed(seed)
ths <- c(0.5, 1, 2, 4)

prune <- function(S, w, R, r2 = 0.5) {   # keep highest-w SNP among candidates in r2 > threshold
  S <- S[order(-w[S])]; keep <- integer(0)
  for (s in S) if (!length(keep) || all(R[s, keep]^2 <= r2)) keep <- c(keep, s)
  sort(keep)
}

one_rep <- function(q = 300, rho = 0.5, ld_r = 0.6, M = 2e5) {
  blocks <- make_ld_ar1(rep(100, 4), rho = ld_r)
  m <- 400; reg1 <- 1:200; reg2 <- 201:400
  w <- make_w(m, prop_nonzero = 0.1, dist = "lnorm")
  nz1 <- intersect(which(w > 0), reg1); nz2 <- intersect(which(w > 0), reg2)
  pr <- cbind(sample(nz1, 3), sample(nz2, 3))
  Rb <- make_Rb(m, "pairs", pairs = pr, pair_rho = c(rho, rho, -rho))
  gc <- make_gcov(q, type = "cluster", rg = 0.6, n_clusters = 10)
  C <- make_intercept(q, overlap = 1, rp = 0.3 * (gc$Rg > 0))
  sim <- simulate_lrcpq(blocks, n = rep(3e5, q), w = w, Rb = Rb, gcov = gc$gcov, intercept = C, M = M)
  R1 <- as.matrix(bdiag(blocks[1:2])); R2 <- as.matrix(bdiag(blocks[3:4]))
  Z1 <- sim$Z[reg1, ]; Z2 <- sim$Z[reg2, ]
  f1 <- lrcq(Z1, ld_from_blocks(blocks[1:2]), sim$n, gc$gcov, M = M, intercept = C, method = "wls", rectify = "C")
  f2 <- lrcq(Z2, ld_from_blocks(blocks[3:4]), sim$n, gc$gcov, M = M, intercept = C, method = "wls", rectify = "C")
  wh1 <- f1$w; wh2 <- f2$w
  rows <- list()
  for (t in ths) for (pr_r2 in c(1, 0.5)) {
    S1 <- screen_snps(wh1, t); S2 <- screen_snps(wh2, t)
    if (pr_r2 < 1) { S1 <- prune(S1, wh1, R1, pr_r2); S2 <- prune(S2, wh2, R2, pr_r2) }
    if (!length(S1) || !length(S2)) next
    fit <- lrcp_distal(Z1, Z2, R1, R2, wh1, wh2, sim$n, gc$gcov, M, C, S1 = S1, S2 = S2, method = "gls")
    ti <- match(pr[, 1], S1); tj <- match(pr[, 2] - 200, S2)
    hit <- !is.na(ti) & !is.na(tj)
    p <- 2 * pnorm(-abs(fit$z)); padj <- matrix(p.adjust(p, "BH"), nrow(p))
    truemask <- matrix(FALSE, length(S1), length(S2)); truemask[cbind(ti[hit], tj[hit])] <- TRUE
    # "near-true": a screened SNP in LD r2>0.5 with a true pair SNP on both sides
    near <- outer(sapply(S1, function(s) any(R1[s, pr[, 1]]^2 > 0.5)), sapply(S2, function(s) any(R2[s, pr[, 2] - 200]^2 > 0.5)))
    nullmask <- !near
    rows[[length(rows) + 1]] <- data.frame(t = t, prune_r2 = pr_r2, n1 = length(S1), n2 = length(S2), n_tests = length(p),
      pairs_screened = sum(hit), mean_rho_hat_true = if (any(hit)) mean(fit$rho[cbind(ti[hit], tj[hit])] * sign(c(rho, rho, -rho)[hit])) else NA,
      typeI_null = mean(p[nullmask] < 0.05), n_disc = sum(padj < 0.05),
      false_disc = sum(padj < 0.05 & nullmask), true_disc = sum(padj < 0.05 & truemask))
  }
  do.call(rbind, rows)
}
res <- do.call(rbind, parallel::mclapply(seq_len(reps), function(r) { set.seed(seed * 1e4 + r); x <- one_rep(); x$rep <- r; x },
                                         mc.cores = as.integer(Sys.getenv("NCORES", "3"))))
saveRDS(res, out)
agg <- aggregate(cbind(n_tests, pairs_screened, mean_rho_hat_true, typeI_null, n_disc, false_disc, true_disc) ~ t + prune_r2, res, mean)
agg$FDP <- aggregate(cbind(false_disc, n_disc) ~ t + prune_r2, res, sum)[, 3] / pmax(aggregate(cbind(false_disc, n_disc) ~ t + prune_r2, res, sum)[, 4], 1)
print(agg)
