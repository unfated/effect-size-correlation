#!/usr/bin/env Rscript
# Calibration of the delete-one-trait-cluster jackknife SE used in the Pan-UKB
# application (30_stage3_lrcq.R), against model SEs (Theorem 3.10), for the
# tag-set OLS estimator. Same generator and windows as run_sim.R.
# Usage: check_jackknife.R <scenario: base|ukb|sparse> <reps> <out_tsv>
suppressPackageStartupMessages(library(lrcpq))
a <- commandArgs(TRUE); scen <- a[1]; reps <- as.integer(a[2]); out <- a[3]
here <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(here, "R", "lrcq_tools.R"))
ld_dir <- "/home/user/data/ld/hm3"
wn <- c("chr22_28000001_31000001", "chr1_2000001_5000001", "chr1_1_3000001")
Rl <- lapply(wn, function(w) read_window(file.path(ld_dir, w))$R)
sizes <- vapply(Rl, nrow, 1L); m <- sum(sizes); off <- c(0, cumsum(sizes))
tags <- lapply(Rl, prune_tags, thr = 0.5)
q <- 300; n <- rep(3e5, q); M <- 1.1e6
prop <- if (scen == "sparse") 0.005 else 0.02
res <- list()
for (r in seq_len(reps)) {
  set.seed(7e5 + r)
  if (scen == "ukb") {
    gc <- make_gcov(q, type = "cluster", rg = 0.5, n_clusters = 5)
    C <- make_intercept(q, overlap = 1, rp = 0.2)
  } else {
    gc <- make_gcov(q, type = "identity"); C <- make_intercept(q)
  }
  w <- make_w(m, idx = sample.int(m, max(1, round(prop * m))), sdlog = 1)
  Z <- simulate_lrcpq(Rl, n, w = w, gcov = gc$gcov, intercept = C, M = M)$Z
  s <- n * gc$h2 / M
  rg <- cov2cor(gc$gcov)
  for (K in c(5, 20, 50)) {
    cl <- cutree(hclust(as.dist(1 - abs(rg)), "average"), k = K)
    if (max(cl) < K) cl <- (seq_len(q) - 1) %% K + 1 + 0 * cl   # identical rg: assign round-robin
    for (b in seq_along(Rl)) {
      ix <- (off[b] + 1):off[b + 1]; P <- tags[[b]]$keep
      f <- lrcq_window(Z[ix, ], Rl[[b]], n, gc$gcov, M, C, method = "ols", tags = P, se = "model")
      D <- Rl[[b]]^2; DT <- D[, P, drop = FALSE]
      Hm <- solve(crossprod(DT), t(DT))
      Yv <- sweep(Z[ix, ]^2, 2, diag(C)) * matrix(s, length(ix), q, byrow = TRUE)
      Uc <- sapply(seq_len(K), function(k) rowSums(Yv[, cl == k, drop = FALSE]))
      oc <- sapply(seq_len(K), function(k) sum(s[cl == k]^2))
      jk <- (Hm %*% (rowSums(Uc) - Uc)) / matrix(sum(oc) - oc, length(P), K, byrow = TRUE)
      se_jk <- sqrt((K - 1) / K * rowSums((jk - rowMeans(jk))^2))
      truth <- as.vector(solve(crossprod(DT), crossprod(DT, D %*% w[ix])))
      res[[length(res) + 1]] <- data.frame(scen = scen, rep = r, K = K, block = b, truth = truth,
                                           est = f$w, se_model = f$se, se_jk = se_jk)
    }
  }
  cat(scen, "rep", r, "\n")
}
d <- do.call(rbind, res)
summ <- do.call(rbind, lapply(split(d, d$K), function(x) {
  null <- x$truth < 0.05
  data.frame(scen = scen, K = x$K[1],
             zsd_model = sd((x$est - x$truth) / x$se_model), zsd_jk = sd((x$est - x$truth) / x$se_jk),
             zsd_model_null = sd((x$est - x$truth)[null] / x$se_model[null]),
             zsd_jk_null = sd((x$est - x$truth)[null] / x$se_jk[null]),
             typeI_model = mean((x$est / x$se_model)[null] > qnorm(0.95)),
             typeI_jk = mean((x$est / x$se_jk)[null] > qt(0.95, x$K[1] - 1)),
             cover_model = mean(abs(x$est - x$truth) < 1.96 * x$se_model),
             cover_jk = mean(abs(x$est - x$truth) < qt(0.975, x$K[1] - 1) * x$se_jk),
             ratio_med = median(x$se_jk / x$se_model))
}))
print(summ)
write.table(summ, out, sep = "\t", quote = FALSE, row.names = FALSE, append = file.exists(out), col.names = !file.exists(out))
