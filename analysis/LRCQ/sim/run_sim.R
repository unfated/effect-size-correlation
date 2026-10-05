#!/usr/bin/env Rscript
# LRCQ simulation study (work plan steps 6-9).
#
# Each replicate: K real UKB in-sample LD windows (HapMap3 SNPs) treated as
# independent blocks; B ~ matrix normal with SNP enrichment w (and optional
# genetic-effect correlation), Z = sqrt(n) R B + noise with intercept matrix C
# (lrcpq::simulate_lrcpq). Per-SNP signal uses a genome-scale M.
# Estimators per window:
#   tag   : tag-set (clump-level) LRCQ, tags greedy-pruned at r^2 < 0.5,
#           responses from all SNPs (primary; theory S3.11)
#   ridge : per-SNP LRCQ with ridge 1e-3 (lrcpq::lrcq_window)
#   + OLS / equal-weight versions of the tag estimator, rectifications A/B/C
#   + comparators: mean chi2, LD-normalised mean chi2, n significant traits,
#     omnibus Wald
#   + annotation-level: category means of tag estimates vs pooled
#     multi-trait s-LDSC style regression with the true categories.
# Usage: run_sim.R <scenario> <rep_from> <rep_to> <outdir>
suppressPackageStartupMessages(library(lrcpq))
args <- commandArgs(TRUE)
scen_name <- args[1]; r0 <- as.integer(args[2]); r1 <- as.integer(args[3]); outdir <- args[4]
here <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(here, "R", "lrcq_tools.R"))

ld_dir <- Sys.getenv("LRCQ_LD_DIR", "/home/user/data/ld/hm3")
win_names <- strsplit(Sys.getenv("LRCQ_SIM_WINDOWS",
  "chr22_28000001_31000001,chr1_2000001_5000001,chr1_1_3000001"), ",")[[1]]

# ---- scenario definitions --------------------------------------------------
base <- list(q = 300, n = 3e5, M = 1.1e6, prop = 0.02, sdlog = 1,
             gtype = "identity", rg = 0, overlap = 0, rp = 0,
             int_true = 1, int_used = 1, rb = "none", rb_rho = 0.5,
             het = FALSE, N_ref = Inf, ridge = 1e-3, tag_r2 = 0.5)
scenarios <- list(
  base = list(),
  q30 = list(q = 30), q100 = list(q = 100), q1000 = list(q = 1000),
  sparse = list(prop = 0.005), dense = list(prop = 0.1), inf = list(prop = 1, sdlog = 0.5),
  ukb = list(gtype = "cluster", rg = 0.5, overlap = 1, rp = 0.2),
  strat = list(int_true = 1.1, int_used = 1.1),
  strat_miss = list(int_true = 1.1, int_used = 1.0),
  rb_local = list(rb = "local"), rb_cross = list(rb = "cross"),
  het = list(het = TRUE),
  ref500 = list(N_ref = 500),
  tag08 = list(tag_r2 = 0.8), tag02 = list(tag_r2 = 0.2)
)
sc <- modifyList(base, scenarios[[scen_name]])

wins <- lapply(win_names, function(w) read_window(file.path(ld_dir, w)))
Rl <- lapply(wins, `[[`, "R")
sizes <- vapply(Rl, nrow, 1L); m <- sum(sizes)
blk <- rep(seq_along(sizes), sizes)
tags <- lapply(Rl, prune_tags, thr = sc$tag_r2)

# noisy reference LD (e.g. a 1000G-sized panel) for the ref500 scenario
ref_ld <- function(R, N) {
  if (!is.finite(N)) return(R)
  X <- matrix(rnorm(N * nrow(R)), N) %*% chol(R + diag(1e-8, nrow(R)))
  Rh <- cor(X)
  e <- eigen(Rh, symmetric = TRUE); v <- pmax(e$values, 1e-4)
  Rh <- e$vectors %*% (v * t(e$vectors)); d <- sqrt(diag(Rh)); Rh / outer(d, d)
}

# synthetic functional annotation: runs of 20 SNPs assigned to 4 categories
make_annot <- function(m) {
  runs <- ceiling(m / 20)
  rep(sample(1:4, runs, TRUE, prob = c(0.7, 0.15, 0.1, 0.05)), each = 20)[1:m]
}
enrich_mult <- c(1, 3, 6, 12)   # relative probability of non-zero w by category

one_rep <- function(rep) {
  set.seed(1e5 * match(scen_name, names(scenarios)) + rep)
  q <- sc$q; n <- rep(sc$n, q)
  annot <- make_annot(m)
  p_nz <- pmin(1, sc$prop * enrich_mult[annot] / mean(enrich_mult[annot]))
  nz <- which(runif(m) < p_nz); if (!length(nz)) nz <- sample.int(m, 1)
  w <- make_w(m, idx = nz, sdlog = sc$sdlog)
  gc <- make_gcov(q, type = sc$gtype, rg = sc$rg, n_clusters = 5)
  C <- make_intercept(q, overlap = sc$overlap, rp = sc$rp, diag_intercept = sc$int_true)
  Rb <- NULL
  if (sc$rb == "local") {      # clusters of 5 non-null SNPs inside the same window
    cl <- unlist(lapply(split(nz, blk[nz]), function(v) split(v, ceiling(seq_along(v) / 5))), recursive = FALSE)
    Rb <- make_Rb(m, "cluster", clusters = cl, rho = sc$rb_rho, signs = TRUE)
  } else if (sc$rb == "cross") { # pairs of non-null SNPs in different windows
    a <- nz[blk[nz] == 1]; b <- nz[blk[nz] == 2]; k <- min(length(a), length(b))
    if (k > 0) Rb <- make_Rb(m, "pairs", pairs = cbind(a[1:k], b[1:k]), pair_rho = sc$rb_rho)
  }
  if (sc$het) {                 # half the traits draw an independent w
    w2 <- make_w(m, idx = sample(nz), sdlog = sc$sdlog)
    half <- seq_len(q) <= q / 2
    s1 <- simulate_lrcpq(Rl, n[half], w = w, gcov = gc$gcov[half, half], intercept = C[half, half], M = sc$M)
    s2 <- simulate_lrcpq(Rl, n[!half], w = w2, gcov = gc$gcov[!half, !half], intercept = C[!half, !half], M = sc$M)
    Z <- cbind(s1$Z, s2$Z)
    sa <- (n * gc$h2 / sc$M)
    w_pw <- (sum(sa[half]^2) * w + sum(sa[!half]^2) * w2) / sum(sa^2)   # power-weighted target
    w_eq <- (w + w2) / 2                                                # equal-weight target
  } else {
    Z <- simulate_lrcpq(Rl, n, w = w, Rb = Rb, gcov = gc$gcov, intercept = C, M = sc$M)$Z
    w_pw <- w_eq <- w
  }
  Cu <- C; diag(Cu) <- sc$int_used
  out <- list(); annot_out <- list()
  off <- c(0, cumsum(sizes))
  for (b in seq_along(Rl)) {
    ix <- (off[b] + 1):off[b + 1]
    Rt <- Rl[[b]]
    Re <- ref_ld(Rt, sc$N_ref)
    tg <- if (is.finite(sc$N_ref)) prune_tags(Re, sc$tag_r2) else tags[[b]]
    P <- tg$keep
    Dt <- Rt^2; DP <- Dt[, P, drop = FALSE]
    proj <- solve(crossprod(DP), crossprod(DP, Dt))      # tag-set estimand map
    t_true <- as.vector(proj %*% w_pw[ix])
    t_true_eq <- as.vector(proj %*% w_eq[ix])
    clump_total <- as.vector(tapply(w_pw[ix], factor(tg$assign, levels = P), sum))
    Zb <- Z[ix, , drop = FALSE]
    fit_tag <- lrcq_clump(Zb, Re, P, n, gc$h2, sc$M, diag(Cu), "wls", N_ref = sc$N_ref)
    fit_tag_ols <- lrcq_clump(Zb, Re, P, n, gc$h2, sc$M, diag(Cu), "ols", N_ref = sc$N_ref)
    # equal-weight tag estimator: average of per-trait OLS solutions
    yb <- sweep(Zb^2, 2, diag(Cu))
    De <- ld_r2(Re, sc$N_ref)[, P, drop = FALSE]
    G <- solve(crossprod(De), t(De))
    fit_tag_eq <- rowMeans(G %*% sweep(yb, 2, n * gc$h2 / sc$M, "/"))
    fr <- lrcq_window(Zb, Re, n, gc$gcov, sc$M, Cu, method = "wls", ridge = sc$ridge,
                      N_ref = sc$N_ref, se = "model")
    cmp <- comparators(Zb, Cu, rowSums(Re^2))
    # per-tag rows
    out[[b]] <- data.frame(block = b, snp = ix[P], w_true = w_pw[ix][P], t_true = t_true,
      t_true_eq = t_true_eq, clump_total = clump_total, tag = fit_tag, tag_ols = fit_tag_ols,
      tag_eq = fit_tag_eq, ridge = fr$w[P], ridge_se = fr$se[P],
      ridge_clump = as.vector(tapply(fr$w, factor(tg$assign, levels = P), sum)),
      cmp[P, ], annot = annot[ix][P])
    annot_out[[b]] <- data.frame(block = b, annot = annot[ix], w = w_pw[ix], ridge = fr$w,
                                 yb = I(yb), De_full = I(ld_r2(Re, sc$N_ref)))
  }
  tab <- do.call(rbind, out)
  # annotation level: (i) mean of ridge per-SNP estimates by category,
  # (ii) pooled multi-trait s-LDSC-style regression on category LD scores
  # (oracle: knows the categories), using all SNPs.
  s <- n * gc$h2 / sc$M
  Lc <- do.call(rbind, lapply(annot_out, function(a) sapply(1:4, function(c) as.matrix(a$De_full) %*% (a$annot == c))))
  Yall <- do.call(rbind, lapply(annot_out, function(a) as.matrix(a$yb)))
  ybar <- as.vector(Yall %*% s) / sum(s^2)
  tau <- as.vector(solve(crossprod(Lc), crossprod(Lc, ybar)))
  an <- do.call(rbind, lapply(annot_out, function(a) a[, c("annot", "w", "ridge")]))
  annot_tab <- data.frame(annot = 1:4, true = tapply(an$w, an$annot, mean),
                          lrcq_ridge = tapply(an$ridge, an$annot, mean), sldsc_pooled = tau)
  list(tab = tab, annot = annot_tab, scen = scen_name, rep = rep, m = m)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
for (r in r0:r1) {
  f <- file.path(outdir, sprintf("%s_rep%03d.rds", scen_name, r))
  if (file.exists(f)) next
  t0 <- Sys.time()
  res <- one_rep(r)
  saveRDS(res, f)
  cat(sprintf("%s rep %d done in %.1fs\n", scen_name, r, as.numeric(Sys.time() - t0, units = "secs")))
}
