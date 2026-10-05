#!/usr/bin/env Rscript
# Stage 3 on Pan-UKB: tag-set (clump-level) LRCQ in every UKB LD window.
#
# Usage: 30_stage3_lrcq.R <zprefix> <ldsc_prefix> <ld_dir> <out_tsv> [method] [N_ref] [trait_ids_file]
#   zprefix     : output of 20_build_zmatrix.py (.Z.f32 + .traits.tsv)
#   ldsc_prefix : output of 21_ldsc_matrix.py (.h2.tsv, .intercept.tsv, .gcov.tsv)
#   method      : wls (default) | equal | ols
#   N_ref       : LD reference size for the r^2 bias correction (337000 for UKB
#                 in-sample windows, 503 for 1000G EUR)
#   trait_ids_file : optional file with one trait_id per line (trait subset; "all" or "-" = all)
#   n_clusters  : trait clusters for the delete-one-cluster jackknife (default 20;
#                 ols/equal only). Clusters: average-linkage clustering of 1 - |r_g|.
#                 Replicates (raw scale, one column per cluster) go to <out>.jk.f32,
#                 rows in the order of <out>.
# Diagonal h2 and intercepts are Pan-UKB's own univariate LDSC values
# (h2_panukb_ldsc, intercept_panukb); off-diagonal c_ab and h_ab come from the
# matrix LDSC fit. Scale: s_a = n_a h2_a / M with M = M_5_50, matching the
# LD scores used for h2; w is renormalised to mean 1 afterwards, which removes
# any common misscaling (theory S1.7).
# Cores: middle 2 Mb of each 3 Mb window ([start+0.5, start+2.5) Mb), extended
# to the chromosome ends for the first and last window. MHC windows skipped.
suppressPackageStartupMessages(library(lrcpq))
a <- commandArgs(TRUE)
zp <- a[1]; lp <- a[2]; ld_dir <- a[3]; out <- a[4]
method <- if (length(a) >= 5) a[5] else "wls"
N_ref <- if (length(a) >= 6) as.numeric(a[6]) else 337000
subset_file <- if (length(a) >= 7 && !(a[7] %in% c("all", "-"))) a[7] else NA
K <- if (length(a) >= 8) as.integer(a[8]) else 20L
here <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(here, "R", "lrcq_tools.R"))

traits_all <- read.delim(paste0(zp, ".traits.tsv"), stringsAsFactors = FALSE)
q_all <- nrow(traits_all)
cols <- seq_len(q_all)
if (!is.na(subset_file)) cols <- which(traits_all$trait_id %in% readLines(subset_file))
traits <- traits_all[cols, ]
q <- nrow(traits)
h2t <- read.delim(paste0(lp, ".h2.tsv"), stringsAsFactors = FALSE)
h2t <- h2t[match(traits$trait_id, h2t$trait_id), ]
stopifnot(identical(h2t$trait_id, traits$trait_id))
rd <- function(f) as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
C <- rd(paste0(lp, ".intercept.tsv"))[traits$trait_id, traits$trait_id]
G <- rd(paste0(lp, ".gcov.tsv"))[traits$trait_id, traits$trait_id]
diag(C) <- h2t$intercept_panukb
# rescale genetic covariances so the diagonal matches the Pan-UKB h2
d_old <- sqrt(pmax(diag(G), 1e-6)); d_new <- sqrt(pmax(h2t$h2_panukb_ldsc, 1e-6))
G <- G / outer(d_old, d_old) * outer(d_new, d_new)
G <- (G + t(G)) / 2
e <- eigen(G, symmetric = TRUE); G <- e$vectors %*% (pmax(e$values, 1e-8) * t(e$vectors))
e <- eigen(C, symmetric = TRUE); C <- e$vectors %*% (pmax(e$values, 1e-4) * t(e$vectors))
n <- traits$N
M <- as.numeric(readLines("/home/user/data/ref/UKBB.EUR.l2.M_5_50"))

uni <- read.delim("/home/user/data/ref/snp_universe.tsv", stringsAsFactors = FALSE)
m_all <- nrow(uni)
zcon <- file(paste0(zp, ".Z.f32"), "rb")
read_rows <- function(ui) {              # universe rows ui (sorted), trait subset cols
  i0 <- min(ui); i1 <- max(ui)
  seek(zcon, (i0 - 1) * q_all * 4)
  X <- matrix(readBin(zcon, "numeric", (i1 - i0 + 1) * q_all, size = 4), ncol = q_all, byrow = TRUE)
  X[ui - i0 + 1, cols, drop = FALSE]
}
Ci <- solve(C)
# trait clusters for the jackknife
RG <- rd(paste0(lp, ".rg.tsv"))[traits$trait_id, traits$trait_id]
RG[!is.finite(RG)] <- 0
RG <- pmin(pmax(RG, -1), 1)
cl <- cutree(hclust(as.dist(1 - abs(RG)), "average"), k = min(K, q))
K <- max(cl)
cat(sprintf("%d traits, %d jackknife clusters (sizes %s)\n", q, K, paste(sort(tabulate(cl), TRUE), collapse = ",")))
write.table(data.frame(trait_id = traits$trait_id, cluster = cl), sub("\\.tsv$", ".clusters.tsv", out),
            sep = "\t", quote = FALSE, row.names = FALSE)
s_all <- n * diag(G) / M
v_all <- if (method == "equal") 1 / s_all else s_all     # u = sum_a v_a y_a, omega = sum_a v_a s_a
jkf <- sub("\\.tsv$", ".jk.f32", out)

files <- list.files(ld_dir, "\\.snps\\.tsv$")
nm <- sub("\\.snps\\.tsv$", "", files)
pr <- do.call(rbind, regmatches(nm, regexec("chr([0-9]+)_([0-9]+)_([0-9]+)", nm)))
win <- data.frame(name = nm, chr = as.integer(pr[, 2]), start = as.numeric(pr[, 3]) - 1,
                  end = as.numeric(pr[, 4]) - 1)
win <- win[order(win$chr, win$start), ]
win$core_lo <- win$start + 0.5e6; win$core_hi <- win$start + 2.5e6
first <- !duplicated(win$chr); last <- !duplicated(win$chr, fromLast = TRUE)
win$core_lo[first] <- 0; win$core_hi[last] <- Inf
win <- win[!(win$chr == 6 & win$end > 25e6 & win$start < 34e6), ]

res <- list()
done <- if (file.exists(out)) unique(read.delim(out)$window) else character(0)
for (i in seq_len(nrow(win))) {
  if (win$name[i] %in% done) next
  t0 <- Sys.time()
  wd <- read_window(file.path(ld_dir, win$name[i]))
  ui <- match(wd$snps$ID, uni$ID)
  stopifnot(!anyNA(ui), all(diff(ui) > 0))
  Z <- read_rows(ui)
  # missing Z: neutral imputation z^2 = c_aa + s_a * l_k (w = 1); flagged in output
  miss <- is.na(Z)
  if (any(miss)) {
    lk <- rowSums(wd$R^2)
    E2 <- outer(lk, n * diag(G) / M) + matrix(diag(C), nrow(Z), q, byrow = TRUE)
    Z[miss] <- sqrt(E2[miss])
  }
  tags <- prune_tags(wd$R, 0.5)
  f <- lrcq_window(Z, wd$R, n, G, M, C, method = method, N_ref = N_ref,
                   tags = tags$keep, se = "model")
  # delete-one-trait-cluster jackknife (closed form for ols / equal weights)
  jk <- NULL
  if (method %in% c("ols", "equal")) {
    D <- lrcpq:::ld_r2(wd$R, N_ref)
    e <- eigen(D, symmetric = TRUE)
    if (min(e$values) < 0) D <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
    DT <- D[, tags$keep, drop = FALSE]
    Hm <- solve(crossprod(DT), t(DT))                       # p x m
    Yv <- (Z^2 - matrix(diag(C), nrow(Z), q, byrow = TRUE)) * matrix(v_all, nrow(Z), q, byrow = TRUE)
    Uc <- sapply(seq_len(K), function(k) rowSums(Yv[, cl == k, drop = FALSE]))   # m x K
    oc <- sapply(seq_len(K), function(k) sum((v_all * s_all)[cl == k]))
    full <- as.vector(Hm %*% rowSums(Uc)) / sum(oc)
    stopifnot(max(abs(full - f$w)) < 1e-6 * max(1, max(abs(f$w))))
    jk <- (Hm %*% (rowSums(Uc) - Uc)) / matrix(sum(oc) - oc, nrow(Hm), K, byrow = TRUE)
  }
  core <- wd$snps$BP >= win$core_lo[i] & wd$snps$BP < win$core_hi[i]
  ct <- core[tags$keep]
  clump_n <- tabulate(match(tags$assign[core], tags$keep), length(tags$keep))
  cmp <- comparators(Z, C, rowSums(wd$R^2))
  P <- tags$keep[ct]
  if (!length(P)) {
    cat(sprintf("%s: %d SNPs, no core tags, skipped\n", win$name[i], nrow(Z)))
    next
  }
  res <- data.frame(window = win$name[i], CHR = wd$snps$CHR[P], BP = wd$snps$BP[P],
                    ID = wd$snps$ID[P], RSID = wd$snps$RSID[P],
                    w_raw = f$w[ct], se = f$se[ct],
                    se_jk = if (is.null(jk)) NA else sqrt((K - 1) / K * rowSums((jk[ct, , drop = FALSE] - rowMeans(jk[ct, , drop = FALSE]))^2)),
                    clump_n = clump_n[ct],
                    l2_window = rowSums(wd$R^2)[P], n_missing = rowSums(miss)[P],
                    cmp[P, , drop = FALSE])
  write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE,
              append = file.exists(out), col.names = !file.exists(out))
  if (!is.null(jk)) {
    jc <- file(jkf, "ab"); writeBin(as.vector(t(jk[ct, , drop = FALSE])), jc, size = 4); close(jc)
  }
  # SNP -> tag map for the core SNPs
  map <- data.frame(ID = wd$snps$ID[core], tag_ID = wd$snps$ID[tags$assign[core]])
  mf <- sub("\\.tsv$", ".snp2tag.tsv", out)
  write.table(map, mf, sep = "\t", quote = FALSE, row.names = FALSE,
              append = file.exists(mf), col.names = !file.exists(mf))
  cat(sprintf("%s: %d SNPs, %d core tags, %.1fs\n", win$name[i], nrow(Z), sum(ct),
              as.numeric(Sys.time() - t0, units = "secs")))
}
close(zcon)
