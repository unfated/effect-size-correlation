#!/usr/bin/env Rscript
# Phenotype-category contrast of annotation enrichment (theory S7.4, S7.4a):
# blood biomarkers ("Biological samples") vs assessment-centre traits.
# lrcpq::lrcq_category with contrast weights and the parametric-bootstrap p
# (the reported test); per-group power-weighted enrichments alongside.
# Windows: a random subset of non-overlapping 1000G windows (compute bound).
# Usage: 47_category_contrast.R <zprefix> <ldsc_prefix> <ld_dir> <annot_npz_ids_tsv> <out_tsv> [n_windows] [n_boot]
suppressPackageStartupMessages(library(lrcpq))
a <- commandArgs(TRUE)
zp <- a[1]; lp <- a[2]; ld_dir <- a[3]; anf <- a[4]; out <- a[5]
nw <- if (length(a) >= 6) as.integer(a[6]) else 100
nb <- if (length(a) >= 7) as.integer(a[7]) else 100
here <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(here, "R", "lrcq_tools.R"))
set.seed(2026)

tr <- read.delim(paste0(zp, ".traits.tsv"), stringsAsFactors = FALSE)
q_all <- nrow(tr)
grp <- ifelse(tr$category_top == "Biological samples", 1, ifelse(tr$category_top == "UK Biobank Assessment Centre", 2, NA))
cols <- which(!is.na(grp)); grp <- grp[cols]
h2t <- read.delim(paste0(lp, ".h2.tsv"), stringsAsFactors = FALSE)
h2t <- h2t[match(tr$trait_id, h2t$trait_id), ][cols, ]
rd <- function(f) as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
ids <- tr$trait_id[cols]
C <- rd(paste0(lp, ".intercept.tsv"))[ids, ids]; G <- rd(paste0(lp, ".gcov.tsv"))[ids, ids]
diag(C) <- h2t$intercept_panukb
d_old <- sqrt(pmax(diag(G), 1e-6)); d_new <- sqrt(pmax(h2t$h2_panukb_ldsc, 1e-6))
G <- G / outer(d_old, d_old) * outer(d_new, d_new); G <- (G + t(G)) / 2
e <- eigen(G, symmetric = TRUE); G <- e$vectors %*% (pmax(e$values, 1e-8) * t(e$vectors))
e <- eigen(C, symmetric = TRUE); C <- e$vectors %*% (pmax(e$values, 1e-4) * t(e$vectors))
n <- tr$N[cols]
M <- as.numeric(readLines("/home/user/data/ref/UKBB.EUR.l2.M_5_50"))
cat(sprintf("group sizes: biomarkers %d, assessment centre %d\n", sum(grp == 1), sum(grp == 2)))

uni <- read.delim("/home/user/data/ref/snp_universe.tsv", stringsAsFactors = FALSE)
an <- read.delim(anf, check.names = FALSE)          # ID + annotation columns (universe order)
stopifnot(identical(an$ID, uni$ID))
zcon <- file(paste0(zp, ".Z.f32"), "rb")
read_rows <- function(ui) {
  i0 <- min(ui); i1 <- max(ui); seek(zcon, (i0 - 1) * q_all * 4)
  X <- matrix(readBin(zcon, "numeric", (i1 - i0 + 1) * q_all, size = 4), ncol = q_all, byrow = TRUE)
  X[ui - i0 + 1, cols, drop = FALSE]
}
files <- sub("\\.snps\\.tsv$", "", list.files(ld_dir, "\\.snps\\.tsv$"))
st <- as.numeric(sub("chr[0-9]+_([0-9]+)_.*", "\\1", files))
chr <- as.integer(sub("chr([0-9]+)_.*", "\\1", files))
keep <- ((st - 1) %% 4e6 == 0) & !(chr == 6 & st > 22e6 & st < 34e6)
pick <- sort(sample(files[keep], min(nw, sum(keep))))
Rl <- list(); Zl <- list(); Al <- list()
for (w in pick) {
  wd <- read_window(file.path(ld_dir, w))
  ui <- match(wd$snps$ID, uni$ID)
  Z <- read_rows(ui)
  miss <- is.na(Z)
  if (any(miss)) {
    lk <- rowSums(wd$R^2)
    E2 <- outer(lk, n * diag(G) / M) + matrix(diag(C), nrow(Z), length(cols), byrow = TRUE)
    Z[miss] <- sqrt(E2[miss])
  }
  Rl[[w]] <- wd$R; Zl[[w]] <- Z; Al[[w]] <- as.matrix(an[ui, -1])
}
close(zcon)
sizes <- vapply(Rl, nrow, 1L)
bid <- rep(seq_along(sizes), sizes); offs <- c(0, cumsum(sizes))
ld <- new_ld(function(idx) {
  b <- unique(bid[idx]); stopifnot(length(b) == 1)
  Rl[[b]][idx - offs[b], idx - offs[b], drop = FALSE]
}, sum(sizes), bid, N = 503)
Z <- do.call(rbind, Zl); A <- do.call(rbind, Al)
cat(sprintf("%d windows, %d SNPs, %d annotations\n", length(pick), nrow(Z), ncol(A)))
res <- list()
for (g in 1:2) {
  u <- numeric(length(cols)); s <- n * diag(G) / M
  u[grp == g] <- s[grp == g] / sum(s[grp == g]^2)
  r <- lrcq_category(Z, ld, A, n, G, M, C, u = u, n_boot = 0)
  r$test <- c("biomarkers", "assessment_centre")[g]; res[[g]] <- r
}
t0 <- Sys.time()
u <- trait_weights(n, G, M, "contrast", groups = grp)
r <- lrcq_category(Z, ld, A, n, G, M, C, u = u, n_boot = nb)
r$test <- "contrast_bio_minus_ac"; res[[3]] <- r
cat(sprintf("contrast with %d bootstrap draws: %.1f min\n", nb, as.numeric(Sys.time() - t0, units = "mins")))
o <- do.call(rbind, lapply(res, function(x) { x$p_norm_diagnostic <- if (is.null(x$p_norm_diagnostic)) NA else x$p_norm_diagnostic; x }))
write.table(o, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(o)
