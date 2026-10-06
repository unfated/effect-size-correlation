# Genome-wide LRCP screen (work plan step 13): conditional test (theory S8.2) for every
# distal pair of LRCQ candidate tags (q_gt0 < 0.05), vectorised with the diagonal
# variance of lrcpq::lrcp_gene and the conservative orientation (larger of the two
# conditional variances). Bootstrap and class-null confirmation of hits is in confirm.R.
# usage: screen.R <candidate_tags.tsv> <tagZ.tsv.gz> <ld_dir> <out_dir>
.libPaths(c("/home/user/rlib2", .libPaths()))
suppressPackageStartupMessages(library(lrcpq))
args <- commandArgs(TRUE); CAND <- args[1]; ZF <- args[2]; LDD <- args[3]; OUT <- args[4]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
M <- 1094844; MIN_DIST <- 5e6
L2 <- "/mnt/project-files/papers/LRCQ/results/real/"
cand <- read.delim(CAND)
Z <- as.matrix(read.delim(gzfile(ZF), row.names = 1, check.names = FALSE))
stopifnot(identical(rownames(Z), cand$ID)); Z[is.na(Z)] <- 0
# robustness (step 14): optional stricter prescreen, e.g. MAXQ=0.01 keeps tags with q_gt0 < 0.01
if (nzchar(Sys.getenv("MAXQ"))) { k <- cand$q_gt0 < as.numeric(Sys.getenv("MAXQ")); cand <- cand[k, ]; Z <- Z[k, ]
  message("prescreen q_gt0 < ", Sys.getenv("MAXQ"), ": ", sum(k), " tags") }
tr <- read.delim(file.path(L2, "panukb.traits.tsv")); stopifnot(identical(tr$trait_id, colnames(Z)))
# robustness (step 14): optional trait subset, TRAITS=indep keeps the 150 'indep' traits of trait-list.tsv
if (Sys.getenv("TRAITS") == "indep") {
  tl <- read.delim(file.path(L2, "trait-list.tsv")); ind <- tl$trait_id[tl$indep %in% c(TRUE, "True", 1)]
  k <- tr$trait_id %in% ind; Z <- Z[, k]; tr <- tr[k, ]; message("trait subset indep: ", sum(k))
}
psd <- function(A, eps = 1e-6) { e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  e$vectors %*% (pmax(e$values, eps * max(e$values)) * t(e$vectors)) }
# plug-ins as in targeted/stage3.R (PLUGIN=ldsc): LRCQ genome-wide cross-trait LDSC
ids <- tr$trait_id
hh <- read.delim(file.path(L2, "panukb_ldsc.h2.tsv"))
rg <- as.matrix(read.delim(file.path(L2, "panukb_ldsc.rg.tsv"), row.names = 1, check.names = FALSE))[ids, ids]
ICm <- as.matrix(read.delim(file.path(L2, "panukb_ldsc.intercept.tsv"), row.names = 1, check.names = FALSE))[ids, ids]
rg[is.na(rg)] <- 0; diag(rg) <- 1; rg <- pmin(pmax(rg, -1), 1); ICm[is.na(ICm)] <- 0
h2 <- pmax(hh$h2_panukb_ldsc[match(ids, hh$trait_id)], 1e-3); ic <- pmax(hh$intercept_panukb[match(ids, hh$trait_id)], 1)
n <- tr$N_eff
gcov <- psd(rg * sqrt(outer(h2, h2))); diag(gcov) <- h2
intercept <- psd(ICm); diag(intercept) <- ic
s <- n * h2 / M; ss <- sum(s^2); Tm <- scale_matrix(n, gcov, M)
# per-window profiles: X = R^-1 Z, tauR = diag R^-1, MG = diag R^-1 G R^-1, v = local variance
m <- nrow(Z); X <- matrix(0, m, ncol(Z)); tauR <- MG <- v <- numeric(m)
for (win in unique(cand$window)) {
  S <- which(cand$window == win)
  sn <- read.delim(file.path(LDD, sprintf("w%s.snps.tsv", win)))
  ix <- match(cand$ID[S], sn$ID); stopifnot(!anyNA(ix))
  R <- matrix(readBin(file.path(LDD, sprintf("w%s.R.f32", win)), "numeric", nrow(sn)^2, size = 4), nrow(sn))[ix, ix, drop = FALSE]
  Zs <- Z[S, , drop = FALSE]; Ri <- solve(R)
  G <- local_moment(Zs, R, n, gcov, M, intercept)
  X[S, ] <- Ri %*% Zs; tauR[S] <- diag(Ri); MG[S] <- diag(Ri %*% G %*% Ri)
  v[S] <- pmax(local_variance(Zs, R, seq_along(S), n, gcov, M, intercept), 0)
}
message("profiles done for ", length(unique(cand$window)), " windows")
sY <- sweep(X, 2, s, "*")
KG <- rowSums((sY %*% intercept) * sY); KT <- rowSums((sY %*% Tm) * sY)
Cm <- tcrossprod(X, sY) / ss
VA <- (outer(tauR, KG) + outer(MG, KT)) / ss^2          # conditional on Z_j (orientation A)
Vb <- t(VA)                                             # conditional on Z_i
zc <- Cm / sqrt(pmax(VA, Vb)); rm(VA, Vb)
zA <- Cm / sqrt((outer(tauR, KG) + outer(MG, KT)) / ss^2)
distal <- outer(cand$CHR, cand$CHR, "!=") | abs(outer(cand$BP, cand$BP, "-")) >= MIN_DIST
ut <- which(upper.tri(zc) & distal, arr.ind = TRUE)
z <- zc[ut]; p <- 2 * pnorm(-abs(z)); q <- p.adjust(p, "BH")
lam <- median(z^2) / qchisq(0.5, 1)
summ <- data.frame(n_tags = m, n_pairs = length(z), lambda_median = lam, sd_z = sd(z),
  frac_p05 = mean(p < 0.05), frac_p1e3 = mean(p < 1e-3), frac_p1e6 = mean(p < 1e-6),
  n_q05 = sum(q < 0.05), n_q01 = sum(q < 0.01), n_bonf = sum(p < 0.05 / length(p)))
print(summ)
write.table(summ, file.path(OUT, "screen_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# QQ data (thinned)
o <- order(p); k <- unique(round(exp(seq(0, log(length(p)), length.out = 4000))))
write.table(data.frame(rank = k, p = p[o][k], exp = k / (length(p) + 1)), file.path(OUT, "screen_qq.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
keep <- which(p < 1e-3); Zs_ <- scale(t(Z))
i <- ut[keep, 1]; j <- ut[keep, 2]
den <- sqrt(v[i] * v[j])
res <- data.frame(ID1 = cand$ID[i], RSID1 = cand$RSID[i], win1 = cand$window[i], ID2 = cand$ID[j],
  RSID2 = cand$RSID[j], win2 = cand$window[j], C = Cm[cbind(i, j)], rho = ifelse(den > 0, Cm[cbind(i, j)] / den, NA),
  z = z[keep], z_A = zA[cbind(i, j)], p = p[keep], q_BH = q[keep],
  r_naive = colSums(Zs_[, i, drop = FALSE] * Zs_[, j, drop = FALSE]) / (ncol(Z) - 1))
res <- res[order(res$p), ]
write.table(res, gzfile(file.path(OUT, "screen_pairs_p1e-3.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(X = X, tauR = tauR, MG = MG, v = v, s = s, Tm = Tm, intercept = intercept, gcov = gcov, n = n,
             KG = KG, KT = KT, ids = cand$ID), file.path(OUT, "profiles.rds"))
print(head(res, 30), digits = 3)
