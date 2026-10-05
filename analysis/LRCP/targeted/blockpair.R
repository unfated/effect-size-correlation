# Block-pair LRCP tests for the six targeted blocks: orientation-free Q over all
# candidate-tag pairs and the named-variant burden covariance, with the conditional
# Gaussian null of theory supplement S8.2 (exact under any dependence among GWAS).
if (nzchar(Sys.getenv("RLIB"))) .libPaths(c(Sys.getenv("RLIB"), .libPaths()))
suppressPackageStartupMessages(library(lrcpq))
args <- commandArgs(TRUE); IN <- args[1]; OUT <- args[2]; z_min <- as.numeric(args[3])
s3 <- readRDS(file.path(OUT, "stage3.rds")); w <- s3$fit$w
Z <- as.matrix(read.delim(gzfile(file.path(IN, "Z.tsv.gz")), row.names = 1, check.names = FALSE))
L <- ld_from_windows(file.path(IN, "ld"), N_ref = 337000)
TARGET <- c(PCSK9 = "1:55505647:G:T", APOB = "2:21263900:G:A", HMGCR = "5:74656539:T:C",
            MDM2 = "12:69216521:T:G", TP53 = "17:7571752:T:G", LDLR = "19:11202306:G:T")
chr <- c(PCSK9 = 1, APOB = 2, HMGCR = 5, MDM2 = 12, TP53 = 17, LDLR = 19)
cand <- lapply(names(TARGET), function(nm) {
  ix <- which(w$chr == chr[[nm]]); t <- match(TARGET[[nm]], w$ID)
  S <- ix[!is.na(w$w_raw[ix]) & w$z[ix] >= z_min & ix != w$tag[t]]
  sort(c(S, t))
}); names(cand) <- names(TARGET)
# trait-cluster jackknife (20 clusters on |C_hat|, delete one cluster) for the named-variant burden
Ch <- as.matrix(read.delim(gzfile(file.path(IN, "C_hat.tsv.gz")), row.names = 1, check.names = FALSE))
cl <- cutree(hclust(as.dist(1 - abs(Ch)), "average"), k = 20)
s_tr <- s3$n * diag(s3$gcov) / s3$M
jack <- function(x, y) {
  est <- function(k) sum((s_tr * x * y)[k]) / sum(s_tr[k]^2)
  jk <- sapply(1:20, function(g) est(cl != g))
  c(C = est(rep(TRUE, length(x))), se = sqrt(19 / 20 * sum((jk - mean(jk))^2)))
}
set.seed(1); out <- list()
for (i in 1:5) for (j in (i + 1):6) {
  a <- names(TARGET)[i]; b <- names(TARGET)[j]; A <- cand[[a]]; B <- cand[[b]]
  RA <- L$ld$get(A); RB <- L$ld$get(B)
  va <- as.numeric(A == match(TARGET[[a]], w$ID)); vb <- as.numeric(B == match(TARGET[[b]], w$ID))
  g <- lrcp_gene(Z[A, ], Z[B, ], RA, RB, VA = cbind(target = va), VB = cbind(target = vb),
                 n = s3$n, gcov = s3$gcov, M = s3$M, intercept = s3$intercept, n_sim = 0)
  t <- lrcp_gene(Z[A, ], Z[B, ], RA, RB, n = s3$n, gcov = s3$gcov, M = s3$M,
                 intercept = s3$intercept, n_sim = 20000)
  jx <- jack(as.vector(crossprod(solve(RA, va), Z[A, ])), as.vector(crossprod(solve(RB, vb), Z[B, ])))
  out[[length(out) + 1]] <- data.frame(block1 = a, block2 = b, n1 = length(A), n2 = length(B),
    z_target_jk = jx[["C"]] / jx[["se"]], r_naive = cor(Z[match(TARGET[[a]], w$ID), ], Z[match(TARGET[[b]], w$ID), ]),
    C_target = g$C[1, 1], z_target = g$z[1, 1], p_target = g$p[1, 1],
    z_gls_gene = if (!is.null(g$z_gls)) g$z_gls[1, 1] else NA, R_AB = if (!is.null(g$R_AB)) as.vector(g$R_AB)[1] else NA, e_B = if (!is.null(g$e_B)) as.vector(g$e_B)[1] else NA,
    Q = t$Q, p_Q = t$p_Q, min_p_tag = min(t$p), p_maxT = min(t$p_maxT))
}
R <- do.call(rbind, out)
R$q_BH_Q <- p.adjust(R$p_Q, "BH"); R$q_BH_maxT <- p.adjust(R$p_maxT, "BH")
write.table(R, file.path(OUT, sprintf("blockpair_z%g.tsv", z_min)), sep = "\t", quote = FALSE, row.names = FALSE)
print(R, digits = 3, row.names = FALSE)
