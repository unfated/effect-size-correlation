# S9: gene-level aggregation. Gene G = 3 causal SNPs in region A, gene H = 3 in region B; all
# 9 cross pairs share correlation rho (coherent sign); within-gene correlation 0.6. Power of single SNP-pair tests (conditional,
# best and Bonferroni-over-9) vs the gene burden test (lrcp_gene with indicator weights), and
# type-I error of the burden test at rho = 0. S1 background (10 clusters, rg 0.6, full overlap).
.libPaths(c("/home/user/rlib2", .libPaths()))
suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
args <- commandArgs(TRUE); out <- args[1]; reps <- as.integer(args[2]); set.seed(as.integer(args[3]))
one <- function(q, rho) {
  blocks <- make_ld_ar1(rep(100, 4), rho = 0.6); m <- 400; reg1 <- 1:200; reg2 <- 201:400
  w <- make_w(m, prop_nonzero = 0.1, dist = "lnorm")
  nz1 <- intersect(which(w > 0), reg1); nz2 <- intersect(which(w > 0), reg2)
  gA <- sample(nz1, 3); gB <- sample(nz2, 3); pr <- as.matrix(expand.grid(gA, gB))
  Rb <- diag(m); Rb[gA, gA] <- 0.6; Rb[gB, gB] <- 0.6; diag(Rb) <- 1   # within-gene correlation 0.6
  Rb[gA, gB] <- rho; Rb[gB, gA] <- rho                                    # coherent cross-gene rho
  gc <- make_gcov(q, type = "cluster", rg = 0.6, n_clusters = 10); C <- make_intercept(q, overlap = 1, rp = 0.3 * (gc$Rg > 0))
  sim <- simulate_lrcpq(blocks, n = rep(3e5, q), w = w, Rb = Rb, gcov = gc$gcov, intercept = C, M = 2e5)
  R1 <- as.matrix(bdiag(blocks[1:2])); R2 <- as.matrix(bdiag(blocks[3:4])); Z1 <- sim$Z[reg1, ]; Z2 <- sim$Z[reg2, ]
  S1 <- which(w[reg1] > 0); S2 <- which(w[reg2] > 0)
  vA <- as.numeric(S1 %in% gA); vB <- as.numeric(S2 %in% (gB - 200))
  tag <- suppressWarnings(lrcp_gene(Z1[S1, ], Z2[S2, ], R1[S1, S1], R2[S2, S2], n = sim$n, gcov = gc$gcov, M = 2e5, intercept = C, n_sim = 0))
  gen <- suppressWarnings(lrcp_gene(Z1[S1, ], Z2[S2, ], R1[S1, S1], R2[S2, S2], VA = cbind(vA), VB = cbind(vB),
                                    n = sim$n, gcov = gc$gcov, M = 2e5, intercept = C, n_sim = 0))
  pp <- tag$p[S1 %in% gA, S2 %in% (gB - 200)]
  # oracle genetic parts including the within-gene correlation (local rho != 0)
  Sw <- sqrt(w) * t(sqrt(w) * Rb)
  GA <- (R1 %*% Sw[reg1, reg1] %*% R1)[S1, S1]; GB <- (R2 %*% Sw[reg2, reg2] %*% R2)[S2, S2]
  geo <- suppressWarnings(lrcp_gene(Z1[S1, ], Z2[S2, ], R1[S1, S1], R2[S2, S2], VA = cbind(vA), VB = cbind(vB),
                                    n = sim$n, gcov = gc$gcov, M = 2e5, intercept = C, GA = GA, GB = GB, n_sim = 0))
  data.frame(q = q, rho = rho, gene_oracleG = as.numeric(geo$p[1, 1] < 0.05), snp_any = mean(pp < 0.05), snp_bonf9 = as.numeric(min(pp) < 0.05 / 9),
             gene = as.numeric(gen$p[1, 1] < 0.05), gene_z = gen$z[1, 1])
}
grid <- expand.grid(q = c(100, 300), rho = c(0, 0.2, 0.3, 0.5))
R <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid) * reps), function(k) { g <- grid[(k - 1) %% nrow(grid) + 1, ]
  tryCatch(one(g$q, g$rho), error = function(e) { message(conditionMessage(e)); NULL }) }, mc.cores = 3))
saveRDS(R, out)
a <- aggregate(cbind(snp_any, snp_bonf9, gene, gene_oracleG) ~ q + rho, R, mean); a$n <- aggregate(gene ~ q + rho, R, length)$gene
print(a, digits = 3)
