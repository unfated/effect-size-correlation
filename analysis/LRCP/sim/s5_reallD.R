# S5: real LD (UKB in-sample, TP53 and MDM2 blocks) and LD-reference mismatch (1000G EUR, n = 503).
# Causal SNPs: 10% of r2 < 0.5 tags (UKB LD), log-normal w; three cross-block pairs with (rho, rho, -rho).
# Z simulated with UKB LD; estimation with UKB LD vs 1000G LD; GLS (estimate, model SE) and conditional test.
.libPaths(c("/home/user/rlib2", .libPaths()))
suppressPackageStartupMessages({library(lrcpq); library(Matrix)})
args <- commandArgs(TRUE); D <- args[1]; out <- args[2]; reps <- as.integer(args[3]); set.seed(as.integer(args[4]))
rd <- function(nm, k) { s <- read.delim(file.path(D, paste0(nm, ".snps.tsv"))); m <- nrow(s)
  matrix(readBin(file.path(D, paste0(nm, ".", k, ".f32")), "numeric", m * m, size = 4), m, m) }
nmA <- "chr17_7317398_8306425"; nmB <- "chr12_67909729_69826542"
RA <- rd(nmA, "Rukb"); RB <- rd(nmB, "Rukb"); RA1 <- rd(nmA, "R1kg"); RB1 <- rd(nmB, "R1kg")
fix <- function(R) { e <- eigen((R + t(R)) / 2, symmetric = TRUE); e$vectors %*% (pmax(e$values, 1e-4) * t(e$vectors)) }
RA <- fix(RA); RB <- fix(RB)
prune <- function(R, r2 = 0.5) { keep <- integer(0); for (s in sample(nrow(R))) if (!length(keep) || all(R[s, keep]^2 < r2)) keep <- c(keep, s); sort(keep) }
mA <- nrow(RA); mB <- nrow(RB); m <- mA + mB
one <- function(q, rho) {
  tA <- prune(RA); tB <- prune(RB)
  cA <- sample(tA, ceiling(0.1 * length(tA))); cB <- sample(tB, ceiling(0.1 * length(tB)))
  w <- numeric(m); w[c(cA, mA + cB)] <- rlnorm(length(cA) + length(cB)); w <- w / mean(w)
  pr <- cbind(sample(cA, 3), mA + sample(cB, 3))
  Rb <- make_Rb(m, "pairs", pairs = pr, pair_rho = c(rho, rho, -rho))
  gc <- make_gcov(q, type = "cluster", rg = 0.6, n_clusters = 10); C <- make_intercept(q, overlap = 1, rp = 0.3 * (gc$Rg > 0))
  sim <- simulate_lrcpq(list(RA, RB), n = rep(3e5, q), w = w, Rb = Rb, gcov = gc$gcov, intercept = C, M = 2e5)
  ZA <- sim$Z[1:mA, ]; ZB <- sim$Z[mA + 1:mB, ]
  SA <- sort(cA); SB <- sort(cB); ti <- match(pr[, 1], SA); tj <- match(pr[, 2] - mA, SB)
  res <- list()
  for (ref in c("ukb", "1kg")) {
    A <- if (ref == "ukb") RA else RA1; B <- if (ref == "ukb") RB else RB1
    f <- suppressWarnings(lrcp_distal(ZA, ZB, A, B, w[1:mA], w[mA + 1:mB], sim$n, gc$gcov, 2e5, C, S1 = SA, S2 = SB,
                                      method = "gls", ridge = if (ref == "1kg") 0.01 else 0))
    g <- suppressWarnings(lrcp_gene(ZA[SA, ], ZB[SB, ], A[SA, SA], B[SB, SB], n = sim$n, gcov = gc$gcov, M = 2e5, intercept = C, n_sim = 0))
    nul <- matrix(TRUE, length(SA), length(SB)); nul[cbind(ti, tj)] <- FALSE
    res[[ref]] <- data.frame(ref = ref, q = q, rho = rho, est = mean(f$rho[cbind(ti, tj)] * sign(c(1, 1, -1))),
      se = mean(f$se[cbind(ti, tj)]), rej_gls_null = mean(abs(f$rho / f$se)[nul] > 1.96, na.rm = TRUE),
      rej_cond_null = mean(abs(g$z)[nul] > 1.96), pow_cond = mean(abs(g$z[cbind(ti, tj)]) > 1.96),
      pow_gls = mean(abs((f$rho / f$se)[cbind(ti, tj)]) > 1.96), nA = length(SA), nB = length(SB))
  }
  do.call(rbind, res)
}
grid <- expand.grid(q = c(100, 300), rho = c(0, 0.5))
R <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid) * reps), function(k) { g <- grid[(k - 1) %% nrow(grid) + 1, ]
  tryCatch(one(g$q, g$rho), error = function(e) { message(conditionMessage(e)); NULL }) }, mc.cores = 3))
saveRDS(R, out)
a <- aggregate(cbind(est, se, rej_gls_null, rej_cond_null, pow_gls, pow_cond, nA, nB) ~ ref + q + rho, R, mean)
a$sd_est <- aggregate(est ~ ref + q + rho, R, sd)$est; a$n <- aggregate(est ~ ref + q + rho, R, length)$est
print(a, digits = 3)
