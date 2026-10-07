# Summarise S1 (distal pairs): bias, SD, SE calibration, coverage, null type-I error.
args <- commandArgs(TRUE); x <- readRDS(args[1]); out <- args[2]
r <- x$res; nl <- x$nulls
key <- c("q", "rho", "traits", "ld_r")
r$atr <- abs(r$truth); s <- sign(ifelse(r$truth == 0, 1, r$truth))
for (m in c("gls", "ols", "naive", "oracle")) r[[paste0(m, "_s")]] <- r[[m]] * s   # sign-aligned estimates
agg <- function(d) {
  t <- d$atr[1]
  data.frame(n = nrow(d),
    gls_bias = mean(d$gls_s) - t, gls_sd = sd(d$gls_s), gls_se = mean(d$gls_se), gls_pse = mean(d$gls_pse),
    gls_cover = mean(abs(d$gls_s - t) <= 1.96 * d$gls_pse),
    ols_bias = mean(d$ols_s) - t, ols_sd = sd(d$ols_s), ols_se = mean(d$ols_se),
    naive_bias = mean(d$naive_s) - t, naive_sd = sd(d$naive_s),
    oracle_bias = mean(d$oracle_s) - t, oracle_sd = sd(d$oracle_s))
}
sp <- split(r, r[, key], drop = TRUE)
S <- do.call(rbind, lapply(sp, function(d) cbind(d[1, key], agg(d))))
S <- S[order(S$traits, S$ld_r, S$rho, S$q), ]
T1 <- aggregate(cbind(gls = abs(z_gls) > 1.96, ols = abs(z_ols) > 1.96, naive = abs(z_naive) > 1.96,
                      gls_1e3 = abs(z_gls) > qnorm(1 - 5e-4)) ~ q + traits + ld_r, data = nl, FUN = mean)
T1$n_null <- aggregate(z_gls ~ q + traits + ld_r, data = nl, FUN = length)$z_gls
T1 <- T1[order(T1$traits, T1$ld_r, T1$q), ]
write.table(format(S, digits = 3), file.path(out, "s1_estimation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(format(T1, digits = 3), file.path(out, "s1_typeI.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(S[S$ld_r == 0.6, c(key, "gls_bias", "gls_sd", "gls_se", "gls_pse", "gls_cover", "ols_bias", "ols_sd", "naive_bias", "oracle_sd")], digits = 2, row.names = FALSE)
print(T1, digits = 2, row.names = FALSE)
