#!/usr/bin/env Rscript
# Summarise LRCQ simulation replicates into tables for the paper.
# Usage: summarise_sim.R <rds_dir> <out_dir>
suppressPackageStartupMessages(library(lrcpq))
args <- commandArgs(TRUE)
indir <- args[1]; outdir <- args[2]
auc <- function(score, truth) {
  truth <- as.logical(truth); if (!any(truth) || all(truth)) return(NA_real_)
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
files <- list.files(indir, "\\.rds$", full.names = TRUE)
per_rep <- list(); annot <- list()
for (f in files) {
  x <- readRDS(f); t <- x$tab
  omni <- t$t_true >= quantile(t$t_true, 0.99)       # top 1% of clumps by true enrichment
  top10 <- t$t_true >= quantile(t$t_true, 0.90)
  nn <- t$t_true >= 0.5                               # clearly non-null clumps
  null <- t$t_true < 0.05
  est <- c("tag", "tag_ols", "tag_eq", "ridge_clump", "mean_chi2", "ldnorm_chi2", "n_sig", "omnibus")
  row <- data.frame(scenario = x$scen, rep = x$rep, m = x$m, n_tags = nrow(t))
  for (e in est) {
    row[[paste0("spearman_", e)]] <- cor(t[[e]], t$t_true, method = "spearman")
    row[[paste0("auc_", e)]] <- auc(t[[e]], omni)
    row[[paste0("auc10_", e)]] <- auc(t[[e]], top10)
    row[[paste0("spearman_nonnull_", e)]] <- if (sum(nn) > 5) cor(t[[e]][nn], t$t_true[nn], method = "spearman") else NA
    row[[paste0("pearson_", e)]] <- cor(t[[e]], t$t_true)
  }
  # SE calibration of the tag-set estimators (standardised residual SD; type I at null tags)
  for (e in c("tag", "tag_ols", "tag_eq")) {
    zr <- (t[[e]] - t$t_true) / t[[paste0(e, "_se")]]
    row[[paste0("zres_sd_", e)]] <- sd(zr)
    row[[paste0("typeI_null_", e)]] <- mean((t[[e]] / t[[paste0(e, "_se")]])[null] > qnorm(0.95))
    row[[paste0("cover95_", e)]] <- mean(abs(zr) < qnorm(0.975))
  }
  # knockoff+ (Method D, theory S3.8) at FDR 10%: realised FDP and power for declaring w > 0
  wv <- t$tag; ts <- sort(unique(abs(wv[wv != 0])))
  thr <- Inf
  for (th in ts) if ((1 + sum(wv <= -th)) / max(sum(wv >= th), 1) <= 0.10) { thr <- th; break }
  disc <- wv >= thr
  row$methodD_n_disc <- sum(disc)
  row$methodD_fdp <- if (sum(disc)) mean(t$t_true[disc] < 0.05) else 0
  row$methodD_power <- mean(disc[t$t_true >= 1])
  for (e in c("tag", "tag_ols", "tag_eq", "ridge_clump")) {
    row[[paste0("bias_", e)]] <- mean(t[[e]] - t$t_true)
    row[[paste0("rmse_", e)]] <- sqrt(mean((t[[e]] - t$t_true)^2))
    row[[paste0("slope_", e)]] <- coef(lm(t[[e]] ~ t$t_true))[2]
  }
  row$pearson_tag_eq_vs_eqtarget <- cor(t$tag_eq, t$t_true_eq)
  row$slope_tag_vs_eqtarget <- coef(lm(t$tag ~ t$t_true_eq))[2]
  row$slope_tag_eq_vs_eqtarget <- coef(lm(t$tag_eq ~ t$t_true_eq))[2]
  # rectification (applied to the tag estimates)
  for (r in c("A", "B", "C", "zero")) {
    wr <- rectify_w(t$tag, r)
    row[[paste0("rmse_rect", r)]] <- sqrt(mean((wr - t$t_true)^2))
    row[[paste0("sum_ratio_rect", r)]] <- sum(wr) / sum(t$t_true)
    row[[paste0("null_zeroed_rect", r)]] <- mean(wr[null] == 0)
    row[[paste0("nonnull_kept_rect", r)]] <- mean(wr[t$t_true >= 1] > 0)
  }
  # ridge SE calibration at null tags
  zz <- t$ridge / t$ridge_se
  row$ridge_z_sd_null <- sd(zz[null]); row$ridge_typeI_null <- mean(zz[null] > qnorm(0.95))
  per_rep[[f]] <- row
  a <- x$annot; a$scenario <- x$scen; a$rep <- x$rep; annot[[f]] <- a
}
pr <- do.call(rbind, per_rep); an <- do.call(rbind, annot)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.table(pr, file.path(outdir, "sim_per_replicate.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
num <- sapply(pr, is.numeric)
agg <- aggregate(pr[, num & !names(pr) %in% c("rep", "m")], list(scenario = pr$scenario), mean, na.rm = TRUE)
nrep <- aggregate(pr$rep, list(scenario = pr$scenario), length); agg$n_rep <- nrep$x[match(agg$scenario, nrep$scenario)]
write.table(agg, file.path(outdir, "sim_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
an_agg <- aggregate(cbind(true, clump_target, lrcq_clump, lrcq_ridge, sldsc_pooled) ~ scenario + annot, an, mean)
an_sd <- aggregate(cbind(lrcq_clump, lrcq_ridge, sldsc_pooled) ~ scenario + annot, an, sd)
names(an_sd)[3:5] <- c("sd_lrcq_clump", "sd_lrcq_ridge", "sd_sldsc_pooled")
write.table(merge(an_agg, an_sd), file.path(outdir, "sim_annotation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(t(agg[, c("scenario", grep("^(pearson|slope|zres_sd|typeI_null|cover95)_tag$|^auc10_|methodD", names(agg), value = TRUE))]))
