# Scale check (lrcpq::check_scale logic, per LD block): observed / expected mean chi-square signal
# under w = 1 for n = total N vs N_eff and M = HapMap3 count vs M_5_50. usage: Rscript check_scale.R <out_dir with stage3.rds>, run in the targeted data dir.
.libPaths(c("/home/user/rlib4", .libPaths())); suppressPackageStartupMessages(library(lrcpq))
OUT <- commandArgs(TRUE)[1]; s3 <- readRDS(file.path(OUT, "stage3.rds"))
Z <- as.matrix(read.delim(gzfile("in/Z.tsv.gz"), row.names = 1, check.names = FALSE))[, s3$keys]
L <- ld_from_windows("in/ld", N_ref = 337000)
l <- numeric(L$ld$m)
for (b in unique(L$ld$block_id)) { ix <- which(L$ld$block_id == b); l[ix] <- rowSums(ld_r2(L$ld$get(ix), L$ld$N)) }
h2 <- diag(s3$gcov); ic <- diag(s3$intercept); n <- s3$n
bl <- L$ld$block_id
for (M in c(1094844, 6805960)) {
  r <- (colMeans(Z^2) - ic) / (n * h2 * mean(l) / M)
  rb <- sapply(unique(bl), function(b) median((colMeans(Z[bl == b, ]^2) - ic) / (n * h2 * mean(l[bl == b]) / M)))
  cat(sprintf("M = %d: median ratio over 444 traits %.2f (traits with h2 > 0.05: %.2f); per block %s\n", M, median(r), median(r[h2 > 0.05]), paste(round(rb, 2), collapse = " ")))
}
nb <- n; nb2 <- read.delim("/mnt/project-files/papers/LRCQ/results/real/trait-list.tsv"); nb2$key <- sub("\\.tsv\\.bgz$", "", basename(nb2$aws_path))
bin <- nb2$binary[match(s3$keys, nb2$key)] %in% c(TRUE, "True")
r <- (colMeans(Z^2) - ic) / (n * h2 * mean(l) / 1094844)
cat(sprintf("M = 1094844, total N: binary median %.2f, quantitative median %.2f\n", median(r[bin]), median(r[!bin])))
neff <- nb2$N_eff[match(s3$keys, nb2$key)]
r2 <- (colMeans(Z^2) - ic) / (neff * h2 * mean(l) / 1094844)
cat(sprintf("M = 1094844, N_eff: binary median %.2f\n", median(r2[bin])))
cat("SNPs", nrow(Z), "blocks", length(unique(bl)), "mean LD score", round(mean(l), 1), "\n")
