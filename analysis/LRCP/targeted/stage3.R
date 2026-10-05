# Targeted LRCP, stage 1-3 on six Berisa blocks (work plan step 12).
# Preliminary plug-ins until the genome-wide LDSC matrices from the LRCQ thread arrive:
#   n = Pan-UKB N_eff; h2 = h2_obs (floored at 1e-3); intercept diag = Pan-UKB ldsc_intercept,
#   off-diagonals = C_hat (Z correlation at near-null SNPs) scaled to those diagonals;
#   gcov = sqrt(h2_a h2_b) * C_hat (Cheverud approximation r_g ~ r_p), projected to PSD.
suppressPackageStartupMessages(library(lrcpq))
args <- commandArgs(TRUE); IN <- args[1]; OUT <- args[2]
M <- 1094844
Z <- as.matrix(read.delim(gzfile(file.path(IN, "Z.tsv.gz")), row.names = 1, check.names = FALSE))
tr <- read.delim(file.path(IN, "traits.tsv"))
stopifnot(identical(colnames(Z), tr$key))
Ch <- as.matrix(read.delim(gzfile(file.path(IN, "C_hat.tsv.gz")), row.names = 1, check.names = FALSE))
psd <- function(A, eps = 1e-6) { e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  e$vectors %*% (pmax(e$values, eps * max(e$values)) * t(e$vectors)) }
n <- tr$N_eff; h2 <- pmax(tr$h2_obs, 1e-3); ic <- pmax(tr$ldsc_intercept, 1)
intercept <- psd(Ch * sqrt(outer(ic, ic))); diag(intercept) <- ic
gcov <- psd(Ch * sqrt(outer(h2, h2))); diag(gcov) <- h2
L <- ld_from_windows(file.path(IN, "ld"), N_ref = 337000)
stopifnot(identical(L$ld$snp$ID, rownames(Z)))
t0 <- Sys.time()
fit <- run_lrcq(Z, L$ld, n, M, gcov = gcov, intercept = intercept, tag_r2 = 0.5,
                method = "wls", rectify = "C", windows = L$windows)
message("lrcq ", round(as.numeric(Sys.time() - t0, units = "secs")), " s")
saveRDS(list(fit = fit, n = n, gcov = gcov, intercept = intercept, M = M), file.path(OUT, "stage3.rds"))
w <- fit$w
print(table(block = w$chr, tag = !is.na(w$w_raw)))
print(summary(w$w_raw))
print(head(w[order(-w$w_raw), c("ID", "RSID", "w_raw", "w", "se", "z")], 25))
print(w[w$ID %in% c("17:7571752:T:G", "12:69216521:T:G", "19:11202306:T:G", "1:55505647:G:T", "5:74656539:T:A", "2:21263900:G:A") |
        w$RSID %in% c("rs78378222", "rs3730556"), ])
