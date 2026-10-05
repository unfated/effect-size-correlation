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
keys <- tr$key
if (Sys.getenv("PLUGIN") == "ldsc") {
  # final plug-ins: LRCQ genome-wide cross-trait LDSC (444 traits); diagonals from Pan-UKB's own
  # LDSC (h2, intercept), off-diagonals from LRCQ's rg and intercept matrices (as in LRCQ stage 3)
  L2 <- "/mnt/project-files/papers/LRCQ/results/real/"
  hh <- read.delim(file.path(L2, "panukb_ldsc.h2.tsv"))
  rg <- as.matrix(read.delim(file.path(L2, "panukb_ldsc.rg.tsv"), row.names = 1, check.names = FALSE))
  ICm <- as.matrix(read.delim(file.path(L2, "panukb_ldsc.intercept.tsv"), row.names = 1, check.names = FALSE))
  tl <- read.delim("/mnt/project-files/papers/LRCQ/results/real/trait-list.tsv")
  tl$key <- sub("\\.tsv\\.bgz$", "", basename(tl$aws_path))
  key_of <- setNames(tl$key, tl$trait_id)
  ids <- hh$trait_id[!is.na(key_of[hh$trait_id]) & key_of[hh$trait_id] %in% colnames(Z)]
  keys <- unname(key_of[ids]); Z <- Z[, keys]; tr <- tr[match(keys, tr$key), ]
  rg <- rg[ids, ids]; rg[is.na(rg)] <- 0; diag(rg) <- 1; rg <- pmin(pmax(rg, -1), 1)
  ICm <- ICm[ids, ids]; ICm[is.na(ICm)] <- 0
  h2 <- pmax(hh$h2_panukb_ldsc[match(ids, hh$trait_id)], 1e-3); ic <- pmax(hh$intercept_panukb[match(ids, hh$trait_id)], 1)
  n <- tr$N_eff
  gcov <- psd(rg * sqrt(outer(h2, h2))); diag(gcov) <- h2
  intercept <- psd(ICm); diag(intercept) <- ic
  message("LDSC plug-ins: ", length(ids), " traits")
} else {
n <- tr$N_eff; h2 <- pmax(tr$h2_obs, 1e-3); ic <- pmax(tr$ldsc_intercept, 1)
intercept <- psd(Ch * sqrt(outer(ic, ic))); diag(intercept) <- ic
gcov <- psd(Ch * sqrt(outer(h2, h2))); diag(gcov) <- h2
}
L <- ld_from_windows(file.path(IN, "ld"), N_ref = 337000)
stopifnot(identical(L$ld$snp$ID, rownames(Z)))
t0 <- Sys.time()
fit <- run_lrcq(Z, L$ld, n, M, gcov = gcov, intercept = intercept, tag_r2 = 0.5,
                method = "wls", rectify = "C", windows = L$windows)
message("lrcq ", round(as.numeric(Sys.time() - t0, units = "secs")), " s")
saveRDS(list(fit = fit, n = n, gcov = gcov, intercept = intercept, M = M, keys = keys), file.path(OUT, "stage3.rds"))
w <- fit$w
print(table(block = w$chr, tag = !is.na(w$w_raw)))
print(summary(w$w_raw))
print(head(w[order(-w$z), c("ID", "w_raw", "w", "se", "z")], 25))
print(w[w$ID %in% c("17:7571752:T:G", "12:69216521:T:G", "1:55505647:G:T", "2:21263900:G:A"), ])
