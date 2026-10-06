# Sensitivity of the TP53 x MDM2 named-variant test to dropping traits (e.g. SHBG, whose
# cis-signal dominates LRCQ's TP53 gene-level enrichment). Same candidates as blockpair.R.
# usage: drop_traits.R <in_dir> <out_dir> <z_min> <regex of trait descriptions to drop>
.libPaths(c("/home/user/rlib2", .libPaths()))
suppressPackageStartupMessages(library(lrcpq))
args <- commandArgs(TRUE); IN <- args[1]; OUT <- args[2]; z_min <- as.numeric(args[3]); rx <- args[4]
s3 <- readRDS(file.path(OUT, "stage3.rds")); w <- s3$fit$w
Z <- as.matrix(read.delim(gzfile(file.path(IN, "Z.tsv.gz")), row.names = 1, check.names = FALSE))[, s3$keys]
tr <- read.delim(file.path(IN, "traits.tsv")); desc <- tr$description[match(s3$keys, tr$key)]
L <- ld_from_windows(file.path(IN, "ld"), N_ref = 337000)
TG <- c(MDM2 = "12:69216521:T:G", TP53 = "17:7571752:T:G"); chr <- c(MDM2 = 12, TP53 = 17)
cand <- lapply(names(TG), function(nm) { ix <- which(w$chr == chr[[nm]]); t <- match(TG[[nm]], w$ID)
  sort(c(ix[!is.na(w$w_raw[ix]) & w$z[ix] >= z_min & ix != w$tag[t]], t)) }); names(cand) <- names(TG)
A <- cand$TP53; B <- cand$MDM2
va <- cbind(as.numeric(A == match(TG[["TP53"]], w$ID))); vb <- cbind(as.numeric(B == match(TG[["MDM2"]], w$ID)))
# inputs for class_null_tp53_mdm2.R (A = MDM2 side, B = TP53 side)
saveRDS(list(ZA = Z[B, ], ZB = Z[A, ], RA = L$ld$get(B), RB = L$ld$get(A), VA = vb, VB = va, n = s3$n, gcov = s3$gcov,
             M = s3$M, intercept = s3$intercept, note = "A = MDM2 block, B = TP53 block; from drop_traits.R"),
        file.path(OUT, "tp53_mdm2_lrcp_gene_inputs.rds"))
run <- function(keep, lab, swap = FALSE) { set.seed(7)
  g <- if (swap) lrcp_gene(Z[B, keep], Z[A, keep], L$ld$get(B), L$ld$get(A), vb, va, n = s3$n[keep],
    gcov = s3$gcov[keep, keep], M = s3$M, intercept = s3$intercept[keep, keep], n_boot = 2000, n_sim = 0) else
    lrcp_gene(Z[A, keep], Z[B, keep], L$ld$get(A), L$ld$get(B), va, vb, n = s3$n[keep],
    gcov = s3$gcov[keep, keep], M = s3$M, intercept = s3$intercept[keep, keep], n_boot = 2000, n_sim = 0)
  data.frame(set = lab, A_side = if (swap) "MDM2" else "TP53", n_traits = sum(keep), C = g$C[1, 1], z = g$z[1, 1], p_norm = g$p_norm_diagnostic[1, 1],
             p_boot = g$p[1, 1], R_AB = g$R_AB[1, 1]) }
drop <- grepl(rx, desc, ignore.case = TRUE)
message("dropping: ", paste(desc[drop], collapse = "; "))
R <- rbind(run(rep(TRUE, length(desc)), "all"), run(!drop, paste("drop", rx)),
           run(rep(TRUE, length(desc)), "all", TRUE), run(!drop, paste("drop", rx), TRUE))
print(R, digits = 3)
write.table(R, file.path(OUT, "tp53_mdm2_drop_traits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
