# Confirm genome-wide screen hits (work plan step 13): for the top tag pairs, rerun
# lrcpq::lrcp_gene on the full windows (conditioning on the other window tags) with the
# conditional parametric bootstrap (both orientations) and with the class null of theory
# S8.5 (profile classes over all candidate tags, class excluding the pair's windows).
# usage: confirm.R <candidate_tags.tsv> <tagZ.tsv.gz> <ld_dir> <out_dir> <n_top> <n_boot>
.libPaths(c(Sys.getenv("RLIB", "/home/user/rlib2"), .libPaths()))
suppressPackageStartupMessages(library(lrcpq))
args <- commandArgs(TRUE); CAND <- args[1]; ZF <- args[2]; LDD <- args[3]; OUT <- args[4]
n_top <- as.integer(args[5]); n_boot <- as.integer(args[6])
cand <- read.delim(CAND)
Z <- as.matrix(read.delim(gzfile(ZF), row.names = 1, check.names = FALSE)); Z[is.na(Z)] <- 0
P <- readRDS(file.path(OUT, "profiles.rds")); stopifnot(identical(P$ids, cand$ID))
M <- 1094844
sc <- read.delim(gzfile(file.path(OUT, "screen_pairs_p1e-3.tsv.gz")))
sel <- sc[sc$q_BH < 0.05 | seq_len(nrow(sc)) <= n_top, ]
cls <- profile_classes(P$X, 0.5)
message("profile classes: ", length(unique(cls)), "; largest ", paste(head(sort(table(cls), TRUE), 5), collapse = " "))
getR <- function(win, ids) { sn <- read.delim(file.path(LDD, sprintf("w%s.snps.tsv", win)))
  ix <- match(ids, sn$ID); stopifnot(!anyNA(ix))   # LD files may hold extra tags (UKB windows overlap)
  matrix(readBin(file.path(LDD, sprintf("w%s.R.f32", win)), "numeric", nrow(sn)^2, size = 4), nrow(sn))[ix, ix, drop = FALSE] }
genes <- read.delim(gzfile("/mnt/project-files/papers/LRCQ/results/real/genes_all.genes.tsv.gz"))
near <- function(id) { ch <- as.integer(sub(":.*", "", id)); bp <- as.integer(strsplit(id, ":")[[1]][2])
  g <- genes[genes$chr == ch, ]; d <- pmax(g$start - bp, bp - g$end, 0); g$gene[which.min(d)] }
set.seed(13); out <- list()
for (r in seq_len(nrow(sel))) {
  a <- match(sel$ID1[r], cand$ID); b <- match(sel$ID2[r], cand$ID)
  A <- which(cand$window == cand$window[a]); B <- which(cand$window == cand$window[b])
  RA <- getR(cand$window[a], cand$ID[A]); RB <- getR(cand$window[b], cand$ID[B])
  va <- cbind(as.numeric(A == a)); vb <- cbind(as.numeric(B == b))
  g <- lrcp_gene(Z[A, , drop = FALSE], Z[B, , drop = FALSE], RA, RB, va, vb, n = P$n, gcov = P$gcov, M = M,
                 intercept = P$intercept, n_boot = n_boot, n_sim = 0)
  vc <- function(t) { use <- which(cls == cls[t] & cand$window != cand$window[a] & cand$window != cand$window[b])
    cc <- if (length(use) >= 5) class_covariance(P$X, P$tauR, P$v, P$intercept, use) else NULL
    if (is.null(cc) || is.null(cc$V)) list(V = P$Tm, q_eff = sum(diag(P$Tm))^2 / sum(P$Tm^2), n_loci = 0L) else cc }
  ca <- vc(a); cb <- vc(b)
  gc <- suppressWarnings(lrcp_gene(Z[A, , drop = FALSE], Z[B, , drop = FALSE], RA, RB, va, vb, n = P$n, gcov = P$gcov,
                 M = M, intercept = P$intercept, VclassA = list(ca$V), VclassB = list(cb$V), n_boot = 0, n_sim = 0))
  out[[r]] <- data.frame(sel[r, c("ID1", "RSID1", "ID2", "RSID2", "rho", "r_naive", "z", "p", "q_BH")],
    gene1 = near(sel$ID1[r]), gene2 = near(sel$ID2[r]),
    z_gene = g$z[1, 1], p_boot = g$p[1, 1], p_norm = g$p_norm_diagnostic[1, 1],
    boot_sd_A = if (is.null(g$boot_sd_A)) NA else g$boot_sd_A[1], boot_sd_B = if (is.null(g$boot_sd_B)) NA else g$boot_sd_B[1], z_gls = g$z_gls[1, 1], R_AB = g$R_AB[1, 1],
    class1 = cls[a], class2 = cls[b], q_eff_class1 = ca$q_eff, q_eff_class2 = cb$q_eff,
    n_class1 = ca$n_loci, n_class2 = cb$n_loci, z_class = gc$z[1, 1], p_class = gc$p[1, 1])
  if (r %% 10 == 0) message(r, " / ", nrow(sel))
}
R <- do.call(rbind, out)
# reported P is the bootstrap P (larger of the two orientations, theory S8.2a); boot_sd_A/B
# are diagnostics only (the self-normalised z* is light-tailed, software diagnostics note)
R$p_report <- R$p_boot
R$testable <- pmin(R$q_eff_class1, R$q_eff_class2) >= 5
write.table(R, file.path(OUT, paste0("confirm_top", Sys.getenv("SUFFIX"), ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(cls, file.path(OUT, "profile_classes.rds"))
print(R[, c("RSID1", "gene1", "RSID2", "gene2", "rho", "r_naive", "z", "p_boot", "boot_sd_A", "boot_sd_B", "p_report", "z_class", "p_class", "q_eff_class1", "q_eff_class2")], digits = 3)
