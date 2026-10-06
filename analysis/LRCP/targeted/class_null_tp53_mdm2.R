# Class null (theory S8.5) for TP53 rs78378222 x MDM2 rs3730556: the class of each variant is
# the set of genome-wide candidate tags (step 13 profiles) whose deconvolved trait profile
# correlates with the variant's at |r| >= 0.5; class covariance from those tags.
.libPaths(c("/home/user/rlib3", .libPaths())); suppressPackageStartupMessages(library(lrcpq))
x <- readRDS(Sys.getenv("INPUTS", "/mnt/project-files/papers/LRCP/results/real/targeted_ldsc/tp53_mdm2_lrcp_gene_inputs.rds"))
P <- readRDS(Sys.getenv("PROFILES", "/home/user/data/gw/out/profiles.rds"))
tr <- read.delim("/mnt/project-files/papers/LRCQ/results/real/panukb.traits.tsv")
tl <- read.delim("/mnt/project-files/papers/LRCQ/results/real/trait-list.tsv")
tl$key <- sub("\\.tsv\\.bgz$", "", basename(tl$aws_path))
ids <- tl$trait_id[match(colnames(x$ZA), tl$key)]; o <- match(ids, tr$trait_id); stopifnot(!anyNA(o))
X <- P$X[, o]; Ic <- P$intercept[o, o]
prof <- function(Z, R, V) as.vector(crossprod(solve(R, V), Z))
for (min_r in c(0.5, 0.3)) {
  cc <- lapply(list(A = prof(x$ZA, x$RA, x$VA), B = prof(x$ZB, x$RB, x$VB)), function(p) {
    use <- which(abs(cor(t(X), p)) >= min_r)
    c(class_covariance(X, P$tauR, P$v, Ic, use), list(n = length(use))) })
  g <- suppressWarnings(lrcp_gene(x$ZA, x$ZB, x$RA, x$RB, x$VA, x$VB, n = x$n, gcov = x$gcov, M = x$M,
         intercept = x$intercept, VclassA = list(cc$A$V), VclassB = list(cc$B$V), n_boot = 0, n_sim = 0))
  cat(sprintf("min|r| %.1f: class sizes MDM2 %d, TP53 %d; q_eff %.1f, %.1f; class-null z %.2f, P %.3f\n",
      min_r, cc$A$n, cc$B$n, cc$A$q_eff, cc$B$q_eff, g$z[1, 1], g$p[1, 1]))
}
