# Locus-specific effective number of traits: trait covariance of the enriched tags'
# Z profiles in each block (minus the null intercept), q_loc = tr(S)^2 / tr(S^2).
args <- commandArgs(TRUE); IN <- args[1]; OUT <- args[2]
s3 <- readRDS(file.path(OUT, "stage3.rds")); w <- s3$fit$w
Z <- as.matrix(read.delim(gzfile(file.path(IN, "Z.tsv.gz")), row.names = 1, check.names = FALSE))
if (!is.null(s3$keys)) Z <- Z[, s3$keys]
nm <- c(`1` = "PCSK9", `2` = "APOB", `5` = "HMGCR", `12` = "MDM2", `17` = "TP53", `19` = "LDLR")
qe <- function(S) sum(diag(S))^2 / sum(S^2)
res <- do.call(rbind, lapply(names(nm), function(c) {
  S <- which(w$chr == as.integer(c) & !is.na(w$w_raw) & w$z >= 3)
  X <- Z[S, , drop = FALSE]
  Sig <- crossprod(X) / length(S) - s3$intercept
  e <- eigen(Sig, symmetric = TRUE); Sig <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
  top <- order(-colMeans(X^2))[1:5]
  data.frame(block = nm[[c]], n_tags = length(S), q_loc = qe(Sig), top_traits = paste(colnames(Z)[top], collapse = "; "))
}))
res$q_global_G <- qe(s3$gcov); res$q_global_P <- qe(s3$intercept)
write.table(res, file.path(OUT, "local_qeff.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res[, c("block", "n_tags", "q_loc", "q_global_G", "q_global_P")], digits = 3)
