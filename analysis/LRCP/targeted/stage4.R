# Targeted LRCP stage 4: distal GLS between all pairs of the six blocks.
# Candidates per block: tags (r2 < 0.5 clumps) with enrichment z >= z_min, plus the
# named motivating/control SNP, which replaces its own tag. Enrichment = max(w_raw, 0)
# on candidates (denominators re-estimated by R-projection, theory S4.7).
suppressPackageStartupMessages(library(lrcpq))
args <- commandArgs(TRUE); IN <- args[1]; OUT <- args[2]; z_min <- as.numeric(args[3])
s3 <- readRDS(file.path(OUT, "stage3.rds")); w <- s3$fit$w
Z <- as.matrix(read.delim(gzfile(file.path(IN, "Z.tsv.gz")), row.names = 1, check.names = FALSE))
if (!is.null(s3$keys)) Z <- Z[, s3$keys]
L <- ld_from_windows(file.path(IN, "ld"), N_ref = 337000)
TARGET <- c(TP53 = "17:7571752:T:G", MDM2 = "12:69216521:T:G", PCSK9 = "1:55505647:G:T",
            LDLR = "19:11202306:G:T", HMGCR = "5:74656539:T:C", APOB = "2:21263900:G:A")
blk <- split(seq_len(nrow(w)), w$chr)
bname <- c(`1` = "PCSK9", `2` = "APOB", `5` = "HMGCR", `12` = "MDM2", `17` = "TP53", `19` = "LDLR")
cand <- list(); wv <- list()
for (b in names(blk)) {
  ix <- blk[[b]]; nm <- bname[[b]]
  t <- match(TARGET[[nm]], w$ID); ttag <- w$tag[t]
  S <- ix[!is.na(w$w_raw[ix]) & w$z[ix] >= z_min & ix != ttag]
  S <- sort(c(S, t))
  wb <- numeric(length(ix)); names(wb) <- ix
  wb[as.character(S)] <- pmax(w$w_raw[S], 0)
  wb[as.character(t)] <- max(w$w_raw[ttag], 1)
  # w of non-candidate tags enters only the null covariance R W R
  tg <- ix[!is.na(w$w_raw[ix]) & !(ix %in% S)]
  wb[as.character(tg)] <- pmax(w$w_raw[tg], 0)
  cand[[nm]] <- match(S, ix); wv[[nm]] <- wb
  message(nm, ": ", length(S), " candidates")
}
res <- list(); nmz <- names(cand)
for (i in 1:5) for (j in (i + 1):6) {
  a <- nmz[i]; b <- nmz[j]; ia <- blk[[names(bname)[bname == a]]]; ib <- blk[[names(bname)[bname == b]]]
  RA <- L$ld$get(ia); RB <- L$ld$get(ib)
  arg <- list(Z1 = Z[ia, ], Z2 = Z[ib, ], R1 = RA, R2 = RB, w1 = wv[[a]], w2 = wv[[b]], n = s3$n,
              gcov = s3$gcov, M = s3$M, intercept = s3$intercept, S1 = cand[[a]], S2 = cand[[b]],
              method = "gls", denominators = "R")
  f0 <- suppressWarnings(do.call(lrcp_distal, arg))
  fp <- suppressWarnings(do.call(lrcp_distal, c(arg, se_type = "plugin")))
  sa <- ia[cand[[a]]]; sb <- ib[cand[[b]]]
  d <- data.frame(block1 = a, block2 = b, id1 = rep(w$ID[sa], length(sb)), id2 = rep(w$ID[sb], each = length(sa)),
                  w1 = rep(w$w_raw[w$tag[sa]], length(sb)), w2 = rep(w$w_raw[w$tag[sb]], each = length(sa)),
                  rho = as.vector(f0$rho), se_null = as.vector(f0$se), se_plugin = as.vector(fp$se))
  d$target_pair <- d$id1 == TARGET[[a]] & d$id2 == TARGET[[b]]
  res[[length(res) + 1]] <- d
}
R <- do.call(rbind, res)
R$z <- R$rho / R$se_null; R$p <- 2 * pnorm(-abs(R$z))
R$q_BH <- p.adjust(R$p, "BH")
R$lo <- R$rho - 1.96 * R$se_plugin; R$hi <- R$rho + 1.96 * R$se_plugin
write.table(R, file.path(OUT, sprintf("stage4_pairs_z%g.tsv", z_min)), sep = "\t", quote = FALSE, row.names = FALSE)
print(R[R$target_pair, c("block1", "block2", "rho", "se_null", "se_plugin", "z", "p", "q_BH")], digits = 3, row.names = FALSE)
cat("\nAll pairs:", nrow(R), " BH<0.05:", sum(R$q_BH < 0.05), "\n")
print(table(R$block1, R$block2, R$q_BH < 0.05))
print(head(R[order(R$p), c("block1", "block2", "id1", "id2", "rho", "se_null", "z", "q_BH")], 20), digits = 3, row.names = FALSE)
