# Gene-level burden tests within genome-wide modules (work plan step 13; theory supplement S8).
# Genes are the nearest genes of module member tags. A gene's tags are all enriched candidate
# tags in the member's LD window with that nearest gene, oriented to the member tag by the sign
# of their deconvolved-profile correlation (S8.1 allele coding). For every distal window pair
# inside a module, lrcpq::lrcp_gene gives the burden C_GH for each gene pair, the conditional
# test given the other window (local-moment G, conditional parametric bootstrap, both
# orientations, larger P; S8.2, S8.2a) and the max-T FWER over gene pairs of the window pair.
# usage: burden.R <candidate_tags.tsv> <tagZ.tsv.gz> <ld_dir> <out_dir> <n_boot> [p_edge]
.libPaths(c(Sys.getenv("RLIB", "/home/user/rlib3"), .libPaths()))
suppressPackageStartupMessages({library(lrcpq); library(parallel)})
args <- commandArgs(TRUE); CAND <- args[1]; ZF <- args[2]; LDD <- args[3]; OUT <- args[4]
n_boot <- as.integer(args[5]); PE <- if (length(args) > 5) as.numeric(args[6]) else 1e-5
M <- 1094844; MIN_DIST <- 5e6
cand <- read.delim(CAND)
Z <- as.matrix(read.delim(gzfile(ZF), row.names = 1, check.names = FALSE)); Z[is.na(Z)] <- 0
P <- readRDS(file.path(OUT, "profiles.rds")); stopifnot(identical(P$ids, cand$ID))
mem <- read.delim(file.path(OUT, "module_members_p1e-05.tsv"))
mods <- read.delim(file.path(OUT, "modules_p1e-05.tsv"))
sc <- read.delim(gzfile(file.path(OUT, "screen_pairs_p1e-3.tsv.gz")))
genes <- read.delim(gzfile("/mnt/project-files/papers/LRCQ/results/real/genes_all.genes.tsv.gz"))
near <- function(id) { ch <- as.integer(sub(":.*", "", id)); bp <- as.integer(strsplit(id, ":")[[1]][2])
  g <- genes[genes$chr == ch, ]; d <- pmax(g$start - bp, bp - g$end, 0); g$gene[which.min(d)] }
getR <- function(win, ids) { sn <- read.delim(file.path(LDD, sprintf("w%s.snps.tsv", win)))
  ix <- match(ids, sn$ID); stopifnot(!anyNA(ix))
  matrix(readBin(file.path(LDD, sprintf("w%s.R.f32", win)), "numeric", nrow(sn)^2, size = 4), nrow(sn))[ix, ix, drop = FALSE] }
# gene units: (module, gene, window) with oriented tag weights
mem$win <- cand$window[match(mem$ID, cand$ID)]
wins <- unique(mem$win)
gmap <- setNames(lapply(wins, function(w) { ix <- which(cand$window == w); setNames(vapply(cand$ID[ix], near, ""), ix) }), wins)
units <- unique(mem[, c("module", "gene", "win")])
units$tags <- lapply(seq_len(nrow(units)), function(u) {
  ix <- as.integer(names(gmap[[units$win[u]]])[gmap[[units$win[u]]] == units$gene[u]])
  lead <- match(mem$ID[mem$module == units$module[u] & mem$gene == units$gene[u] & mem$win == units$win[u]][1], cand$ID)
  s <- sign(cor(t(P$X[ix, , drop = FALSE]), P$X[lead, ]))[, 1]; s[s == 0] <- 1
  setNames(s, ix) })
units$n_tags <- lengths(units$tags)
units$chr <- cand$CHR[match(units$win, cand$window)]
units$pos <- vapply(seq_len(nrow(units)), function(u) median(cand$BP[as.integer(names(units$tags[[u]]))]), 0)
message("gene units: ", nrow(units), "; tags per gene: ", paste(range(units$n_tags), collapse = "-"))
# distal window pairs within each module
wp <- do.call(rbind, lapply(split(seq_len(nrow(units)), units$module), function(ix) {
  w <- unique(units$win[ix]); if (length(w) < 2) return(NULL)
  cb <- t(combn(w, 2)); data.frame(module = units$module[ix[1]], winA = cb[, 1], winB = cb[, 2]) }))
far <- function(a, b) { ua <- units[units$win == a, ][1, ]; ub <- units[units$win == b, ][1, ]
  ua$chr != ub$chr || abs(ua$pos - ub$pos) >= MIN_DIST }
wp <- wp[mapply(far, wp$winA, wp$winB), ]
message("distal window pairs: ", nrow(wp))
edge <- sc[sc$p < PE, ]
run_wp <- function(k) {
  set.seed(1000 + k); m <- wp$module[k]
  uA <- which(units$module == m & units$win == wp$winA[k]); uB <- which(units$module == m & units$win == wp$winB[k])
  A <- which(cand$window == wp$winA[k]); B <- which(cand$window == wp$winB[k])
  V <- function(us, S) sapply(us, function(u) { v <- numeric(length(S)); t <- units$tags[[u]]
    v[match(as.integer(names(t)), S)] <- t; v })
  VA <- matrix(V(uA, A), length(A)); VB <- matrix(V(uB, B), length(B))
  g <- lrcp_gene(Z[A, , drop = FALSE], Z[B, , drop = FALSE], getR(wp$winA[k], cand$ID[A]), getR(wp$winB[k], cand$ID[B]),
                 VA, VB, n = P$n, gcov = P$gcov, M = M, intercept = P$intercept, n_boot = n_boot, n_sim = 0)
  do.call(rbind, lapply(seq_along(uA), function(i) do.call(rbind, lapply(seq_along(uB), function(j) {
    tA <- cand$ID[as.integer(names(units$tags[[uA[i]]]))]; tB <- cand$ID[as.integer(names(units$tags[[uB[j]]]))]
    e <- edge[(edge$ID1 %in% tA & edge$ID2 %in% tB) | (edge$ID1 %in% tB & edge$ID2 %in% tA), ]
    s <- sc[(sc$ID1 %in% tA & sc$ID2 %in% tB) | (sc$ID1 %in% tB & sc$ID2 %in% tA), ]
    data.frame(module = m, geneA = units$gene[uA[i]], geneB = units$gene[uB[j]], winA = wp$winA[k], winB = wp$winB[k],
               n_tagsA = units$n_tags[uA[i]], n_tagsB = units$n_tags[uB[j]], linked = nrow(e) > 0,
               best_tag_p = if (nrow(s)) min(s$p) else NA, best_tag_z = if (nrow(s)) s$z[which.min(s$p)] else NA,
               C = g$C[i, j], z = g$z[i, j], p_boot = g$p[i, j], p_maxT_wp = g$p_maxT[i, j],
               p_norm = g$p_norm_diagnostic[i, j], z_gls = g$z_gls[i, j], R_AB = g$R_AB[i, j]) }))))
}
res <- do.call(rbind, mclapply(seq_len(nrow(wp)), run_wp, mc.cores = as.integer(Sys.getenv("NCORES", "4")),
                               mc.preschedule = FALSE))
res$q_BH <- p.adjust(res$p_boot, "BH")
res$narrow_module <- (mods$min_q_eff_class[match(res$module, mods$module)] < 5) %in% TRUE  # NA (no class null run) = not narrow, as in fig4.py
res <- res[order(res$module, res$p_boot), ]
write.table(res, file.path(OUT, "module_gene_burden.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res[, c("module", "geneA", "geneB", "n_tagsA", "n_tagsB", "linked", "z", "best_tag_z", "p_boot", "p_maxT_wp", "q_BH")], digits = 3)
