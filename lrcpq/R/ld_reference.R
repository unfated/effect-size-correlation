# Real LD references: PLINK reader, LD blocks, allele harmonisation.

#' Berisa and Pickrell (2016) EUR LD blocks (GRCh37)
#'
#' Approximately independent LD blocks from LDetect, shipped with the
#' package (1,703 blocks).
#' @return Data frame with \code{chr} (integer), \code{start}, \code{stop}.
#' @export
ld_blocks_eur <- function() {
  f <- system.file("extdata", "ldetect_EUR_hg19.bed", package = "lrcpq")
  b <- utils::read.table(f, header = TRUE, strip.white = TRUE,
                         stringsAsFactors = FALSE)
  data.frame(chr = as.integer(sub("chr", "", b$chr)), start = b$start,
             stop = b$stop)
}

#' Read a PLINK .bim file
#' @param prefix PLINK file prefix.
#' @export
read_bim <- function(prefix) {
  b <- utils::read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE,
                         colClasses = c("integer", "character", "numeric",
                                        "integer", "character", "character"))
  names(b) <- c("chr", "rsid", "cm", "pos", "a1", "a2")
  b
}

#' Read genotypes from a PLINK .bed file
#'
#' Pure-R reader for SNP-major PLINK 1 binary files.
#' @param prefix PLINK file prefix.
#' @param snps Integer indices of SNPs (rows of the .bim) to read; default all.
#' @param n Number of individuals (read from the .fam when \code{NULL}).
#' @return Individuals-by-SNPs matrix of A1 allele counts (\code{NA} missing).
#' @export
read_bed <- function(prefix, snps = NULL, n = NULL) {
  if (is.null(n)) n <- length(readLines(paste0(prefix, ".fam")))
  nbytes <- ceiling(n / 4)
  if (is.null(snps)) {
    nsnp <- length(readLines(paste0(prefix, ".bim")))
    snps <- seq_len(nsnp)
  }
  con <- file(paste0(prefix, ".bed"), "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 3)
  if (!identical(magic, as.raw(c(0x6c, 0x1b, 0x01))))
    stop("not a SNP-major PLINK .bed file")
  # read the contiguous span once, then decode
  lo <- min(snps); hi <- max(snps)
  seek(con, 3 + (lo - 1) * nbytes)
  raw <- readBin(con, "raw", (hi - lo + 1) * nbytes)
  raw <- matrix(raw, nbytes, hi - lo + 1)[, snps - lo + 1, drop = FALSE]
  v <- as.integer(raw)
  codes <- rbind(bitwAnd(v, 3L), bitwAnd(bitwShiftR(v, 2L), 3L),
                 bitwAnd(bitwShiftR(v, 4L), 3L), bitwShiftR(v, 6L))
  # 00 hom A1 -> 2, 01 missing, 10 het -> 1, 11 hom A2 -> 0
  lut <- c(2L, NA_integer_, 1L, 0L)
  G <- matrix(lut[codes + 1L], 4 * nbytes, length(snps))
  G[seq_len(n), , drop = FALSE]
}

#' Harmonise summary-statistic alleles to a reference
#'
#' Matches by chromosome and position and returns the sign that aligns each
#' Z-score to the reference A1 allele. Ambiguous (A/T, C/G) SNPs and allele
#' mismatches are dropped.
#'
#' @param chr,pos,effect,other Summary statistics: chromosome, position,
#'   effect allele (e.g. Pan-UKB \code{alt}) and other allele (\code{ref}).
#' @param ref A data frame with \code{chr}, \code{pos}, \code{a1}, \code{a2}
#'   (e.g. from [read_bim()]).
#' @param drop_ambiguous Drop strand-ambiguous SNPs.
#' @return Data frame with \code{i} (row in the summary statistics),
#'   \code{j} (row in \code{ref}) and \code{sign} (+1 / -1).
#' @export
harmonise_alleles <- function(chr, pos, effect, other, ref,
                              drop_ambiguous = TRUE) {
  key_s <- paste(chr, pos, sep = ":")
  key_r <- paste(ref$chr, ref$pos, sep = ":")
  j <- match(key_s, key_r)
  i <- which(!is.na(j))
  j <- j[i]
  e <- toupper(effect[i]); o <- toupper(other[i])
  a1 <- toupper(ref$a1[j]); a2 <- toupper(ref$a2[j])
  sign <- ifelse(e == a1 & o == a2, 1, ifelse(e == a2 & o == a1, -1, NA))
  if (drop_ambiguous) {
    amb <- paste0(e, o) %in% c("AT", "TA", "CG", "GC")
    sign[amb] <- NA
  }
  keep <- !is.na(sign)
  data.frame(i = i[keep], j = j[keep], sign = sign[keep])
}

#' LD reference from per-chromosome PLINK files
#'
#' Builds an \code{lrcpq_ld} object whose \code{get()} computes LD on demand
#' from the reference genotypes (e.g. 1000 Genomes EUR, see
#' \code{inst/scripts/build_1000g_eur_hm3.sh}). Genotypes of one chromosome
#' are cached at a time.
#'
#' @param prefixes Named character vector of PLINK prefixes, names are
#'   chromosomes (e.g. \code{c("22" = "eur_hm3_chr22")}).
#' @param snps Optional data frame of the SNPs to use (columns \code{chr},
#'   \code{pos}), in genomic order; defaults to all reference SNPs. Only
#'   reference SNPs matching \code{snps} are kept.
#' @param blocks LD blocks (default [ld_blocks_eur()]).
#' @param maf_min Drop reference SNPs with lower minor allele frequency.
#' @return An \code{lrcpq_ld} object; \code{$snp} holds \code{chr},
#'   \code{pos}, \code{rsid}, \code{a1}, \code{a2}, \code{maf}, \code{block}.
#' @export
ld_from_plink <- function(prefixes, snps = NULL, blocks = ld_blocks_eur(),
                          maf_min = 0.01) {
  bims <- lapply(names(prefixes), function(cn) {
    b <- read_bim(prefixes[[cn]])
    b$file <- prefixes[[cn]]
    b$row <- seq_len(nrow(b))
    b
  })
  bim <- do.call(rbind, bims)
  if (!is.null(snps)) {
    keep <- paste(bim$chr, bim$pos) %in% paste(snps$chr, snps$pos)
    bim <- bim[keep, ]
  }
  # MAF filter, computed chromosome by chromosome
  bim$maf <- NA_real_
  for (f in unique(bim$file)) {
    sel <- which(bim$file == f)
    G <- read_bed(f, bim$row[sel])
    p <- colMeans(G, na.rm = TRUE) / 2
    bim$maf[sel] <- pmin(p, 1 - p)
  }
  bim <- bim[bim$maf >= maf_min, ]
  bim <- bim[order(bim$chr, bim$pos), ]
  bim$block <- assign_blocks(bim$chr, bim$pos, blocks)
  N <- length(readLines(paste0(bim$file[1], ".fam")))
  cache <- new.env()
  cache$file <- ""
  get <- function(idx) {
    files <- unique(bim$file[idx])
    Gs <- NULL
    for (f in files) {
      if (cache$file != f) {
        sel <- which(bim$file == f)
        G <- read_bed(f, bim$row[sel])
        G <- scale(G)
        G[is.na(G)] <- 0
        cache$G <- G
        cache$map <- sel
        cache$file <- f
      }
      cols <- match(idx[bim$file[idx] == f], cache$map)
      Gs <- cbind(Gs, cache$G[, cols, drop = FALSE])
    }
    R <- crossprod(Gs) / (N - 1)
    diag(R) <- 1
    R
  }
  snp <- bim[, c("chr", "pos", "rsid", "a1", "a2", "maf", "block")]
  rownames(snp) <- NULL
  new_ld(get, nrow(bim), bim$block, snp, N)
}
