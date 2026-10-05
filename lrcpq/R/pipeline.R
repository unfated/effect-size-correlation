# End-to-end pipelines for the four-stage LRCQ/LRCP estimator.

#' Run LRCQ end to end (stages 1-3)
#'
#' Stage 1-2 values are taken from \code{gcov}/\code{intercept} when given
#' (e.g. genome-wide LDSC, the Pan-UKB h2 manifest, or the truth in a
#' simulation); otherwise they are estimated with [ldsc_matrix()] from
#' \code{Z} and \code{ldscore}, which needs genome-scale data. Stage 3 then
#' estimates enrichment per tag clump (default) or per SNP.
#'
#' @param Z m-by-q Z-scores aligned to \code{ld}.
#' @param ld An \code{lrcpq_ld} object.
#' @param n Sample sizes.
#' @param M Number of SNPs heritability refers to (the size of the SNP
#'   universe behind \code{D}, e.g. 1,094,844 HapMap3 SNPs).
#' @param gcov,intercept Stage 1-2 inputs, or \code{NULL} to estimate.
#' @param ldscore LD scores for stage 1-2 estimation.
#' @param tag_r2 Clump threshold (\code{NULL} for per-SNP estimates).
#' @param method,rectify Passed to [lrcq()].
#' @param windows Optional precomputed windows (see [ld_from_windows()]).
#' @param ... Passed to [lrcq()].
#' @return List with \code{stage12} (gcov, intercept) and \code{w}
#'   (data frame from [lrcq()], with SNP metadata when available).
#' @export
run_lrcq <- function(Z, ld, n, M, gcov = NULL, intercept = NULL, ldscore = NULL,
                     tag_r2 = 0.5, method = "wls", rectify = "C",
                     windows = NULL, ...) {
  if (is.null(gcov)) {
    if (is.null(ldscore)) stop("give gcov/intercept or ldscore for stage 1-2")
    st <- ldsc_matrix(Z, ldscore, n, M)
    gcov <- st$gcov
    if (is.null(intercept)) intercept <- st$intercept
  }
  if (is.null(intercept)) intercept <- diag(ncol(Z))
  w <- lrcq(Z, ld, n, gcov, M, intercept, method = method, rectify = rectify,
            tag_r2 = tag_r2, windows = windows, ...)
  if (!is.null(ld$snp)) w <- cbind(ld$snp, w)
  list(stage12 = list(gcov = gcov, intercept = intercept), w = w)
}

#' Run LRCP for distal window pairs (stage 4)
#'
#' Screens tags on stage-3 enrichment, then estimates genetic-effect
#' correlation for every pair of screened tags in different windows that
#' share no LD (different chromosomes, or further apart than
#' \code{min_distance}), using [lrcp_distal()] and optionally [lrcp_mle()].
#' P-values are adjusted by Benjamini-Hochberg (primary) and
#' Benjamini-Yekutieli (conservative check) across all pairs.
#'
#' @param Z,ld,n,M As in [run_lrcq()].
#' @param gcov,intercept Stage 1-2 values.
#' @param lrcq_fit Output of [run_lrcq()] (its \code{w} table).
#' @param threshold Screening cutoff on enrichment.
#' @param z_min Optional screening cutoff on the enrichment z-score.
#' @param min_distance Minimum distance (bp) between windows on the same
#'   chromosome.
#' @param method LRCP estimator for [lrcp_distal()].
#' @param denominators See [lrcp_distal()].
#' @param mle Also run [lrcp_mle()] per window pair.
#' @return Data frame with one row per tag pair.
#' @export
run_lrcp <- function(Z, ld, n, M, gcov, intercept = NULL, lrcq_fit,
                     threshold = 10, z_min = -Inf, min_distance = 5e6,
                     method = "gls", denominators = "R", mle = FALSE) {
  q <- ncol(Z)
  if (is.null(intercept)) intercept <- diag(q)
  wt <- lrcq_fit$w
  cand <- which(!is.na(wt$w_raw) & wt$w_raw >= threshold & wt$z >= z_min)
  if (length(cand) < 2) return(data.frame())
  blk <- ld$block_id[cand]
  groups <- split(cand, blk)
  chr <- if (!is.null(ld$snp$chr)) ld$snp$chr else rep(1, ld$m)
  pos <- if (!is.null(ld$snp$pos)) ld$snp$pos else seq_len(ld$m)
  out <- list()
  gnames <- names(groups)
  for (i in seq_along(groups)) for (j in seq_along(groups)) {
    if (j <= i) next
    A <- groups[[i]]; B <- groups[[j]]
    same_chr <- chr[A[1]] == chr[B[1]]
    if (same_chr && min(abs(outer(pos[A], pos[B], "-"))) < min_distance) next
    # each candidate's own block as its region
    ia <- which(ld$block_id == ld$block_id[A[1]])
    ib <- which(ld$block_id == ld$block_id[B[1]])
    RA <- ld$get(ia); RB <- ld$get(ib)
    wa <- numeric(length(ia)); wa[match(A, ia)] <- wt$w_raw[A]
    wb <- numeric(length(ib)); wb[match(B, ib)] <- wt$w_raw[B]
    f <- lrcp_distal(Z[ia, , drop = FALSE], Z[ib, , drop = FALSE], RA, RB, wa, wb,
                     n, gcov, M, intercept, S1 = match(A, ia), S2 = match(B, ib),
                     method = method, denominators = denominators)
    res <- data.frame(snp1 = rep(A, length(B)), snp2 = rep(B, each = length(A)),
                      rho = as.vector(f$rho), se = as.vector(f$se))
    if (mle) {
      g <- lrcp_mle(Z[ia, , drop = FALSE], Z[ib, , drop = FALSE], RA, RB, wa, wb,
                    n, gcov, M, intercept, S1 = match(A, ia), S2 = match(B, ib))
      res$rho_mle <- as.vector(g$rho); res$se_mle <- as.vector(g$se)
      res$lrt_window_pair <- g$lrt
    }
    out[[length(out) + 1]] <- res
  }
  res <- do.call(rbind, out)
  if (is.null(res)) return(data.frame())
  res$z <- res$rho / res$se
  res$p <- 2 * stats::pnorm(-abs(res$z))
  res$q_BH <- stats::p.adjust(res$p, "BH")
  res$q_BY <- stats::p.adjust(res$p, "BY")
  if (!is.null(ld$snp)) {
    id <- if (!is.null(ld$snp$ID)) ld$snp$ID else paste(chr, pos, sep = ":")
    res$id1 <- id[res$snp1]; res$id2 <- id[res$snp2]
  }
  res
}
