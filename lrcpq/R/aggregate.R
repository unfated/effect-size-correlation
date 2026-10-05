# Aggregated inference (theory supplement S7, S8):
#   category-level LRCQ (annotation enrichment, annotation regression,
#   phenotype-category contrasts) and gene/pathway-level LRCP with the
#   conditional Gaussian null.

#' Exact variance of a weighted quadratic form in Z-scores
#'
#' Theorem 3.10 / eq. (7.1) of the theory supplement: for Gaussian Z with
#' \eqn{Cov(vec Z) = \Gamma \otimes R + T \otimes G},
#' \deqn{Var(\sum_a u_a z_a' A z_a) = 2[\tau_{RR} u'(\Gamma\circ\Gamma)u +
#'   2\tau_{RG} u'(\Gamma\circ T)u + \tau_{GG} u'(T\circ T)u],}
#' \eqn{\tau_{XY} = tr(AXAY)}.
#'
#' @param A Symmetric m-by-m matrix of the quadratic form.
#' @param u Length-q trait weights.
#' @param R LD matrix.
#' @param G \eqn{R S_w R} (genetic part of the Z covariance, per unit of
#'   \eqn{t_{ab}}).
#' @param Gamma q-by-q intercept matrix.
#' @param Tm q-by-q scale matrix \eqn{t_{ab} = \sqrt{n_a n_b} h_{ab}/M}.
#' @export
quadform_var <- function(A, u, R, G, Gamma, Tm) {
  AR <- A %*% R
  AG <- A %*% G
  tRR <- sum(AR * t(AR)); tRG <- sum(AR * t(AG)); tGG <- sum(AG * t(AG))
  2 * (tRR * sum(u * ((Gamma^2) %*% u)) + 2 * tRG * sum(u * ((Gamma * Tm) %*% u)) +
         tGG * sum(u * ((Tm^2) %*% u)))
}

#' Pseudo-inverse of a symmetric PSD matrix
#' @keywords internal
psd_inverse <- function(D, tol = 1e-8) {
  e <- eigen((D + t(D)) / 2, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  e$vectors[, keep, drop = FALSE] %*% (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
}

#' Trait weights for LRCQ aggregates
#'
#' @param n,gcov,M As in [lrcq()].
#' @param type \code{"power"}: \eqn{u = s/\|s\|^2} (OLS, Proposition 3.2);
#'   \code{"equal"}: \eqn{u_a = 1/(q s_a)}; \code{"contrast"}: \eqn{u_a =
#'   s_a/\|s_{g_1}\|^2} on group 1 and \eqn{-s_a/\|s_{g_2}\|^2} on group 2.
#' @param groups For \code{"contrast"}: length-q vector with values 1, 2 or
#'   \code{NA} (trait unused).
#' @export
trait_weights <- function(n, gcov, M, type = c("power", "equal", "contrast"),
                          groups = NULL) {
  type <- match.arg(type)
  h2 <- if (is.null(dim(gcov))) gcov else diag(gcov)
  s <- n * h2 / M
  switch(type,
    power = s / sum(s^2),
    equal = 1 / (length(s) * s),
    contrast = {
      u <- numeric(length(s))
      g1 <- which(groups == 1); g2 <- which(groups == 2)
      u[g1] <- s[g1] / sum(s[g1]^2)
      u[g2] <- -s[g2] / sum(s[g2]^2)
      u
    })
}

#' Category-level LRCQ: fold-enrichment of a SNP annotation
#'
#' Clump-level, model-free category enrichment (theory supplement S7.2a):
#' \deqn{\hat E_c = \sum_k f_k \hat W_k / \sum_k f_k g_k,}
#' where \eqn{\hat W_k} are tag-set clump totals ([lrcq_window()] with
#' \code{tags}), \eqn{g_k} the clump sizes and \eqn{f_k = |c \cap k|/g_k}.
#' Unlike the per-SNP average \eqn{1_c'\hat w/|c|}, whose variance is
#' unbounded when a category splits tight LD clumps, this stays well
#' conditioned. Its estimand is the overlap-weighted enrichment of the clumps
#' that \eqn{c} touches; it equals \eqn{E_c} when \eqn{c} is a union of
#' clumps or \eqn{w} is flat within clumps, and is otherwise attenuated
#' toward the clump mean. Use [lrcq_annot_regression()] for the
#' annotation-model estimand. A significant result means "clumps overlapping
#' \eqn{c} are enriched"; attributing it to the SNPs in \eqn{c} rather than
#' their proxies needs the annotation model or fine-mapping.
#'
#' The variance is eq. (7.1) with \eqn{A = Diag(L'f)/\sum f_k g_k}, exact
#' under cross-trait correlation and sample overlap. Its genetic part needs a
#' plug-in \eqn{\hat G = R\,Diag(\hat w_+)R}; following S7.4a the plug-in is
#' always clipped at 0, and by default debiased: the quadratic term
#' \eqn{\tau_{GG}} is reduced by \eqn{\sum_{st} Cov(\hat W_s, \hat W_t)
#' [(RAR)_{st}]^2} (floored at 0). With
#' \code{u = trait_weights(..., "contrast", groups)} it returns the
#' phenotype-category contrast \eqn{\hat E_c^{(g_1)} - \hat E_c^{(g_2)}}
#' (S7.4).
#'
#' @param Z m-by-q Z-scores aligned to \code{ld}.
#' @param ld An \code{lrcpq_ld} object (blocks are treated as independent).
#' @param annot Logical or 0/1 vector (length m) marking the annotation, or a
#'   matrix with one column per annotation.
#' @param n,gcov,M,intercept As in [lrcq()].
#' @param u Trait weights (default power weights, matching OLS).
#' @param tag_r2 Pruning threshold defining clumps.
#' @param plugin Rectification applied to the power-weighted clump totals
#'   that enter the genetic part of the variance (then clipped at 0).
#' @param w_plugin Optional per-SNP enrichment for the variance instead (e.g.
#'   the truth in a simulation, or estimates from independent traits); used
#'   as given, without debiasing.
#' @param variance \code{"debiased"} (default) or \code{"naive"} plug-in.
#' @param n_boot For contrasts: number of parametric-bootstrap draws for the
#'   studentised contrast under \eqn{H_0} (S7.4a; default 200 for contrasts,
#'   0 to skip). When it runs, \code{p} is the bootstrap p-value (the test
#'   to report) and the normal-reference value is kept as
#'   \code{p_norm_diagnostic}, a diagnostic only: the debiased plug-in is
#'   not a safe substitute for the bootstrap (S7.4a). Draws
#'   \eqn{Z^* \sim N(0, \Gamma\otimes R + T\otimes \hat G_+)} from the pooled,
#'   clipped fit and recomputes the whole statistic, plug-in included.
#'
#' @section Phenotype-category contrasts: with sparse, heavy-tailed
#'   enrichment a same-data plug-in makes the contrast z-statistic
#'   self-normalising, so the normal-reference test is conservative (in a
#'   4-cluster, full-overlap simulation the naive plug-in gave null SD 0.81
#'   and no rejections at 5\%; the debiased plug-in is closer). The
#'   debiased plug-in is used only as the bootstrap's studentising
#'   denominator. The bootstrap p-value is calibrated (theory supplement check C3; in the
#'   package's own check with q = 40 traits in 4 clusters and 6 dominant
#'   loci, type I error 0.040 at 5\% against 0.030 naive and 0.073 debiased)
#'   and is the default contrast test.
#' @return Data frame with one row per annotation: \code{size}, \code{E},
#'   \code{se}, \code{z} (for \eqn{E = 1}, or 0 for a contrast), \code{p},
#'   and \code{p_norm_diagnostic} when the bootstrap ran.
#' @export
lrcq_category <- function(Z, ld, annot, n, gcov, M, intercept = NULL,
                          u = NULL, tag_r2 = 0.5, plugin = "A",
                          w_plugin = NULL, variance = c("debiased", "naive"),
                          n_boot = NULL) {
  variance <- match.arg(variance)
  Z <- as.matrix(Z)
  q <- ncol(Z)
  annot <- as.matrix(annot) * 1
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  if (is.null(u)) u <- trait_weights(n, gcov, M, "power")
  upow <- trait_weights(n, gcov, M, "power")
  Tm <- scale_matrix(n, gcov, M)
  Cd <- diag(intercept)
  K <- ncol(annot)
  contrast <- any(u < 0)
  if (is.null(n_boot)) n_boot <- if (contrast) 200 else 0
  # trait-weight constants of eq. (7.1)
  qf <- function(x, X) sum(x * (X %*% x))
  cu <- c(qf(u, intercept^2), qf(u, intercept * Tm), qf(u, Tm^2))
  cp <- c(qf(upow, intercept^2), qf(upow, intercept * Tm), qf(upow, Tm^2))
  # per-block setup that does not depend on Z
  blocks <- list()
  for (b in unique(ld$block_id)) {
    idx <- which(ld$block_id == b)
    if (!any(annot[idx, ] != 0)) next
    R <- ld$get(idx)
    D <- ld_r2(R, ld$N)
    e <- eigen(D, symmetric = TRUE)
    if (min(e$values) < 0) D <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
    tg <- prune_ld(R, tag_r2)
    DT <- D[, tg, drop = FALSE]
    L <- solve(crossprod(DT), t(DT))                       # p x m
    clump <- tg[apply(abs(R[, tg, drop = FALSE]), 1, which.max)]
    gk <- as.vector(table(factor(clump, levels = tg)))
    lf <- matrix(0, length(idx), K); den <- numeric(K)
    for (k in seq_len(K)) {
      fk <- as.vector(tapply(annot[idx, k], factor(clump, levels = tg), sum)) / gk
      fk[is.na(fk)] <- 0
      lf[, k] <- as.vector(crossprod(L, fk))
      den[k] <- sum(fk * gk)
    }
    blocks[[length(blocks) + 1]] <- list(idx = idx, R = R, L = L, tg = tg, lf = lf,
                                         den = den, Rh = psd_factor(R))
  }
  # statistic for one set of Z-scores: numerators and variances summed over blocks
  stat <- function(Zs) {
    num <- vr <- den <- numeric(K)
    for (bl in blocks) {
      Y <- Zs[bl$idx, , drop = FALSE]^2 - matrix(Cd, length(bl$idx), q, byrow = TRUE)
      ybar <- as.vector(Y %*% u)
      R <- bl$R
      if (!is.null(w_plugin)) {
        wexp <- pmax(w_plugin[bl$idx], 0)
      } else {
        Wpow <- pmax(rectify_w(as.vector(bl$L %*% (Y %*% upow)), plugin), 0)
        wexp <- numeric(length(bl$idx)); wexp[bl$tg] <- Wpow
      }
      G <- R %*% (wexp * R)
      if (variance == "debiased" && is.null(w_plugin)) {
        # Cov of the power-weighted clump totals (S7.4a, item 1)
        Vy <- 2 * (cp[1] * R^2 + 2 * cp[2] * R * G + cp[3] * G^2)
        LT <- bl$L
        CovW <- LT %*% Vy %*% t(LT)
        Rt <- R[, bl$tg, drop = FALSE]
      }
      for (k in seq_len(K)) {
        a <- bl$lf[, k]
        if (!any(a != 0)) next
        num[k] <- num[k] + sum(a * ybar)
        den[k] <- den[k] + bl$den[k]
        AR <- a * R; AG <- a * G
        tRR <- sum(AR * t(AR)); tRG <- sum(AR * t(AG)); tGG <- sum(AG * t(AG))
        if (variance == "debiased" && is.null(w_plugin)) {
          RAR <- crossprod(Rt, a * Rt)
          tGG <- max(tGG - sum(CovW * RAR^2), 0)
        }
        vr[k] <- vr[k] + 2 * (tRR * cu[1] + 2 * tRG * cu[2] + tGG * cu[3])
      }
    }
    E <- num / den
    se <- sqrt(vr) / den
    list(E = E, se = se, z = if (contrast) E / se else (E - 1) / se)
  }
  st <- stat(Z)
  out <- data.frame(annotation = colnames(annot) %||% seq_len(K),
                    size = colSums(annot), E = st$E, se = st$se, z = st$z,
                    p = 2 * stats::pnorm(-abs(st$z)))
  if (contrast && n_boot > 0) {
    # H0 fit: pooled, clipped power-weighted clump totals per block
    Th <- psd_factor(Tm); Ch <- psd_factor(intercept)
    Gh <- lapply(blocks, function(bl) {
      Y <- Z[bl$idx, , drop = FALSE]^2 - matrix(Cd, length(bl$idx), q, byrow = TRUE)
      W <- pmax(rectify_w(as.vector(bl$L %*% (Y %*% upow)), plugin), 0)
      wexp <- numeric(length(bl$idx)); wexp[bl$tg] <- W
      psd_factor(bl$R %*% (wexp * bl$R))
    })
    zb <- replicate(n_boot, {
      Zs <- Z
      for (i in seq_along(blocks)) {
        bl <- blocks[[i]]; mb <- length(bl$idx)
        Zs[bl$idx, ] <- bl$Rh %*% matrix(stats::rnorm(mb * q), mb) %*% t(Ch) +
          Gh[[i]] %*% matrix(stats::rnorm(mb * q), mb) %*% t(Th)
      }
      stat(Zs)$z
    })
    zb <- matrix(zb, nrow = K)
    out$p_norm_diagnostic <- out$p
    out$p <- (1 + rowSums(abs(zb) >= abs(st$z))) / (n_boot + 1)
  }
  out
}

#' Annotation regression for LRCQ
#'
#' \eqn{\hat\tau = (\sum_b A_b' D_b^2 A_b)^{-1} \sum_b A_b' D_b \bar y_b}
#' (theory supplement S7.3): conditional effects of annotations on \eqn{w},
#' stable under tight LD because no \eqn{D^{-1}} is needed. The estimand is
#' the \eqn{D^2}-weighted projection of \eqn{w} on the annotation span.
#' Standard errors are exact (eq. 7.1 with
#' \eqn{A_b = Diag(D_b A_b S^{-1} e_k)}), with a delete-one-block jackknife
#' alongside.
#'
#' @param annot m-by-K annotation matrix; a column of 1s is added if absent.
#' @inheritParams lrcq_category
#' @return Data frame with \code{tau}, \code{se} (exact), \code{se_jk},
#'   \code{z}.
#' @export
lrcq_annot_regression <- function(Z, ld, annot, n, gcov, M, intercept = NULL,
                                  u = NULL) {
  Z <- as.matrix(Z)
  q <- ncol(Z)
  A <- as.matrix(annot)
  if (!any(apply(A, 2, function(x) all(x == 1)))) A <- cbind(base = 1, A)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  if (is.null(u)) u <- trait_weights(n, gcov, M, "power")
  Tm <- scale_matrix(n, gcov, M)
  Cd <- diag(intercept)
  K <- ncol(A)
  blocks <- unique(ld$block_id)
  XtX <- array(0, c(K, K, length(blocks)))
  Xty <- matrix(0, K, length(blocks))
  cache <- vector("list", length(blocks))
  upow <- trait_weights(n, gcov, M, "power")
  ybar_pow <- function(idx) as.vector((Z[idx, , drop = FALSE]^2 -
    matrix(Cd, length(idx), q, byrow = TRUE)) %*% upow)
  for (bi in seq_along(blocks)) {
    idx <- which(ld$block_id == blocks[bi])
    R <- ld$get(idx)
    D <- ld_r2(R, ld$N)
    ybar <- as.vector((Z[idx, , drop = FALSE]^2 -
                         matrix(Cd, length(idx), q, byrow = TRUE)) %*% u)
    DA <- D %*% A[idx, , drop = FALSE]
    XtX[, , bi] <- crossprod(DA)
    Xty[, bi] <- crossprod(DA, ybar)
    # genetic part of the variance from rectified clump totals (model-free)
    e <- eigen(D, symmetric = TRUE)
    Dp <- if (min(e$values) < 0) e$vectors %*% (pmax(e$values, 0) * t(e$vectors)) else D
    tg <- prune_ld(R, 0.5)
    DT <- Dp[, tg, drop = FALSE]
    Wk <- rectify_w(as.vector(solve(crossprod(DT), crossprod(DT, ybar_pow(idx)))), "A")
    wexp <- numeric(length(idx)); wexp[tg] <- pmax(Wk, 0)
    cache[[bi]] <- list(idx = idx, R = R, DA = DA, G = R %*% (wexp * R))
  }
  S <- apply(XtX, c(1, 2), sum); sv <- rowSums(Xty)
  tau <- solve(S, sv)
  Si <- solve(S)
  vr <- numeric(K)
  for (bi in seq_along(blocks)) {
    cb <- cache[[bi]]
    G <- cb$G
    for (k in seq_len(K)) {
      d <- as.vector(cb$DA %*% Si[, k])
      vr[k] <- vr[k] + quadform_var(diag(d, length(d)), u, cb$R, G, intercept, Tm)
    }
  }
  nb <- length(blocks)
  jk <- t(vapply(seq_len(nb), function(b) solve(S - XtX[, , b], sv - Xty[, b]), numeric(K)))
  se_jk <- sqrt((nb - 1) / nb * colSums(sweep(jk, 2, colMeans(jk))^2))
  data.frame(annotation = colnames(A) %||% seq_len(K), tau = tau,
             se = sqrt(vr), se_jk = se_jk, z = tau / sqrt(vr))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Gene-level LRCP between two distal windows
#'
#' Burden covariance \eqn{\hat C_{GH} = \|s\|^{-2}\sum_a s_a (g_G'z_{Aa})
#' (h_H'z_{Ba})}, \eqn{g_G = R_A[T,T]^{-1} v_G} (theory supplement eq. 8.1),
#' with the conditional Gaussian null of eq. (8.2): given window B, the
#' statistic is linear in window A's Z-scores, so its null variance needs
#' only window A's trait covariance \eqn{\Sigma_x = \tau_R \Gamma + \tau_G T}.
#' This is the primary per-pair test (S8.5): unlike the GLS null SE it does
#' not depend on window B's estimated enrichment. Within the window pair, a
#' max-|T| adjustment controls the family-wise error rate.
#'
#' @section Non-separable loci (S8.5): loci acting through a narrow set of
#'   mediators (e.g. lipid loci through LDL) have a class-specific trait
#'   covariance with few effective dimensions. Two independent loci of such
#'   a class have nearly collinear trait profiles, which a test built on the
#'   genome-wide \eqn{T} reads as correlation. Supply \code{VclassA}
#'   (and \code{VclassB}) from [class_covariance()] on held-out loci of the
#'   same profile class ([profile_classes()]) to test against the class null
#'   \eqn{\Sigma_x = \tau_R\Gamma + \kappa_A \hat V_{class}}, with
#'   \eqn{\kappa_A} trace-matched to the realised signal of \eqn{x}. When
#'   \eqn{q_{eff}(\hat V_{class}) < 5} alignment within the class is not
#'   testable (a warning is given); report \eqn{\hat\rho} descriptively as
#'   shared-mediator alignment.
#'
#' @param ZA,ZB Tag Z-scores (tags-by-q) of the two windows.
#' @param RA,RB Tag LD matrices.
#' @param VA,VB Tag-by-gene weight matrices (e.g. indicators of each gene's
#'   tags); identity gives tag-level results.
#' @param n,gcov,M,intercept As in [lrcq()]. \code{gcov} must hold genome-wide
#'   genetic covariances (LDSC rg), never a diagonal matrix.
#' @param GA,GB Genetic part \eqn{R S_w R} on each window's tags (default from
#'   [local_variance()] assuming no local correlation, clipped at 0).
#' @param VclassA,VclassB Optional class trait covariances (q-by-q), one
#'   matrix for all genes or a list with one per gene.
#' @param n_sim Gaussian draws for the max-|T| adjustment (0 to skip).
#' @return List with \code{C} (genes A by genes B), \code{se}, \code{z},
#'   \code{p} (conditional test), \code{p_maxT} (FWER-adjusted within the
#'   window pair), the GLS-type null \code{se_gls} and \code{z_gls} (model on
#'   both sides), the diagnostics \code{R_AB}
#'   (\eqn{(Sy)'\Sigma_x(Sy)/tr(S\Sigma_x S\Sigma_y)};
#'   \eqn{z/z_{gls} = R_{AB}^{-1/2}}; about \eqn{1 \pm \sqrt{2/q_{eff}}}
#'   under the model) and \code{e_B} (realised over model-predicted signal
#'   energy of side B; \eqn{e_B \ll 1} means B's enrichment is overstated),
#'   \code{q_eff_class} (when a class null is used), \code{Q}
#'   (orientation-free Frobenius statistic of the tag-level matrix) and
#'   \code{p_Q} (Liu-Satterthwaite; model null).
#' @export
lrcp_gene <- function(ZA, ZB, RA, RB, VA = diag(nrow(ZA)), VB = diag(nrow(ZB)),
                      n, gcov, M, intercept = NULL, GA = NULL, GB = NULL,
                      VclassA = NULL, VclassB = NULL, n_sim = 10000) {
  ZA <- as.matrix(ZA); ZB <- as.matrix(ZB)
  q <- ncol(ZA)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  s <- n * diag(gcov) / M
  ss <- sum(s^2)
  Tm <- scale_matrix(n, gcov, M)
  RA <- as.matrix(RA); RB <- as.matrix(RB)
  if (is.null(GA)) {
    vA <- pmax(local_variance(ZA, RA, seq_len(nrow(ZA)), n, gcov, M, intercept), 0)
    GA <- RA %*% (vA * RA)
  }
  if (is.null(GB)) {
    vB <- pmax(local_variance(ZB, RB, seq_len(nrow(ZB)), n, gcov, M, intercept), 0)
    GB <- RB %*% (vB * RB)
  }
  as_list <- function(V, k) if (is.null(V)) NULL else if (is.list(V)) V else rep(list(V), k)
  stat <- function(VA, VB, VcA = NULL, VcB = NULL) {
    g <- solve(RA, VA); h <- solve(RB, VB)
    nA <- ncol(g); nB <- ncol(h)
    X <- crossprod(g, ZA)          # genesA x q
    Y <- crossprod(h, ZB)          # genesB x q
    C <- X %*% (s * t(Y)) / ss
    sY <- t(Y) * s                 # q x genesB, s_a y_aH
    MR <- crossprod(g, RA %*% g); MG <- crossprod(g, GA %*% g)
    tRB <- diag(crossprod(h, RB %*% h)); tGB <- diag(crossprod(h, GB %*% h))
    # trait covariance of each gene's projected profile, genetic part
    gen_x <- function(i) Tm * MG[i, i]
    kap <- NULL
    if (!is.null(VcA)) {
      kap <- vapply(seq_len(nA), function(i)
        max(sum(X[i, ]^2) - MR[i, i] * sum(diag(intercept)), 0) / sum(diag(VcA[[i]])), 0)
      gen_x <- function(i) kap[i] * VcA[[i]]
    }
    gen_y <- function(j) Tm * tGB[j]
    if (!is.null(VcB)) {
      kb <- vapply(seq_len(nB), function(j)
        max(sum(Y[j, ]^2) - tRB[j] * sum(diag(intercept)), 0) / sum(diag(VcB[[j]])), 0)
      gen_y <- function(j) kb[j] * VcB[[j]]
    }
    # conditional covariance of vec(C) (A fastest), eq. (8.2)
    KG <- crossprod(sY, intercept %*% sY)
    if (is.null(VcA)) {
      KT <- crossprod(sY, Tm %*% sY)
      V <- (kronecker(KG, MR) + kronecker(KT, MG)) / ss^2
    } else {
      # class null: genes' genetic parts are kappa_i V_i; cross-gene terms use
      # the LD-implied correlation of the genes with the averaged class covariance
      V <- kronecker(KG, MR)
      dg <- sqrt(pmax(diag(MG), 1e-300))
      for (i in seq_len(nA)) for (k in seq_len(nA)) {
        Vik <- if (i == k) gen_x(i) else
          MG[i, k] / (dg[i] * dg[k]) * sqrt(kap[i] * kap[k]) * (VcA[[i]] + VcA[[k]]) / 2
        Bik <- crossprod(sY, Vik %*% sY)
        V[i + (seq_len(nB) - 1) * nA, k + (seq_len(nB) - 1) * nA] <-
          V[i + (seq_len(nB) - 1) * nA, k + (seq_len(nB) - 1) * nA] + Bik
      }
      V <- V / ss^2
    }
    # GLS-type null variance and diagnostics per gene pair
    vg <- RAB <- eB <- matrix(0, nA, nB)
    for (i in seq_len(nA)) {
      Sx <- MR[i, i] * intercept + gen_x(i)
      for (j in seq_len(nB)) {
        Sy <- tRB[j] * intercept + gen_y(j)
        vg[i, j] <- sum(outer(s, s) * Sx * Sy) / ss^2
        RAB[i, j] <- sum(sY[, j] * (Sx %*% sY[, j])) / ss^2 / vg[i, j]
      }
    }
    for (j in seq_len(nB))
      eB[, j] <- (sum(s * Y[j, ]^2) - tRB[j] * sum(s * diag(intercept))) /
        (tGB[j] * sum(s * diag(Tm)))
    list(C = C, V = V, vg = vg, RAB = RAB, eB = eB)
  }
  nA <- ncol(as.matrix(VA)); nB <- ncol(as.matrix(VB))
  VcA <- as_list(VclassA, nA); VcB <- as_list(VclassB, nB)
  st <- stat(as.matrix(VA), as.matrix(VB), VcA, VcB)
  se <- matrix(sqrt(pmax(diag(st$V), 0)), nrow(st$C))
  z <- st$C / se
  p <- 2 * stats::pnorm(-abs(z))
  se_gls <- sqrt(pmax(st$vg, 0))
  qe <- NULL
  if (!is.null(VcA)) {
    qe <- vapply(VcA, function(V) { e <- pmax(eigen(V, TRUE, TRUE)$values, 0); sum(e)^2 / sum(e^2) }, 0)
    if (any(qe < 5)) warning("class effective trait count below 5: alignment within the class ",
                             "is not testable; report rho descriptively (theory S8.5)", call. = FALSE)
  }
  p_maxT <- NULL
  if (n_sim > 0) {
    sd <- sqrt(pmax(diag(st$V), 1e-300))
    Cor <- st$V / outer(sd, sd)
    L <- psd_factor(Cor)
    sims <- L %*% matrix(stats::rnorm(nrow(L) * n_sim), nrow(L))
    mx <- apply(abs(sims), 2, max)
    p_maxT <- matrix(vapply(abs(as.vector(z)), function(t) mean(mx >= t), 0), nrow(st$C))
  }
  # orientation-free Q on tag level (model null)
  tl <- stat(diag(nrow(ZA)), diag(nrow(ZB)))
  Q <- sum(tl$C^2)
  ev <- eigen(tl$V, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > 0]
  # Liu-Satterthwaite: Q ~ a chi2_d, a = sum ev^2 / sum ev, d = (sum ev)^2 / sum ev^2
  a <- sum(ev^2) / sum(ev); d <- sum(ev)^2 / sum(ev^2)
  list(C = st$C, se = se, z = z, p = p, p_maxT = p_maxT,
       se_gls = se_gls, z_gls = st$C / se_gls, R_AB = st$RAB, e_B = st$eB,
       q_eff_class = qe, Q = Q, p_Q = stats::pchisq(Q / a, d, lower.tail = FALSE))
}

#' Cluster loci by their trait profiles
#'
#' Hierarchical clustering (average linkage) of the sign-free correlations
#' \eqn{|corr(x_t, x_{t'})|} between projected trait profiles, used to define
#' locus classes for the class null of [lrcp_gene()] (theory supplement
#' S8.5).
#'
#' @param X Loci-by-q matrix of projected Z-score profiles
#'   (\eqn{R_{TT}^{-1} Z_T} rows for enriched tags, genome-wide).
#' @param min_abs_cor Profiles are joined while their average absolute
#'   correlation exceeds this value.
#' @return Integer class labels.
#' @export
profile_classes <- function(X, min_abs_cor = 0.5) {
  X <- as.matrix(X)
  if (nrow(X) < 2) return(rep(1L, nrow(X)))
  d <- stats::as.dist(1 - abs(stats::cor(t(X))))
  stats::cutree(stats::hclust(d, method = "average"), h = 1 - min_abs_cor)
}

#' Class trait covariance from held-out loci
#'
#' \eqn{\hat V_{class} = avg_t (x_t x_t' - \tau_{R,t}\Gamma)/\tau_{G,t}} over
#' loci of one profile class, excluding the loci under test and their LD
#' blocks (theory supplement S8.5).
#'
#' @param X Loci-by-q projected profiles.
#' @param tau_R,tau_G Per-locus noise and signal scales: for tag \eqn{t} of a
#'   tag set \eqn{T}, \eqn{\tau_R = (R_{TT}^{-1})_{tt}} and \eqn{\tau_G} its
#'   R-projected variance ([local_variance()]).
#' @param intercept q-by-q intercept matrix \eqn{\Gamma}.
#' @param use Logical or index vector of the held-out loci to average.
#' @return List with \code{V} (projected to the PSD cone), \code{q_eff}
#'   (\eqn{(tr V)^2 / tr V^2}) and \code{n_loci}.
#' @export
class_covariance <- function(X, tau_R, tau_G, intercept, use = seq_len(nrow(X))) {
  X <- as.matrix(X)
  if (is.logical(use)) use <- which(use)
  use <- use[tau_G[use] > 0]
  if (!length(use)) return(list(V = NULL, q_eff = NA_real_, n_loci = 0L))
  V <- Reduce(`+`, lapply(use, function(t)
    (tcrossprod(X[t, ]) - tau_R[t] * intercept) / tau_G[t])) / length(use)
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  ev <- pmax(e$values, 0)
  V <- e$vectors %*% (ev * t(e$vectors))
  list(V = V, q_eff = sum(ev)^2 / sum(ev^2), n_loci = length(use))
}
