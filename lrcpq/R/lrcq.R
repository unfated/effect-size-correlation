# Stage 3: LRCQ estimation of the genetic-effect enrichment w.
#
# Regression model (chapter 6, LRCQ key equation, R_beta = 0):
#   y_ka = z_ka^2 - C_aa = s_a (D w)_k + e_ka,   s_a = n_a h2_a / M,  D = R o R.
# Stacking traits gives the design X = s (x) D, so (fast OLS/WLS derivation)
#   X' Omega X = D' diag(omega) D,   omega_k = sum_a s_a^2 / v_ka
#   X' Omega y = D' u,               u_k     = sum_a s_a y_ka / v_ka
# and the q-trait regression collapses to one m-by-m solve.

#' Rectify negative enrichment estimates
#'
#' Chapter 6.3.2 rectification-adjustment methods, applied only when the
#' vector contains negative values.
#' \describe{
#'   \item{A}{zero every element \eqn{\le |\min(\hat w)|}.}
#'   \item{B}{zero all negatives plus the same number of smallest positives.}
#'   \item{C}{zero all elements up to the cutoff \eqn{c} at which the sum of
#'     positives in \eqn{(0, c]} is closest to the absolute sum of negatives
#'     (ties to the lower cutoff).}
#' }
#' @param w Estimated enrichment.
#' @param method One of \code{"none"}, \code{"A"}, \code{"B"}, \code{"C"},
#'   \code{"zero"} (just truncate negatives at 0).
#' @export
rectify_w <- function(w, method = c("none", "A", "B", "C", "zero")) {
  method <- match.arg(method)
  if (method == "none" || !any(w < 0, na.rm = TRUE)) return(w)
  neg <- which(w < 0)
  out <- w
  if (method == "zero") {
    out[neg] <- 0
  } else if (method == "A") {
    out[w <= abs(min(w, na.rm = TRUE))] <- 0
  } else if (method == "B") {
    pos <- which(w > 0)
    k <- min(length(neg), length(pos))
    small <- pos[order(w[pos])][seq_len(k)]
    out[c(neg, small)] <- 0
  } else if (method == "C") {
    total_neg <- -sum(w[neg])
    pos <- which(w > 0)
    sp <- sort(w[pos])
    cs <- cumsum(sp)
    # cutoff whose cumulative positive mass is closest to the negative mass
    # (ties go to the lower cutoff)
    j <- which.min(abs(cs - total_neg))
    cutoff <- sp[j]
    out[w <= cutoff] <- 0
  }
  out
}

#' LRCQ enrichment estimate within one LD window
#'
#' @param Z m-by-q Z-scores of the window SNPs.
#' @param R m-by-m LD matrix of the window.
#' @param n Length-q sample sizes.
#' @param gcov q-by-q genetic covariance (or length-q heritabilities).
#' @param M Number of SNPs heritability refers to (genome-wide).
#' @param intercept q-by-q intercept matrix (or length-q univariate
#'   intercepts). Responses are \eqn{z^2 - C_{aa}}.
#' @param method \code{"ols"}, \code{"wls"} (one reweighting step from OLS),
#'   \code{"irls"} (iterate weights to convergence), or \code{"equal"}.
#'   WLS weights are \eqn{1/\mathrm{Var}(z^2) = 1/(2\sigma^4)}. OLS/WLS
#'   estimate the power-weighted mean enrichment
#'   \eqn{\sum_a s_a^2 w_a / \sum_a s_a^2} when enrichment differs by trait;
#'   \code{"equal"} gives \eqn{\hat w = q^{-1}\sum_a D^{-1}(z_a^2 - C_{aa})/s_a},
#'   unbiased for the simple mean over traits but noisier.
#' @param ridge Ridge penalty as a fraction of the mean diagonal of
#'   \eqn{D'\Omega D}; 0 gives the exact fast OLS/WLS solution.
#' @param N_ref Reference panel size for the r-squared bias correction.
#' @param Rb Optional genetic-effect correlation for the window. When given,
#'   the chapter's cross terms are included as an offset evaluated at the
#'   current \eqn{\hat w} (\code{cross_terms = TRUE} in [lrcq()]).
#' @param pair_sum See [sigma_w()].
#' @param se \code{"model"}: sandwich SE under Gaussian Z with the full
#'   Isserlis covariance across SNPs and traits; \code{"none"}.
#' @param max_iter,tol IRLS controls.
#' @param w_plugin Enrichment used to evaluate the SE formula; default is
#'   \code{w_hat} rectified by method A, which keeps the SE calibrated
#'   (raw plug-in overstates SEs by roughly 20\% when most SNPs are null).
#' @return List with \code{w}, \code{se}, \code{iter}.
#' @export
lrcq_window <- function(Z, R, n, gcov, M, intercept = NULL,
                        method = c("wls", "ols", "irls", "equal"), ridge = 0,
                        N_ref = Inf, Rb = NULL, pair_sum = c("once", "chapter"),
                        se = c("model", "none"), max_iter = 20, tol = 1e-6,
                        w_plugin = NULL) {
  method <- match.arg(method)
  se <- match.arg(se)
  pair_sum <- match.arg(pair_sum)
  Z <- as.matrix(Z)
  m <- nrow(Z)
  q <- ncol(Z)
  if (is.null(dim(gcov))) gcov <- diag(gcov, q)
  if (is.null(intercept)) intercept <- diag(q)
  if (is.null(dim(intercept))) intercept <- diag(intercept, q)
  h2 <- diag(gcov)
  s <- n * h2 / M
  Cd <- diag(intercept)
  R <- as.matrix(R)
  D <- ld_r2(R, N_ref)
  Y <- Z^2 - matrix(Cd, m, q, byrow = TRUE)

  offset_fn <- function(w) {
    if (is.null(Rb)) return(rep(0, m))
    S <- sigma_w(pmax(w, 0), Rb, pair_sum)
    diag(S) <- 0
    rowSums((R %*% S) * R)
  }
  solve_w <- function(V, off) {
    Omega <- matrix(s^2, m, q, byrow = TRUE) / V
    omega <- rowSums(Omega)
    u <- rowSums((Y - outer(off, s)) * matrix(s, m, q, byrow = TRUE) / V)
    A <- crossprod(D, D * omega)
    if (ridge > 0) A <- A + diag(ridge * mean(diag(A)), m)
    list(w = as.vector(solve(A, crossprod(D, u))), A = A)
  }
  var_fn <- function(w, off) {
    mu <- outer(pmax(as.vector(D %*% w) + off, 0), s) +
      matrix(Cd, m, q, byrow = TRUE)
    2 * mu^2
  }

  V0 <- if (method == "equal") matrix(s^2, m, q, byrow = TRUE) else matrix(1, m, q)
  off <- rep(0, m)
  fit <- solve_w(V0, off)
  iter <- 1
  if (!is.null(Rb)) {
    off <- offset_fn(fit$w)
    fit <- solve_w(V0, off)
  }
  if (method %in% c("wls", "irls")) {
    n_it <- if (method == "wls") 1 else max_iter
    for (it in seq_len(n_it)) {
      w_old <- fit$w
      V <- var_fn(fit$w, off)
      fit <- solve_w(V, off)
      if (!is.null(Rb)) off <- offset_fn(fit$w)
      iter <- iter + 1
      if (max(abs(fit$w - w_old)) < tol * max(1, max(abs(w_old)))) break
    }
  }
  w_hat <- fit$w
  se_w <- rep(NA_real_, m)
  if (se == "model") {
    V <- if (method %in% c("ols", "equal")) V0 else var_fn(w_hat, off)
    Fm <- matrix(s, m, q, byrow = TRUE) / V
    cmat <- scale_matrix(n, gcov, M)
    # Plugging raw w_hat into the variance inflates it (noise in null SNPs
    # is clipped at zero), so the default plug-in is the rectified estimate.
    w_se <- if (is.null(w_plugin)) rectify_w(w_hat, "A") else w_plugin
    Ak <- R %*% sigma_w(pmax(w_se, 0), Rb, pair_sum) %*% R
    Vu <- 2 * (Ak^2 * (Fm %*% (cmat^2) %*% t(Fm)) +
                 2 * Ak * R * (Fm %*% (cmat * intercept) %*% t(Fm)) +
                 R^2 * (Fm %*% (intercept^2) %*% t(Fm)))
    Ainv <- solve(fit$A)
    G <- Ainv %*% t(D)
    se_w <- sqrt(pmax(rowSums((G %*% Vu) * G), 0))
  }
  list(w = w_hat, se = se_w, iter = iter)
}

#' Genome-wide LRCQ: genetic-effect enrichment for every SNP
#'
#' Intercepts \eqn{C_{aa}} must come from stage 2 (genome-wide LDSC) and are
#' held fixed: with free per-trait intercepts inside a window the design
#' \eqn{[s \otimes D, I_q \otimes 1]} is one rank short and \eqn{w} is not
#' identifiable. From squared Z-scores alone only \eqn{diag(R S_w R)} is
#' identified, so \eqn{\hat w} estimates \eqn{w} under the working assumption
#' that local genetic-effect correlation is zero; otherwise its bias is
#' \eqn{D^{-1} diag(R\, C\, R)} with \eqn{C} the off-diagonal part of
#' \eqn{S_w} (theory supplement S1.5, S1.6).
#'
#' Runs [lrcq_window()] over sliding windows of LD blocks (chapter 6 sliding
#' window strategy) and keeps the core-block estimates.
#'
#' @param Z m-by-q Z-score matrix aligned to the LD reference SNPs.
#' @param ld An \code{lrcpq_ld} object.
#' @param window_size Blocks per window (odd).
#' @param rectify Rectification method for negative estimates, applied
#'   genome-wide after estimation (see [rectify_w()]).
#' @param cross_terms Include the chapter's \eqn{r_\beta} cross terms using
#'   \code{Rb} (default \code{FALSE}: chi-square-only LRCQ as in the
#'   abstract and Table 6.1).
#' @param Rb m-by-m genetic-effect correlation (needed if
#'   \code{cross_terms = TRUE}).
#' @param chr Optional chromosome per SNP (windows do not cross chromosomes).
#' @param ... Passed to [lrcq_window()].
#' @inheritParams lrcq_window
#' @return Data frame with \code{w_raw}, \code{w} (rectified), \code{se},
#'   \code{z} (\code{w_raw / se}) and \code{block}.
#' @export
lrcq <- function(Z, ld, n, gcov, M, intercept = NULL, window_size = 3,
                 method = c("wls", "ols", "irls", "equal"), rectify = "none",
                 cross_terms = FALSE, Rb = NULL, chr = NULL, ...) {
  method <- match.arg(method)
  Z <- as.matrix(Z)
  stopifnot(nrow(Z) == ld$m)
  m <- nrow(Z)
  w <- se <- rep(NA_real_, m)
  for (win in make_windows(ld$block_id, window_size, chr)) {
    Rw <- ld$get(win$idx)
    Rbw <- if (cross_terms) as.matrix(Rb[win$idx, win$idx]) else NULL
    f <- lrcq_window(Z[win$idx, , drop = FALSE], Rw, n, gcov, M, intercept,
                     method = method, N_ref = ld$N, Rb = Rbw, ...)
    core <- win$idx[win$core]
    w[core] <- f$w[win$core]
    se[core] <- f$se[win$core]
  }
  data.frame(snp = seq_len(m), block = ld$block_id, w_raw = w,
             w = rectify_w(w, rectify), se = se, z = w / se)
}

#' Brute-force stacked regression for checking the fast LRCQ solution
#'
#' Builds the \eqn{mq \times m} design \eqn{X = s \otimes D} explicitly and
#' calls \code{lm.wfit}. Only for small problems and tests.
#' @inheritParams lrcq_window
#' @param V Optional m-by-q matrix of response variances (WLS weights 1/V).
#' @export
lrcq_bruteforce <- function(Z, R, n, gcov, M, intercept = NULL, V = NULL) {
  Z <- as.matrix(Z)
  m <- nrow(Z)
  q <- ncol(Z)
  h2 <- if (is.null(dim(gcov))) gcov else diag(gcov)
  Cd <- if (is.null(intercept)) rep(1, q) else if (is.null(dim(intercept))) intercept else diag(intercept)
  s <- n * h2 / M
  D <- as.matrix(R)^2
  X <- kronecker(matrix(s, q, 1), D)
  y <- as.vector(Z^2 - matrix(Cd, m, q, byrow = TRUE))
  wt <- if (is.null(V)) rep(1, m * q) else 1 / as.vector(V)
  stats::lm.wfit(X, y, wt)$coefficients
}
