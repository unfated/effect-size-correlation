sim_two_regions <- function(q = 30, rho = c(0.6, -0.4, 0), gcov = NULL, C = NULL) {
  blocks <- make_ld_ar1(c(30, 30), 0.5); m <- 60
  w <- numeric(m); idx1 <- c(5, 12, 20); idx2 <- 30 + c(8, 15, 25); w[c(idx1, idx2)] <- 10
  Rb <- as.matrix(make_Rb(m, "pairs", pairs = cbind(idx1, idx2), pair_rho = rho))
  if (is.null(gcov)) gcov <- make_gcov(q, h2 = rep(0.3, q))$gcov
  if (is.null(C)) C <- diag(q)
  list(blocks = blocks, w = w, Rb = Rb, gcov = gcov, C = C, n = rep(3e5, q), M = 2e4)
}

test_that("distal LRCP regression is unbiased, calibrated, and halves under the chapter form", {
  set.seed(11)
  s <- sim_two_regions(q = 30)
  est <- se <- chap <- NULL
  for (r in 1:120) {
    sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, Rb = s$Rb, gcov = s$gcov, M = s$M)
    a <- list(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]], s$w[1:30], s$w[31:60], s$n, s$gcov, s$M)
    f <- do.call(lrcp_distal, a)
    g <- do.call(lrcp_distal, c(a, list(pair_sum = "chapter")))
    est <- rbind(est, diag(f$rho)); se <- rbind(se, diag(f$se)); chap <- rbind(chap, diag(g$rho))
  }
  sdv <- apply(est, 2, sd)
  expect_true(all(abs(colMeans(est) - c(0.6, -0.4, 0)) < 4 * sdv / sqrt(120)))
  expect_equal(colMeans(chap), colMeans(est) / 2, tolerance = 1e-8)
  expect_true(all(abs(colMeans(se) / sdv - 1) < 0.25))
})

test_that("LRCP SEs account for correlated traits and sample overlap", {
  set.seed(12)
  q <- 30
  gc <- make_gcov(q, h2 = rep(0.3, q), type = "cluster", rg = 0.7, n_clusters = 2)
  s <- sim_two_regions(q, rho = c(0, 0, 0), gcov = gc$gcov, C = make_intercept(q, 1, 0.4))
  zz <- NULL
  for (r in 1:120) {
    sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, gcov = s$gcov, intercept = s$C, M = s$M)
    f <- lrcp_distal(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]],
                     s$w[1:30], s$w[31:60], s$n, s$gcov, s$M, s$C)
    zz <- c(zz, as.vector(f$z))
  }
  expect_lt(abs(sd(zz) - 1), 0.15)
})

test_that("LRCP MLE recovers rho and the gradient is exact", {
  set.seed(13)
  s <- sim_two_regions(q = 30)
  est <- NULL
  for (r in 1:15) {
    sim <- simulate_lrcpq(s$blocks, n = s$n, w = s$w, Rb = s$Rb, gcov = s$gcov, M = s$M)
    f <- lrcp_mle(sim$Z[1:30, ], sim$Z[31:60, ], s$blocks[[1]], s$blocks[[2]],
                  s$w[1:30], s$w[31:60], s$n, s$gcov, s$M)
    expect_equal(f$convergence, 0)
    expect_gt(f$lrt, -1e-6)
    est <- rbind(est, diag(f$rho))
  }
  expect_true(all(abs(colMeans(est) - c(0.6, -0.4, 0)) < 0.12))
})

test_that("trait covariance helper matches the moment equation", {
  R <- make_ld_ar1(5, 0.5)[[1]]; w <- c(2, 0, 1, 1, 1)
  g <- make_gcov(3, h2 = c(0.2, 0.3, 0.4), type = "cs", rg = 0.5)$gcov
  n <- c(1e4, 2e4, 3e4); C <- make_intercept(3, 0.5, 0.2)
  Mk <- expected_trait_cov(2, R, w, NULL, n, g, 100, C)
  for (a in 1:3) for (b in 1:3)
    expect_equal(Mk[a, b], expected_zz(R, w, NULL, n, g, 100, C, a, b)[2, 2])
  expect_equal(q_eff(diag(3)), 3)
})
