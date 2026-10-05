test_that("quadform_var matches Monte Carlo", {
  set.seed(51)
  R <- make_ld_ar1(8, 0.5)[[1]]; w <- c(5, 0, 0, 3, 0, 0, 0, 2) * 8 / 10
  q <- 6; gc <- make_gcov(q, h2 = rep(0.3, q), type = "cs", rg = 0.5)
  C <- make_intercept(q, 1, 0.3); n <- rep(1e4, q); M <- 100
  A <- diag(runif(8)); u <- runif(q)
  Tm <- scale_matrix(n, gc$gcov, M); G <- R %*% (w * R)
  v <- quadform_var(A, u, R, G, C, Tm)
  x <- replicate(4000, { Z <- simulate_lrcpq(R, n = n, w = w, gcov = gc$gcov, intercept = C, M = M)$Z
    sum(u * colSums(Z * (A %*% Z))) })
  expect_lt(abs(var(x) / v - 1), 0.1)
})

test_that("category enrichment and annotation regression are calibrated", {
  set.seed(52)
  blocks <- make_ld_ar1(rep(30, 8), 0.6); ld <- ld_from_blocks(blocks); m <- ld$m
  annot <- rbinom(m, 1, 0.3)
  w <- rep(1, m); w[annot == 1] <- 2.5; w <- w * m / sum(w)   # dense, no clumping effects
  q <- 20; gc <- make_gcov(q, h2 = rep(0.3, q), type = "cluster", rg = 0.6, n_clusters = 4)
  C <- make_intercept(q, 1, 0.3); n <- rep(3e5, q); M <- 1e4
  r <- replicate(80, {
    Z <- simulate_lrcpq(blocks, n = n, w = w, gcov = gc$gcov, intercept = C, M = M)$Z
    ar <- lrcq_annot_regression(Z, ld, annot, n, gc$gcov, M, C)
    c(ar$tau, ar$se)
  })
  sdv <- apply(r[1:2, ], 1, sd)
  expect_true(all(abs(rowMeans(r[3:4, ]) / sdv - 1) < 0.3))
  # tau estimand: w = tau0 + tau1 * annot exactly here
  truth <- c(min(w), max(w) - min(w))
  expect_true(all(abs(rowMeans(r[1:2, ]) - truth) < 4 * sdv / sqrt(80)))
})

test_that("gene-level LRCP null is calibrated under correlated, overlapping traits", {
  set.seed(53)
  blocks <- make_ld_ar1(c(30, 30), 0.5); m <- 60
  TA <- c(3, 6, 10, 14); TB <- 30 + TA; w <- numeric(m); w[c(TA, TB)] <- 10
  q <- 30; gc <- make_gcov(q, h2 = rep(0.3, q), type = "cluster", rg = 0.6, n_clusters = 3)
  C <- make_intercept(q, 1, 0.3); n <- rep(3e5, q); M <- 2e4
  V <- kronecker(diag(2), matrix(1, 2, 1))
  z <- replicate(150, {
    Z <- simulate_lrcpq(blocks, n = n, w = w, gcov = gc$gcov, intercept = C, M = M)$Z
    R1 <- blocks[[1]]
    lrcp_gene(Z[TA, ], Z[TB, ], R1[TA, TA], R1[TA, TA], V, V, n, gc$gcov, M, C,
              GA = (R1 %*% (w[1:30] * R1))[TA, TA], n_sim = 0)$z
  })
  expect_lt(abs(sd(as.vector(z)) - 1), 0.12)
})

test_that("power helpers reproduce the S9 tables", {
  expect_equal(round(se_lrcq(0, 4e5, 0.1, qe = 100), 1), 5.0)
  expect_equal(round(power_lrcq(c(10, 20, 50), 4e5, 0.1, qe = 100), 2), c(0.01, 0.29, 0.97))
  expect_equal(round(se_lrcp(20, 20, 0, 4e5, 0.1, qe = 300) / 20, 2), 0.14)
  qe <- effective_traits(diag(3), diag(3))
  expect_equal(unname(qe), c(3, 3, 3))
})

test_that("category contrast: debiased plug-in and bootstrap are calibrated", {
  set.seed(54)
  blocks <- make_ld_ar1(rep(30, 8), 0.6); ld <- ld_from_blocks(blocks); m <- ld$m
  w <- rep(0.02, m); loc <- c(5, 40, 75, 130, 170, 220); w[loc] <- runif(6, 40, 80)
  w <- w * m / sum(w)
  annot <- rbinom(m, 1, 0.3); annot[loc[1:3]] <- 1
  q <- 40; gc <- make_gcov(q, h2 = rep(0.3, q), type = "cluster", rg = 0.6, n_clusters = 4)
  C <- make_intercept(q, 1, 0.3); n <- rep(3e5, q); M <- 2e4
  u <- trait_weights(n, gc$gcov, M, "contrast", rep(c(1, 1, 2, 2), each = 10))
  r <- t(replicate(100, {
    Z <- simulate_lrcpq(blocks, n = n, w = w, gcov = gc$gcov, intercept = C, M = M)$Z
    b <- lrcq_category(Z, ld, annot, n, gc$gcov, M, C, u = u, n_boot = 40)
    c(b$z, b$p)
  }))
  expect_lt(abs(sd(r[, 1]) - 1), 0.2)
  expect_lt(abs(mean(r[, 2] <= 0.1) - 0.1), 0.08)
})
