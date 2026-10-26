test_that("decompose_uq returns calibrated intervals on clean input", {
  set.seed(11)
  n_g <- 400; k <- 4
  mu <- matrix(runif(n_g * k, 1, 50), n_g, k,
               dimnames = list(paste0("g", 1:n_g), c("A","B","C","D")))
  mu <- sweep(mu, 2, colSums(mu), "/") * 1e4
  mu_se <- mu * 0.05
  f1 <- c(A=.5, B=.2, C=.2, D=.1); f2 <- c(A=.1, B=.4, C=.3, D=.2)
  obs <- log2(as.vector(f2 %*% t(mu)) / as.vector(f1 %*% t(mu)))
  names(obs) <- rownames(mu)

  out <- decompose_uq(obs, f1, f2, mu, mu_se,
                      f_se_ctrl = rep(.005, 4), f_se_case = rep(.005, 4),
                      n_boot = 300, seed = 2)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), n_g)
  expect_true(all(c("gene","obs","comp_med","resid_med","resid_lo",
                    "resid_hi","p_same_sign","direction","call") %in% names(out)))
  ## With no within-cell change the residual should be near zero
  expect_lt(median(abs(out$resid_med)), 0.15)
  ## Intervals must be ordered
  expect_true(all(out$resid_lo <= out$resid_hi))
  ## Calls must be one of the three allowed values
  expect_true(all(out$call %in% c("interpretable","over_predicted","undetermined")))
})

test_that("decompose_uq recovers a known within-cell change", {
  set.seed(12)
  n_g <- 300; k <- 3
  mu <- matrix(runif(n_g * k, 1, 50), n_g, k,
               dimnames = list(paste0("g", 1:n_g), c("A","B","C")))
  mu <- sweep(mu, 2, colSums(mu), "/") * 1e4
  mu_se <- mu * 0.02
  f1 <- c(A=.6, B=.3, C=.1)
  f2 <- f1                                  # no composition change
  delta <- numeric(n_g); delta[1:30] <- 1.5
  mu_case <- mu * (2^delta)
  for (i in 1:k) mu_case[, i] <- mu_case[, i] / sum(mu_case[, i]) * sum(mu[, i])
  obs <- log2(as.vector(f2 %*% t(mu_case)) / as.vector(f1 %*% t(mu)))
  names(obs) <- rownames(mu)

  out <- decompose_uq(obs, f1, f2, mu_case * 0 + mu, mu_se,
                      f_se_ctrl = rep(.002, 3), f_se_case = rep(.002, 3),
                      n_boot = 300, seed = 3)
  ## Residuals for the changed genes should sit near their true value
  expect_lt(median(abs(out$resid_med[1:30] - delta[1:30])), 0.5)
})

test_that("decompose_uq validates its input", {
  mu <- matrix(1:12, 4, 3, dimnames = list(paste0("g", 1:4), c("A","B","C")))
  obs <- stats::setNames(1:4, paste0("g", 1:4))
  ## no shared cell-type names at all
  expect_error(decompose_uq(obs, c(X=.5,Y=.5), c(X=.2,Y=.8), mu, mu),
               "cannot determine orientation")
  ## a single shared cell type is also insufficient
  expect_error(decompose_uq(obs, c(A=1), c(A=1), mu, mu),
               "cannot determine orientation")
  ## missing gene names
  mu2 <- mu; rownames(mu2) <- NULL
  expect_error(decompose_uq(obs, c(A=.5,B=.2,C=.3), c(A=.2,B=.3,C=.5), mu2, mu2),
               "gene names")
  ## mismatched shapes
  expect_error(decompose_uq(obs, c(A=.5,B=.2,C=.3), c(A=.2,B=.3,C=.5),
                            mu, mu[1:2, ]),
               "same shape")
})

test_that("null_calibrate pools rather than averages", {
  set.seed(21)
  n_g <- 200
  ## Complementary pairs: column 1 and 2 are exact negatives
  a <- rnorm(n_g, 0, .3)
  null_mat <- cbind(a, -a)
  rownames(null_mat) <- paste0("g", 1:n_g)
  observed <- setNames(rnorm(n_g, 0, .1), paste0("g", 1:n_g))
  out <- null_calibrate(null_mat, observed)
  ## Pooling keeps the spread; averaging would collapse it to zero
  expect_gt(out$null_quantiles[["95%"]], 0.3)
  expect_equal(length(out$null_pooled), 2 * n_g)
})

test_that("null_calibrate flags large residuals", {
  set.seed(22)
  n_g <- 1000
  null_mat <- matrix(rnorm(n_g * 30, 0, .2), n_g, 30,
                     dimnames = list(paste0("g", 1:n_g), NULL))
  observed <- setNames(rnorm(n_g, 0, .2), paste0("g", 1:n_g))
  observed[1:20] <- 3
  out <- null_calibrate(null_mat, observed)
  expect_true(all(out$genes$exceeds_99[1:20]))
  expect_lt(sum(out$genes$exceeds_95[21:n_g]) / (n_g - 20), 0.10)
  expect_true(all(out$genes$pct >= 0 & out$genes$pct <= 1))
})

test_that("null_calibrate validates its input", {
  ## no gene names on the null matrix
  m <- matrix(stats::rnorm(100), 10, 10)
  expect_error(null_calibrate(m, stats::rnorm(10)), "rownames")
  ## too few null residuals to calibrate
  m2 <- matrix(stats::rnorm(20), 10, 2,
               dimnames = list(paste0("g", 1:10), NULL))
  observed <- stats::setNames(stats::rnorm(10), paste0("g", 1:10))
  expect_error(null_calibrate(m2, observed), "too few null residuals")
})

test_that("scale_null follows the variance formula", {
  expect_equal(scale_null(4, 4, 4), 1)
  expect_lt(scale_null(76, 13, 4), 0.5)
  expect_gt(scale_null(2, 2, 4), 1)
  expect_error(scale_null(1, 5), "at least 2")
})

test_that("reference_qc flags unreliable cell types", {
  set.seed(31)
  mu <- matrix(runif(200 * 3, 1, 100), 200, 3,
               dimnames = list(NULL, c("good","ok","poor")))
  se <- cbind(good = mu[,1] * 0.02, ok = mu[,2] * 0.15, poor = mu[,3] * 0.8)
  qc <- reference_qc(mu, se, n_cells = c(good = 5000, ok = 1200, poor = 60))
  expect_equal(nrow(qc), 3)
  expect_false(qc$flagged[qc$celltype == "good"])
  expect_true(qc$flagged[qc$celltype == "poor"])
  expect_equal(qc$n_cells[qc$celltype == "poor"], 60)
})
