test_that("group composition is the MEAN, not the median or the first sample", {
  ## Reviewer mutation test: switching group composition from the mean to the
  ## median, or to the first sample of the group, left every existing test
  ## green, because every fixture had identical rows within a group.
  ## Here the rows differ, so the three statistics give different answers.
  ##
  ## UPDATE (v6.0.0): compose_prediction() gained a `convention` argument and
  ## the default changed to "mean_of_log". The mutation this test guards
  ## against -- replacing the MEAN by the MEDIAN or by the FIRST ROW -- applies
  ## to BOTH conventions, so both are now tested explicitly.
  mu <- matrix(c(10, 1,
                 1, 10), nrow = 2, byrow = TRUE,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  ## four control samples with DIFFERENT composition
  f_ctrl <- rbind(c(0.9, 0.1), c(0.7, 0.3), c(0.5, 0.5), c(0.1, 0.9))
  ## two case samples, also different
  f_case <- rbind(c(0.2, 0.8), c(0.4, 0.6))
  colnames(f_ctrl) <- colnames(f_case) <- c("A", "B")
  frac <- rbind(f_ctrl, f_case)
  is_case <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)

  ## ---- log_of_mean: mean the PROPORTIONS first, then take the log ratio ----
  got_lom <- compose_prediction(frac, mu, is_case, convention = "log_of_mean")
  mean_c <- colMeans(f_ctrl); mean_a <- colMeans(f_case)
  want_lom <- log2(as.vector(mean_a %*% mu) / as.vector(mean_c %*% mu))
  expect_equal(as.numeric(got_lom), as.numeric(want_lom), tolerance = 1e-12)

  ## ---- mean_of_log: log each SAMPLE first, then average within group -------
  got_mol <- compose_prediction(frac, mu, is_case, convention = "mean_of_log")
  want_mol <- as.numeric(colMeans(log2(f_case %*% mu)) -
                         colMeans(log2(f_ctrl %*% mu)))
  expect_equal(as.numeric(got_mol), want_mol, tolerance = 1e-12)

  ## ---- the two conventions must NOT coincide on this fixture --------------
  ## They differ here because the within-group rows differ. If they agreed,
  ## the `convention` argument would be doing nothing.
  expect_false(isTRUE(all.equal(as.numeric(got_lom), as.numeric(got_mol))))

  ## ---- the mutation guard: median and first row must differ from the mean --
  for (conv in c("log_of_mean", "mean_of_log")) {
    got <- compose_prediction(frac, mu, is_case, convention = conv)
    med_c <- apply(f_ctrl, 2, stats::median); med_a <- apply(f_case, 2, stats::median)
    if (conv == "log_of_mean") {
      by_median <- log2(as.vector(med_a %*% mu) / as.vector(med_c %*% mu))
    } else {
      by_median <- as.numeric(log2(med_a %*% mu) - log2(med_c %*% mu))
    }
    expect_false(isTRUE(all.equal(as.numeric(got), as.numeric(by_median))),
                 info = paste("median must differ under", conv))

    first_c <- f_ctrl[1, ]; first_a <- f_case[1, ]
    if (conv == "log_of_mean") {
      by_first <- log2(as.vector(first_a %*% mu) / as.vector(first_c %*% mu))
    } else {
      by_first <- as.numeric(log2(first_a %*% mu) - log2(first_c %*% mu))
    }
    expect_false(isTRUE(all.equal(as.numeric(got), as.numeric(by_first))),
                 info = paste("first row must differ under", conv))
  }
})

test_that("the default convention is mean_of_log and it is recorded", {
  mu <- matrix(c(10, 1, 1, 10), nrow = 2, byrow = TRUE,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  frac <- matrix(c(0.9, 0.1, 0.7, 0.3, 0.2, 0.8, 0.4, 0.6), nrow = 4,
                 byrow = TRUE, dimnames = list(NULL, c("A", "B")))
  ic <- c(FALSE, FALSE, TRUE, TRUE)

  d_def <- compose_prediction(frac, mu, ic)
  d_mol <- compose_prediction(frac, mu, ic, convention = "mean_of_log")
  expect_equal(as.numeric(d_def), as.numeric(d_mol), tolerance = 1e-15)
  expect_identical(attr(d_def, "convention"), "mean_of_log")

  ## log_of_mean must be reachable and must differ on this fixture
  d_lom <- compose_prediction(frac, mu, ic, convention = "log_of_mean")
  expect_identical(attr(d_lom, "convention"), "log_of_mean")
  expect_false(isTRUE(all.equal(as.numeric(d_def), as.numeric(d_lom))))

  ## an unknown convention must be rejected, not silently defaulted
  expect_error(compose_prediction(frac, mu, ic, convention = "nonsense"))
})

test_that("pct is resid / obs, not obs / resid", {
  ## Reviewer mutation: pct = obs / resid passed every test in the suite.
  mu <- matrix(c(10, 1, 5, 2, 8, 3), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2", "g3")))
  frac <- rbind(c(0.9, 0.1), c(0.1, 0.9), c(0.3, 0.7))
  colnames(frac) <- c("A", "B")
  obs <- c(g1 = 2, g2 = -1, g3 = 0.5)
  d <- decompose_composition(obs, frac, mu, c(FALSE, TRUE, TRUE),
                             validate = FALSE)
  want <- d$resid / d$obs
  expect_equal(d$pct, want, tolerance = 1e-12)
  ## and it must NOT be the reciprocal
  expect_false(isTRUE(all.equal(d$pct, d$obs / d$resid)))
})

test_that("mean_expr is the composition-weighted mean, not a plain mean", {
  ## Reviewer mutation: an unweighted mean of mu passed every test.
  mu <- matrix(c(100, 1, 1), nrow = 3,
               dimnames = list(c("A", "B", "C"), "g1"))
  ## strongly unequal composition, so weighted != unweighted
  frac <- rbind(c(0.98, 0.01, 0.01), c(0.98, 0.01, 0.01),
                c(0.01, 0.98, 0.01), c(0.01, 0.98, 0.01))
  colnames(frac) <- c("A", "B", "C")
  obs <- c(g1 = 1)
  d <- decompose_composition(obs, frac, mu, c(FALSE, FALSE, TRUE, TRUE),
                             validate = FALSE)
  fc <- (0.98 + 0.01) / 2; fa <- (0.01 + 0.98) / 2
  fb <- fc
  pooled <- c(A = (fc + fa) / 2, B = fb, C = fb)
  ## recompute exactly as the function documents
  fm <- colMeans(frac)
  want <- as.numeric(sum(colMeans(frac) * mu[, 1]))
  expect_equal(d$mean_expr, want, tolerance = 1e-12)
  ## unweighted mean would be 34; the weighted value must not be 34
  expect_false(abs(d$mean_expr - mean(mu[, 1])) < 1e-9)
})

test_that("the zero-denominator guard returns NA rather than Inf", {
  ## Reviewer mutation: removing the `num <= 0 | den <= 0` guard passed every
  ## existing test, because the only fixture produced 0/0.
  ## row 1 = A = (0, 0), row 2 = B = (5, 5)
  mu <- matrix(c(0, 0, 5, 5), nrow = 2, byrow = TRUE,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  ## cell type A has zero expression for both genes
  frac <- rbind(c(1, 0), c(1, 0), c(0, 1), c(0, 1))
  colnames(frac) <- c("A", "B")
  ## case group has only type B, control group only type A: denominator is 0
  is_case <- c(FALSE, FALSE, TRUE, TRUE)
  out <- compose_prediction(frac, mu, is_case)
  expect_true(all(is.na(out)))
  expect_false(any(is.infinite(out)))
})

test_that("decompose_composition rejects duplicated frac cell types", {
  mu <- matrix(c(10, 1, 5, 2), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  frac <- matrix(c(0.5, 0.3, 0.2,
                   0.5, 0.3, 0.2,
                   0.2, 0.3, 0.5,
                   0.2, 0.3, 0.5), nrow = 4, byrow = TRUE,
                 dimnames = list(paste0("s", 1:4), c("A", "A", "B")))
  expect_error(
    decompose_composition(c(g1 = 1, g2 = 2), frac, mu,
                          c(FALSE, FALSE, TRUE, TRUE), validate = FALSE),
    "duplicated cell-type column names")
})

test_that("decompose_composition rejects duplicated mu cell types", {
  mu <- matrix(c(10, 1, 5, 2), nrow = 2,
               dimnames = list(c("A", "A"), c("g1", "g2")))
  frac <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4, byrow = TRUE,
                 dimnames = list(paste0("s", 1:4), c("A", "B")))
  expect_error(
    decompose_composition(c(g1 = 1, g2 = 2), frac, mu,
                          c(FALSE, FALSE, TRUE, TRUE), validate = FALSE),
    "duplicated cell-type row names")
})

test_that("decompose_composition rejects an unnamed mean_expr", {
  mu <- matrix(c(10, 1, 5, 2), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  frac <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4, byrow = TRUE,
                 dimnames = list(paste0("s", 1:4), c("A", "B")))
  expect_error(
    decompose_composition(c(g1 = 1, g2 = 2), frac, mu,
                          c(FALSE, FALSE, TRUE, TRUE),
                          mean_expr = c(10, 20), validate = FALSE),
    "must be named")
})

test_that("decompose_uq accepts both mu orientations and gives the same answer", {
  ## cell types x genes (the decompose_composition convention)
  mu_ctg <- matrix(c(10, 1, 5, 2, 8, 3, 7, 4), nrow = 2, byrow = TRUE,
                   dimnames = list(c("A", "B"), c("g1", "g2", "g3", "g4")))
  ## genes x cell types (the historical decompose_uq convention)
  mu_gtc <- t(mu_ctg)
  od <- c(g1 = 1, g2 = -1, g3 = 0.5, g4 = 0.2)
  f1 <- c(A = 0.6, B = 0.4); f2 <- c(A = 0.4, B = 0.6)
  r1 <- decompose_uq(od, f1, f2, mu_ctg, mu_ctg,
                     f_se_ctrl = c(0.01, 0.01), f_se_case = c(0.01, 0.01),
                     n_boot = 20, seed = 3)
  r2 <- decompose_uq(od, f1, f2, mu_gtc, mu_gtc,
                     f_se_ctrl = c(0.01, 0.01), f_se_case = c(0.01, 0.01),
                     n_boot = 20, seed = 3)
  expect_equal(r1$gene, r2$gene)
  expect_equal(r1$resid_med, r2$resid_med, tolerance = 1e-12)
})

test_that("decompose_uq rejects an ambiguous mu", {
  ## cell-type names on BOTH dimensions
  mu <- matrix(1, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(
    decompose_uq(c(A = 1, B = 2), c(A = 0.6, B = 0.4), c(A = 0.4, B = 0.6),
                 mu, mu, n_boot = 5),
    "ambiguous")
})

test_that("decompose_uq rejects a f_se of the wrong length", {
  mu <- matrix(c(10, 1, 5, 2, 8, 3, 7, 4), nrow = 2, byrow = TRUE,
               dimnames = list(c("A", "B"), c("g1", "g2", "g3", "g4")))
  od <- c(g1 = 1, g2 = -1, g3 = 0.5, g4 = 0.2)
  f1 <- c(A = 0.6, B = 0.4); f2 <- c(A = 0.4, B = 0.6)
  ## two cell types in use, but f_se has three elements
  expect_error(
    decompose_uq(od, f1, f2, mu, mu, f_se_ctrl = c(0.01, 0.01, 0.01),
                 f_se_case = c(0.01, 0.01), n_boot = 10),
    "has length")
})

test_that("decompose_uq rejects a f_case missing a cell-type name", {
  mu <- matrix(c(10, 1, 5, 2, 8, 3, 7, 4), nrow = 2, byrow = TRUE,
               dimnames = list(c("A", "B"), c("g1", "g2", "g3", "g4")))
  od <- c(g1 = 1, g2 = -1, g3 = 0.5, g4 = 0.2)
  expect_error(
    decompose_uq(od, c(A = 0.6, B = 0.4), c(A = 0.5), mu, mu,
                 f_se_ctrl = c(0.01, 0.01), f_se_case = c(0.01, 0.01),
                 n_boot = 10),
    "missing name")
})

test_that("null_calibrate labels the 95th percentile correctly", {
  ## Reviewer mutation: with probs lacking 0.95, exceeds_95 was NA or silently
  ## received the 99th percentile.
  set.seed(71)
  n_g <- 400
  nm <- matrix(stats::rnorm(n_g * 20, 0, 0.2), n_g, 20,
               dimnames = list(paste0("g", 1:n_g), NULL))
  obs <- stats::setNames(stats::rnorm(n_g, 0, 0.2), paste0("g", 1:n_g))
  ## probs WITHOUT 0.95
  out <- null_calibrate(nm, obs, probs = c(0.5, 0.9, 0.99))
  expect_false(any(is.na(out$genes$exceeds_95)))
  q95 <- stats::quantile(out$null_pooled, 0.95, names = FALSE)
  q99 <- stats::quantile(out$null_pooled, 0.99, names = FALSE)
  expect_equal(out$genes$exceeds_95, out$genes$abs_resid > q95)
  expect_equal(out$genes$exceeds_99, out$genes$abs_resid > q99)
  ## the two flags must not be identical in general
  expect_false(identical(out$genes$exceeds_95, out$genes$exceeds_99))
})

test_that("marker_consistency rejects a malformed pass_range", {
  g <- paste0("g", 1:50)
  set.seed(81)
  comp <- stats::setNames(stats::rnorm(50, 0, .5), g)
  obs <- comp * 1.0
  mk <- list(A = g[1:6])
  expect_error(marker_consistency(obs, comp, mk, pass_range = 1.0),
               "two finite numbers")
  expect_error(marker_consistency(obs, comp, mk, pass_range = c(2, 1)),
               "increasing")
})

test_that("cross_cohort_agreement rejects mislabelled comp_list", {
  g <- paste0("g", 1:200)
  set.seed(91)
  o <- list(A = stats::setNames(stats::rnorm(200), g),
            B = stats::setNames(stats::rnorm(200), g))
  c1 <- list(X = o$A, Y = o$B)   # wrong names
  expect_error(cross_cohort_agreement(o, c1), "names must match")
})

test_that("cross_cohort_agreement rejects duplicated cohort names", {
  g <- paste0("g", 1:200)
  set.seed(92)
  o <- list(A = stats::setNames(stats::rnorm(200), g),
            A = stats::setNames(stats::rnorm(200), g))
  expect_error(cross_cohort_agreement(o), "unique")
})
