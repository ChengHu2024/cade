# ============================================================================
# Regression tests for the defects found by the adversarial audit
#
# Each test names the audit finding it pins down, so that a future regression
# points straight back at the reason the guard exists.
# ============================================================================

## ---- P3: single-gene input --------------------------------------------------
test_that("P3: attribute_sensitivity returns ONE row for one gene", {
  d <- data.frame(gene = "a", obs = 1, pred = 0.4, resid = 0.6,
                  pct = 0.6, mean_expr = 1)
  s <- attribute_sensitivity(d)
  expect_equal(nrow(s), 1)
  expect_equal(s$gene, "a")
  expect_equal(s$n_thresholds, 6)              # seq(0.2, 0.45, 0.05)
  expect_true(s$frac_residual >= 0 && s$frac_residual <= 1)
})

test_that("P3: a single threshold does not collapse the matrix", {
  d <- data.frame(gene = c("a", "b"), obs = c(1, 1),
                  pred = c(0.4, 0.1), resid = c(0.6, 0.9),
                  pct = c(0.6, 0.9), mean_expr = 1)
  s <- attribute_sensitivity(d, thresholds = 0.3)
  expect_equal(nrow(s), 2)
  expect_true(all(s$n_thresholds == 1))
})

test_that("P3: stability of a single gene matches the multi-gene case", {
  ## build a 3-gene frame, take each gene alone, and check that the stability
  ## verdict is identical in both runs
  d <- data.frame(gene = c("a", "b", "c"), obs = c(1, 1, 1),
                  pred = c(0.4, 0.35, 0.05), resid = c(0.6, 0.65, 0.95),
                  pct = c(0.6, 0.65, 0.95), mean_expr = 1)
  all3 <- attribute_sensitivity(d)
  for (g in d$gene) {
    one <- attribute_sensitivity(d[d$gene == g, ])
    expect_equal(one$stable, all3$stable[all3$gene == g],
                 label = sprintf("stable for gene %s", g))
    expect_equal(one$frac_residual, all3$frac_residual[all3$gene == g],
                 label = sprintf("frac_residual for gene %s", g))
  }
})


## ---- P7: signed pct with an overshooting composition term -------------------
test_that("P7: an overshooting composition term is not called composition-driven", {
  ## obs = +1, pred = +3 (overshoot), resid = -2  -> pct = -2
  d <- data.frame(gene = "g", obs = 1, pred = 3, resid = -2,
                  pct = -2, mean_expr = 1)
  out <- attribute_genes(d)
  expect_equal(out$class, "over-predicted")
  expect_match(out$class_reason, "opposite in sign")
})

test_that("P7: a large residual with an overshoot is not hidden", {
  ## the case the audit found: residual is 95% of the change but pct is signed
  d <- data.frame(gene = "g", obs = 1, pred = 1.95, resid = -0.95,
                  pct = -0.95, mean_expr = 1)
  expect_equal(attribute_genes(d)$class, "over-predicted")
})

test_that("P7: same-sign genes are still classified by pct", {
  d <- data.frame(gene = c("a", "b", "c"), obs = c(1, 1, 1),
                  pred = c(0.95, 0.5, 0.05), resid = c(0.05, 0.5, 0.95),
                  pct = c(0.05, 0.5, 0.95), mean_expr = 1)
  expect_equal(attribute_genes(d)$class,
               c("composition-driven", "mixed", "residual"))
})

test_that("P7: a negative observed change works the same way", {
  d <- data.frame(gene = "g", obs = -1, pred = -0.1, resid = -0.9,
                  pct = 0.9, mean_expr = 1)
  expect_equal(attribute_genes(d)$class, "residual")
})


## ---- P6: is_case containing NA must stop, not silently produce all-NA -------
test_that("P6: is_case with NA is rejected by validation", {
  G <- 30
  set.seed(2)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(0.6, 0.4), 6), ncol = 2, byrow = TRUE,
               dimnames = list(paste0("s", 1:6), c("A", "B")))
  isc <- c(FALSE, FALSE, NA, TRUE, TRUE, TRUE)
  chk <- check_inputs(stats::setNames(stats::rnorm(G), colnames(mu)),
                      fr, mu, isc)
  expect_false(chk$ok)
  expect_true(any(grepl("NA", chk$errors)))
})

test_that("P6: a legitimate call never yields all-NA classification", {
  G <- 30
  set.seed(3)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(0.6, 0.4, 0.3, 0.7), 3), ncol = 2, byrow = TRUE,
               dimnames = list(paste0("s", 1:6), c("A", "B")))
  isc <- c(rep(FALSE, 3), rep(TRUE, 3))
  obs <- stats::setNames(stats::rnorm(G), colnames(mu))
  d <- decompose_composition(obs, fr, mu, isc)
  a <- attribute_genes(d)
  expect_false(all(a$class == "unchanged"))
})


## ---- P4 / P5: duplicate names must not be silently merged -------------------
test_that("P4: duplicated cell-type names are rejected before any arithmetic", {
  G <- 30
  set.seed(4)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(0.6, 0.4, 0.3, 0.7), 3), ncol = 2, byrow = TRUE)
  colnames(fr) <- c("A", "A")                  # duplicate on purpose
  obs <- stats::setNames(stats::rnorm(G), colnames(mu))
  expect_error(compose_prediction(fr, mu, c(rep(FALSE, 3), rep(TRUE, 3))),
               "cell type\\(s\\) in `mu` but not in `frac`|duplicated")
})

test_that("P5: duplicated mu row names are rejected", {
  ## Regression: this test used to pass because `frac` carried a cell type
  ## that `mu` lacked, not because of the duplication. Here `frac` has exactly
  ## the cell types `mu` declares, so duplication is the only defect.
  G <- 30
  set.seed(5)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "A"), paste0("g", seq_len(G))))
  fr <- matrix(rep(1, 6), ncol = 1,
               dimnames = list(paste0("s", 1:6), "A"))
  expect_error(compose_prediction(fr, mu, c(rep(FALSE, 3), rep(TRUE, 3))),
               "duplicated cell-type row names")
})


## ---- P8: a log-scaled mu must stop ----------------------------------------
test_that("P8: a log-scaled mu stops decompose_composition", {
  G <- 40
  set.seed(6)
  mu <- matrix(stats::runif(2 * G, 100, 600), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(0.6, 0.4, 0.3, 0.7), 3), ncol = 2, byrow = TRUE,
               dimnames = list(paste0("s", 1:6), c("A", "B")))
  obs <- stats::setNames(stats::rnorm(G), colnames(mu))
  isc <- c(rep(FALSE, 3), rep(TRUE, 3))

  expect_silent(decompose_composition(obs, fr, mu, isc))       # linear: fine
  expect_error(decompose_composition(obs, fr, log2(mu + 1), isc),
               "log-scaled|compressed")
})
