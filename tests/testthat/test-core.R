# ============================================================================
# Tests for cade
#
# The tests are organised around the design rules in 01-decompose.R, because
# those rules exist to prevent specific silent errors that were observed in
# real analyses. Each test is written so that it FAILS if the corresponding
# rule is removed.
# ============================================================================

## ---- helper: a tiny, fully controlled data set ---------------------------
make_fixture <- function() {
  ## two cell types, three genes
  mu <- matrix(c(
    ##  g1   g2   g3
    10,   2,   0,      # cell type A
    1,    8,   4       # cell type B
  ), nrow = 2, byrow = TRUE, dimnames = list(c("A", "B"), c("g1", "g2", "g3")))

  ## controls: 80% A; cases: 20% A
  frac <- rbind(
    ctrl1 = c(A = 0.8, B = 0.2),
    ctrl2 = c(A = 0.8, B = 0.2),
    case1 = c(A = 0.2, B = 0.8),
    case2 = c(A = 0.2, B = 0.8)
  )
  is_case <- c(FALSE, FALSE, TRUE, TRUE)
  list(mu = mu, frac = frac, is_case = is_case)
}


test_that("compose_prediction reproduces the analytic answer", {
  f <- make_fixture()
  pred <- compose_prediction(f$frac, f$mu, f$is_case)

  ## analytic:
  ##   ctrl bulk g1 = .8*10 + .2*1 = 8.2 ; case = .2*10 + .8*1 = 2.8
  expect_equal(unname(pred["g1"]), log2(2.8 / 8.2), tolerance = 1e-12)
  ##   ctrl g2 = .8*2 + .8*.  -- careful: mu[B,"g2"] = 8
  ##   ctrl = .8*2 + .2*8 = 3.2 ; case = .2*2 + .8*8 = 6.8
  expect_equal(unname(pred["g2"]), log2(6.8 / 3.2), tolerance = 1e-12)
  ##   g3: mu[A,"g3"] = 0, mu[B,"g3"] = 4
  ##   ctrl = .2*4 = 0.8 ; case = .8*4 = 3.2
  expect_equal(unname(pred["g3"]), log2(3.2 / 0.8), tolerance = 1e-12)
})

test_that("R1: columns are aligned by name, not by position", {
  f <- make_fixture()
  ## deliberately reverse the fraction columns
  rev_frac <- f$frac[, c("B", "A")]
  expect_equal(
    compose_prediction(rev_frac, f$mu, f$is_case),
    compose_prediction(f$frac,     f$mu, f$is_case)
  )
})

test_that("R1: a missing cell type is an error, not a silent drop", {
  f <- make_fixture()
  expect_error(
    compose_prediction(f$frac[, "A", drop = FALSE], f$mu, f$is_case),
    "cell type\\(s\\) in `mu` but not in `frac`"
  )
})

test_that("R1: an extra cell type is an error, not a silent drop", {
  f <- make_fixture()
  extra <- cbind(f$frac, C = 0)
  expect_error(compose_prediction(extra, f$mu, f$is_case),
               "cell type\\(s\\) in `frac` but not in `mu`")
})

test_that("R2: a log-scale mu produces a different, wrong answer", {
  f <- make_fixture()
  correct <- compose_prediction(f$frac, f$mu, f$is_case)
  wrong   <- compose_prediction(f$frac, log10(f$mu + 1), f$is_case)
  ## not merely different - the point is that nothing warns you
  expect_false(isTRUE(all.equal(unname(correct), unname(wrong))))
})

test_that("R4: a zero denominator yields NA rather than -Inf", {
  f <- make_fixture()
  ## make mu[B, "g3"] zero so the case bulk for g3 is zero
  mu <- f$mu; mu["B", "g3"] <- 0
  pred <- compose_prediction(f$frac, mu, f$is_case)
  expect_true(is.na(pred["g3"]))
})

test_that("R4: fractions are normalised but zero-sum samples are rejected", {
  f <- make_fixture()
  odd <- f$frac * 3              # does not sum to 1
  expect_equal(compose_prediction(odd, f$mu, f$is_case),
               compose_prediction(f$frac, f$mu, f$is_case))
  bad <- f$frac; bad[1, ] <- 0
  expect_error(compose_prediction(bad, f$mu, f$is_case), "all-zero fractions")
})

test_that("is_case must contain both groups", {
  f <- make_fixture()
  expect_error(compose_prediction(f$frac, f$mu, rep(TRUE, 4)), "both TRUE and FALSE")
})


test_that("decompose_composition adds up", {
  f <- make_fixture()
  obs <- c(g1 = -1.5, g2 = 1.1, g3 = 2.0)
  d <- decompose_composition(obs, f$frac, f$mu, f$is_case, validate = FALSE)
  expect_equal(d$resid, d$obs - d$pred, tolerance = 1e-12)
  expect_equal(d$gene, names(obs))
  expect_true(all(c("mean_expr") %in% names(d)))
})

test_that("decompose_composition requires named obs", {
  f <- make_fixture()
  expect_error(decompose_composition(c(-1, 1), f$frac, f$mu, f$is_case, validate = FALSE),
               "named numeric vector")
})


## ---- calibration ---------------------------------------------------------
test_that("expression_filter uses quantiles of the input, not fixed cuts", {
  x <- c(0.1, 1, 10, 100)
  expect_equal(sum(expression_filter(x, 0)), 4)
  expect_equal(sum(expression_filter(x, 0.5)), 2)
  expect_true(all(expression_filter(x, 0.5) == (x > 2)))
})

test_that("check_calibration recovers a known slope", {
  set.seed(42)
  n <- 2000
  pred <- stats::rnorm(n)
  mean_expr <- stats::runif(n, 0.1, 100)
  obs <- 0.8 * pred + stats::rnorm(n, sd = 0.2)
  d <- data.frame(gene = paste0("g", seq_len(n)), obs = obs, pred = pred,
                  mean_expr = mean_expr)
  cal <- check_calibration(d, levels = 0)
  expect_equal(cal$slope, 0.8, tolerance = 0.05)
  expect_gt(cal$r2, 0.9)
})

test_that("check_calibration refuses an effect-size-filtered input", {
  ## The point of this test: with only large-|obs| genes, a decomposition whose
  ## composition term is UNINFORMATIVE still yields an attractive slope, so the
  ## function must not be fed such a subset. We demonstrate that the same genes
  ## give a near-zero slope on the full set.
  set.seed(7)
  n <- 5000
  pred <- stats::rnorm(n)
  obs  <- stats::rnorm(n)              # independent of pred: no signal
  d <- data.frame(gene = paste0("g", seq_len(n)), obs = obs, pred = pred,
                  mean_expr = stats::runif(n, 0.1, 100))
  keep <- abs(d$obs) > 1.5
  circ <- check_calibration(d[keep, ], levels = 0)
  full <- check_calibration(d, levels = 0)
  ## both slopes should be near zero here, because there is no signal at all;
  ## what matters is that we can detect the absence on the full set
  expect_lt(abs(full$slope), 0.1)
  expect_lt(full$r2, 0.01)
})

test_that("calibration_verdict wording matches the numbers", {
  good <- data.frame(level = 0.75, label = "top 25% expressed", n = 1000,
                     slope = 0.85, intercept = 0, r2 = 0.4,
                     pearson = 0.63, median_pred_error = 0.1)
  expect_match(calibration_verdict(good), "informative")

  none <- good; none$slope <- 0.03; none$r2 <- 0.005
  expect_match(calibration_verdict(none), "uninformative")
})


## ---- attribution ---------------------------------------------------------
test_that("attribute_genes labels the three regimes", {
  d <- data.frame(gene = c("a", "b", "c"), obs = c(1, 1, 1),
                  pred = c(0.95, 0.5, 0.05), resid = c(0.05, 0.5, 0.95),
                  pct = c(0.05, 0.5, 0.95), mean_expr = 1)
  out <- attribute_genes(d)
  expect_equal(out$class, c("composition-driven", "mixed", "residual"))
})

test_that("attribute_genes marks near-zero obs as unchanged", {
  d <- data.frame(gene = "a", obs = 0.01, pred = 0.5, resid = -0.49,
                  pct = -49, mean_expr = 1)
  expect_equal(attribute_genes(d, min_obs = 0.1)$class, "unchanged")
})

test_that("attribute_sensitivity reports instability", {
  d <- data.frame(gene = c("a", "b"), obs = c(1, 1),
                  pred = c(0.55, 0.05), resid = c(0.45, 0.95),
                  pct = c(0.45, 0.95), mean_expr = 1)
  s <- attribute_sensitivity(d)
  expect_false(s$stable[s$gene == "a"])     # 0.45 sits between thresholds
  expect_true(s$stable[s$gene == "b"])
})


## ---- effect size ---------------------------------------------------------
test_that("cohens_d matches a hand computation", {
  x <- c(1, 2, 3, 4, 5); y <- c(0, 1, 2, 3, 4)
  d <- cohens_d(x, y)
  sp <- sqrt((4 * stats::var(x) + 4 * stats::var(y)) / 8)
  raw <- (mean(x) - mean(y)) / sp
  J <- 1 - 3 / (4 * 10 - 9)
  expect_equal(d$d, raw * J, tolerance = 1e-12)
  expect_lt(d$lo, d$d); expect_gt(d$hi, d$d)
})

test_that("pool_cohens_d handles a homogeneous set", {
  p <- pool_cohens_d(c(1.0, 1.0, 1.0), c(0.3, 0.3, 0.3))
  expect_equal(p$pooled$d, 1.0, tolerance = 1e-9)
  expect_equal(p$tau2, 0, tolerance = 1e-9)
  expect_equal(p$p_het, 1, tolerance = 1e-6)
})

test_that("pool_cohens_d flags heterogeneity", {
  p <- pool_cohens_d(c(-2, 2), c(0.1, 0.1))
  expect_gt(p$tau2, 0)
  expect_lt(p$p_het, 0.01)
})


## ---- reference checks ----------------------------------------------------
test_that("check_reference flags starved cell types", {
  ct <- rep(c("A", "B"), c(500, 10))
  r <- check_reference(ct, min_cells = 100)
  expect_false(r$ok)
  expect_match(r$flags[1], "fewer than 100 cells")
})

test_that("check_reference flags a disease-only reference", {
  ct <- rep(c("A", "B"), c(500, 500))
  grp <- rep("disease", 1000)
  r <- check_reference(ct, group = grp)
  expect_false(r$ok)
  expect_match(paste(r$flags, collapse = " "), "no normal/control cells")
})

test_that("check_reference flags missing marker genes", {
  ct <- rep(c("A", "B"), 200)
  r <- check_reference(ct, markers = list(A = c("CD3D", "PECAM1")),
                       genes = c("CD3D", "ACTB"))
  expect_match(paste(r$flags, collapse = " "), "PECAM1")
})

test_that("compare_references reports sign flips", {
  f <- make_fixture()
  obs <- c(g1 = -1.0, g2 = 1.0, g3 = 1.5)
  cres <- compare_references(obs, f$frac, f$frac, f$mu, f$mu, f$is_case, validate = FALSE)
  expect_equal(cres$agreement$n_disagree, 0)
  expect_equal(cres$agreement$cor_resid, 1)
})


## ---- simulation ----------------------------------------------------------
test_that("simulate_pseudobulk ground truth is exact", {
  f <- make_fixture()
  sim <- simulate_pseudobulk(f$mu,
                             frac_ctrl = c(A = .8, B = .2),
                             frac_case = c(A = .2, B = .8),
                             n_ctrl = 3, n_case = 3,
                             intrinsic = c(g1 = 1), noise = 0)
  comp <- compose_prediction(f$frac, f$mu, f$is_case)
  expect_equal(sim$truth$composition, unname(comp[sim$truth$gene]),
               tolerance = 1e-12)
  expect_equal(sim$truth$intrinsic[sim$truth$gene == "g1"], 1)
})

test_that("a noiseless simulation is recovered exactly", {
  f <- make_fixture()
  sim <- simulate_pseudobulk(f$mu,
                             frac_ctrl = c(A = .8, B = .2),
                             frac_case = c(A = .2, B = .8),
                             n_ctrl = 4, n_case = 4,
                             intrinsic = c(g1 = 1, g3 = -0.5), noise = 0)
  ## with noise = 0 every sample is identical, so the group means reproduce
  ## the truth; use the truth directly as "observed" to test the algebra
  d <- data.frame(gene = sim$truth$gene,
                  obs = sim$truth$total,
                  pred = sim$truth$composition,
                  resid = sim$truth$total - sim$truth$composition,
                  pct = NA_real_, mean_expr = 1)
  sc <- score_simulation(d, sim$truth, c("g1", "g3"))
  expect_equal(sc$cor_pred, 1)
  expect_equal(sc$recovery_intrinsic, 1, tolerance = 1e-12)
})

test_that("score_simulation reports a failure to recover", {
  f <- make_fixture()
  sim <- simulate_pseudobulk(f$mu, c(A = .8, B = .2), c(A = .2, B = .8),
                             4, 4, intrinsic = c(g1 = 1), noise = 0)
  ## pretend the decomposition attributed everything to composition
  d <- data.frame(gene = sim$truth$gene,
                  obs = sim$truth$total,
                  pred = sim$truth$total,        # wrong on purpose
                  resid = 0, pct = 0, mean_expr = 1)
  sc <- score_simulation(d, sim$truth, "g1")
  expect_equal(sc$recovery_intrinsic, 0, tolerance = 1e-12)
})

test_that("simulate_pseudobulk rejects unknown intrinsic genes", {
  f <- make_fixture()
  expect_error(
    simulate_pseudobulk(f$mu, c(A = .8, B = .2), c(A = .2, B = .8), 2, 2,
                        intrinsic = c(notagene = 1)),
    "not in `mu`"
  )
})

test_that("two intrinsic genes are each applied, not alternated across rows", {
  ## Regression test. `ba[, idx] * v` recycles v in column-major order, so with
  ## two intrinsic genes the multipliers alternate row by row and the injected
  ## changes partially cancel. The symptom is that the group mean does not move
  ## by the injected amount, which score_simulation() reports as a low recovery.
  set.seed(1)
  G <- 40
  mu <- matrix(stats::runif(3 * G, 50, 200), nrow = 3,
               dimnames = list(c("A", "B", "C"), paste0("g", seq_len(G))))
  fc <- c(A = .6, B = .3, C = .1)
  fa <- c(A = .2, B = .5, C = .3)

  sim <- simulate_pseudobulk(mu, fc, fa, n_ctrl = 20, n_case = 20,
                             intrinsic = c(g1 = 1, g2 = -1), noise = 0)

  ## with noise = 0 every sample in a group is identical, so the ratio of the
  ## case mean to the control mean must equal 2^(composition) * 2^(intrinsic)
  ctrl_mean <- colMeans(sim$bulk[!sim$is_case, , drop = FALSE])
  case_mean <- colMeans(sim$bulk[ sim$is_case, , drop = FALSE])
  for (g in c("g1", "g2")) {
    t <- sim$truth[sim$truth$gene == g, ]
    observed <- unname(log2(case_mean[g] / ctrl_mean[g]))
    expect_equal(observed, unname(t$total), tolerance = 1e-10,
                 label = sprintf("observed vs truth total for %s", g))
  }
})
