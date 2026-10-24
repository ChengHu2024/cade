# ============================================================================
# Tests for automatic validation and self-test
#
# The theme: a wrong input must STOP the computation. A number that looks
# ordinary but is wrong is the failure this layer exists to prevent.
# ============================================================================

make_valid <- function() {
  G <- 60
  set.seed(11)
  mu <- matrix(stats::runif(3 * G, 50, 300), nrow = 3,
               dimnames = list(c("A", "B", "C"), paste0("g", seq_len(G))))
  fr <- rbind(matrix(rep(c(.6, .3, .1), 10), ncol = 3, byrow = TRUE),
              matrix(rep(c(.2, .5, .3), 10), ncol = 3, byrow = TRUE))
  colnames(fr) <- c("A", "B", "C")
  list(obs = stats::setNames(stats::rnorm(G), colnames(mu)),
       frac = fr, mu = mu,
       is_case = c(rep(FALSE, 10), rep(TRUE, 10)))
}


test_that("a valid input passes validation", {
  v <- make_valid()
  chk <- check_inputs(v$obs, v$frac, v$mu, v$is_case)
  expect_true(chk$ok)
  expect_length(chk$errors, 0)
})

## ---- the semantic trap that motivates the whole layer --------------------
test_that("a log-scaled mu STOPS the computation", {
  v <- make_valid()
  mu_log <- log2(v$mu + 1)                     # wrong on purpose
  chk <- check_inputs(v$obs, v$frac, mu_log, v$is_case)
  expect_false(chk$ok)
  expect_true(any(grepl("log-scaled|compressed", chk$errors)))
})

test_that("decompose_composition refuses a log-scaled mu rather than returning a number", {
  v <- make_valid()
  mu_log <- log2(v$mu + 1)
  expect_error(
    decompose_composition(v$obs, v$frac, mu_log, v$is_case),
    "input validation failed"
  )
})

test_that("decompose_composition still accepts a correct linear mu", {
  v <- make_valid()
  d <- decompose_composition(v$obs, v$frac, v$mu, v$is_case)
  expect_s3_class(d, "cade_decomp")
  expect_equal(d$resid, d$obs - d$pred, tolerance = 1e-12)
})

test_that("validate = FALSE bypasses the check", {
  v <- make_valid()
  mu_log <- log2(v$mu + 1)
  d <- decompose_composition(v$obs, v$frac, mu_log, v$is_case, validate = FALSE)
  expect_s3_class(d, "cade_decomp")
})

## ---- structural traps ----------------------------------------------------
test_that("duplicated cell-type names are rejected", {
  v <- make_valid()
  fr <- v$frac; colnames(fr) <- c("A", "A", "C")
  expect_false(check_inputs(v$obs, fr, v$mu, v$is_case)$ok)
})

test_that("duplicated mu row names are rejected", {
  v <- make_valid()
  mu <- v$mu; rownames(mu) <- c("A", "A", "C")
  expect_false(check_inputs(v$obs, v$frac, mu, v$is_case)$ok)
})

test_that("a case/control name mismatch is rejected", {
  v <- make_valid()
  fr <- v$frac; colnames(fr) <- c("A", "b", "C")   # lower-case b
  chk <- check_inputs(v$obs, fr, v$mu, v$is_case)
  expect_false(chk$ok)
  expect_true(any(grepl("not in `mu`", chk$errors)))
})

test_that("an all-zero cell type in mu is rejected", {
  v <- make_valid()
  mu <- v$mu; mu["C", ] <- 0
  chk <- check_inputs(v$obs, v$frac, mu, v$is_case)
  expect_false(chk$ok)
  expect_true(any(grepl("all-zero", chk$errors)))
})

test_that("wrong is_case length is rejected", {
  v <- make_valid()
  expect_false(check_inputs(v$obs, v$frac, v$mu, v$is_case[1:5])$ok)
})

test_that("is_case with only one group is rejected", {
  v <- make_valid()
  expect_false(check_inputs(v$obs, v$frac, v$mu, rep(TRUE, 20))$ok)
  expect_false(check_inputs(v$obs, v$frac, v$mu, rep(FALSE, 20))$ok)
})

test_that("NA in obs is rejected", {
  v <- make_valid()
  v$obs[3] <- NA_real_
  expect_false(check_inputs(v$obs, v$frac, v$mu, v$is_case)$ok)
})

test_that("integer mu above 1000 warns about raw counts", {
  v <- make_valid()
  mu <- round(v$mu * 100)
  chk <- check_inputs(v$obs, v$frac, mu, v$is_case)
  expect_true(chk$ok)                       # warning, not error
  expect_true(any(grepl("raw counts", chk$warnings)))
})

test_that("genes absent from mu produce a warning, not an error", {
  v <- make_valid()
  obs <- v$obs
  names(obs)[1:2] <- c("not_in_mu_1", "not_in_mu_2")
  chk <- check_inputs(obs, v$frac, v$mu, v$is_case)
  expect_true(chk$ok)
  expect_true(any(grepl("absent from `mu`", chk$warnings)))
})

test_that("strict mode promotes warnings to errors", {
  v <- make_valid()
  obs <- v$obs
  names(obs)[1] <- "not_in_mu"
  expect_true(check_inputs(obs, v$frac, v$mu, v$is_case, strict = FALSE)$ok)
  expect_false(check_inputs(obs, v$frac, v$mu, v$is_case, strict = TRUE)$ok)
})


## ---- self test -----------------------------------------------------------
test_that("self_test passes on a healthy installation", {
  st <- self_test(n_genes = 300, n_ctrl = 15, n_case = 15, verbose = FALSE)
  expect_true(st$pass)
  expect_gt(st$recovery, 0.90)
  expect_gt(st$cor_pred, 0.99)
  expect_true(all(st$checks$pass))
})

test_that("self_test reports the criteria it evaluated", {
  st <- self_test(n_genes = 200, n_ctrl = 10, n_case = 10, verbose = FALSE)
  ## 四条：三条算术判据 + 一条"报告层不得改动数字"
  expect_equal(nrow(st$checks), 4)
  expect_true(all(c("criterion", "measured", "threshold", "pass") %in%
                  names(st$checks)))
  expect_true(all(st$checks$pass))
})

test_that("self_test includes the reporting-invariance criterion", {
  st <- self_test(n_genes = 200, n_ctrl = 10, n_case = 10, verbose = FALSE)
  expect_true(any(grepl("reporting does not alter", st$checks$criterion)))
  i <- grep("reporting does not alter", st$checks$criterion)
  expect_equal(st$checks$measured[i], 1)      # bit-identical
})

## ---- print method --------------------------------------------------------
test_that("print.cade_decomp summarises without error", {
  v <- make_valid()
  d <- decompose_composition(v$obs, v$frac, v$mu, v$is_case)
  expect_output(print(d), "cade decomposition")
  expect_output(print(d), "composition term")
})
