# ============================================================================
# Tests for verify_deposit()
#
# The theme: a row permutation that leaves every sample-level statistic intact
# can still destroy every gene-level conclusion. Every scenario below is a
# known-answer case built from one synthetic 500 x 10 matrix.
#
#   A / D  rows in the correct order
#            -> A1, A2, A3 all pass; verdict "ok"
#
#   B      intensity and detection p moved together, probe IDs left behind
#            -> this is the A1 blind spot: A1 does not change. A2 (with an
#               order-aligned reference) fails and the permutation is recovered,
#               so the overall verdict is "fail".
#
#   B2     intensity and detection p moved with different permutations
#            -> A1 fails: the within-column pairs are broken.
#
#   C      only intensity moved; probe IDs and detection p stay in order
#            -> C = 只有 intensity 被行置换，probe_ids 与 detection_p 保持原序
#               → 这对应 GSE47472：39,426 行强度块错位，探针 ID 与 P 值块不错位
#               A1 fails on every column; pi is recovered and
#               intensity[pi, ] restores the correct order.
# ============================================================================

## One synthetic deposit with known answers. I > 0 so that A3 can use log2().
## `sign` selects the encoding of the detection block: raw Illumina p-values
## fall with signal ("negative"), a detection score or -log10(p) rises with it
## ("positive"). The original case measured the positive convention.
deposit_fixture <- function(n = 500L, m = 10L, seed = 20261025,
                            sign = c("negative", "positive")) {
  sign <- match.arg(sign)
  set.seed(seed)
  row_effect <- stats::rnorm(n, 2, 0.6)
  I <- exp(row_effect + matrix(stats::rnorm(n * m, 0, 0.15), n, m))
  dimnames(I) <- list(paste0("probe", seq_len(n)),
                      paste0("sample", seq_len(m)))
  z <- (I - mean(I)) / stats::sd(I)
  P <- if (sign == "negative") stats::plogis(-z) else stats::plogis(z)
  r  <- stats::rnorm(n)
  cc <- stats::rnorm(m)
  ## S = log2(I) + r_i + c_j exactly: A3 has a zero residual by construction.
  S <- log2(I) + r + rep(cc, each = n)
  list(I = I, P = P, S = S)
}


## ---------------------------------------------------------------------------
## A / D: the correct row order
## ---------------------------------------------------------------------------

test_that("D: a faithful deposit passes all three checks", {
  d <- deposit_fixture()
  res <- verify_deposit(d$I, d$P, probe_ids = rownames(d$I),
                        reference = d$S, reference_intensity = d$I)

  # the returned object is the documented class
  expect_s3_class(res, "cade_deposit")
  # nothing is out of order, so no check may fail
  expect_equal(res$verdict, "ok")
  # A1: intensity and detection p agree in every column
  expect_true(res$a1$pass)
  expect_equal(res$a1$n_fail, 0L)
  # A2: the row-mean fingerprint is exact here
  expect_true(res$a2$pass)
  expect_equal(res$a2$rho, 1, tolerance = 1e-12)
  # A3: the additive identity holds to floating-point precision
  expect_true(res$a3$pass)
  expect_equal(res$a3$slope, 1, tolerance = 1e-8)
  expect_lt(res$a3$residual_sd, 1e-10)
  expect_lt(res$a3$max_resid, 1e-10)
  expect_lt(res$a3$unexplained_var, 1e-20)
  # a perfect fit must not be hiding duplicated reference columns
  expect_length(res$a3$value_identical_columns, 0)
  # dimensions are reported back
  expect_equal(res$n_probes, nrow(d$I))
  expect_equal(res$n_samples, ncol(d$I))
  # the row order matches the reference exactly, so pi is the identity
  expect_true(res$permutation_recoverable)
  expect_true(res$permutation$is_identity)
  # the report carries the no-inference convention of the original case
  expect_match(res$summary, "不对成因作任何推断")
})

test_that("A: correct row order gives a high per-column |rho|", {
  d <- deposit_fixture()
  res <- verify_deposit(d$I, d$P, probe_ids = rownames(d$I), checks = "A1")

  # one correlation per sample column, named after the columns
  expect_equal(length(res$a1$rho_per_column), ncol(d$I))
  expect_equal(names(res$a1$rho_per_column), colnames(d$P))
  # every column is far above the default 0.7 threshold
  expect_true(all(abs(res$a1$rho_per_column) > 0.99))
  # and the summary counts no failures
  expect_equal(res$a1$n_fail, 0L)
  # only A1 was requested, so the other checks are recorded as not run
  expect_true(is.na(res$a2$pass))
  expect_true(is.na(res$a3$pass))
  expect_equal(res$checks_run, "A1")
  # A1 alone is enough here
  expect_equal(res$verdict, "ok")
})

test_that("A1 accepts either sign convention and reports which one it saw", {
  neg <- deposit_fixture(sign = "negative")
  pos <- deposit_fixture(sign = "positive")

  r_neg <- verify_deposit(neg$I, neg$P, checks = "A1")
  r_pos <- verify_deposit(pos$I, pos$P, checks = "A1")

  # raw p-value encoding: intensity and p-value move in opposite directions
  expect_true(r_neg$a1$pass)
  expect_equal(r_neg$a1$sign_convention, "negative")
  expect_true(all(r_neg$a1$rho_per_column < -0.99))
  # detection-score encoding: both rise together, and this is the convention
  # actually measured in the original case
  expect_true(r_pos$a1$pass)
  expect_equal(r_pos$a1$sign_convention, "positive")
  expect_true(all(r_pos$a1$rho_per_column > 0.99))

  # a file with both conventions in different columns is labelled "mixed",
  # and the |rho| criterion still passes
  z <- (neg$I - mean(neg$I)) / stats::sd(neg$I)
  P_mixed <- neg$P
  P_mixed[, 1:5]  <- stats::plogis(z[, 1:5])
  P_mixed[, 6:10] <- stats::plogis(-z[, 6:10])
  r_mixed <- verify_deposit(neg$I, P_mixed, checks = "A1")
  expect_true(r_mixed$a1$pass)
  expect_equal(r_mixed$a1$sign_convention, "mixed")
})


test_that("a passing measurement near its threshold is reported as suspect", {
  d <- deposit_fixture()

  # |rho| is 1.0; with min_rho = 0.99 the column passes but sits inside the
  # marginal band (threshold + 0.05), so the verdict is not a clean "ok"
  near <- verify_deposit(d$I, d$P, checks = "A1", min_rho = 0.99)
  expect_true(near$a1$pass)
  expect_equal(near$verdict, "suspect")

  # with the default threshold the same data gives a clean pass
  far <- verify_deposit(d$I, d$P, checks = "A1")
  expect_equal(far$verdict, "ok")
})


## ---------------------------------------------------------------------------
## B: both blocks moved together -- the A1 blind spot
## ---------------------------------------------------------------------------

test_that("B: moving both blocks together is invisible to A1 but caught by A2 and pi", {
  d <- deposit_fixture()
  set.seed(101)
  perm <- sample.int(nrow(d$I))
  I_b <- d$I[perm, , drop = FALSE]
  P_b <- d$P[perm, , drop = FALSE]
  base <- verify_deposit(d$I, d$P, checks = "A1")

  res <- verify_deposit(I_b, P_b, probe_ids = rownames(d$I),
                        reference = d$S, reference_intensity = d$I)

  # A1 is invariant to a common row permutation: the correlations are the same
  # as for the correctly ordered matrix (the blind spot, asserted explicitly)
  expect_true(res$a1$pass)
  expect_equal(res$a1$rho_per_column, base$a1$rho_per_column,
               tolerance = 1e-12)
  # A2 sees the permutation because row means have moved to other positions
  expect_false(res$a2$pass)
  expect_lt(res$a2$rho, 0.9)
  # the row tuples are intact, so pi is recoverable and non-identity
  expect_true(res$permutation_recoverable)
  expect_false(res$permutation$is_identity)
  # direction convention: intensity[pi, ] is in the reference row order
  expect_equal(res$permutation$pi, order(perm))
  expect_identical(I_b[res$permutation$pi, , drop = FALSE], d$I)
  # the same pi also restores the detection block, which moved with intensity
  expect_identical(P_b[res$permutation$pi, , drop = FALSE], d$P)
  # overall verdict is fail even though A1 passed
  expect_equal(res$verdict, "fail")
})


## ---------------------------------------------------------------------------
## C: the original case -- only intensity is out of order
## ---------------------------------------------------------------------------
## C = 只有 intensity 被行置换，probe_ids 与 detection_p 保持原序
##    → 这对应 GSE47472：39,426 行强度块错位，探针 ID 与 P 值块不错位
test_that("C: the original case, only intensity out of order, fails A1 and is repairable", {
  d <- deposit_fixture()
  set.seed(102)
  perm <- sample.int(nrow(d$I))
  I_c <- d$I[perm, , drop = FALSE]   # probe IDs and P stay in the original order
  res <- verify_deposit(I_c, d$P, probe_ids = rownames(d$I),
                        reference_intensity = d$I)

  # A1 fails on every column: the within-column pairs are broken
  expect_false(res$a1$pass)
  expect_equal(res$a1$n_fail, ncol(d$I))
  expect_lt(max(abs(res$a1$rho_per_column)), 0.1)
  # and the overall verdict is fail
  expect_equal(res$verdict, "fail")
  # the permutation is recovered with the documented direction
  expect_true(res$permutation_recoverable)
  expect_equal(res$permutation$pi, order(perm))

  # applying pi repairs the matrix: intensity[pi, ] equals the correct order
  I_fix <- I_c[res$permutation$pi, , drop = FALSE]
  expect_identical(I_fix, d$I)

  # after the repair, A1 recovers: this checks both detection and correctness
  # of pi, not merely that the un-repaired data failed
  fixed <- verify_deposit(I_fix, d$P, reference_intensity = d$I)
  expect_true(fixed$a1$pass)
  expect_gt(min(abs(fixed$a1$rho_per_column)), 0.9)
  expect_equal(fixed$verdict, "ok")
})

test_that("B2: different permutations of intensity and detection p fail A1", {
  d <- deposit_fixture()
  set.seed(103)
  perm_i <- sample.int(nrow(d$I))
  perm_p <- sample.int(nrow(d$I))

  res <- verify_deposit(d$I[perm_i, , drop = FALSE],
                        d$P[perm_p, , drop = FALSE],
                        probe_ids = rownames(d$I))

  # unlike B, the two blocks moved differently, so A1 sees the misalignment
  expect_false(res$a1$pass)
  expect_equal(res$a1$n_fail, ncol(d$I))
  expect_lt(max(abs(res$a1$rho_per_column)), 0.2)
  expect_equal(res$verdict, "fail")
})

test_that("pi direction: reference[i, ] == intensity[pi[i], ], not the inverse", {
  d <- deposit_fixture()
  set.seed(104)
  pi_true <- sample.int(nrow(d$I))
  # build the deposit with a known pi_true: I_dep[pi_true[i], ] == I[i, ]
  I_dep <- d$I[order(pi_true), , drop = FALSE]

  res <- verify_deposit(I_dep, d$P, reference_intensity = d$I)

  # the solver must return pi_true itself, not its inverse
  expect_identical(res$permutation$pi, pi_true)
  expect_false(identical(res$permutation$pi, order(pi_true)))
  # the defining equation holds row for row
  expect_identical(I_dep[res$permutation$pi, , drop = FALSE], d$I)
})


## ---------------------------------------------------------------------------
## A2: the row-mean fingerprint
## ---------------------------------------------------------------------------

test_that("A2: row means agree across halves and collapse when rows are moved", {
  d <- deposit_fixture()
  n <- nrow(d$I)
  half2 <- d$I[, 6:10, drop = FALSE]

  # row means from one half of the samples predict row means from the other
  # half; this is the fingerprint that a row permutation destroys
  ok <- verify_deposit(d$I, d$P, reference_intensity = half2, checks = "A2")
  expect_true(ok$a2$pass)
  expect_gt(ok$a2$rho, 0.9)
  expect_equal(ok$verdict, "ok")

  # permuting the rows of that half breaks the per-probe correspondence
  set.seed(201)
  half2_perm <- half2[sample.int(n), , drop = FALSE]
  bad <- verify_deposit(d$I, d$P, reference_intensity = half2_perm,
                        checks = "A2")
  expect_false(bad$a2$pass)
  expect_lt(bad$a2$rho, 0.5)
  expect_equal(bad$verdict, "fail")

  # a full row permutation of the reference collapses the fingerprint further
  set.seed(202)
  whole <- d$I[sample.int(n), , drop = FALSE]
  bad2 <- verify_deposit(d$I, d$P, reference_intensity = whole, checks = "A2")
  expect_false(bad2$a2$pass)
  expect_lt(abs(bad2$a2$rho), 0.2)
  expect_equal(bad2$verdict, "fail")
})

test_that("A2 without a reference is skipped and recorded as NA", {
  d <- deposit_fixture()
  res <- verify_deposit(d$I, d$P, checks = "A2")

  # the specification requires a silent skip to be recorded as NA, not as pass
  expect_true(is.na(res$a2$pass))
  expect_true(is.na(res$a2$rho))
  expect_match(res$a2$reference_note, "reference_intensity")
  # no check could run, so the verdict must not claim "ok"
  expect_equal(res$verdict, "suspect")
  expect_false(res$permutation_recoverable)
  expect_match(res$permutation$reason, "reference_intensity")
})


## ---------------------------------------------------------------------------
## A3: the additive identity
## ---------------------------------------------------------------------------

test_that("A3: the exact additive identity has residual SD at rounding scale", {
  d <- deposit_fixture()
  res <- verify_deposit(d$I, d$P, reference = d$S, checks = "A3")

  # S = log2(I) + r_i + c_j by construction, so the fit is exact
  expect_true(res$a3$pass)
  expect_equal(res$a3$slope, 1, tolerance = 1e-8)
  expect_lt(res$a3$residual_sd, 1e-10)
  expect_lt(res$a3$max_resid, 1e-10)
  expect_lt(res$a3$unexplained_var, 1e-20)
  # a perfect fit is only accepted when no reference columns are duplicates
  expect_length(res$a3$value_identical_columns, 0)
  expect_equal(res$verdict, "ok")
})

test_that("A3: 0.01 noise is detected and fails", {
  d <- deposit_fixture()
  n <- nrow(d$I)
  m <- ncol(d$I)
  set.seed(2026)
  S_noise <- d$S + matrix(stats::rnorm(n * m, 0, 0.01), n, m)

  res <- verify_deposit(d$I, d$P, reference = S_noise, checks = "A3")

  # residual SD is at the injected noise scale, far above 1e-4
  expect_false(res$a3$pass)
  expect_gt(res$a3$residual_sd, 1e-3)
  expect_lt(res$a3$residual_sd, 0.03)
  expect_gt(res$a3$max_resid, 0.01)
  expect_equal(res$verdict, "fail")
})

test_that("A3: value-identical reference columns force failure even with a perfect fit", {
  d <- deposit_fixture()
  n <- nrow(d$I)
  m <- ncol(d$I)
  ## Duplicate a sample column in both the intensity and the reference matrix,
  ## with its additive column effect made equal, so the additive identity is
  ## still exact while the reference has two value-identical columns.
  I_dup <- d$I
  I_dup[, 5] <- I_dup[, 3]
  set.seed(9)
  cc <- stats::rnorm(m)
  cc[5] <- cc[3]
  S_dup <- log2(I_dup) + stats::rnorm(n) + rep(cc, each = n)

  res <- verify_deposit(I_dup, d$P, reference = S_dup, checks = "A3")

  # the additive fit is exact ...
  expect_lt(res$a3$residual_sd, 1e-10)
  # ... but two reference columns are identical, so the check must fail
  expect_false(res$a3$pass)
  expect_length(res$a3$value_identical_columns, 1)
  expect_equal(res$a3$value_identical_columns[[1]], c(3L, 5L))
  expect_equal(res$verdict, "fail")
})

test_that("A3: a row permutation in intensity breaks the identity", {
  d <- deposit_fixture()
  set.seed(301)
  perm <- sample.int(nrow(d$I))

  res <- verify_deposit(d$I[perm, , drop = FALSE], d$P,
                        reference = d$S, checks = "A3")

  # log2(I) is no longer aligned with S, so the residual is large
  expect_false(res$a3$pass)
  expect_gt(res$a3$residual_sd, 0.1)
  expect_equal(res$verdict, "fail")
})

test_that("A3: tolerance_digits controls the residual threshold", {
  d <- deposit_fixture()
  n <- nrow(d$I)
  m <- ncol(d$I)
  set.seed(302)
  S5 <- d$S + matrix(stats::rnorm(n * m, 0, 0.05), n, m)

  # 0.05 noise is inside a 0.1 threshold
  loose <- verify_deposit(d$I, d$P, reference = S5, checks = "A3",
                          tolerance_digits = 1)
  expect_true(loose$a3$pass)
  expect_lt(loose$a3$residual_sd, 0.1)

  # the default 1e-4 threshold rejects the same matrix
  strict <- verify_deposit(d$I, d$P, reference = S5, checks = "A3")
  expect_false(strict$a3$pass)
})


## ---------------------------------------------------------------------------
## Permutation recovery failure modes
## ---------------------------------------------------------------------------

test_that("permutation recovery reports why pi cannot be solved", {
  d <- deposit_fixture()

  # duplicated row tuples make the permutation non-identifiable
  dup <- d$I
  dup[2, ] <- dup[1, ]
  r1 <- verify_deposit(dup, d$P, reference_intensity = d$I, checks = "A2")
  expect_false(r1$permutation_recoverable)
  expect_match(r1$permutation$reason, "duplicate")

  # row tuples of different widths cannot be compared
  r2 <- verify_deposit(d$I, d$P, reference_intensity = d$I[, 1:5, drop = FALSE],
                       checks = "A2")
  expect_false(r2$permutation_recoverable)
  expect_match(r2$permutation$reason, "column")
  # A2 itself still works because it only needs the row means
  expect_true(r2$a2$pass)
})


## ---------------------------------------------------------------------------
## Input validation and the print method
## ---------------------------------------------------------------------------

test_that("input validation stops with explicit messages", {
  d <- deposit_fixture()

  # intensity is mandatory
  expect_error(verify_deposit(), "intensity")
  # detection p-values and reference must have the intensity dimensions
  expect_error(verify_deposit(d$I, d$P[1:100, , drop = FALSE]),
               "same dimensions")
  expect_error(verify_deposit(d$I, d$P, reference = d$S[, 1:5, drop = FALSE]),
               "same dimensions")
  # the fingerprint reference must have one row per probe
  expect_error(
    verify_deposit(d$I, d$P, reference_intensity = d$I[1:100, , drop = FALSE]),
    "same number of rows")
  # probe IDs must describe every row
  expect_error(verify_deposit(d$I, d$P, probe_ids = rep("p", 10)),
               "probe_ids")
  # A1 cannot run without detection p-values, and must say so
  expect_error(verify_deposit(d$I), "A1")
  # unknown check names and out-of-range thresholds are rejected
  expect_error(verify_deposit(d$I, d$P, checks = "A9"))
  expect_error(verify_deposit(d$I, d$P, min_rho = 1.5), "min_rho")
  # non-finite values are not silently propagated
  I_na <- d$I
  I_na[1, 1] <- NA_real_
  expect_error(verify_deposit(I_na, d$P), "non-finite")
  # A3 needs a positive matrix because it uses log2()
  I_zero <- d$I
  I_zero[1, 1] <- 0
  expect_error(verify_deposit(I_zero, d$P, reference = d$S, checks = "A3"),
               "positive")
})

test_that("print.cade_deposit prints a table-style report", {
  d <- deposit_fixture()
  res <- verify_deposit(d$I, d$P, reference = d$S, reference_intensity = d$I)

  out <- capture.output(print(res))

  # the header names the object and its dimensions
  expect_true(any(grepl("cade deposit verification", out)))
  expect_true(any(grepl("500 probes x 10 samples", out)))
  # the verdict is stated in the first lines
  expect_true(any(grepl("verdict: OK", out)))
  # every check has a row in the table
  expect_true(any(grepl("A1", out)))
  expect_true(any(grepl("A2", out)))
  expect_true(any(grepl("A3", out)))
  # and the printable summary follows the table
  expect_true(any(grepl("沉积矩阵核验报告", out)))
})
