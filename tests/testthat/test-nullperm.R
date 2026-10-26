# ============================================================================
# Tests for null_permutation() and print.cade_nullperm().
#
# Every case is synthetic and small: the full 20000-permutation reproduction on
# GSE232911 lives in the analysis repository, not in the package test suite.
#
# A note on K. With very few cell types the permutation null has only K! distinct
# atoms, so p is coarse even when n_perm is large. Tests that need a smooth null
# use K = 5; the rest use the default K = 3.
# ============================================================================

## Deterministic toy problem: K cell types, G genes, n samples, both groups.
## `obs` is built from the true composition term plus noise so that r_obs is
## positive and well away from 0 and 1.
make_toy <- function(seed = 1, G = 40L, K = 3L, n = 12L, noise = 0.2) {
  set.seed(seed)
  mu <- matrix(stats::runif(K * G, 1, 50), K, G,
               dimnames = list(paste0("ct", seq_len(K)),
                               paste0("g", seq_len(G))))
  fr <- matrix(stats::runif(n * K), n, K,
               dimnames = list(paste0("s", seq_len(n)),
                               paste0("ct", seq_len(K))))
  fr <- fr / rowSums(fr)
  is_case <- rep(c(FALSE, TRUE), length.out = n)
  comp <- comp_mean_of_log(mu, fr, is_case)
  obs <- stats::setNames(comp + stats::rnorm(G, 0, noise),
                         paste0("g", seq_len(G)))
  list(mu = mu, frac = fr, is_case = is_case, obs = obs)
}

## The documented mean-of-log composition term, written independently of the
## implementation: log2 per sample, mean within group, difference of the means.
comp_mean_of_log <- function(mu, fr, is_case) {
  cs <- t(fr[is_case, , drop = FALSE])
  ct <- t(fr[!is_case, , drop = FALSE])
  rowMeans(log2(pmax(t(mu) %*% cs, 1e-12))) -
    rowMeans(log2(pmax(t(mu) %*% ct, 1e-12)))
}

# ---------------------------------------------------------------------------
# 1. r_obs is the correlation of the composition term and obs
# ---------------------------------------------------------------------------

test_that("r_obs equals the direct correlation of comp and obs", {
  d <- make_toy()
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case,
                          n_perm = 5L, seed = 1L, progress = FALSE)
  comp <- comp_mean_of_log(d$mu, d$frac, d$is_case)

  expect_equal(out$r_obs, stats::cor(comp, d$obs), tolerance = 1e-12)
  expect_identical(out$convention, "mean_of_log")
  expect_identical(out$stat, "pearson")
  expect_identical(out$mu_orientation, "celltypes_x_genes")
  expect_equal(out$n, length(d$obs))
  expect_equal(out$K, nrow(d$mu))
})

test_that("r_obs under log_of_mean follows the group-mean definition", {
  d <- make_toy()
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 2L,
                          convention = "log_of_mean", progress = FALSE)
  ## package convention: mu is K x G, so the linear predictor is t(mu) %*% f
  fc <- colMeans(d$frac[d$is_case, , drop = FALSE])
  fn <- colMeans(d$frac[!d$is_case, , drop = FALSE])
  comp <- log2(pmax(as.numeric(t(d$mu) %*% fc), 1e-12)) -
    log2(pmax(as.numeric(t(d$mu) %*% fn), 1e-12))

  expect_equal(out$r_obs, stats::cor(comp, d$obs), tolerance = 1e-12)
  expect_identical(out$convention, "log_of_mean")
})

# ---------------------------------------------------------------------------
# 2. identity permutation, and the direction of a supplied permutation
# ---------------------------------------------------------------------------

test_that("the identity permutation puts r_obs inside the null", {
  d <- make_toy()
  K <- nrow(d$mu)
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, progress = FALSE,
                          perm = matrix(seq_len(K), nrow = 1L))

  expect_equal(out$n_perm, 1L)
  expect_equal(out$null, out$r_obs)
  expect_equal(out$p, 1)
  expect_equal(out$n_ge, 1L)
})

test_that("perm reorders the cell-type axis of mu, not the genes", {
  d <- make_toy(K = 4L)
  pi <- c(2L, 4L, 1L, 3L)                 # 4-cycle: not self-inverse
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, progress = FALSE,
                          perm = matrix(pi, nrow = 1L))
  comp <- comp_mean_of_log(d$mu[pi, , drop = FALSE], d$frac, d$is_case)

  expect_equal(out$null, stats::cor(comp, d$obs), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. the two conventions really differ
# ---------------------------------------------------------------------------

test_that("the two conventions give different answers", {
  ## Random proportions make the linear predictor vary across samples, so the
  ## Jensen gap between mean-of-log and log-of-mean is non-zero (~5e-3 here).
  d <- make_toy()
  a <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 2L,
                        progress = FALSE, convention = "mean_of_log")
  b <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 2L,
                        progress = FALSE, convention = "log_of_mean")

  expect_gt(abs(a$r_obs - b$r_obs), 1e-3)
})

# ---------------------------------------------------------------------------
# 4. shape of the returned object and the Monte Carlo arithmetic
# ---------------------------------------------------------------------------

test_that("the returned object has the documented fields", {
  d <- make_toy()
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 7L,
                          seed = 42L, progress = FALSE)

  expect_s3_class(out, "cade_nullperm")
  expect_named(out, c("r_obs", "null", "null_mean", "null_sd", "null_lo",
                      "null_hi", "z", "n_ge", "p", "p_se", "n_perm", "n",
                      "K", "convention", "stat", "seed", "elapsed",
                      "mu_orientation"))
  expect_length(out$null, 7L)
  expect_equal(out$n_perm, 7L)
  expect_equal(out$seed, 42L)
  expect_true(is.numeric(out$elapsed) && out$elapsed >= 0)

  expect_equal(out$p, (1 + out$n_ge) / (1 + out$n_perm))
  expect_gt(out$p, 0)
  expect_lte(out$p, 1)
  expect_equal(out$p_se, sqrt(out$p * (1 - out$p) / out$n_perm))

  expect_equal(out$null_mean, mean(out$null))
  expect_equal(out$null_sd, stats::sd(out$null))
  expect_equal(out$null_lo, stats::quantile(out$null, 0.025, names = FALSE))
  expect_equal(out$null_hi, stats::quantile(out$null, 0.975, names = FALSE))
  expect_equal(out$z, (out$r_obs - out$null_mean) / out$null_sd)
  expect_equal(out$n_ge, sum(out$null >= out$r_obs))
})

test_that("the same seed gives the same null distribution", {
  d <- make_toy()
  a <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 9L,
                        seed = 99L, progress = FALSE)
  b <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 9L,
                        seed = 99L, progress = FALSE)
  expect_identical(a$null, b$null)
})

test_that("the Monte Carlo near-0.05 warning fires when it should", {
  ## K = 5 so the null is smooth, noise tuned so p lands at 13/201 = 0.0647.
  ## n_perm = 200 < 10000 and p in [0.04, 0.09], so the warning must fire.
  set.seed(2026)
  P <- t(replicate(200L, sample.int(5L, 5L)))
  d <- make_toy(seed = 1L, G = 40L, K = 5L, n = 12L, noise = 0.80)

  expect_warning(
    out <- null_permutation(d$obs, d$frac, d$mu, d$is_case,
                            progress = FALSE, perm = P),
    "p is close to 0.05")
  expect_gte(out$p, 0.04)
  expect_lte(out$p, 0.09)
})

test_that("no near-0.05 warning when the p value is far from 0.05", {
  d <- make_toy()                       # obs tracks comp closely: p is small
  expect_silent(
    null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 50L,
                     seed = 3L, progress = FALSE))
})

# ---------------------------------------------------------------------------
# 5. print method
# ---------------------------------------------------------------------------

test_that("print() works and reports the result", {
  d <- make_toy()
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 5L,
                          seed = 7L, progress = FALSE)

  expect_output(print(out), "cade null permutation")
  expect_output(print(out), "convention mean_of_log")
  expect_output(print(out), "genes")
  expect_output(print(out), "p \\(Monte Carlo SE\\)")
  ## print() returns its argument invisibly (captured so the suite stays quiet)
  txt <- utils::capture.output(v <- withVisible(print(out)))
  expect_false(v$visible)
  expect_identical(v$value, out)
})

# ---------------------------------------------------------------------------
# 6. small hand-worked example: 3 genes, 2 cell types
# ---------------------------------------------------------------------------

test_that("a 3-gene, 2-cell-type example matches a hand calculation", {
  mu <- rbind(A = c(10, 20, 30), B = c(40, 50, 60))
  colnames(mu) <- c("g1", "g2", "g3")
  fr <- rbind(s1 = c(A = 0.2, B = 0.8),
              s2 = c(A = 0.4, B = 0.6),
              s3 = c(A = 0.7, B = 0.3),
              s4 = c(A = 0.9, B = 0.1))
  obs <- c(g1 = 1, g2 = 2, g3 = 3)
  is_case <- c(FALSE, FALSE, TRUE, TRUE)

  ## mu %*% f per sample, written out by hand:
  ##   s1: 0.2*A + 0.8*B = (34, 44, 54)     s2: 0.4*A + 0.6*B = (28, 38, 48)
  ##   s3: 0.7*A + 0.3*B = (19, 29, 39)     s4: 0.9*A + 0.1*B = (13, 23, 33)
  comp <- (log2(c(19, 29, 39)) + log2(c(13, 23, 33))) / 2 -
    (log2(c(34, 44, 54)) + log2(c(28, 38, 48))) / 2

  ## swapping the two cell types puts B on top: mu[c(2, 1), ] %*% f
  ##   s1: 0.2*B + 0.8*A = (16, 26, 36)     s2: (22, 32, 42)
  ##   s3: 0.7*B + 0.3*A = (31, 41, 51)     s4: (37, 47, 57)
  comp_swapped <- (log2(c(31, 41, 51)) + log2(c(37, 47, 57))) / 2 -
    (log2(c(16, 26, 36)) + log2(c(22, 32, 42))) / 2

  out <- null_permutation(obs, fr, mu, is_case, progress = FALSE,
                          perm = rbind(c(1L, 2L), c(2L, 1L)))

  expect_equal(out$r_obs, stats::cor(comp, obs), tolerance = 1e-12)
  expect_equal(out$null, c(stats::cor(comp, obs),
                           stats::cor(comp_swapped, obs)),
               tolerance = 1e-12)
  ## identity is >= r_obs, the swap is not: n_ge = 1 of 2, p = (1 + 1) / 3
  expect_equal(out$n_ge, 1L)
  expect_equal(out$p, 2 / 3)
})

# ---------------------------------------------------------------------------
# 7. input validation and orientation detection
# ---------------------------------------------------------------------------

test_that("misaligned cell types are rejected with a specific message", {
  d <- make_toy()
  bad <- d$mu
  rownames(bad) <- paste0("x", seq_len(nrow(bad)))
  expect_error(
    null_permutation(d$obs, d$frac, bad, d$is_case, n_perm = 2L, progress = FALSE),
    "same set of cell types")

  ## dropping a proportion column leaves rows not summing to 1: the warning is
  ## expected here, only the error matters
  expect_error(
    suppressWarnings(
      null_permutation(d$obs, d$frac[, 1:2, drop = FALSE], d$mu, d$is_case,
                       n_perm = 2L, progress = FALSE)),
    "same set of cell types")
})

test_that("obs must be a named numeric vector", {
  d <- make_toy()
  expect_error(
    null_permutation(unname(d$obs), d$frac, d$mu, d$is_case,
                     n_perm = 2L, progress = FALSE),
    "named")
  expect_error(
    null_permutation(as.character(d$obs), d$frac, d$mu, d$is_case,
                     n_perm = 2L, progress = FALSE),
    "numeric")
})

test_that("a proportions table with a grouping column is rejected with a hint", {
  d <- make_toy()
  df <- as.data.frame(d$frac)
  df$group <- rep("case", nrow(df))
  expect_error(
    null_permutation(d$obs, df, d$mu, d$is_case, n_perm = 2L, progress = FALSE),
    "non-numeric")
})

test_that("is_case must be logical, aligned and two-sided", {
  d <- make_toy()
  expect_error(
    null_permutation(d$obs, d$frac, d$mu, as.integer(d$is_case),
                     n_perm = 2L, progress = FALSE),
    "logical")
  expect_error(
    null_permutation(d$obs, d$frac, d$mu, d$is_case[-1L],
                     n_perm = 2L, progress = FALSE),
    "length")
  expect_error(
    null_permutation(d$obs, d$frac, d$mu, rep(TRUE, nrow(d$frac)),
                     n_perm = 2L, progress = FALSE),
    "control")
})

test_that("mu in genes x cell types layout is detected and transposed", {
  d <- make_toy()
  expect_warning(
    out_t <- null_permutation(d$obs, d$frac, t(d$mu), d$is_case,
                              n_perm = 4L, progress = FALSE),
    "genes x cell types")

  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case,
                          n_perm = 4L, progress = FALSE)
  expect_equal(out_t$r_obs, out$r_obs, tolerance = 1e-12)
  expect_identical(out_t$mu_orientation, "celltypes_x_genes")
})

test_that("genes restricts the test and is validated", {
  d <- make_toy()
  sel <- paste0("g", 1:15)
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 3L,
                          progress = FALSE, genes = sel)
  comp <- comp_mean_of_log(d$mu[, sel, drop = FALSE], d$frac, d$is_case)

  expect_equal(out$n, 15L)
  expect_equal(out$r_obs, stats::cor(comp, d$obs[sel]), tolerance = 1e-12)

  ## fewer than 3 genes cannot give a correlation
  expect_error(
    null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 2L,
                     progress = FALSE, genes = c("g1", "g2")),
    "at least 3")

  ## a gene missing from obs
  expect_error(
    null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 2L,
                     progress = FALSE, genes = c("g1", "g2", "not_a_gene")),
    "names\\(obs\\)")

  ## a gene present in obs but missing from mu
  obs_extra <- c(d$obs, not_in_mu = 1)
  expect_error(
    null_permutation(obs_extra, d$frac, d$mu, d$is_case, n_perm = 2L,
                     progress = FALSE, genes = c("g1", "g2", "not_in_mu")),
    "colnames\\(mu\\)")
})

test_that("bad permutations are rejected", {
  d <- make_toy()                                  # K = 3
  expect_error(
    null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 2L,
                     progress = FALSE, perm = matrix(c(1, 1, 2), nrow = 1L)),
    "not permutations")
  expect_error(
    null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 2L,
                     progress = FALSE, perm = matrix(c(1, 2), nrow = 1L)),
    "has 2 column")
})

# ---------------------------------------------------------------------------
# 8. statistic, cores and progress options
# ---------------------------------------------------------------------------

test_that("spearman ranks before correlating", {
  d <- make_toy()
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 2L,
                          stat = "spearman", progress = FALSE)
  comp <- comp_mean_of_log(d$mu, d$frac, d$is_case)

  expect_identical(out$stat, "spearman")
  expect_equal(out$r_obs, stats::cor(comp, d$obs, method = "spearman"),
               tolerance = 1e-12)
})

test_that("n_cores > 1 falls back with a warning on Windows", {
  skip_if_not(.Platform$OS.type == "windows", "Windows-specific fallback")
  d <- make_toy()
  expect_warning(
    out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 4L,
                            n_cores = 2L, progress = FALSE),
    "mclapply")
  expect_equal(out$n_perm, 4L)
  expect_length(out$null, 4L)
})

test_that("n_cores = 1 (default) runs the same arithmetic", {
  d <- make_toy()
  out <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 4L,
                          n_cores = 1L, progress = FALSE)
  expect_equal(out$n_perm, 4L)
})

test_that("progress printing does not change the result", {
  d <- make_toy()
  quiet <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 6L,
                            seed = 11L, progress = FALSE)
  expect_output(
    loud <- null_permutation(d$obs, d$frac, d$mu, d$is_case, n_perm = 6L,
                             seed = 11L, progress = TRUE),
    "null_permutation")
  expect_identical(quiet$null, loud$null)
})
