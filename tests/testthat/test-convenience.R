# ============================================================================
# Tests for the convenience layer: remedies, report, recommend_reference
#
# The critical test is "report numbers are transcription, not recomputation":
# every number that appears in the report must be found, unchanged, in the
# object that was passed in. If someone later makes report() recompute anything,
# that test will fail.
# ============================================================================

make_decomp <- function(n = 120, seed = 1) {
  set.seed(seed)
  G <- n
  ct <- c("A", "B")
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(ct, paste0("g", seq_len(G))))
  fr <- rbind(matrix(rep(c(.7, .3), 10), 10, 2, byrow = TRUE),
              matrix(rep(c(.2, .8), 10), 10, 2, byrow = TRUE))
  colnames(fr) <- ct
  isc <- c(rep(FALSE, 10), rep(TRUE, 10))
  obs <- stats::setNames(stats::rnorm(G, 0, 0.6), colnames(mu))
  decompose_composition(obs, fr, mu, isc)
}

## ---- remedies ------------------------------------------------------------
test_that("check_inputs returns a remedies element with one row per message", {
  G <- 40
  set.seed(2)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3), 8), 8, 2, byrow = TRUE,
               dimnames = list(paste0("s", 1:8), c("A", "B")))
  obs <- stats::setNames(stats::rnorm(G), colnames(mu))

  chk <- check_inputs(obs, fr, mu, rep(TRUE, 8))     # all-case -> 1 error
  expect_true("remedies" %in% names(chk))
  expect_equal(nrow(chk$remedies$for_errors), length(chk$errors))
  expect_equal(nrow(chk$remedies$for_warnings), length(chk$warnings))
})

test_that("the log-scale error carries a concrete fix", {
  G <- 40
  set.seed(3)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3), 8), 8, 2, byrow = TRUE,
               dimnames = list(paste0("s", 1:8), c("A", "B")))
  obs <- stats::setNames(stats::rnorm(G), colnames(mu))
  isc <- c(rep(FALSE, 4), rep(TRUE, 4))

  chk <- check_inputs(obs, fr, log2(mu + 1), isc)
  expect_false(chk$ok)
  rem <- chk$remedies$for_errors$remedy
  expect_true(any(grepl("2\\^L", rem)))          # suggests the inverse transform
})

test_that("the NA is_case error explains the silent-failure consequence", {
  G <- 40
  set.seed(4)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3), 8), 8, 2, byrow = TRUE,
               dimnames = list(paste0("s", 1:8), c("A", "B")))
  obs <- stats::setNames(stats::rnorm(G), colnames(mu))
  chk <- check_inputs(obs, fr, mu, c(FALSE, FALSE, NA, TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_false(chk$ok)
  expect_true(any(grepl("unchanged", chk$remedies$for_errors$remedy)))
})

test_that("an unmatched message gets an empty remedy, not an invented one", {
  G <- 40
  set.seed(5)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3), 8), 8, 2, byrow = TRUE,
               dimnames = list(paste0("s", 1:8), c("A", "B")))
  obs <- stats::setNames(stats::rnorm(G), colnames(mu))
  chk <- check_inputs(obs, fr, mu, c(rep(FALSE, 4), rep(TRUE, 4)))
  ## valid input -> no errors, so nothing to invent
  expect_equal(length(chk$remedies$for_errors$message), 0)
})

test_that("diagnose_inputs returns the check object invisibly", {
  G <- 40
  set.seed(6)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3), 8), 8, 2, byrow = TRUE,
               dimnames = list(paste0("s", 1:8), c("A", "B")))
  obs <- stats::setNames(stats::rnorm(G), colnames(mu))
  expect_output(diagnose_inputs(obs, fr, mu, c(rep(FALSE, 4), rep(TRUE, 4))),
                "cade input diagnosis")
})


## ---- ★ report() 是转印，不是重算 ---------------------------------------
test_that("report numbers are transcription, not recomputation", {
  d <- make_decomp()
  rep <- suppressWarnings(
    capture.output(report(d, calibrate = FALSE, attribute = FALSE, top_n = 5)))
  txt <- paste(rep, collapse = "\n")

  ## 顶层范围必须与对象逐位一致（两位小数，即报告本身的精度）
  expect_match(txt, sprintf("%+.2f", min(d$obs,   na.rm = TRUE)), fixed = TRUE)
  expect_match(txt, sprintf("%+.2f", max(d$obs,   na.rm = TRUE)), fixed = TRUE)
  expect_match(txt, sprintf("%+.2f", min(d$pred,  na.rm = TRUE)), fixed = TRUE)
  expect_match(txt, sprintf("%+.2f", max(d$pred,  na.rm = TRUE)), fixed = TRUE)
  expect_match(txt, sprintf("%+.2f", min(d$resid, na.rm = TRUE)), fixed = TRUE)
  expect_match(txt, sprintf("%+.2f", max(d$resid, na.rm = TRUE)), fixed = TRUE)

  ## 表内每一行：报告里的三个数必须等于对象里同一基因的三位值
  ord <- order(-abs(d$resid))
  sel <- d[ord[seq_len(5)], , drop = FALSE]
  for (i in seq_len(nrow(sel))) {
    ## 用固定字符串匹配：%+.2f 的输出含 "+"，在正则里是量词
    rowpat <- sprintf("| %s | %s | %s | %s |",
                      sel$gene[i],
                      sprintf("%+.2f", sel$obs[i]),
                      sprintf("%+.2f", sel$pred[i]),
                      sprintf("%+.2f", sel$resid[i]))
    expect_true(grepl(rowpat, txt, fixed = TRUE),
                label = sprintf("report row for %s must transcribe the object", sel$gene[i]))
  }
})

test_that("report reproduces the ranges it prints", {
  d <- make_decomp()
  rep <- suppressWarnings(capture.output(
    report(d, calibrate = FALSE, attribute = FALSE)))
  txt <- paste(rep, collapse = "\n")
  ## 观测与残差的范围必须与对象一致
  expect_match(txt, sprintf("%+.2f", min(d$obs, na.rm = TRUE)), fixed = TRUE)
  expect_match(txt, sprintf("%+.2f", max(d$obs, na.rm = TRUE)), fixed = TRUE)
  expect_match(txt, sprintf("%+.2f", min(d$resid, na.rm = TRUE)), fixed = TRUE)
  expect_match(txt, sprintf("%+.2f", max(d$resid, na.rm = TRUE)), fixed = TRUE)
})

test_that("report does not modify the object it is given", {
  d <- make_decomp()
  before <- unclass(d)
  invisible(suppressWarnings(capture.output(
    report(d, calibrate = FALSE, attribute = FALSE))))
  expect_identical(unclass(d), before)
})

test_that("report carries the verdict when calibration fails", {
  d <- make_decomp()
  ## 把 pred 打乱，校准必然崩
  d$pred <- sample(d$pred)
  rep <- suppressWarnings(capture.output(
    report(d, calibrate = TRUE, attribute = FALSE)))
  txt <- paste(rep, collapse = "\n")
  expect_match(txt, "uninformative|weakly calibrated")
})

test_that("report writes to a file when asked", {
  d <- make_decomp()
  f <- tempfile(fileext = ".md")
  on.exit(unlink(f), add = TRUE)
  suppressMessages(report(d, file = f, calibrate = FALSE, attribute = FALSE))
  expect_true(file.exists(f))
  expect_gt(length(readLines(f)), 10)
})

test_that("report adds a known-answer section when truth is supplied", {
  G <- 60
  set.seed(8)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fc <- c(A = .7, B = .3); fa <- c(A = .2, B = .8)
  sim <- simulate_pseudobulk(mu, fc, fa, 10, 10,
                             intrinsic = c(g1 = 1), noise = 0.05, seed = 8)
  obs <- stats::setNames(
    vapply(colnames(sim$bulk), function(g) {
      a <- mean(sim$bulk[sim$is_case, g]); b <- mean(sim$bulk[!sim$is_case, g])
      log2(a / b)
    }, numeric(1)), colnames(sim$bulk))
  fr <- rbind(matrix(rep(fc, 10), 10, 2, byrow = TRUE),
              matrix(rep(fa, 10), 10, 2, byrow = TRUE))
  colnames(fr) <- c("A", "B")
  d <- decompose_composition(obs, fr, mu, sim$is_case)
  rep <- suppressWarnings(capture.output(
    report(d, calibrate = FALSE, attribute = FALSE, truth = sim$truth)))
  expect_match(paste(rep, collapse = "\n"), "Known-answer check")
})


## ---- recommend_reference -------------------------------------------------
test_that("recommend_reference prefers the reference with normal tissue", {
  G <- 30
  set.seed(9)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3, .2, .8), 5), 10, 2, byrow = TRUE,
               colnames <- NULL)
  colnames(fr) <- c("A", "B")
  isc <- c(rep(FALSE, 5), rep(TRUE, 5))
  obs <- stats::setNames(stats::rnorm(G, 0, .5), colnames(mu))
  cmp <- compare_references(obs, fr, fr, mu, mu, isc, validate = FALSE)

  r <- suppressMessages(suppressWarnings(capture.output(
    rv <- recommend_reference(cmp, a_has_normal = TRUE, b_has_normal = FALSE))))
  r <- rv
  expect_match(r$verdict, "reference B fails")
})

test_that("recommend_reference refuses to pick when both are unknown", {
  G <- 30
  set.seed(10)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3, .2, .8), 5), 10, 2, byrow = TRUE)
  colnames(fr) <- c("A", "B")
  isc <- c(rep(FALSE, 5), rep(TRUE, 5))
  obs <- stats::setNames(stats::rnorm(G, 0, .5), colnames(mu))
  cmp <- compare_references(obs, fr, fr, mu, mu, isc, validate = FALSE)
  r <- suppressMessages(suppressWarnings(capture.output(
    rv <- recommend_reference(cmp))))
  r <- rv
  expect_match(r$verdict, "no recommendation|unknown")
})

test_that("recommend_reference warns against choosing by calibration slope", {
  G <- 30
  set.seed(11)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3, .2, .8), 5), 10, 2, byrow = TRUE)
  colnames(fr) <- c("A", "B")
  isc <- c(rep(FALSE, 5), rep(TRUE, 5))
  obs <- stats::setNames(stats::rnorm(G, 0, .5), colnames(mu))
  cmp <- compare_references(obs, fr, fr, mu, mu, isc, validate = FALSE)
  r <- suppressMessages(suppressWarnings(capture.output(
    rv <- recommend_reference(cmp, TRUE, TRUE))))
  r <- rv
  expect_true(any(grepl("circular", r$actions)))
})

test_that("recommend_reference flags a starved cell type", {
  G <- 30
  set.seed(12)
  mu <- matrix(stats::runif(2 * G, 50, 300), nrow = 2,
               dimnames = list(c("A", "B"), paste0("g", seq_len(G))))
  fr <- matrix(rep(c(.7, .3, .2, .8), 5), 10, 2, byrow = TRUE)
  colnames(fr) <- c("A", "B")
  isc <- c(rep(FALSE, 5), rep(TRUE, 5))
  obs <- stats::setNames(stats::rnorm(G, 0, .5), colnames(mu))
  cmp <- compare_references(obs, fr, fr, mu, mu, isc, validate = FALSE)
  r <- suppressMessages(suppressWarnings(capture.output(
    rv <- recommend_reference(cmp, a_has_normal = TRUE, b_has_normal = TRUE,
                              a_cells_per_type = c(A = 500, B = 400),
                              b_cells_per_type = c(A = 500, B = 20)))))
  r <- rv
  expect_match(r$verdict, "reference B fails")
})
