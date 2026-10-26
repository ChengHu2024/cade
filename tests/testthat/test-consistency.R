test_that("marker_consistency recovers a known ratio", {
  set.seed(41)
  g <- paste0("g", 1:400)
  comp <- setNames(rnorm(400, 0, 0.3), g)
  mk <- list(A = g[1:10], B = g[11:20], C = g[21:30])
  comp[mk$A] <- rep(1.2, 10)
  comp[mk$B] <- rep(-0.8, 10)
  comp[mk$C] <- rep(2.0, 10)
  obs <- comp
  obs[mk$A] <- comp[mk$A] * 1.0     # ratio 1.0
  obs[mk$B] <- comp[mk$B] * 0.5     # ratio 0.5
  obs[mk$C] <- comp[mk$C] * 3.0     # ratio 3.0

  r <- marker_consistency(obs, comp, mk)
  expect_equal(nrow(r), 3)
  expect_equal(r$ratio_median[r$celltype == "A"], 1.0)
  expect_equal(r$ratio_median[r$celltype == "B"], 0.5)
  expect_equal(r$ratio_median[r$celltype == "C"], 3.0)
  expect_equal(r$verdict[r$celltype == "A"], "pass")
  expect_equal(r$verdict[r$celltype == "B"], "pass")   # 0.5 is the boundary
  expect_equal(r$verdict[r$celltype == "C"], "off")
})

test_that("marker_consistency flags sign flips", {
  set.seed(42)
  g <- paste0("g", 1:200)
  comp <- setNames(rnorm(200, 0, 0.3), g)
  mk <- list(A = g[1:8])
  comp[mk$A] <- rep(1.5, 8)
  obs <- comp
  obs[mk$A] <- -0.4 * comp[mk$A]      # opposite sign -> ratio negative
  r <- marker_consistency(obs, comp, mk)
  expect_equal(r$verdict, "sign_flipped")
  expect_lt(r$ratio_median, 0)
})

test_that("marker_consistency skips cell types with too few markers", {
  set.seed(43)
  g <- paste0("g", 1:100)
  comp <- setNames(rnorm(100, 0, .3), g)
  mk <- list(A = g[1:10], Tiny = c("g1"))
  r <- marker_consistency(comp, comp, mk)
  expect_equal(nrow(r), 1)
  expect_equal(r$celltype, "A")
})

test_that("marker_consistency excludes unstable small comp terms", {
  set.seed(44)
  g <- paste0("g", 1:100)
  ## six markers, of which five have a composition term near zero
  comp <- setNames(c(rep(0.001, 5), runif(95, 0.4, 1)), g)
  obs  <- comp * 1.0
  mk <- list(A = g[1:6])
  ## with a 0.1 floor only one marker survives, which is below min_markers
  expect_error(marker_consistency(obs, comp, mk, min_abs_comp = 0.1),
               "enough usable markers")
  ## without the floor all six are usable
  r <- marker_consistency(obs, comp, mk, min_abs_comp = 0)
  expect_equal(nrow(r), 1)
  expect_equal(r$n_markers, 6)
})

test_that("marker_consistency validates input", {
  expect_error(marker_consistency(c(1, 2), c(1, 2)),
               "named vectors")
  expect_error(marker_consistency(setNames(1:5, letters[1:5]),
                                  setNames(1:5, letters[1:5]),
                                  unname(list(a = "a"))),
               "named list")
})

test_that("cross_cohort_agreement recovers expected correlation", {
  set.seed(51)
  g <- paste0("g", 1:800)
  comp <- setNames(rnorm(800, 0, .4), g)
  o1 <- setNames(comp + rnorm(800, 0, .05), g)
  o2 <- setNames(comp + rnorm(800, 0, .05), g)
  res <- cross_cohort_agreement(list(A = o1, B = o2), list(A = comp, B = comp))
  expect_equal(nrow(res$obs), 1)
  expect_gt(res$obs$pearson, 0.9)
  expect_gt(res$obs$sign_agree, 90)
  ## the composition term is identical in both cohorts
  expect_equal(res$comp$pearson, 1)
})

test_that("cross_cohort_agreement handles more than two cohorts", {
  set.seed(52)
  g <- paste0("g", 1:300)
  o <- list(A = setNames(rnorm(300), g),
            B = setNames(rnorm(300), g),
            C = setNames(rnorm(300), g))
  res <- cross_cohort_agreement(o)
  expect_equal(nrow(res$obs), 3)       # three pairs
  expect_null(res$comp)
})

test_that("cross_cohort_agreement validates input", {
  expect_error(cross_cohort_agreement(list(A = 1:10)), "at least two")
  expect_error(cross_cohort_agreement(list(1:10, 1:10)), "must be named")
  expect_error(cross_cohort_agreement(list(A = setNames(1:100, paste0("g", 1:100)),
                                           B = setNames(1:100, paste0("g", 1:100))),
                                      list(A = 1)), "same length")
})

test_that("CADE_markers is a named list of gene vectors", {
  expect_type(CADE_markers, "list")
  expect_true(all(nzchar(names(CADE_markers))))
  expect_true(all(vapply(CADE_markers, is.character, logical(1))))
  expect_gte(min(vapply(CADE_markers, length, integer(1))), 4)
})
