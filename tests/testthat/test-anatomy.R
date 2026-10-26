test_that("anatomical gate catches the failures it was built for", {
  mk <- function(mat) { m <- as.matrix(mat); m / rowSums(m) }

  ## 真实失败案例：肺组织 0% 内皮、0% 淋巴
  f_bad <- mk(rbind(
    c(epithelial = .31, endothelial = .000005, fibroblast = .24, T_cell = 0, myeloid = .10, other = .35)
  ))
  expect_error(anatomical_plausibility(f_bad, tissue = "lung"),
               "anatomically implausible")
  expect_error(anatomical_plausibility(f_bad, tissue = "lung"), "endothelium")

  ## 真实失败案例：正常乳腺 0% 脂肪
  f_br <- mk(rbind(
    c(epithelial = .39, adipocyte = 0, endothelial = .05, fibroblast = .43, T_cell = .02, other = .11)
  ))
  expect_error(anatomical_plausibility(f_br, tissue = "breast"), "adipocyte")

  ## 通过的案例
  f_ok <- mk(rbind(
    c(epithelial = .30, endothelial = .12, fibroblast = .25, T_cell = .10, myeloid = .10, other = .13)
  ))
  expect_silent(anatomical_plausibility(f_ok, tissue = "lung"))
})

test_that("action controls the failure mode", {
  f_bad <- cbind(endothelial = 0, fibroblast = 1)
  rules <- list(endo = list(types = "endothelial", min = 0.01))
  expect_error(anatomical_plausibility(f_bad, rules = rules, action = "error"))
  expect_warning(anatomical_plausibility(f_bad, rules = rules, action = "warn"),
                 "anatomically implausible")
  r <- anatomical_plausibility(f_bad, rules = rules, action = "report")
  expect_false(r$pass)
})

test_that("the gate refuses to guess", {
  f <- cbind(a = .5, b = .5)
  expect_error(anatomical_plausibility(f), "supply rules")
  expect_error(anatomical_plausibility(f, rules = list(x = list(types = "a"))), "min")
  expect_error(anatomical_plausibility(f, rules = list(x = list(types = "zzz", min = .1))),
               "none of the named cell types")
})

test_that("rules validate their sample selector", {
  f <- matrix(1 / 3, 3, 3, dimnames = list(c("s1", "s2", "s3"), c("a", "b", "c")))
  expect_error(anatomical_plausibility(f, rules = list(x = list(types = "a", min = .1,
                                                               samples = c(TRUE, FALSE)))),
               "length")
  expect_error(anatomical_plausibility(f, rules = list(x = list(types = "a", min = .1,
                                                               samples = c("nope")))),
               "selects no rows")
  r <- anatomical_plausibility(f, rules = list(x = list(types = "a", min = .1, samples = "s2")),
                               action = "report")
  expect_equal(r$n_samples, 1L)
})

test_that("NA and negative fractions are refused", {
  expect_error(anatomical_plausibility(cbind(a = c(NA, .5), b = c(.5, .5)),
                                       rules = list(x = list(types = "a", min = 0))),
               "NA")
  expect_error(anatomical_plausibility(cbind(a = -.1, b = 1.1),
                                       rules = list(x = list(types = "a", min = 0))),
               "non-negative")
})

test_that("regex type selectors match multiple columns", {
  f <- matrix(1 / 4, 2, 4,
              dimnames = list(NULL, c("endothelial cell of artery",
                                      "vein endothelial cell", "fibroblast", "T cell")))
  r <- anatomical_plausibility(f, rules = list(e = list(types = "endothelial", min = .1)))
  expect_equal(r$min_fraction, .5)
})

test_that("built-in rule sets expose their floors", {
  for (t in c("lung", "breast", "aorta", "blood")) {
    rr <- anatomical_rules(t)
    expect_true(length(rr) >= 2)
    for (x in rr) expect_true(!is.null(x$min) && x$min > 0 && x$min <= 1)
  }
  expect_error(anatomical_rules("nope"), "arg")
})
