test_that("decompose_composition rejects genes absent from mu", {
  mu <- matrix(c(10, 1, 5, 2, 8, 3), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2", "g3")))
  fr <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4, byrow = TRUE,
               dimnames = list(paste0("s", 1:4), c("A", "B")))
  obs <- c(g1 = 1, g2 = -1, g3 = 0.5, gX = 2, gY = -2, gZ = 1)

  ## Regression: `obs - pred` used to align by POSITION, so extra genes in
  ## `obs` were silently recycled and every residual came out wrong.
  expect_error(
    decompose_composition(obs, fr, mu, c(FALSE, FALSE, TRUE, TRUE),
                          validate = FALSE),
    "absent from `mu`")
})

test_that("decompose_composition aligns obs and pred by name", {
  mu <- matrix(c(10, 1, 5, 2, 8, 3), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2", "g3")))
  fr <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4, byrow = TRUE,
               dimnames = list(paste0("s", 1:4), c("A", "B")))
  ## order deliberately reversed relative to mu
  obs <- c(g3 = 0.5, g1 = 1, g2 = -1)
  d <- decompose_composition(obs, fr, mu, c(FALSE, FALSE, TRUE, TRUE),
                             validate = FALSE)
  ## residual must equal obs - pred, gene by gene
  want <- obs[d$gene] - d$pred
  expect_equal(d$resid, unname(want), tolerance = 1e-12)
  expect_equal(d$gene, names(obs))
})

test_that("decompose_composition rejects is_case with NA", {
  mu <- matrix(c(10, 1, 5, 2), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  fr <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4, byrow = TRUE,
               dimnames = list(paste0("s", 1:4), c("A", "B")))
  expect_error(
    decompose_composition(c(g1 = 1, g2 = 2), fr, mu,
                          c(TRUE, NA, FALSE, TRUE), validate = FALSE),
    "contains NA")
})

test_that("decompose_composition rejects an all-NA obs", {
  mu <- matrix(c(10, 1, 5, 2), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  fr <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4, byrow = TRUE,
               dimnames = list(paste0("s", 1:4), c("A", "B")))
  expect_error(
    decompose_composition(c(g1 = NA, g2 = NA), fr, mu,
                          c(FALSE, FALSE, TRUE, TRUE), validate = FALSE),
    "entirely NA")
})

test_that("decompose_composition rejects duplicated gene names", {
  mu <- matrix(c(10, 1, 5, 2), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  fr <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4, byrow = TRUE,
               dimnames = list(paste0("s", 1:4), c("A", "B")))
  expect_error(
    decompose_composition(c(g1 = 1, g1 = 2), fr, mu,
                          c(FALSE, FALSE, TRUE, TRUE), validate = FALSE),
    "duplicated")
})

test_that("decompose_composition warns on a cell type absent from a group", {
  mu <- matrix(c(10, 1, 5, 2), nrow = 2,
               dimnames = list(c("A", "B"), c("g1", "g2")))
  fr <- matrix(c(1, 0, 1, 0, 1, 0, 1, 0), nrow = 4, byrow = TRUE,
               dimnames = list(paste0("s", 1:4), c("A", "B")))
  expect_warning(
    decompose_composition(c(g1 = 1, g2 = 2), fr, mu,
                          c(FALSE, FALSE, TRUE, TRUE), validate = FALSE),
    "zero mean fraction")
})