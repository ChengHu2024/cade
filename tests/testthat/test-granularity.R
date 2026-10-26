# ============================================================================
# Tests for the crossed granularity experiment
#
# The fixture is deliberately small (10 cell types, 600 genes, 20 samples) so
# the whole file runs in a couple of seconds. The marker blocks make the
# aggregation weights matter: merging Fibroblast and VSMC is the coarsening
# that carries the largest composition change in vascular tissue, so the four
# weight conventions cannot give the same correlation there.
# ============================================================================

gran_fixture <- function(n_gene = 600, n_sample = 20) {
  set.seed(11)
  types <- c("B_cell", "Endothelial", "Fibroblast", "Macrophage", "Mast",
             "Monocyte", "NK_cell", "Plasma", "T_cell", "VSMC")
  mu <- matrix(1, length(types), n_gene,
               dimnames = list(types, paste0("g", seq_len(n_gene))))

  ## 40 marker genes per cell type, strongly expressed, plus a shared block.
  for (i in seq_along(types)) {
    idx <- ((i - 1) * 40 + 1):(i * 40)
    mu[i, idx] <- 40 + 10 * i
  }
  ## Different total RNA per type (so pool_rna is not equal to the others).
  mu <- mu * seq(0.6, 2.4, length.out = length(types))
  attr(mu, "cell_counts") <- stats::setNames(
    c(120, 3000, 900, 2200, 400, 700, 150, 500, 1600, 5000), types)
  attr(mu, "n_cells") <- 300L

  frac <- matrix(stats::runif(n_sample * length(types), 0.01, 1),
                 n_sample, length(types),
                 dimnames = list(paste0("s", seq_len(n_sample)), types))
  is_case <- rep(c(FALSE, TRUE), each = n_sample / 2)
  ## Vascular contrast: VSMC collapses and Fibroblast expands in cases.
  frac[is_case, "VSMC"] <- frac[is_case, "VSMC"] * 0.15
  frac[is_case, "Fibroblast"] <- frac[is_case, "Fibroblast"] * 4
  frac <- frac / rowSums(frac)
  list(mu = mu, frac = frac, is_case = is_case, types = types)
}

test_that("the crossed design runs on synthetic data and returns its contract", {
  fx <- gran_fixture()
  g <- granularity_crossed(fx$mu, fx$frac, fx$is_case,
                           levels = c(6, 10), seeds = 1:2)

  expect_s3_class(g, "cade_granularity")
  expect_true(all(c("raw", "matrix_r", "sd_by_ref", "identifiable",
                    "summary", "level_groups", "level_maps",
                    "merge_rules_used", "r_definition", "proxies") %in%
                    names(g)))
  expect_identical(g$r_definition,
                   "cor(comp, obs) without deconvolution")
  expect_type(g$merge_rules_used, "list")
  expect_setequal(names(g$merge_rules_used), c("6", "10"))
  expect_s3_class(g$raw, "data.frame")
  expect_identical(names(g$raw),
                   c("true_level", "ref_level", "weight", "seed", "r"))
  expect_true(all(c("true_level", "ref_level", "weight", "seed", "r") %in%
                    names(g$raw)))
  expect_true(is.numeric(g$raw$r))
  expect_true(mean(is.finite(g$raw$r)) > 0.9)

  expect_true(is.matrix(g$matrix_r))
  expect_identical(dim(g$matrix_r), c(2L, 2L))
  expect_identical(rownames(g$matrix_r), c("6", "10"))
  expect_identical(colnames(g$matrix_r), c("6", "10"))

  expect_named(g$sd_by_ref, c("6", "10"))
  expect_true(all(is.finite(g$sd_by_ref)))
  expect_named(g$identifiable, c("6", "10"))
  expect_type(g$identifiable, "logical")
  expect_true(is.character(g$summary) && length(g$summary) > 0)
})

test_that("raw has exactly one row per fit", {
  fx <- gran_fixture()
  levels <- c(6, 10)
  weights <- c("atlas_cells", "pool_cells", "pool_rna", "equal")
  seeds <- 1:3
  g <- granularity_crossed(fx$mu, fx$frac, fx$is_case, levels = levels,
                           weights = weights, seeds = seeds)

  expect_equal(nrow(g$raw),
               length(levels)^2 * length(weights) * length(seeds))
  expect_identical(nrow(g$raw), 2L * 2L * 4L * 3L)
  expect_identical(sort(unique(g$raw$weight)), sort(weights))
  expect_identical(sort(unique(g$raw$seed)), seeds)
  expect_equal(sort(unique(g$raw$true_level)), sort(levels))
  expect_equal(sort(unique(g$raw$ref_level)), sort(levels))
  expect_identical(nrow(unique(g$raw[, c("true_level", "ref_level",
                                        "weight", "seed")])),
                   nrow(g$raw))
})

test_that("the four aggregation weights are not interchangeable", {
  fx <- gran_fixture()
  g <- granularity_crossed(fx$mu, fx$frac, fx$is_case,
                           levels = c(6, 10), seeds = 1L)

  ## Same truth, same reference level, same seed: only the weight changes.
  sub <- g$raw[g$raw$true_level == 10 & g$raw$ref_level == 6, ]
  expect_identical(nrow(sub), 4L)
  expect_identical(length(unique(round(sub$r, 10))), 4L)
  expect_gt(diff(range(sub$r)), 1e-3)

  ## At the finest reference level no rows are merged, so the weight cannot
  ## change the answer for a fixed truth level; this guards against a weight
  ## accidentally leaking into the composition rather than the reference.
  for (tl in unique(g$raw$true_level)) {
    vals <- g$raw$r[g$raw$ref_level == 10 & g$raw$true_level == tl]
    expect_identical(length(unique(round(vals, 10))), 1L)
  }
})

test_that("the merged reference has the requested number of compartments", {
  fx <- gran_fixture()
  g <- granularity_crossed(fx$mu, fx$frac, fx$is_case,
                           levels = c(6, 10), seeds = 1L)

  expect_identical(length(g$level_groups[["6"]]), 6L)
  expect_identical(length(g$level_groups[["10"]]), 10L)
  expect_identical(length(unique(g$level_groups[["6"]])), 6L)
  m6 <- g$level_maps[["6"]]
  expect_true("Fibroblast" %in% names(m6))
  expect_true("VSMC" %in% names(m6))
  ## Fibroblast and VSMC form the coarse Stromal compartment at K = 6.
  expect_identical(unname(m6[["Fibroblast"]]), unname(m6[["VSMC"]]))
  mr6 <- g$merge_rules_used[["6"]]
  expect_identical(sort(mr6[["Stromal"]]), sort(c("Fibroblast", "VSMC")))
  expect_identical(mr6[["B_plasma"]], c("B_cell", "Plasma"))
  ## The 10-type input is explicitly flagged as not the 15-type ladder.
  expect_true(any(grepl("merge ladder adapted for 10 input cell type",
                        g$summary)))

  ## Three levels: 15 -> 10 -> 6 is exercised with a synthetic 15-type atlas.
  types15 <- c(fx$types, "Alveolar_epi", "Airway_epi", "Dendritic",
               "Neutrophil", "Pericyte")
  set.seed(3)
  mu15 <- matrix(stats::rgamma(length(types15) * 200, shape = 2, rate = 0.4),
                 length(types15), 200,
                 dimnames = list(types15, paste0("g", 1:200)))
  attr(mu15, "cell_counts") <- stats::setNames(
    c(120, 3000, 900, 2200, 400, 700, 150, 500, 1600, 5000,
      800, 700, 600, 300, 250), types15)
  frac15 <- matrix(stats::runif(10 * length(types15), 0.01, 1),
                   10, length(types15),
                   dimnames = list(paste0("s", 1:10), types15))
  frac15 <- frac15 / rowSums(frac15)
  g15 <- granularity_crossed(mu15, frac15, rep(c(FALSE, TRUE), each = 5),
                             levels = c(15, 10, 6), seeds = 1L)
  expect_identical(vapply(g15$level_groups, length, integer(1)),
                   c(`6` = 6L, `10` = 10L, `15` = 15L))
  expect_identical(dim(g15$matrix_r), c(3L, 3L))
  expect_identical(nrow(g15$raw), 3L * 3L * 4L * 1L)
  ## Fifteen input types: no "adapted ladder" warning.
  expect_false(any(grepl("merge ladder adapted", g15$summary)))
})

test_that("unknown cell-type names stop with a clear error", {
  fx <- gran_fixture()
  bad_mu <- rbind(fx$mu, `Weird_type` = 1)
  bad_frac <- cbind(fx$frac, Weird_type = 0.01)
  bad_frac <- bad_frac / rowSums(bad_frac)

  expect_error(
    granularity_crossed(bad_mu, bad_frac, fx$is_case,
                        levels = c(6, 10), seeds = 1L),
    "not covered by the built-in merge ladder"
  )
  expect_error(
    granularity_crossed(bad_mu, bad_frac, fx$is_case,
                        levels = c(6, 10), seeds = 1L),
    "Weird_type"
  )

  ## The identity level does not merge anything, so it is allowed.
  expect_error(
    granularity_crossed(bad_mu, bad_frac, fx$is_case,
                        levels = 11, seeds = 1L),
    NA
  )

  ## A custom map removes the ambiguity.
  custom <- stats::setNames(rep(paste0("C", 1:6), length.out = 11),
                            rownames(bad_mu))
  expect_error(
    granularity_crossed(bad_mu, bad_frac, fx$is_case, levels = 6,
                        seeds = 1L, level_map = custom),
    NA
  )
})

test_that("the print method reports the design and the verdict", {
  fx <- gran_fixture()
  g <- granularity_crossed(fx$mu, fx$frac, fx$is_case,
                           levels = c(6, 10), seeds = 1L)
  expect_error(print(g), NA)
  expect_output(print(g), "CADE crossed granularity experiment")
  expect_output(print(g), "reference-level spread")
  expect_invisible(print(g))
})

test_that("the experiment is reproducible from the seeds", {
  fx <- gran_fixture()
  g1 <- granularity_crossed(fx$mu, fx$frac, fx$is_case,
                            levels = c(6, 10), seeds = 1:2)
  g2 <- granularity_crossed(fx$mu, fx$frac, fx$is_case,
                            levels = c(6, 10), seeds = 1:2)
  expect_equal(g1$raw$r, g2$raw$r)
})

test_that("input validation rejects structural mistakes", {
  fx <- gran_fixture()

  ## one cell type removed from `frac`
  expect_error(
    granularity_crossed(fx$mu, fx$frac[, -1], fx$is_case,
                        levels = c(6, 10), seeds = 1L),
    "cell types"
  )
  ## a single group
  expect_error(
    granularity_crossed(fx$mu, fx$frac, rep(TRUE, nrow(fx$frac)),
                        levels = c(6, 10), seeds = 1L),
    "TRUE and FALSE"
  )
  ## more compartments than cell types
  expect_error(
    granularity_crossed(fx$mu, fx$frac, fx$is_case, levels = 99,
                        seeds = 1L),
    "exceed"
  )
  ## unknown weight
  expect_error(
    granularity_crossed(fx$mu, fx$frac, fx$is_case, levels = c(6, 10),
                        weights = "magic", seeds = 1L),
    "unknown weight"
  )
  ## non-positive sd threshold
  expect_error(
    granularity_crossed(fx$mu, fx$frac, fx$is_case, levels = c(6, 10),
                        seeds = 1L, sd_threshold = 0),
    "sd_threshold"
  )
})

test_that("n_cores > 1 falls back to serial on Windows", {
  skip_if(.Platform$OS.type != "windows")
  fx <- gran_fixture()
  expect_warning(
    granularity_crossed(fx$mu, fx$frac, fx$is_case, levels = 6,
                        seeds = 1L, n_cores = 2L),
    "Windows"
  )
})
