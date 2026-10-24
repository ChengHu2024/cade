# ============================================================================
# cade: automatic input validation and self-test
#
# Why this file exists
# --------------------
# Every guard-rail in 01-04 catches a *specific* arithmetic mistake. They do
# not catch a wrong *input*. A CP10k matrix passed where raw counts were
# expected, or a reference that is already log-scaled, produces a number that
# looks entirely ordinary and is wrong.
#
# The three layers here are ordered by how early they can act:
#
#   check_inputs()  - structural and semantic checks on obs/frac/mu. Called
#                     automatically by decompose_composition(), so a caller
#                     cannot forget it.
#   self_test()     - a known-answer simulation. Runs in about a second. If the
#                     package cannot recover a change it injected itself, no
#                     result computed with it is trustworthy.
#
# Nothing here changes the arithmetic. It only refuses, warns, or attests.
# ============================================================================

#' Validate the inputs to a composition attribution
#'
#' Performs structural and semantic checks on the three inputs before any
#' arithmetic is done. Called automatically by [decompose_composition()].
#'
#' @param obs Named numeric vector of observed log2 fold changes.
#' @param frac Fraction matrix (samples x cell types).
#' @param mu Per-cell mean expression (cell types x genes), linear scale.
#' @param is_case Logical vector over rows of `frac`.
#' @param min_genes Minimum number of genes to accept, default 10.
#' @param strict If `TRUE`, warnings are promoted to errors. Default `FALSE`.
#'
#' @return A list with `ok` (logical), `errors` (character), `warnings`
#'   (character) and `details` (a list of measured quantities). Never silently
#'   returns: if `ok` is `FALSE` the caller should stop.
#'
#' @details
#' The semantic checks that matter most:
#'
#' * **Compressed `mu`.** If the largest per-cell mean is below 30, `mu` is
#'   probably log-scaled or otherwise transformed. A linear CP10k matrix
#'   normally reaches several hundred for at least one gene (mitochondrial or
#'   ribosomal transcripts are always present). A log-scaled `mu` makes the
#'   ratio inside `log2()` meaningless while producing a plausible number.
#'
#' * **Integer `mu` above 1000.** Raw counts rather than per-cell means. Using
#'   counts instead of CP10k shifts the ratio between cell types by their
#'   library-size ratio, which is not a property of the biology.
#'
#' * **Duplicate names.** `match()` silently takes the first, so a duplicated
#'   cell-type label merges two populations without notice.
#'
#' * **Non-normalised fractions.** Harmless in practice because
#'   [compose_prediction()] renormalises, but a row summing to 0.3 or 3.0
#'   usually means the wrong columns were selected.
#'
#' @examples
#' mu <- matrix(c(500, 20, 300, 40), nrow = 2,
#'              dimnames = list(c("A", "B"), c("g1", "g2")))
#' fr <- matrix(c(0.7, 0.3, 0.7, 0.3, 0.3, 0.7, 0.3, 0.7), nrow = 4,
#'              byrow = TRUE, dimnames = list(paste0("s", 1:4), c("A", "B")))
#' ok <- check_inputs(c(g1 = 0.2, g2 = -0.1), fr, mu,
#'                    is_case = c(FALSE, FALSE, TRUE, TRUE))
#' ok$ok
#' @export
check_inputs <- function(obs, frac, mu, is_case, min_genes = 10,
                         strict = FALSE) {

  errors <- character(0)
  warns  <- character(0)
  add_err <- function(...) errors <<- c(errors, sprintf(...))
  add_wrn <- function(...) warns  <<- c(warns,  sprintf(...))

  ## ---- obs -------------------------------------------------------------
  if (!is.numeric(obs)) add_err("`obs` is not numeric.")
  if (is.null(names(obs))) add_err("`obs` has no names; gene identity is required.")
  if (length(obs) < min_genes)
    add_err("`obs` has %d genes; at least %d are required.", length(obs), min_genes)
  if (any(!is.finite(obs)))
    add_err("`obs` contains %d non-finite value(s).", sum(!is.finite(obs)))
  if (any(duplicated(names(obs))))
    add_err("`obs` has duplicated gene names: %s",
            paste(utils::head(unique(names(obs)[duplicated(names(obs))]), 5), collapse = ", "))
  obs_range <- if (length(obs) && any(is.finite(obs))) range(obs, na.rm = TRUE) else c(NA, NA)
  if (all(is.finite(obs_range)) && all(obs_range >= 0))
    add_wrn("every `obs` value is >= 0; a two-group comparison normally has both signs.")

  ## ---- frac ------------------------------------------------------------
  frac <- as.matrix(frac)
  if (!is.numeric(frac)) add_err("`frac` is not numeric.")
  if (is.null(colnames(frac))) add_err("`frac` has no cell-type column names.")
  if (any(!is.finite(frac))) add_err("`frac` contains non-finite values.")
  if (any(frac < 0, na.rm = TRUE)) add_err("`frac` contains negative values.")
  if (any(duplicated(colnames(frac))))
    add_err("`frac` has duplicated cell-type names: %s",
            paste(unique(colnames(frac)[duplicated(colnames(frac))]), collapse = ", "))
  rs <- rowSums(frac)
  if (any(rs <= 0)) add_err("%d row(s) of `frac` sum to zero.", sum(rs <= 0))
  frac_rowsum_range <- if (length(rs)) range(rs[rs > 0], na.rm = TRUE) else c(NA, NA)
  if (all(is.finite(frac_rowsum_range)) && any(abs(frac_rowsum_range - 1) > 0.05))
    add_wrn("some `frac` rows do not sum to 1 (range %.2f to %.2f); renormalised internally.",
            frac_rowsum_range[1], frac_rowsum_range[2])

  ## ---- mu --------------------------------------------------------------
  mu <- as.matrix(mu)
  if (!is.numeric(mu)) add_err("`mu` is not numeric.")
  if (is.null(rownames(mu))) add_err("`mu` has no cell-type row names.")
  if (any(!is.finite(mu))) add_err("`mu` contains %d non-finite value(s).", sum(!is.finite(mu)))
  if (any(mu < 0, na.rm = TRUE)) add_err("`mu` contains negative values.")
  if (any(duplicated(rownames(mu))))
    add_err("`mu` has duplicated cell-type row names: %s",
            paste(unique(rownames(mu)[duplicated(rownames(mu))]), collapse = ", "))

  ## ---- alignment -------------------------------------------------------
  if (!is.null(colnames(frac)) && !is.null(rownames(mu))) {
    only_mu   <- setdiff(rownames(mu), colnames(frac))
    only_frac <- setdiff(colnames(frac), rownames(mu))
    if (length(only_mu))
      add_err("cell type(s) present in `mu` but not in `frac`: %s",
              paste(only_mu, collapse = ", "))
    if (length(only_frac))
      add_err("cell type(s) present in `frac` but not in `mu`: %s",
              paste(only_frac, collapse = ", "))
    if (!length(only_mu) && !length(only_frac) &&
        !setequal(colnames(frac), rownames(mu)))
      add_err("cell-type names differ in case, space or encoding between `frac` and `mu`.")
  }

  ## ---- genes -----------------------------------------------------------
  common <- intersect(names(obs), colnames(mu))
  if (length(common) < min_genes)
    add_err("only %d gene(s) shared between `obs` and `mu`; at least %d required.",
            length(common), min_genes)
  n_dropped <- length(names(obs)) - length(common)
  if (n_dropped > 0)
    add_wrn("%d gene(s) in `obs` are absent from `mu` and will be returned as NA.",
            n_dropped)

  ## ---- is_case ---------------------------------------------------------
  ## Guard each subsequent check on the earlier ones having passed: `any()`,
  ## `all()` and `sum()` all return NA once `is_case` contains an NA, and an
  ## NA in `if (...)` is an error rather than a message. The audit found this
  ## by passing c(FALSE, FALSE, NA, TRUE, TRUE, TRUE).
  is_case_clean <- is.logical(is_case) && !anyNA(is_case)
  if (!is.logical(is_case))
    add_err("`is_case` is not logical (class: %s).", class(is_case)[1])
  else if (anyNA(is_case))
    add_err(paste0(
      "`is_case` contains %d NA. Every downstream step would propagate NA and ",
      "every gene would be labelled 'unchanged' without an error. ",
      "Supply a complete logical vector."), sum(is.na(is_case)))
  if (length(is_case) != nrow(frac))
    add_err("`is_case` has length %d but `frac` has %d rows.", length(is_case), nrow(frac))
  else if (is_case_clean) {
    n_case <- sum(is_case)
    n_ctrl <- sum(!is_case)
    if (n_case == 0)  add_err("`is_case` has no TRUE (no case samples).")
    if (n_ctrl == 0)  add_err("`is_case` has no FALSE (no control samples).")
    if (n_case > 0 && n_ctrl > 0 && min(n_case, n_ctrl) < 2)
      add_wrn("only %d sample(s) in the smaller group; group means are unstable.",
              min(n_case, n_ctrl))
  }

  ## ---- semantic checks on mu (the silent ones) -------------------------
  mu_max <- suppressWarnings(max(mu, na.rm = TRUE))
  mu_mean <- suppressWarnings(mean(mu, na.rm = TRUE))
  mu_is_integer <- isTRUE(all.equal(mu, round(mu)))

  if (is.finite(mu_max) && mu_max < 30)
    add_err(paste0(
      "the largest per-cell mean in `mu` is ", format(mu_max, digits = 3),
      ". A linear CP10k reference normally exceeds 30 for at least one gene. ",
      "`mu` looks log-scaled (or otherwise compressed); the ratio inside log2() ",
      "would be meaningless. Pass a LINEAR matrix."))

  if (is.finite(mu_max) && mu_max > 1000 && mu_is_integer)
    add_wrn(paste0(
      "`mu` holds integers with a maximum of ", format(mu_max, digits = 3),
      ", which looks like raw counts rather than per-cell means. Counts and CP10k ",
      "differ by the library size of each cell type, so the composition term would ",
      "be scaled by that ratio. Consider dividing by the library size first."))

  ## cell types with no expression at all
  dead <- rownames(mu)[apply(mu, 1, function(r) all(r == 0, na.rm = TRUE))]
  if (length(dead))
    add_err("cell type(s) with all-zero expression in `mu`: %s",
            paste(dead, collapse = ", "))

  ## ---- verdict ---------------------------------------------------------
  if (strict) { errors <- c(errors, warns); warns <- character(0) }

  list(
    ok       = length(errors) == 0,
    errors   = errors,
    warnings = warns,
    remedies = .input_remedies(errors, warns, list(
      cell_types_frac = colnames(frac), cell_types_mu = rownames(mu),
      mu_max = mu_max, mu_all_integer = mu_is_integer,
      n_genes_shared = length(common), n_gene_dropped = n_dropped,
      n_case = if (is.logical(is_case)) sum(is_case == TRUE) else NA_integer_,
      n_control = if (is.logical(is_case)) sum(is_case == FALSE) else NA_integer_)),
    details  = list(
      n_genes_obs        = length(obs),
      n_genes_shared     = length(common),
      n_gene_dropped     = n_dropped,
      genes_dropped      = utils::head(setdiff(names(obs), common), 20),
      cell_types_frac    = colnames(frac),
      cell_types_mu      = rownames(mu),
      n_cell_types       = length(rownames(mu)),
      obs_range          = obs_range,
      frac_rowsum_range  = frac_rowsum_range,
      mu_max             = mu_max,
      mu_mean            = mu_mean,
      mu_all_integer     = mu_is_integer,
      n_case             = if (is.logical(is_case)) sum(is_case == TRUE) else NA_integer_,
      n_control          = if (is.logical(is_case)) sum(is_case == FALSE) else NA_integer_
    )
  )
}


#' Self-test: can this installation recover a change it injected itself?
#'
#' Builds a pseudo-bulk data set from a synthetic reference with a known
#' composition change and a known within-cell change, decomposes it, and
#' reports whether the injected change was recovered.
#'
#' @param n_genes Number of synthetic genes, default 400.
#' @param n_ctrl,n_case Samples per group, default 20 each.
#' @param noise Log-normal noise, default 0.05.
#' @param seed Integer seed passed to [simulate_pseudobulk()].
#' @param verbose Print a short report.
#'
#' @return A list with `pass` (logical), `recovery` (numeric), `cor_pred`
#'   (numeric) and `checks` (a data frame of individual criteria).
#'
#' @details
#' Run this after installing, and again whenever the R or package versions
#' change. It takes about a second and does not need any data.
#'
#' The criteria are deliberately loose (`recovery > 0.9`, `cor_pred > 0.99`)
#' so that ordinary floating-point and platform variation cannot cause a false
#' failure, while a genuine regression in the arithmetic will.
#'
#' @examples
#' st <- self_test(verbose = FALSE)
#' st$pass
#' @export
self_test <- function(n_genes = 400, n_ctrl = 20, n_case = 20,
                      noise = 0.05, seed = 20261023, verbose = TRUE) {
  set.seed(seed)
  G <- as.integer(n_genes)
  mu <- matrix(stats::runif(3 * G, 50, 300), nrow = 3,
               dimnames = list(c("A", "B", "C"), paste0("g", seq_len(G))))
  ## cell-type specific genes so that composition carries real signal
  mu["A", 1] <- 900; mu["B", 1] <- 3; mu["C", 1] <- 2
  mu["B", 2] <- 900; mu["A", 2] <- 3; mu["C", 2] <- 2

  fc <- c(A = 0.60, B = 0.30, C = 0.10)
  fa <- c(A = 0.20, B = 0.50, C = 0.30)

  sim <- simulate_pseudobulk(mu, fc, fa, n_ctrl = n_ctrl, n_case = n_case,
                             intrinsic = c(g1 = 1, g2 = -1), noise = noise,
                             seed = seed)

  ## observed log2 fold change from the simulated bulk
  obs <- vapply(colnames(sim$bulk), function(g) {
    a <- mean(sim$bulk[sim$is_case, g]); b <- mean(sim$bulk[!sim$is_case, g])
    if (a <= 0 || b <= 0) NA_real_ else log2(a / b)
  }, numeric(1))

  frac <- rbind(matrix(rep(fc, n_ctrl), ncol = 3, byrow = TRUE),
                matrix(rep(fa, n_case), ncol = 3, byrow = TRUE))
  colnames(frac) <- c("A", "B", "C")

  d  <- decompose_composition(obs, frac, mu, sim$is_case, validate = FALSE)
  sc <- score_simulation(d, sim$truth, intrinsic_genes = c("g1", "g2"))

  ## ---- 第 4 条：报告层不得改动任何数字 ----
  ## 这一条是"加入 report() 之后真实性是否还成立"的直接证据：
  ## 渲染前后把对象逐字段比对，任何改动都会让它失败。
  before <- as.data.frame(d)[, c("gene", "obs", "pred", "resid")]
  invisible(suppressWarnings(suppressMessages(
    capture.output(report(d, calibrate = FALSE, attribute = FALSE, top_n = 5)))))
  after <- as.data.frame(d)[, c("gene", "obs", "pred", "resid")]
  identical_after <- isTRUE(all.equal(before, after, tolerance = 0))

  checks <- data.frame(
    criterion = c("composition term recovers the truth",
                  "injected within-cell change recovered",
                  "no false positives above tolerance",
                  "reporting does not alter the decomposition"),
    measured  = c(sc$cor_pred, sc$recovery_intrinsic, sc$false_positive_rate,
                  as.numeric(identical_after)),
    threshold = c("> 0.99", "in (0.90, 1.10)", "< 0.02", "== 1 (bit-identical)"),
    pass      = c(isTRUE(sc$cor_pred > 0.99),
                  isTRUE(sc$recovery_intrinsic > 0.90 &&
                         sc$recovery_intrinsic < 1.10),
                  isTRUE(sc$false_positive_rate < 0.02),
                  identical_after),
    stringsAsFactors = FALSE
  )

  pass <- all(checks$pass)

  if (verbose) {
    cat("cade self-test\n")
    cat(sprintf("  R %s.%s\n", R.version$major, R.version$minor))
    cat(sprintf("  cade %s\n",
                tryCatch(as.character(utils::packageVersion("cade")),
                         error = function(e) "development")))
    cat(sprintf("  %d genes, %d vs %d samples, noise %.2f\n\n",
                G, n_ctrl, n_case, noise))
    for (i in seq_len(nrow(checks)))
      cat(sprintf("  [%s] %-42s measured %.4f  (%s)\n",
                  if (checks$pass[i]) "PASS" else "FAIL",
                  checks$criterion[i], checks$measured[i], checks$threshold[i]))
    cat(sprintf("\n  %s\n", if (pass) "self-test PASSED" else "self-test FAILED"))
  }

  list(pass = pass, recovery = sc$recovery_intrinsic,
       cor_pred = sc$cor_pred, checks = checks)
}


#' Print method for a decomposition
#'
#' @param x Output of [decompose_composition()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.cade_decomp <- function(x, ...) {
  n <- nrow(x)
  n_na <- sum(!is.finite(x$pred))
  cat(sprintf("<cade decomposition> %d genes", n))
  if (n_na) cat(sprintf(" (%d with no usable composition term)", n_na))
  cat("\n")
  if (n) {
    rng <- function(v) sprintf("[%+.2f, %+.2f]", min(v, na.rm = TRUE), max(v, na.rm = TRUE))
    cat(sprintf("  observed log2FC   %s\n", rng(x$obs)))
    cat(sprintf("  composition term  %s\n", rng(x$pred)))
    cat(sprintf("  residual          %s\n", rng(x$resid)))
  }
  invisible(x)
}
