# ============================================================================
# cade: uncertainty of the comp-obs correlation
#
# WHY THIS FUNCTION EXISTS
#
# decompose_uq() propagates PARAMETER uncertainty: it perturbs `mu` and the
# proportions and reports how the residual moves. That is not what a reader
# needs in order to judge a reported comp-obs correlation.
#
# The correlation itself has sampling uncertainty, and the resampling UNIT
# decides which uncertainty is being reported. Papers routinely say
# "bootstrap 95% CI" without saying what was resampled, and the two units
# give intervals that differ by a factor of three to four:
#
#   unit = "gene"    resample genes; the interval reflects which genes
#                    happened to be selected. This is the WRONG unit for a
#                    claim about patients.
#
#   unit = "sample"  resample patients within each group; the interval
#                    reflects sampling of the cohort. This is the unit that
#                    a "95% CI of the correlation" normally means.
#
# `resample` then decides WHICH quantities move with the resample:
#
#   "proportions"        re-estimate proportions from the resampled samples,
#                        keep the observed log fold changes fixed. The
#                        interval thus carries proportion-estimation noise
#                        only and UNDERSTATES true between-cohort variability.
#
#   "proportions_and_obs" also resample the observed log fold changes across
#                        samples. This is the only variant that can approach
#                        the true sampling variability of the correlation.
#
# Neither argument has a default. A caller who does not state them cannot
# get an interval, because the number would be uninterpretable.
# ============================================================================

#' Sampling uncertainty of the comp-obs correlation
#'
#' @param frac A matrix or data frame of cell-type fractions, one row per
#'   sample, one column per cell type.
#' @param obs Named numeric vector of observed log2 fold changes, one per gene.
#' @param mu Per-cell mean expression (cell types x genes), linear scale.
#' @param is_case Logical vector over rows of `frac`.
#' @param resample One of `"proportions"` or `"proportions_and_obs"`. No
#'   default: the choice changes what the interval means.
#' @param unit One of `"sample"` or `"gene"`. No default. `"sample"` resamples
#'   patients within each group; `"gene"` resamples genes. For any statement
#'   about a cohort, `"sample"` is the correct unit.
#' @param n_boot Number of bootstrap replicates.
#' @param seed Integer seed.
#' @param ci Confidence level.
#' @param convention Sample-aggregation convention for the composition term;
#'   see [compose_prediction()].
#' @param expr Optional matrix of per-sample expression, samples x genes, on a
#'   LINEAR scale. Required when `resample = "proportions_and_obs"`: without
#'   per-sample data the observed log fold changes cannot be resampled, and a
#'   function that silently returned the same interval as
#'   `"proportions"` would be lying. When supplied, each bootstrap replicate
#'   recomputes `obs` from the resampled samples as well as re-estimating the
#'   proportions.
#'
#' @return A list of class `cade_corrci` with `r`, `lo`, `hi`, `width`,
#'   `n_boot`, `frac_boot`, and a `caveat` string that states in words what
#'   the interval does and does not cover.
#'
#' @details
#' With `resample = "proportions"` only the proportions are re-estimated, so
#' the interval reflects proportion noise alone. It does **not** represent the
#' true between-cohort variability of the correlation, which is wider. The
#' returned `caveat` says so, so that the number cannot be quoted without the
#' caveat attached.
#'
#' With `unit = "gene"` the interval is about gene selection, not about
#' patients. It is provided because some pipelines report it, but it should
#' not be described as a confidence interval for the cohort-level
#' correlation.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' fr <- matrix(runif(40), 20, 2, dimnames = list(NULL, c("A", "B")))
#' fr <- fr / rowSums(fr)
#' mu <- matrix(runif(200, 10, 100), 2, 100,
#'              dimnames = list(c("A", "B"), paste0("g", 1:100)))
#' ob <- stats::setNames(stats::rnorm(100), colnames(mu))
#' ic <- rep(c(FALSE, TRUE), each = 10)
#' cor_ci(fr, ob, mu, ic, resample = "proportions", unit = "sample")
#' }
#' @export
cor_ci <- function(frac, obs, mu, is_case,
                   resample = c("proportions", "proportions_and_obs"),
                   unit = c("sample", "gene"),
                   n_boot = 1000L, seed = 1L, ci = 0.95,
                   convention = c("mean_of_log", "log_of_mean"),
                   expr = NULL) {

  if (missing(resample))
    stop("`resample` must be stated explicitly: \"proportions\" or ",
         "\"proportions_and_obs\". The default was removed because the choice ",
         "changes what the interval means.", call. = FALSE)
  if (missing(unit))
    stop("`unit` must be stated explicitly: \"sample\" or \"gene\". The ",
         "default was removed because the two units give intervals that ",
         "differ by a factor of three to four.", call. = FALSE)
  resample <- match.arg(resample)
  unit     <- match.arg(unit)
  convention <- match.arg(convention)

  ## `proportions_and_obs` needs per-sample expression. Without it the
  ## observed log fold changes cannot move, and returning the same interval
  ## as `proportions` would misrepresent what was done.
  if (resample == "proportions_and_obs" && is.null(expr))
    stop("resample = \"proportions_and_obs\" requires `expr`, a matrix of ",
         "per-sample expression (samples x genes, linear scale). Without it ",
         "the observed log fold changes cannot be resampled, and the result ",
         "would be identical to resample = \"proportions\".", call. = FALSE)
  if (!is.null(expr)) {
    expr <- as.matrix(expr)
    if (nrow(expr) != nrow(frac))
      stop(sprintf("`expr` has %d rows but `frac` has %d samples.",
                   nrow(expr), nrow(frac)), call. = FALSE)
    if (any(!is.finite(expr))) stop("`expr` contains non-finite values.")
    if (any(expr < 0)) stop("`expr` contains negative values; must be linear.")
  }

  frac <- as.matrix(frac)
  mu   <- as.matrix(mu)
  if (!is.logical(is_case)) stop("`is_case` must be logical.")
  if (anyNA(is_case)) stop("`is_case` contains NA.")
  if (length(is_case) != nrow(frac))
    stop("`is_case` length does not match nrow(frac).")
  if (all(is_case) || !any(is_case))
    stop("`is_case` must contain both TRUE and FALSE.")
  if (is.null(names(obs))) stop("`obs` must be named.")
  if (is.null(colnames(mu))) stop("`mu` must have gene column names.")
  if (is.null(rownames(mu))) stop("`mu` must have cell-type row names.")

  ct <- rownames(mu)
  if (!setequal(colnames(frac), ct))
    stop(sprintf("cell-type names differ between `frac` and `mu`.\n  only in frac: %s\n  only in mu: %s",
                 paste(setdiff(colnames(frac), ct), collapse = ", "),
                 paste(setdiff(ct, colnames(frac)), collapse = ", ")), call. = FALSE)
  frac <- frac[, ct, drop = FALSE]
  rs <- rowSums(frac)
  if (any(rs <= 0)) stop("some samples have all-zero fractions.")
  frac <- frac / rs

  genes <- intersect(names(obs), colnames(mu))
  if (length(genes) < 3)
    stop(sprintf("only %d gene(s) shared between `obs` and `mu`.", length(genes)))
  obs  <- obs[genes]
  MUg  <- mu[, genes, drop = FALSE]

  comp_of <- function(fr_, conv) {
    FA <- fr_[is_case, , drop = FALSE]
    FB <- fr_[!is_case, , drop = FALSE]
    if (conv == "log_of_mean") {
      fa <- colMeans(FA); fb <- colMeans(FB)
      as.numeric(log2(pmax(as.vector(fa %*% MUg), 1e-12) /
                      pmax(as.vector(fb %*% MUg), 1e-12)))
    } else {
      Lc <- log2(pmax(FA %*% MUg, 1e-12))
      Ln <- log2(pmax(FB %*% MUg, 1e-12))
      as.numeric(colMeans(Lc) - colMeans(Ln))
    }
  }

  c0 <- comp_of(frac, convention)
  keep <- is.finite(c0) & is.finite(obs)
  if (sum(keep) < 3) stop("fewer than 3 finite gene pairs.")
  c0k <- c0[keep]; ok <- obs[keep]
  r_obs <- stats::cor(c0k, ok)
  gkeep <- genes[keep]

  set.seed(seed)
  idx_case <- which(is_case)
  idx_ctrl <- which(!is_case)
  n_g <- sum(keep)
  rs_boot <- numeric(n_boot)

  ## Pre-compute the per-sample log expression once when observations are to
  ## be resampled, so that each replicate only subsets rows.
  if (!is.null(expr)) {
    genes_all <- colnames(mu)
    if (is.null(colnames(expr)))
      stop("`expr` must have gene column names matching `mu`.", call. = FALSE)
    ex_genes <- intersect(genes_all, colnames(expr))
    if (length(ex_genes) < 3)
      stop("fewer than 3 genes shared between `expr` and `mu`.", call. = FALSE)
    EX <- log2(pmax(expr[, ex_genes, drop = FALSE], 1e-12))
    ## keep only the genes that survived the finite mask on the observed side
    ex_keep <- ex_genes[ex_genes %in% gkeep]
    if (length(ex_keep) < 3)
      stop("fewer than 3 genes survive the finite mask in `expr`.", call. = FALSE)
    EX <- EX[, ex_keep, drop = FALSE]
    obs_genes_b <- ex_keep
  }

  for (b in seq_len(n_boot)) {
    if (unit == "gene") {
      gi <- sample(n_g, n_g, replace = TRUE)
      cm <- c0k[gi]
      ok_b <- ok[gi]
      rs_boot[b] <- if (stats::sd(cm) > 0) stats::cor(cm, ok_b) else NA_real_
      next
    }

    ## unit == "sample"
    si_c <- sample(idx_case, length(idx_case), replace = TRUE)
    si_n <- sample(idx_ctrl, length(idx_ctrl), replace = TRUE)
    row_idx <- c(si_c, si_n)

    fr_b <- frac[row_idx, , drop = FALSE]
    cm <- comp_of(fr_b, convention)
    cm <- cm[keep]

    if (resample == "proportions_and_obs") {
      ## Recompute obs from the resampled samples as well.
      if (convention == "log_of_mean") {
        m_c <- colMeans(EX[si_c, , drop = FALSE])
        m_n <- colMeans(EX[si_n, , drop = FALSE])
        ok_b <- as.numeric(m_c - m_n)
      } else {
        ok_b <- as.numeric(log2(pmax(colMeans(exp(EX[si_c, , drop = FALSE] * log(2))), 1e-12) /
                                pmax(colMeans(exp(EX[si_n, , drop = FALSE] * log(2))), 1e-12)))
      }
      ## align: both cm and ok_b are in the order of obs_genes_b
      idx_cm <- match(obs_genes_b, gkeep)
      cm <- cm[idx_cm]
    } else {
      ok_b <- ok
    }

    rs_boot[b] <- if (stats::sd(cm, na.rm = TRUE) > 0)
      stats::cor(cm, ok_b, use = "complete.obs") else NA_real_
  }

  a <- (1 - ci) / 2
  lo <- stats::quantile(rs_boot, a,     names = FALSE, na.rm = TRUE)
  hi <- stats::quantile(rs_boot, 1 - a, names = FALSE, na.rm = TRUE)

  caveat <- if (resample == "proportions")
    paste0("This interval resamples ", unit, "s and re-estimates proportions only; ",
           "the observed log fold changes were held fixed. It therefore reflects ",
           "noise in proportion estimation and DOES NOT represent the true ",
           "between-cohort variability of the correlation, which is wider.")
  else
    paste0("This interval resamples ", unit, "s and recomputes both the proportions ",
           "and the observed log fold changes from `expr`. It still does not ",
           "propagate uncertainty in the reference profile; see decompose_uq() ",
           "for that.")

  if (unit == "gene")
    caveat <- paste(caveat,
                    "The resampling unit is the gene, so this interval describes gene selection, not cohort sampling.")

  out <- list(r = r_obs, lo = lo, hi = hi, width = hi - lo,
              n_boot = n_boot, frac_boot = rs_boot,
              resample = resample, unit = unit,
              convention = convention, n_genes = n_g,
              used_expr = !is.null(expr),
              caveat = caveat)
  class(out) <- "cade_corrci"
  out
}

#' @export
print.cade_corrci <- function(x, ...) {
  cat("<cade_corrci>\n")
  cat(sprintf("  r = %.4f   95%% CI [%.4f, %.4f]   width %.4f\n",
              x$r, x$lo, x$hi, x$width))
  cat(sprintf("  n_boot = %d   unit = %s   resample = %s   genes = %d\n",
              x$n_boot, x$unit, x$resample, x$n_genes))
  cat(sprintf("  convention = %s   expr supplied = %s\n",
              x$convention, isTRUE(x$used_expr)))
  cat("  caveat: ", x$caveat, "\n", sep = "")
  invisible(x)
}
