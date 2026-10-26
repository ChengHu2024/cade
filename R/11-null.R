# ============================================================================
# Empirical null calibration and reference diagnostics
#
# The residual produced by decompose_composition() has no scale: a value of
# 0.5 is not intrinsically large or small. null_calibrate() builds an
# empirical null from comparisons among normal samples of the same study, and
# reports each gene's position in that null.
#
# IMPORTANT (documented limitation): residuals from normal-versus-normal
# comparisons are compared by POOLING, not by averaging. Averaging over
# complementary splits cancels the residual exactly, because a 4-vs-4 split and
# its complement give opposite residuals. Only non-complementary splits are
# used, and their absolute residuals are pooled.
# ============================================================================

#' Empirical null distribution for the within-cell residual
#'
#' Builds an empirical null for the residual by decomposing comparisons among
#' normal samples of the same study, then reports each gene's position in that
#' null.
#'
#' @param resid_normal Matrix of residuals from normal-versus-normal
#'   decompositions, one column per comparison, genes in rows. Column names are
#'   optional.
#' @param resid_observed Named numeric vector of observed residuals whose
#'   position is to be evaluated. Names must match \code{rownames(resid_normal)}.
#' @param probs Quantile levels to report for the null.
#'
#' @return A list with
#'   \item{null_quantiles}{quantiles of the pooled absolute null residuals;}
#'   \item{null_pooled}{the pooled absolute null residuals;}
#'   \item{genes}{a data frame with \code{gene}, \code{resid},
#'     \code{abs_resid}, \code{pct} (the empirical percentile within the null),
#'     \code{exceeds_95} and \code{exceeds_99}.}
#'
#' @details
#' The null pooled here describes a single group-level comparison of the same
#' size as the normal-versus-normal splits used to build it. A case-control
#' comparison with larger groups has less noise, so its null is narrower. Match
#' group sizes before comparing across studies: for groups of size
#' \eqn{n_A, n_B} and a reference null built at \eqn{n_0} per side, the
#' approximate scaling factor is
#' \deqn{\sqrt{((1/n_A + 1/n_B) / (2/n_0))}}
#'
#' The null does NOT correct for a systematic difference between the reference
#' and the study cohort, which affects cell-type-specific genes in a consistent
#' direction. It controls donor-to-donor and batch variation only.
#'
#' @examples
#' set.seed(1)
#' n_g <- 500
#' null_mat <- matrix(rnorm(n_g * 20, 0, 0.2), n_g, 20,
#'                    dimnames = list(paste0("g", 1:n_g), NULL))
#' observed <- setNames(rnorm(n_g, 0, 0.2), paste0("g", 1:n_g))
#' observed[1:10] <- observed[1:10] * 5
#' out <- null_calibrate(null_mat, observed)
#' head(out$genes)
#' table(out$genes$exceeds_99)
#' @export
null_calibrate <- function(resid_normal, resid_observed,
                           probs = c(.5, .75, .90, .95, .99)) {
  resid_normal <- as.matrix(resid_normal)
  if (is.null(rownames(resid_normal)))
    stop("resid_normal must carry gene names in rownames")
  if (!is.null(names(resid_observed))) {
    g <- intersect(rownames(resid_normal), names(resid_observed))
    if (length(g) < 10) stop("fewer than 10 genes in common")
    resid_normal <- resid_normal[g, , drop = FALSE]
    resid_observed <- resid_observed[g]
  } else {
    if (length(resid_observed) != nrow(resid_normal))
      stop("resid_observed length must match nrow(resid_normal)")
    g <- rownames(resid_normal)
  }
  pooled <- abs(as.vector(resid_normal))
  pooled <- pooled[!is.na(pooled)]
  if (length(pooled) < 100) stop("too few null residuals to calibrate")
  ## A9: `exceeds_95`/`exceeds_99` used to index into `nl` by matching probs.
  ## When `probs` did not contain .95 the first index returned NA, and a probs
  ## vector without .95 but with .99 silently attached the 99th percentile to
  ## the 95% flag. Compute both thresholds directly instead.
  if (any(diff(probs) <= 0))
    stop("probs must be strictly increasing")
  if (any(probs <= 0 | probs >= 1))
    stop("probs must lie strictly between 0 and 1")
  nl  <- stats::quantile(pooled, probs, names = FALSE)
  a   <- abs(resid_observed)
  pct <- vapply(a, function(x) mean(pooled < x), numeric(1))
  q95 <- stats::quantile(pooled, 0.95, names = FALSE)
  q99 <- stats::quantile(pooled, 0.99, names = FALSE)
  list(
    null_quantiles = stats::setNames(nl, paste0(probs * 100, "%")),
    null_pooled = pooled,
    genes = data.frame(gene = g, resid = as.numeric(resid_observed),
                       abs_resid = a, pct = pct,
                       exceeds_95 = a > q95, exceeds_99 = a > q99,
                       stringsAsFactors = FALSE))
}

#' @rdname null_calibrate
#' @param n_a,n_b Group sizes of the case-control comparison.
#' @param n_0 Side size of the normal-versus-normal splits used for the null.
#' @param probs Level to scale.
#' @return \code{scale_null} returns the multiplicative factor to apply to the
#'   null quantiles.
#' @export
scale_null <- function(n_a, n_b, n_0 = 4, probs = 0.99) {
  if (n_a < 2 || n_b < 2 || n_0 < 2) stop("group sizes must be at least 2")
  sqrt((1 / n_a + 1 / n_b) / (2 / n_0))
}

#' Diagnostic for a single-cell reference matrix
#'
#' Reports, for each cell type, the number of cells behind the reference mean and
#' the median relative standard error across genes.
#'
#' @param mu Numeric matrix, genes in rows, cell types in columns.
#' @param mu_se Numeric matrix of standard errors, same shape.
#' @param n_cells Named integer vector of cell counts per cell type. Optional.
#' @param rel_se_soft Median relative standard error above which a cell type is
#'   flagged.
#'
#' @return A data frame with one row per cell type.
#' @examples
#' set.seed(1)
#' mu <- matrix(runif(4 * 200, 1, 10), 200, 4,
#'              dimnames = list(NULL, c("A", "B", "C", "D")))
#' se <- mu * 0.2
#' reference_qc(mu, se, n_cells = c(A = 5000, B = 2000, C = 800, D = 300))
#' @export
reference_qc <- function(mu, mu_se, n_cells = NULL, rel_se_soft = 0.25) {
  mu <- as.matrix(mu); mu_se <- as.matrix(mu_se)
  if (!identical(dim(mu), dim(mu_se))) stop("mu and mu_se must have the same shape")
  cts <- colnames(mu)
  out <- lapply(seq_along(cts), function(k) {
    m <- mu[, k]; s <- mu_se[, k]
    ok <- m > 1e-6 & !is.na(m) & !is.na(s)
    rel <- if (sum(ok) > 10) s[ok] / m[ok] else NA_real_
    data.frame(celltype = cts[k],
               n_cells = if (!is.null(n_cells)) unname(n_cells[cts[k]]) else NA_integer_,
               median_rel_se = stats::median(rel, na.rm = TRUE),
               q75_rel_se = stats::quantile(rel, .75, names = FALSE, na.rm = TRUE),
               pct_rel_se_gt_50 = 100 * mean(rel > .5, na.rm = TRUE),
               flagged = is.na(stats::median(rel, na.rm = TRUE)) ||
                         stats::median(rel, na.rm = TRUE) > rel_se_soft,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
