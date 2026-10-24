# ============================================================================
# cade: calibration diagnostics
#
# Regress observed log2 fold change on the composition-implied log2 fold
# change. With a good reference and a correct pipeline the slope is near 1 and
# the intercept near 0: if composition explains the data, then obs ~ pred.
#
# THE ONE THING THAT MATTERS: HOW GENES ARE FILTERED.
#
#   Filtering by |obs| (effect size) is CIRCULAR. Selecting genes with large
#   observed changes and then reporting how well composition predicts those
#   changes guarantees a favourable slope, because pred and obs are correlated
#   by construction whenever composition carries any signal.
#
#   Filtering by EXPRESSION LEVEL is NOT circular. Low-expressed genes have
#   noisy log2 fold changes; including them dilutes the slope toward zero.
#   Removing them by an expression rank is independent of the outcome.
#
# Hence every function here filters by expression, never by effect size.
# ============================================================================

#' Expression-rank filter
#'
#' Selects genes whose mean expression is above a quantile of the distribution.
#'
#' @param mean_expr Numeric vector of mean expression (linear scale).
#' @param keep Quantile in `[0, 1)`. `0` keeps everything; `0.5` keeps the
#'   upper half; `0.75` keeps the upper quarter.
#' @param method `"quantile"` (default) uses the empirical quantile of
#'   `mean_expr`; `"absolute"` uses a fixed cut on `mean_expr` instead, which is
#'   only meaningful when `mean_expr` is on a known scale (CP10k, FPKM, ...)
#'   and therefore not comparable across platforms.
#'
#' @return Logical vector, `TRUE` for kept genes.
#' @examples
#' expression_filter(c(0.1, 1, 10, 100), keep = 0.5)
#' @export
expression_filter <- function(mean_expr, keep = 0, method = c("quantile", "absolute")) {
  method <- match.arg(method)
  mean_expr <- as.numeric(mean_expr)
  if (any(!is.finite(mean_expr)))
    stop("`mean_expr` contains non-finite values.")
  if (keep < 0 || keep >= 1)
    stop("`keep` must be in [0, 1).")
  if (keep == 0) return(rep(TRUE, length(mean_expr)))
  if (method == "quantile") {
    cut <- stats::quantile(mean_expr, keep, na.rm = TRUE)
    return(mean_expr > cut)
  }
  mean_expr > keep
}


#' Calibration of a composition-attribution decomposition
#'
#' @param decomp A data frame from [decompose_composition()], or any data frame
#'   with columns `gene`, `obs`, `pred` and `mean_expr`.
#' @param levels Numeric vector of expression quantiles at which to report the
#'   calibration. Default reports all genes plus the upper half and upper
#'   quarter.
#' @param min_n Minimum number of genes for a stratum to be reported; strata
#'   below this are dropped rather than reported with an unstable slope.
#'
#' @return A data frame with one row per stratum and columns `level`, `n`,
#'   `slope`, `intercept`, `r2`, `pearson`, `median_pred_error`.
#'
#' @details
#' The slope is an ordinary least squares slope of `obs` on `pred`. A slope
#' well below 1 with a good R-squared means the composition term carries the
#' right direction but is scaled down - typical of a reference whose cell
#' states do not match the tissue. A slope near 0 with near-zero R-squared
#' means the composition prediction is uninformative for this dataset: the
#' decomposition should not be interpreted at all.
#'
#' **Never** call this with genes pre-filtered on `|obs|`.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(gene = paste0("g", 1:100),
#'                 obs = rnorm(100),
#'                 mean_expr = runif(100, 0.1, 100))
#' d$pred <- d$obs * 0.8 + rnorm(100, sd = 0.3)
#' check_calibration(d)
#' @export
check_calibration <- function(decomp, levels = c(0, 0.5, 0.75), min_n = 30) {

  need <- c("gene", "obs", "pred", "mean_expr")
  miss <- setdiff(need, names(decomp))
  if (length(miss))
    stop(sprintf("`decomp` is missing column(s): %s", paste(miss, collapse = ", ")))

  d <- decomp[is.finite(decomp$obs) & is.finite(decomp$pred) &
              is.finite(decomp$mean_expr), , drop = FALSE]
  if (nrow(d) < min_n)
    stop(sprintf("only %d usable genes; need at least %d.", nrow(d), min_n))

  rows <- lapply(levels, function(lv) {
    keep <- expression_filter(d$mean_expr, lv)
    n <- sum(keep)
    if (n < min_n) return(NULL)
    fit <- stats::lm(d$obs[keep] ~ d$pred[keep])
    data.frame(
      level             = lv,
      label             = if (lv == 0) "all genes" else
                          sprintf("top %d%% expressed", round((1 - lv) * 100)),
      n                 = n,
      slope             = unname(stats::coef(fit)[2]),
      intercept         = unname(stats::coef(fit)[1]),
      r2                = summary(fit)$r.squared,
      pearson           = stats::cor(d$obs[keep], d$pred[keep]),
      median_pred_error = stats::median(abs(d$obs[keep] - d$pred[keep])),
      stringsAsFactors  = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out))
    stop("no stratum reached `min_n`.")
  rownames(out) <- NULL
  out
}


#' Verdict on a calibration table
#'
#' Turns the numeric output of [check_calibration()] into a short, defensible
#' statement, so that the wording in a manuscript does not drift away from what
#' the numbers support.
#'
#' @param cal Output of [check_calibration()].
#' @param good Lower and upper bounds on an acceptable slope at the most
#'   stringent stratum, default `c(0.7, 1.5)`.
#'
#' @return A one-element character vector.
#' @examples
#' set.seed(1)
#' d <- data.frame(gene = paste0("g", 1:200),
#'                 pred = rnorm(200), mean_expr = runif(200, 0.1, 100))
#' d$obs <- d$pred * 0.9 + rnorm(200, sd = 0.3)
#' calibration_verdict(check_calibration(d))
#' @export
calibration_verdict <- function(cal, good = c(0.7, 1.5)) {
  last <- cal[nrow(cal), ]
  ok <- last$slope >= good[1] && last$slope <= good[2]
  if (ok) {
    sprintf("composition prediction is informative: slope %.2f (R2 %.2f) on the %s stratum (n = %d)",
            last$slope, last$r2, last$label, last$n)
  } else if (last$slope < good[1] && last$r2 < 0.05) {
    sprintf("composition prediction is uninformative: slope %.2f, R2 %.3f on the %s stratum; do not interpret the decomposition",
            last$slope, last$r2, last$label)
  } else {
    sprintf("composition prediction is weakly calibrated: slope %.2f, R2 %.2f on the %s stratum; interpret direction only, not magnitude",
            last$slope, last$r2, last$label)
  }
}
