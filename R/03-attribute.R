# ============================================================================
# cade: gene-level attribution and effect sizes
#
# The decomposition gives a continuous residual. Turning it into discrete
# labels ("composition-driven" / "mixed" / "residual") is convenient for
# a table but introduces thresholds that are not implied by the algebra. The
# manuscript's position is therefore:
#
#   report the continuous residual spectrum; if a classification is required,
#   report how sensitive it is to the threshold.
#
# `attribute_genes()` therefore always returns the continuous values together
# with the label, and `attribute_sensitivity()` quantifies label stability.
# ============================================================================

#' Label genes by where their change comes from
#'
#' @param decomp Output of [decompose_composition()].
#' @param lower,upper Thresholds on `pct = resid / obs`.
#'   `pct <= lower` -> "composition-driven" (the residual explains little of
#'   the observed change); `pct >= upper` -> "residual"; otherwise
#'   "mixed".
#' @param min_obs Genes with `|obs|` below this are labelled "unchanged"
#'   regardless of `pct`, because `pct` is unstable when `obs` is near zero.
#'
#' @return `decomp` with two added columns, `class` and `class_reason`.
#'
#' @details
#' The defaults `lower = 0.3`, `upper = 0.7` are conventions, not results.
#' They are chosen so that "composition-driven" means the residual is under a
#' third of the observed change. Pass explicitly chosen values if a different
#' convention is needed, and always report the sensitivity.
#'
#' `pct` is **signed**, because `resid / obs` carries the sign of `resid`
#' relative to `obs`. When the composition term overshoots - `pred` larger in
#' magnitude than `obs` and opposite in sign - `pct` can be close to 1 while the
#' residual still accounts for most of the observed change. Those genes are
#' labelled `"over-predicted"` rather than `"composition-driven"`, because a
#' signed ratio is not a meaningful share in that regime. An adversarial audit
#' found this mislabelling; it is the reason the fourth class exists.
#'
#' @examples
#' d <- data.frame(gene = c("a", "b", "c"), obs = c(1, 1, 1),
#'                 pred = c(0.95, 0.5, 0.05), resid = c(0.05, 0.5, 0.95),
#'                 pct = c(0.05, 0.5, 0.95), mean_expr = c(1, 1, 1))
#' attribute_genes(d)
#' @export
attribute_genes <- function(decomp, lower = 0.3, upper = 0.7, min_obs = 0) {
  if (!all(c("obs", "resid", "pct") %in% names(decomp)))
    stop("`decomp` must contain `obs`, `resid` and `pct`.")
  if (lower >= upper) stop("`lower` must be < `upper`.")

  out <- decomp
  out$class <- NA_character_
  out$class_reason <- NA_character_

  unchanged <- is.na(out$pct) | abs(out$obs) < min_obs
  out$class[unchanged]        <- "unchanged"
  out$class_reason[unchanged] <- "obs near zero; pct undefined or unstable"

  ## P7: opposite signs mean the composition term overshoots the observation.
  overshoot <- !unchanged & !is.na(out$resid) & !is.na(out$obs) &
               sign(out$resid) != sign(out$obs)

  out$class[overshoot]        <- "over-predicted"
  out$class_reason[overshoot] <- paste0(
    "residual (", round(out$resid[overshoot], 2), ") is opposite in sign to the ",
    "observed change (", round(out$obs[overshoot], 2),
    "); signed pct = ", round(out$pct[overshoot], 2),
    " is not a meaningful share")

  ok <- !unchanged & !overshoot
  out$class[ok] <- ifelse(
    out$pct[ok] <= lower, "composition-driven",
    ifelse(out$pct[ok] >= upper, "residual", "mixed"))
  out$class_reason[ok] <- sprintf("pct = %.2f (lower = %.2f, upper = %.2f)",
                                  out$pct[ok], lower, upper)
  out
}


#' How stable is a gene's label across thresholds?
#'
#' @param decomp Output of [decompose_composition()].
#' @param thresholds Numeric vector of candidate `lower` values, with
#'   `upper = 1 - lower` throughout.
#' @param genes Optional subset of genes.
#' @param min_obs As in [attribute_genes()].
#'
#' @return A data frame with one row per gene and the fraction of thresholds
#'   at which the gene is labelled `residual`.
#'
#' @details
#' A gene labelled intrusive at only one of several thresholds is not a stable
#' finding. The manuscript reports label stability rather than a single label
#' for exactly this reason.
#'
#' @examples
#' d <- data.frame(gene = c("a", "b"), obs = c(1, 1),
#'                 pred = c(0.5, 0.2), resid = c(0.5, 0.8),
#'                 pct = c(0.5, 0.8), mean_expr = c(1, 1))
#' attribute_sensitivity(d)
#' @export
attribute_sensitivity <- function(decomp, thresholds = seq(0.2, 0.45, by = 0.05),
                                  genes = NULL, min_obs = 0) {
  if (!is.null(genes)) decomp <- decomp[decomp$gene %in% genes, , drop = FALSE]
  if (!nrow(decomp)) stop("no genes selected.")

  thresholds <- thresholds[thresholds > 0 & thresholds < 0.5]
  if (!length(thresholds))
    stop("`thresholds` must contain at least one value in (0, 0.5).")

  ## Always build an n_genes x n_thresholds character matrix. sapply() would
  ## simplify to a vector when n_genes == 1 or when n_thresholds == 1, and the
  ## subsequent rowMeans() would then recycle silently. vapply() with an
  ## explicit template cannot collapse.
  cls <- vapply(thresholds, function(th) {
    as.character(attribute_genes(decomp, lower = th, upper = 1 - th,
                                 min_obs = min_obs)$class)
  }, FUN.VALUE = character(nrow(decomp)))

  ## vapply() returns a plain vector whenever FUN.VALUE has length 1, which
  ## happens as soon as there is a SINGLE gene - the very case this guard
  ## exists for. Reshape with the gene count as the row count: the elements are
  ## stored column-major, so matrix(x, nrow = n_genes) is correct for all four
  ## (n_genes, n_thresholds) combinations, including 1 x 1.
  if (is.null(dim(cls))) cls <- matrix(cls, nrow = nrow(decomp))
  stopifnot(nrow(cls) == nrow(decomp), ncol(cls) == length(thresholds))

  data.frame(
    gene              = decomp$gene,
    obs               = decomp$obs,
    pred              = decomp$pred,
    resid             = decomp$resid,
    pct               = decomp$pct,
    mean_obs_effect   = decomp$obs,
    n_thresholds      = length(thresholds),
    frac_residual     = rowMeans(cls == "residual"),
    frac_composition  = rowMeans(cls == "composition-driven"),
    frac_mixed        = rowMeans(cls == "mixed"),
    stable            = apply(cls, 1, function(r) length(unique(r)) == 1),
    stringsAsFactors  = FALSE
  )
}


#' Cohen's d with a confidence interval
#'
#' @param case,control Numeric vectors.
#' @param conf Confidence level, default 0.95.
#'
#' @return A one-row data frame with `n_case`, `n_control`, `d`, `se`, `lo`,
#'   `hi`.
#'
#' @details
#' Uses the pooled-SD form with the Hedges small-sample correction
#' `J = 1 - 3/(4*(n1+n2)-9)`, and a normal-approximation interval. For small
#' samples this interval is slightly optimistic; the effect size itself is
#' unaffected.
#'
#' A narrow interval on a large `d` means the two groups really do differ on
#' that scale, so an inconsistency between cohorts is heterogeneity rather than
#' imprecision. Say so explicitly instead of describing the interval as "wide".
#'
#' @examples
#' cohens_d(rnorm(30, 1), rnorm(30, 0))
#' @export
cohens_d <- function(case, control, conf = 0.95) {
  x <- as.numeric(case); y <- as.numeric(control)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) stop("need at least 2 finite values per group.")

  sp <- sqrt(((n1 - 1) * stats::var(x) + (n2 - 1) * stats::var(y)) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0)
    return(data.frame(n_case = n1, n_control = n2, d = NA_real_, se = NA_real_,
                      lo = NA_real_, hi = NA_real_))
  d  <- (mean(x) - mean(y)) / sp
  J  <- if (n1 + n2 > 3) 1 - 3 / (4 * (n1 + n2) - 9) else 1
  d  <- d * J
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2 - 2)))
  z  <- stats::qnorm(1 - (1 - conf) / 2)
  data.frame(n_case = n1, n_control = n2, d = d, se = se,
             lo = d - z * se, hi = d + z * se)
}


#' Pool Cohen's d across cohorts (fixed-effect with heterogeneity test)
#'
#' @param d Numeric vector of per-cohort effect sizes.
#' @param se Numeric vector of standard errors.
#' @param study Optional character vector of cohort labels.
#'
#' @return A list with `pooled` (a one-row data frame), `by_study`
#'   (per-cohort rows), `tau2` (between-study variance, DerSimonian-Laird),
#'   `Q`, `df`, `p_het`, and `i2`.
#'
#' @details
#' A near-zero `tau2` with a large pooled `d` supports a consistent effect.
#' A large `tau2` with a significant `Q` means the pooled number should not be
#' reported as a single effect; report the cohorts side by side instead.
#'
#' @examples
#' pool_cohens_d(c(0.9, 1.1, 1.0), c(0.3, 0.3, 0.35))
#' @export
pool_cohens_d <- function(d, se, study = NULL) {
  ok <- is.finite(d) & is.finite(se) & se > 0
  d <- d[ok]; se <- se[ok]
  if (is.null(study)) study <- paste0("cohort", seq_along(d)) else study <- study[ok]
  k <- length(d)
  if (k < 2) stop("need at least 2 cohorts.")

  w  <- 1 / se^2
  fe <- sum(w * d) / sum(w)
  Q  <- sum(w * (d - fe)^2)
  df <- k - 1
  p_het <- stats::pchisq(Q, df = df, lower.tail = FALSE)
  C  <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / C)
  i2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0

  wr <- 1 / (se^2 + tau2)
  pooled <- sum(wr * d) / sum(wr)
  pooled_se <- sqrt(1 / sum(wr))
  z <- stats::qnorm(0.975)

  list(
    pooled = data.frame(k = k, d = pooled, se = pooled_se,
                        lo = pooled - z * pooled_se, hi = pooled + z * pooled_se,
                        tau2 = tau2, Q = Q, df = df, p_het = p_het, i2 = i2),
    by_study = data.frame(study = study, d = d, se = se,
                          lo = d - z * se, hi = d + z * se),
    tau2 = tau2, Q = Q, df = df, p_het = p_het, i2 = i2
  )
}
