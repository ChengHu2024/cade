# ============================================================================
# cade: core decomposition
#
# The central quantity:
#
#   pred(g)   = log2( sum_k f_case[k] * mu[k, g] / sum_k f_ctrl[k] * mu[k, g] )
#   resid(g)  = obs(g) - pred(g)
#
# where obs(g) is the observed log2 fold change of gene g, f* are the mean
# cell-type fractions in each group, and mu is the per-cell mean expression
# matrix (cell type x gene) on a LINEAR scale (e.g. CP10k).
#
# DESIGN RULES (learned the hard way; each one prevents a silent error):
#
#   R1. The column order of `frac` and the row order of `mu` must agree.
#       They usually do NOT, because deconvolution tools and reference tables
#       are produced by different code. `decompose_composition()` reorders and
#       verifies; it never trusts the input order.
#
#   R2. `mu` must be linear. If it is already log-transformed the ratio inside
#       the log2 becomes meaningless and the sign can flip.
#
#   R3. `obs` must be a log2 fold change, not a difference of logs with an
#       added pseudocount on already-log data. `check_scale()` flags the
#       classic double-log mistake: 2^E -> log2(2^E + 1).
#
#   R4. Guard against zero denominators. A cell type can have zero mean
#       expression for a gene; the ratio is then undefined, not infinity.
# ============================================================================

#' Compose a composition-implied log2 fold change
#'
#' For each gene, computes the log2 ratio of the bulk expression predicted by
#' the case composition to that predicted by the control composition.
#'
#' @param frac A matrix or data frame of cell-type fractions, one row per
#'   sample, one column per cell type. Columns are matched to `mu` by NAME, not
#'   by position.
#' @param mu A matrix of per-cell mean expression with cell types in rows and
#'   genes in columns, on a LINEAR scale. Row names must be cell-type labels.
#' @param is_case Logical vector, one element per row of `frac`; `TRUE` marks
#'   case samples.
#' @param genes Optional character vector restricting the output.
#'
#' @return A named numeric vector of log2 fold changes, one per gene. `NA` where
#'   the denominator is zero or a cell type is missing.
#'
#' @details
#' Cell types present in only one of `frac` and `mu` raise an error rather than
#' being silently dropped: dropping a cell type changes the composition
#' prediction for every gene.
#'
#' @examples
#' mu <- matrix(c(10, 1, 5, 2), nrow = 2,
#'              dimnames = list(c("A", "B"), c("g1", "g2")))
#' fr <- matrix(c(0.8, 0.2,
#'                0.8, 0.2,
#'                0.2, 0.8,
#'                0.2, 0.8), nrow = 4, byrow = TRUE,
#'              dimnames = list(paste0("s", 1:4), c("A", "B")))
#' compose_prediction(fr, mu, is_case = c(FALSE, FALSE, TRUE, TRUE))
#' @export
compose_prediction <- function(frac, mu, is_case, genes = NULL) {

  frac <- as.matrix(frac)
  mu   <- as.matrix(mu)

  if (!is.logical(is_case)) stop("`is_case` must be logical.")
  if (length(is_case) != nrow(frac))
    stop(sprintf("`is_case` has length %d but `frac` has %d rows.",
                 length(is_case), nrow(frac)))
  if (all(is_case) || !any(is_case))
    stop("`is_case` must contain both TRUE and FALSE.")
  if (any(!is.finite(frac)))
    stop("`frac` contains non-finite values.")
  if (any(frac < 0, na.rm = TRUE))
    stop("`frac` contains negative values; fractions must be >= 0.")
  if (is.null(rownames(mu)))
    stop("`mu` must have cell-type row names.")
  if (any(!is.finite(mu)))
    stop("`mu` contains non-finite values (NA/Inf).")

  ## ---- R1: align by name, then prove it -------------------------------
  ct <- rownames(mu)
  missing_in_frac <- setdiff(ct, colnames(frac))
  extra_in_frac   <- setdiff(colnames(frac), ct)
  if (length(missing_in_frac))
    stop(sprintf("cell type(s) in `mu` but not in `frac`: %s",
                 paste(missing_in_frac, collapse = ", ")))
  if (length(extra_in_frac))
    stop(sprintf("cell type(s) in `frac` but not in `mu`: %s",
                 paste(extra_in_frac, collapse = ", ")))

  frac <- frac[, ct, drop = FALSE]
  if (!identical(colnames(frac), rownames(mu)))
    stop("internal: column alignment failed")   # can only fire if names are
                                                # duplicated or NA

  ## ---- R4: proportions need not sum to 1, but should be close ---------
  rs <- rowSums(frac)
  if (any(rs <= 0))
    stop("some samples have all-zero fractions.")
  frac <- frac / rs                              # normalise; harmless if
                                                 # already normalised

  f_ctrl <- colMeans(frac[!is_case, , drop = FALSE])
  f_case <- colMeans(frac[ is_case, , drop = FALSE])

  genes <- if (is.null(genes)) colnames(mu) else intersect(genes, colnames(mu))
  if (!length(genes)) stop("no genes in common.")

  num <- as.numeric(crossprod(f_case, mu[, genes, drop = FALSE]))
  den <- as.numeric(crossprod(f_ctrl, mu[, genes, drop = FALSE]))

  ## ---- R4: guard the ratio --------------------------------------------
  bad <- !is.finite(num) | !is.finite(den) | num <= 0 | den <= 0
  out <- rep(NA_real_, length(genes))
  out[!bad] <- log2(num[!bad] / den[!bad])
  names(out) <- genes

  attr(out, "f_ctrl") <- setNames(f_ctrl, ct)
  attr(out, "f_case") <- setNames(f_case, ct)
  out
}


#' Decompose observed log2 fold changes into composition and cell-intrinsic parts
#'
#' @param obs Named numeric vector of observed log2 fold changes, one per gene.
#' @param frac Fraction matrix (samples x cell types).
#' @param mu Per-cell mean expression (cell types x genes), linear scale.
#' @param is_case Logical vector over rows of `frac`.
#' @param mean_expr Optional named numeric vector of mean expression used for
#'   calibration filtering. If `NULL`, it is taken from the case+control average
#'   of `mu` weighted by the pooled composition.
#' @param validate Run [check_inputs()] before decomposing. Default `TRUE`;
#'   set `FALSE` only inside functions that have already validated.
#'
#' @return A data frame of class `cade_decomp` with columns `gene`, `obs`,
#'   `pred`, `resid`, `pct`, `mean_expr`. `pct` is `resid / obs` and is `NA`
#'   when `obs` is near zero. Attributes `f_ctrl` and `f_case` hold the mean
#'   composition of each group.
#'
#' @details
#' `compose_prediction()` is the only place where the composition term is
#' computed. This function merely joins it to `obs` and subtracts.
#'
#' With `validate = TRUE` the function stops if [check_inputs()] returns any
#' error. The checks that stop it are structural (misalignment, wrong length,
#' non-finite values) and semantic (a `mu` that is log-scaled rather than
#' linear). Warnings are printed but do not stop the computation.
#'
#' @examples
#' set.seed(1)
#' mu <- matrix(stats::runif(2 * 12, 30, 400), nrow = 2,
#'              dimnames = list(c("A", "B"), paste0("g", 1:12)))
#' fr <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4,
#'              byrow = TRUE, dimnames = list(paste0("s", 1:4), c("A", "B")))
#' ob <- stats::setNames(stats::rnorm(12), colnames(mu))
#' d <- decompose_composition(ob, fr, mu, is_case = c(FALSE, FALSE, TRUE, TRUE))
#' head(d)
#' @export
decompose_composition <- function(obs, frac, mu, is_case, mean_expr = NULL,
                                  validate = TRUE) {

  if (validate) {
    chk <- check_inputs(obs, frac, mu, is_case)
    if (length(chk$warnings)) {
      for (w in chk$warnings) warning(w, call. = FALSE)
    }
    if (!chk$ok) {
      stop(paste0("input validation failed:\n",
                  paste0("  * ", chk$errors, collapse = "\n"),
                  "\nRun check_inputs() for the full report."),
           call. = FALSE)
    }
  }

  if (is.null(names(obs))) stop("`obs` must be a named numeric vector.")
  pred <- compose_prediction(frac, mu, is_case, genes = names(obs))

  resid <- obs - pred
  pct <- ifelse(is.na(obs) | abs(obs) < 1e-12, NA_real_, resid / obs)

  if (is.null(mean_expr)) {
    fc <- attr(pred, "f_ctrl"); fa <- attr(pred, "f_case")
    pooled <- (fc + fa) / 2
    mean_expr <- as.numeric(crossprod(pooled, mu[names(fc), names(pred),
                                                 drop = FALSE]))
    names(mean_expr) <- names(pred)
  } else {
    mean_expr <- mean_expr[names(pred)]
  }

  out <- data.frame(gene = names(pred),
                    obs   = as.numeric(obs),
                    pred  = as.numeric(pred),
                    resid = as.numeric(resid),
                    pct   = as.numeric(pct),
                    mean_expr = as.numeric(mean_expr),
                    stringsAsFactors = FALSE)
  attr(out, "f_ctrl") <- attr(pred, "f_ctrl")
  attr(out, "f_case") <- attr(pred, "f_case")
  class(out) <- c("cade_decomp", "data.frame")
  out
}
