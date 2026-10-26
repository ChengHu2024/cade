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
#' @param convention How to aggregate over samples. `"mean_of_log"` (default)
#'   takes the log2 of each sample's predicted bulk and then averages within
#'   group; this matches limma's logFC. `"log_of_mean"` averages the
#'   proportions within group first and then takes the log2 ratio, which is
#'   what most deconvolution pipelines report. **The two are not
#'   interchangeable**: on GSE232911 they differ by 0.14 in the resulting
#'   comp-obs correlation. Always report which one was used.
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
compose_prediction <- function(frac, mu, is_case, genes = NULL,
                              convention = c("mean_of_log", "log_of_mean")) {

  convention <- match.arg(convention)

  frac <- as.matrix(frac)
  mu   <- as.matrix(mu)

  if (!is.logical(is_case)) stop("`is_case` must be logical.")
  if (anyNA(is_case)) stop("`is_case` contains NA; every sample must be case or control.")
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
  if (is.null(colnames(frac)))
    stop("`frac` must have cell-type column names.")
  if (any(!is.finite(mu)))
    stop("`mu` contains non-finite values (NA/Inf).")
  if (any(mu < 0))
    stop("`mu` contains negative values; reference means must be >= 0.")

  ## A1: duplicated `frac` column names silently collapse a cell type, because
  ##     `frac[, ct]` returns only the first match for each name.
  if (anyDuplicated(colnames(frac)))
    stop(sprintf("`frac` has duplicated cell-type column names: %s",
                 paste(unique(colnames(frac)[duplicated(colnames(frac))]),
                       collapse = ", ")), call. = FALSE)

  ## A2: duplicated `mu` row names silently double one cell type's weight in
  ##     the crossprod.
  if (anyDuplicated(rownames(mu)))
    stop(sprintf("`mu` has duplicated cell-type row names: %s",
                 paste(unique(rownames(mu)[duplicated(rownames(mu))]),
                       collapse = ", ")), call. = FALSE)

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

  ## A cell type that is absent from an entire group contributes nothing to
  ## the composition term but still occupies a column; warn rather than
  ## silently reducing the reference to the remaining types.
  absent <- ct[f_ctrl == 0 | f_case == 0]
  if (length(absent))
    warning(sprintf(paste0("cell type(s) with zero mean fraction in at least ",
                           "one group: %s. Genes expressed only in these ",
                           "types will return NA."),
                    paste(absent, collapse = ", ")), call. = FALSE)

  if (is.null(genes)) {
    genes <- colnames(mu)
  } else {
    dropped <- setdiff(genes, colnames(mu))
    if (length(dropped))
      stop(sprintf(paste0("%d gene(s) in `obs` are absent from `mu` ",
                          "(e.g. %s).\n  Dropping them would silently ",
                          "misalign `obs` and the composition term.\n  ",
                          "Intersect `names(obs)` with `colnames(mu)` first."),
                   length(dropped), paste(head(dropped, 5), collapse = ", ")),
           call. = FALSE)
    genes <- genes[genes %in% colnames(mu)]
  }
  if (!length(genes)) stop("no genes in common.")

  MUg <- mu[, genes, drop = FALSE]

  if (convention == "log_of_mean") {
    ## ---------------------------------------------------------------
    ## log-of-mean: average the proportions within each group first, then
    ## take the log2 ratio. This is the historical behaviour and is what
    ## most deconvolution pipelines report.
    ## ---------------------------------------------------------------
    num <- as.numeric(crossprod(f_case, MUg))
    den <- as.numeric(crossprod(f_ctrl, MUg))
    bad <- !is.finite(num) | !is.finite(den) | num <= 0 | den <= 0
    out <- rep(NA_real_, length(genes))
    out[!bad] <- log2(num[!bad] / den[!bad])
  } else {
    ## ---------------------------------------------------------------
    ## mean-of-log: take the log of each SAMPLE's predicted bulk, then
    ## average within group. This matches limma's logFC and is the
    ## convention used in the companion paper.
    ##
    ## The two conventions are not interchangeable. On GSE232911 they
    ## differ by 0.14 in the resulting comp-obs correlation
    ## (0.2846 mean-of-log versus 0.4251 log-of-mean), so the convention
    ## must always be reported.
    ## ---------------------------------------------------------------
    FA <- frac[is_case, , drop = FALSE]
    FB <- frac[!is_case, , drop = FALSE]
    ## (n x K) %*% (K x G) -> (n x G); no transpose here.
    ## MUg is already cell-types x genes, so FA %*% MUg gives one row per
    ## sample and one column per gene.
    ##
    ## DO NOT mask zeros with pmax(): if a gene is not expressed in any cell
    ## type present in a group, the predicted bulk is exactly 0 and the log2
    ## is -Inf. Masking it to 1e-12 turns an UNDEFINED ratio into a large
    ## finite number, which silently defeats the zero-denominator guard that
    ## log_of_mean has. Let -Inf through and convert non-finite results to NA
    ## at the end, so both conventions share the same NA semantics.
    Lc <- log2(FA %*% MUg)
    Ln <- log2(FB %*% MUg)
    num <- colMeans(Lc)
    den <- colMeans(Ln)
    out <- as.numeric(num - den)
    ## -Inf - (-Inf) is NaN; -Inf - finite is -Inf. Both mean "undefined".
    out[!is.finite(out)] <- NA_real_
  }
  names(out) <- genes

  attr(out, "f_ctrl")     <- setNames(f_ctrl, ct)
  attr(out, "f_case")     <- setNames(f_case, ct)
  attr(out, "convention") <- convention
  out
}


#' Decompose observed log2 fold changes into composition and a residual part
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
#' @param convention Sample-aggregation convention; see
#'   [compose_prediction()]. Stored on the result as
#'   `attr(x, "convention")`. **The convention changes the answer**: on
#'   GSE232911, `mean_of_log` gives a comp-obs correlation of 0.2846 and
#'   `log_of_mean` gives 0.4251. It must be reported alongside any
#'   `comp`-`obs` agreement statistic.
#'
#' @return A data frame of class `cade_decomp` with columns `gene`, `obs`,
#'   `pred`, `resid`, `pct`, `mean_expr`. `pct` is `resid / obs` and is `NA`
#'   when `obs` is near zero. Attributes `f_ctrl` and `f_case` hold the mean
#'   composition of each group.
#'
#' @section What `resid` is and is not:
#' `resid` is `obs - pred`. Writing `d_case` and `d_ctrl` for the two groups'
#' deviations from the reference, the residual is
#' `resid = d_case - d_ctrl`.
#'
#' This is an algebraic identity. It holds for every input, including a
#' reference matrix that is pure noise (verified: max error 9e-16 with a random
#' reference), so it carries no evidence about the reference and no evidence
#' that the residual measures regulation within cells.
#'
#' `d_ctrl` is never estimated. Whenever the control samples differ from the
#' reference, `resid` is shifted by that difference. With zero within-cell
#' change in both groups, a reference whose cell-type profiles are swapped
#' produces residuals as large as -5.5. Treat `resid` as an unexplained
#' remainder, not as a within-cell effect, unless the control group can be
#' shown to match the reference.
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
                                  validate = TRUE,
                                  convention = c("mean_of_log", "log_of_mean")) {

  convention <- match.arg(convention)

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
  if (anyNA(names(obs)) || any(!nzchar(names(obs))))
    stop("`obs` has missing or empty gene names.")
  if (anyDuplicated(names(obs)))
    stop("`obs` has duplicated gene names; aggregate them first.")
  if (all(is.na(obs)))
    stop("`obs` is entirely NA; nothing to decompose.")

  pred <- compose_prediction(frac, mu, is_case, genes = names(obs),
                            convention = convention)

  ## Align by NAME. `obs - pred` would align by POSITION and, when the
  ## lengths differ by an integer multiple, recycle silently and give
  ## wrong residuals with no error.
  stopifnot(identical(names(pred), names(obs)))
  resid <- obs[names(pred)] - pred
  pct <- ifelse(is.na(obs) | abs(obs) < 1e-12, NA_real_, resid / obs)

  if (is.null(mean_expr)) {
    fc <- attr(pred, "f_ctrl"); fa <- attr(pred, "f_case")
    pooled <- (fc + fa) / 2
    mean_expr <- as.numeric(crossprod(pooled, mu[names(fc), names(pred),
                                                 drop = FALSE]))
    names(mean_expr) <- names(pred)
  } else {
    if (is.null(names(mean_expr)))
      stop("`mean_expr` must be named; an unnamed vector would silently ",
           "produce NA for every gene.", call. = FALSE)
    miss <- setdiff(names(pred), names(mean_expr))
    if (length(miss))
      stop(sprintf("`mean_expr` is missing %d gene(s) (e.g. %s).",
                   length(miss), paste(head(miss, 5), collapse = ", ")),
           call. = FALSE)
    mean_expr <- mean_expr[names(pred)]
  }

  out <- data.frame(gene = names(pred),
                    obs   = as.numeric(obs),
                    pred  = as.numeric(pred),
                    resid = as.numeric(resid),
                    pct   = as.numeric(pct),
                    mean_expr = as.numeric(mean_expr),
                    stringsAsFactors = FALSE)
  attr(out, "f_ctrl")     <- attr(pred, "f_ctrl")
  attr(out, "f_case")     <- attr(pred, "f_case")
  attr(out, "convention") <- convention
  class(out) <- c("cade_decomp", "data.frame")
  out
}
