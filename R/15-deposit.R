# ============================================================================
# cade: verification of a deposited expression matrix
#
# Why this file exists
# --------------------
# A deposit is a file; an analysis is a program. The two can disagree without
# either one containing a step that reports an error. If the rows of an
# intensity block are not in the same order as the probe-ID column, then every
# number is attached to the wrong gene, and every gene-level conclusion is
# wrong, while the matrix still looks perfectly ordinary: the sample
# correlation matrix, the clustering, and the mean-variance relationship are
# all unchanged by a row permutation, because they do not use row identity.
#
# The three checks below are deliberately orthogonal:
#
#   A1  within-column monotone association of intensity and detection
#       p-values. A strong association is expected whenever both blocks come
#       from the same array. This check is invariant to any *common* row
#       permutation of the two blocks, so it does not see a whole-block
#       permutation; it sees a *relative* permutation of the two blocks.
#
#   A2  row-mean fingerprint against an order-aligned reference. A row
#       permutation moves each row mean to a different position, so this
#       check does see a whole-block permutation. It needs an input whose row
#       order is known to be correct.
#
#   A3  additive identity against a second representation of the same data.
#       Residuals at file-rounding scale mean that the two matrices agree
#       element by element up to an additive row effect plus an additive
#       column effect. A row permutation of one block breaks that identity.
#
# Reporting convention
# --------------------
# This file reports numerical verdicts only. It does not infer why a deposit
# is misaligned. Following the wording used for the original case, the
# conclusion is at most "the deposited file does not faithfully represent the
# analysed data"; nothing here states how that happened or who caused it.
# ============================================================================


## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

## Coerce one input to a finite numeric matrix, or stop with a message that
## names the argument.
.deposit_matrix <- function(x, what) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x))
    stop("`", what, "` is not numeric; a numeric matrix is required.",
         call. = FALSE)
  storage.mode(x) <- "double"
  if (nrow(x) < 2L)
    stop("`", what, "` has ", nrow(x), " row(s); at least 2 are required.",
         call. = FALSE)
  if (ncol(x) < 1L)
    stop("`", what, "` has no columns.", call. = FALSE)
  if (any(!is.finite(x)))
    stop("`", what, "` contains ", sum(!is.finite(x)),
         " non-finite value(s); handle missing values before verification.",
         call. = FALSE)
  x
}

## Two-way centring: subtract the row means and the column means, add back the
## grand mean. The result is orthogonal to every row-constant and every
## column-constant matrix, which is exactly the nuisance space of the additive
## model used by A3.
.deposit_double_center <- function(x) {
  x - rowMeans(x) - rep(colMeans(x), each = nrow(x)) + mean(x)
}

## Column labels for reports; synthetic labels when the input has none.
.deposit_colnames <- function(x) {
  cn <- colnames(x)
  if (is.null(cn)) paste0("V", seq_len(ncol(x))) else cn
}

## Pairs of columns that are identical value for value. Two columns are the
## same when max(abs(x - y)) == 0, i.e. exact equality on every element.
## The full comparison is O(m^2) in the number of columns and is used for
## m <= 100. Above that the columns are hashed first: a column is represented
## by its %.17g strings, which round-trip to the same double, so only columns
## that are truly identical share a hash. The hash is O(n * m).
.deposit_identical_column_pairs <- function(x, full_scan_max = 100L) {
  x <- as.matrix(x)
  m <- ncol(x)
  if (m < 2L) return(list())

  if (m <= full_scan_max) {
    out <- list()
    for (j in seq_len(m - 1L)) {
      for (k in seq.int(j + 1L, m)) {
        if (all(x[, j] == x[, k])) out[[length(out) + 1L]] <- c(j, k)
      }
    }
    return(out)
  }

  keys <- apply(x, 2, function(v) paste(sprintf("%.17g", v), collapse = "\r"))
  groups <- split(seq_len(m), keys)
  out <- list()
  for (g in groups) {
    if (length(g) < 2L) next
    for (j in seq_len(length(g) - 1L)) {
      for (k in seq.int(j + 1L, length(g))) out[[length(out) + 1L]] <- g[c(j, k)]
    }
  }
  out
}

## One key per row, built from the row's values in up to `max_cols` columns.
## More columns make keys more specific but cost memory; ten columns are enough
## to identify a row in practice and keep a 40,000-row deposit small.
## %.17g is the shortest format guaranteed to round-trip a double, so equality
## of keys is equality of the underlying numbers.
.deposit_row_keys <- function(x, max_cols = 10L) {
  x <- as.matrix(x)
  k <- min(ncol(x), max_cols)
  parts <- lapply(seq_len(k), function(j) sprintf("%.17g", x[, j]))
  do.call(paste, c(parts, list(sep = "\r")))
}

## Recover the permutation that maps the rows of `reference_intensity` onto the
## rows of `intensity`. Returns pi such that intensity[pi, ] is in the reference
## row order. The mapping must be unique (no duplicate row keys), complete
## (every reference row has a counterpart) and bijective.
##
## This is a repair-oriented diagnostic: when it succeeds, the files are known
## to contain the same row tuples and the misalignment is fully described by
## pi. When it fails, the reason is reported instead of guessing.
.deposit_recover_permutation <- function(intensity, reference_intensity) {
  n <- nrow(intensity)
  if (ncol(reference_intensity) != ncol(intensity)) {
    return(list(ok = FALSE, pi = NULL, is_identity = NA,
                source = "reference_intensity",
                reason = sprintf(paste0(
                  "`reference_intensity` has %d column(s) but `intensity` has ",
                  "%d; row tuples are not comparable"),
                  ncol(reference_intensity), ncol(intensity))))
  }

  ki <- .deposit_row_keys(intensity)
  kr <- .deposit_row_keys(reference_intensity)

  if (anyDuplicated(ki)) {
    return(list(ok = FALSE, pi = NULL, is_identity = NA,
                source = "reference_intensity",
                reason = sprintf(paste0(
                  "`intensity` has %d duplicate row tuple(s); the permutation ",
                  "is not identifiable"), sum(duplicated(ki)))))
  }
  if (anyDuplicated(kr)) {
    return(list(ok = FALSE, pi = NULL, is_identity = NA,
                source = "reference_intensity",
                reason = sprintf(paste0(
                  "`reference_intensity` has %d duplicate row tuple(s); the ",
                  "permutation is not identifiable"), sum(duplicated(kr)))))
  }

  pi <- match(kr, ki)
  if (anyNA(pi)) {
    return(list(ok = FALSE, pi = NULL, is_identity = NA,
                source = "reference_intensity",
                reason = sprintf(paste0(
                  "%d row(s) of `reference_intensity` have no exact counterpart ",
                  "in `intensity`"), sum(is.na(pi)))))
  }
  if (anyDuplicated(pi)) {
    return(list(ok = FALSE, pi = NULL, is_identity = NA,
                source = "reference_intensity",
                reason = "the row mapping is not a bijection"))
  }

  list(ok = TRUE, pi = pi,
       is_identity = all(pi == seq_len(n)),
       source = "reference_intensity",
       reason = sprintf("complete bijection over %d row(s); %d row(s) moved",
                        n, sum(pi != seq_len(n))))
}


## ---------------------------------------------------------------------------
## The verifier
## ---------------------------------------------------------------------------

#' Verify that a deposited matrix is internally consistent
#'
#' @description
#' Runs three orthogonal checks on a deposited expression matrix and returns a
#' verdict with a printable report. Check A1 tests the within-column monotone
#' association between intensity and detection p-values; A2 tests a row-mean
#' fingerprint against an order-aligned reference; A3 tests an additive
#' identity against a second representation of the same data. When a
#' reference-intensity matrix is supplied, the function also attempts to
#' recover the row permutation that would bring the two matrices into the same
#' order. The function reports numerical verdicts only and makes no inference
#' about the cause of a misalignment.
#'
#' @param intensity Numeric matrix of expression intensities, rows = probes,
#'   columns = samples. Required.
#' @param detection_p Numeric matrix of detection p-values with the same
#'   dimensions as \code{intensity}. Required only when check \code{"A1"} is
#'   requested.
#' @param probe_ids Optional atomic vector of probe identifiers, one per row of
#'   \code{intensity}. Used for length validation and bookkeeping; row identity
#'   is not inferred from it.
#' @param reference Optional numeric matrix of the same dimensions as
#'   \code{intensity}: a second representation of the same data (for example a
#'   series matrix). Required for check \code{"A3"}.
#' @param reference_intensity Optional numeric matrix whose row means can be
#'   compared with those of \code{intensity}; it must have the same number of
#'   rows and is assumed to be in the correct row order. Required for check
#'   \code{"A2"} and for permutation recovery.
#' @param checks Character vector, any subset of \code{"A1"}, \code{"A2"},
#'   \code{"A3"}. Default all three.
#' @param min_rho Minimum \strong{absolute} Spearman correlation between
#'   intensity and detection p-values for each column, default 0.7. The
#'   absolute value is used because the sign depends on whether the file stores
#'   p-values (low = detected, negative association) or detection scores
#'   (positive association).
#' @param min_fingerprint Minimum Spearman correlation between the row means of
#'   \code{intensity} and of \code{reference_intensity}, default 0.9.
#' @param tolerance_digits Residual SD of the additive identity must not exceed
#'   \code{10^(-tolerance_digits)}, default 4, i.e. \code{1e-4}.
#'
#' @return A list of class \code{"cade_deposit"} with elements:
#'   \itemize{
#'     \item \code{a1}: \code{pass}, \code{rho_per_column} (signed Spearman
#'       correlations), \code{n_fail}, \code{worst_column},
#'       \code{sign_convention} (\code{"positive"}, \code{"negative"} or
#'       \code{"mixed"}, from the median sign of the column correlations).
#'     \item \code{a2}: \code{pass}, \code{rho}, \code{reference_note}.
#'     \item \code{a3}: \code{pass}, \code{slope}, \code{residual_sd},
#'       \code{max_resid}, \code{unexplained_var}, \code{value_identical_columns}.
#'     \item \code{verdict}: \code{"ok"}, \code{"suspect"} or \code{"fail"}.
#'     \item \code{summary}: a character string ready to be printed.
#'     \item \code{n_probes}, \code{n_samples}: dimensions of \code{intensity}.
#'     \item \code{permutation_recoverable}: \code{TRUE} when a complete
#'       bijection between the rows of \code{intensity} and
#'       \code{reference_intensity} was found (including the identity, when
#'       no rows moved).
#'     \item \code{permutation}: a list with \code{pi} (the recovered
#'       permutation, or \code{NULL}), \code{is_identity}, \code{source} and
#'       \code{reason}. The direction convention is fixed:
#'       \code{reference[i, ] == intensity[pi[i], ]} for every row \code{i},
#'       so the repaired matrix is \code{intensity[pi, ]}, not
#'       \code{intensity[order(pi), ]}.
#'   }
#'   Additional bookkeeping fields (\code{checks_run},
#'   \code{checks_requested}, \code{thresholds}) are included.
#'
#' @details
#' \strong{What each check can and cannot see.}
#'
#' A1 is a per-column Spearman correlation between \code{intensity} and
#' \code{detection_p}. A real array deposit shows a strong association: the
#' detection p-value of a probe is determined by the negative-control
#' distribution and therefore moves with the signal. The correlation is
#' \strong{invariant to a common row permutation of the two blocks}: permuting
#' both matrices by the same vector reorders the pairs within every column and
#' leaves each rank correlation unchanged. A1 therefore detects a
#' \emph{relative} row permutation between the two blocks, but it is blind by
#' construction to a permutation that moves both blocks together.
#'
#' The \emph{sign} of that association is an encoding convention, not a
#' property of the data. Raw Illumina detection p-values fall as signal rises
#' (negative association), while a file that stores a detection score or
#' \code{-log10(p)} rises with signal (positive). In the original case the
#' measured correlations were positive, which is consistent with the second
#' encoding; the sign must not be assumed for a new file. The criterion is
#' therefore \code{|rho| >= min_rho}, \code{rho_per_column} reports the signed
#' values, and \code{a1$sign_convention} reports which convention the file
#' appears to use.
#'
#' A2 compares the row means of \code{intensity} with those of an
#' \code{reference_intensity} whose row order is known to be correct. A row
#' permutation moves each row mean to a different position, so A2 detects a
#' whole-block permutation. If \code{reference_intensity} is not supplied, A2
#' is skipped and recorded as \code{NA}.
#'
#' A3 fits \code{S_ij = a * log2(I_ij) + r_i + c_j} by two-way centring,
#' where \code{S} is \code{reference} and \code{I} is \code{intensity}. With
#' \code{S_c} and \code{L_c} the two-way-centred versions of \code{S} and
#' \code{log2(I)}, the slope is \code{sum(S_c * L_c) / sum(L_c^2)} and the
#' residual is \code{S_c - a * L_c}. The check passes when the residual SD is
#' at file-rounding scale (at most \code{10^(-tolerance_digits)}) \strong{and}
#' \code{reference} contains no two value-identical columns. Identical columns
#' mean that the second representation duplicates a sample column, so a
#' perfect fit proves nothing. If \code{reference} is not supplied, A3 is
#' skipped and recorded as \code{NA}.
#'
#' \strong{The structure of the original case.} In the original deposit, 39,426
#' intensity rows were out of the row order used by the probe-ID column and by
#' the detection-p block: the intensity block alone was out of order, while the
#' probe IDs and the p-values stayed in the original order. That is the case A1
#' is designed to catch, and it is the case reproduced by check C in the tests.
#' The opposite situation, in which intensity and detection p-values move
#' together while the probe IDs stay behind, leaves every within-column pair
#' intact and is invisible to A1 by construction; it is caught by A2 (with an
#' order-aligned reference) and described by the recovered permutation.
#'
#' \strong{Permutation recovery.} When \code{reference_intensity} is supplied,
#' the function builds one key per row from the row's values (using
#' \code{\%.17g}, which round-trips a double exactly) and looks for a complete,
#' unique, bijective match between the two row sets. If one exists, \code{pi}
#' is returned with a fixed direction convention:
#' \code{reference_intensity[i, ] == intensity[pi[i], ]} for every row
#' \code{i}, so the repaired matrix is \code{intensity[pi, ]}. This is the
#' opposite of the convention in which \code{new[pi[i]] <- old[i]}; the
#' direction is stated explicitly and is asserted in the tests. Recovery fails,
#' with a stated reason, when rows are duplicated (the permutation is not
#' identifiable), when any reference row has no exact counterpart, or when the
#' two matrices have different numbers of columns. Recovery does not run
#' without \code{reference_intensity}.
#'
#' \strong{Verdict.} \code{"fail"} if any executed check fails, or if a
#' non-identity permutation is recovered (the rows are demonstrably not in the
#' reference order). \code{"suspect"} if no check fails but an executed
#' measurement sits close to its threshold, or if no requested check could be
#' executed. \code{"ok"} otherwise. A check that was requested but could not
#' run is recorded as \code{NA} and listed in the summary; \code{"ok"} is a
#' statement about the checks that ran, not about the checks that were skipped.
#'
#' \strong{No inference about cause.} The function states numerical facts about
#' the files supplied. It does not infer how a deposit became misaligned.
#' Following the wording used for the original case, the strongest conclusion
#' available is that the deposited file does not faithfully represent the
#' analysed data.
#'
#' Inputs must be finite; missing values have to be handled before calling.
#' For A3, \code{intensity} must be strictly positive because the check uses
#' \code{log2()}.
#'
#' @seealso \code{\link{check_inputs}}, \code{\link{self_test}}
#'
#' @examples
#' set.seed(1)
#' I <- matrix(exp(rnorm(80 * 4, 2, 0.5)), 80, 4,
#'             dimnames = list(paste0("p", 1:80), paste0("s", 1:4)))
#' P <- plogis(-(I - mean(I)) / sd(I))
#' res <- verify_deposit(I, P, probe_ids = rownames(I))
#' res$verdict
#' print(res)
#' @export
verify_deposit <- function(intensity = NULL, detection_p = NULL,
                           probe_ids = NULL, reference = NULL,
                           reference_intensity = NULL,
                           checks = c("A1", "A2", "A3"),
                           min_rho = 0.7, min_fingerprint = 0.9,
                           tolerance_digits = 4) {

  ## ---- argument validation ----------------------------------------------
  if (is.null(intensity))
    stop("`intensity` is required (rows = probes, columns = samples).",
         call. = FALSE)

  if (!is.character(checks) || length(checks) < 1L)
    stop("`checks` must contain at least one of \"A1\", \"A2\", \"A3\".",
         call. = FALSE)
  checks <- match.arg(toupper(checks), c("A1", "A2", "A3"), several.ok = TRUE)
  checks <- unique(checks)

  .check_unit <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0 || v > 1)
      stop("`", nm, "` must be a single number in [0, 1].", call. = FALSE)
    as.numeric(v)
  }
  min_rho         <- .check_unit(min_rho, "min_rho")
  min_fingerprint <- .check_unit(min_fingerprint, "min_fingerprint")

  if (!is.numeric(tolerance_digits) || length(tolerance_digits) != 1L ||
      !is.finite(tolerance_digits) || tolerance_digits < 0 ||
      tolerance_digits > 15)
    stop("`tolerance_digits` must be a single number between 0 and 15.",
         call. = FALSE)

  ## ---- inputs ------------------------------------------------------------
  I  <- .deposit_matrix(intensity, "intensity")
  P  <- if (!is.null(detection_p)) .deposit_matrix(detection_p, "detection_p")
  R  <- if (!is.null(reference)) .deposit_matrix(reference, "reference")
  RI <- if (!is.null(reference_intensity))
    .deposit_matrix(reference_intensity, "reference_intensity")

  if (!is.null(P) && !identical(dim(P), dim(I)))
    stop("`detection_p` must have the same dimensions as `intensity` (",
         nrow(I), " x ", ncol(I), "), got ", nrow(P), " x ", ncol(P), ".",
         call. = FALSE)
  if (!is.null(R) && !identical(dim(R), dim(I)))
    stop("`reference` must have the same dimensions as `intensity` (",
         nrow(I), " x ", ncol(I), "), got ", nrow(R), " x ", ncol(R), ".",
         call. = FALSE)
  if (!is.null(RI) && nrow(RI) != nrow(I))
    stop("`reference_intensity` must have the same number of rows as ",
         "`intensity` (", nrow(I), "), got ", nrow(RI), ".", call. = FALSE)

  if (!is.null(probe_ids)) {
    if (is.data.frame(probe_ids)) {
      if (ncol(probe_ids) != 1L)
        stop("`probe_ids` must be a single column or an atomic vector.",
             call. = FALSE)
      probe_ids <- probe_ids[[1L]]
    }
    if (!is.atomic(probe_ids))
      stop("`probe_ids` must be an atomic vector.", call. = FALSE)
    if (length(probe_ids) != nrow(I))
      stop("`probe_ids` has length ", length(probe_ids),
           " but `intensity` has ", nrow(I), " rows.", call. = FALSE)
  }

  if ("A1" %in% checks && is.null(P))
    stop("check A1 requires `detection_p`; supply it or drop \"A1\" from ",
         "`checks`.", call. = FALSE)

  if (!is.null(P) && !is.null(colnames(P)) && !is.null(colnames(I)) &&
      !identical(colnames(P), colnames(I)))
    warning("column names of `detection_p` and `intensity` differ; columns ",
            "are paired by position.", call. = FALSE)

  ## =========================================================================
  ## A1: intensity vs detection p-value, per column
  ## =========================================================================
  ## The signed correlation is reported, but the criterion uses |rho|. The
  ## sign is an encoding convention, not a property of the data: Illumina
  ## detection p-values fall as signal rises (negative association), while a
  ## file that stores a detection score rises with signal (positive). A
  ## strong monotone association in either direction is what a faithful
  ## deposit guarantees.
  a1 <- list(pass = NA, rho_per_column = numeric(0), n_fail = NA_integer_,
             worst_column = NA_character_, worst_rho = NA_real_,
             sign_convention = NA_character_, threshold = min_rho,
             note = "A1 not requested")
  if ("A1" %in% checks) {
    cn <- if (!is.null(colnames(P))) colnames(P) else .deposit_colnames(I)
    rho <- vapply(seq_len(ncol(I)), function(j) {
      suppressWarnings(stats::cor(I[, j], P[, j], method = "spearman"))
    }, numeric(1))
    names(rho) <- cn
    abs_rho <- abs(rho)
    n_fail  <- sum(!is.finite(abs_rho) | abs_rho < min_rho)
    worst   <- which.min(ifelse(is.finite(abs_rho), abs_rho, -Inf))
    ## The median sign describes which encoding the file appears to use. It is
    ## reported, never used as a pass criterion; see the |rho| rationale above.
    med_rho <- suppressWarnings(stats::median(rho, na.rm = TRUE))
    sign_convention <- if (!is.finite(med_rho)) {
      NA_character_
    } else if (med_rho > 0) {
      "positive"
    } else if (med_rho < 0) {
      "negative"
    } else {
      "mixed"
    }
    a1 <- list(
      pass = n_fail == 0L,
      rho_per_column = rho,
      n_fail = as.integer(n_fail),
      worst_column = unname(cn[worst]),
      worst_rho = unname(abs_rho[worst]),
      sign_convention = sign_convention,
      threshold = min_rho,
      note = if (n_fail == 0L) {
        "all columns passed"
      } else {
        sprintf("%d column(s) below |rho| = %.2f", n_fail, min_rho)
      }
    )
    if (any(!is.finite(rho)))
      a1$note <- paste0(a1$note, "; undefined rho (constant or fully tied ",
                        "column) counts as a failure")
  }

  ## =========================================================================
  ## A2: row-mean fingerprint against an order-aligned reference
  ## =========================================================================
  a2 <- list(pass = NA, rho = NA_real_, reference_note = "A2 not requested",
             threshold = min_fingerprint)
  if ("A2" %in% checks) {
    if (is.null(RI)) {
      a2$reference_note <- paste0("no `reference_intensity` supplied; A2 ",
                                  "skipped and recorded as NA")
    } else {
      rho2 <- suppressWarnings(stats::cor(rowMeans(I), rowMeans(RI),
                                          method = "spearman"))
      a2 <- list(
        pass = is.finite(rho2) && rho2 >= min_fingerprint,
        rho = unname(rho2),
        reference_note = sprintf(
          "`reference_intensity` is %d x %d; rows are paired with `intensity` by position",
          nrow(RI), ncol(RI)),
        threshold = min_fingerprint
      )
    }
  }

  ## =========================================================================
  ## A3: additive identity S = a * log2(I) + r_i + c_j
  ## =========================================================================
  a3 <- list(pass = NA, slope = NA_real_, residual_sd = NA_real_,
             max_resid = NA_real_, unexplained_var = NA_real_,
             value_identical_columns = list(),
             threshold = 10^(-tolerance_digits), note = "A3 not requested")
  if ("A3" %in% checks) {
    if (is.null(R)) {
      a3$note <- paste0("no `reference` supplied; A3 skipped and recorded ",
                        "as NA")
    } else if (any(I <= 0)) {
      stop("check A3 requires positive `intensity` because it uses log2(); ",
           sum(I <= 0), " value(s) are <= 0.", call. = FALSE)
    } else {
      Sc <- .deposit_double_center(R)
      Lc <- .deposit_double_center(log2(I))
      denom <- sum(Lc * Lc)
      tot   <- sum(Sc * Sc)

      if (!is.finite(denom) || denom <= 0) {
        a3$pass <- FALSE
        a3$note <- paste0("log2(intensity) has no variation after two-way ",
                          "centring; the additive identity is vacuous")
      } else if (!is.finite(tot) || tot <= 0) {
        a3$pass <- FALSE
        a3$note <- paste0("`reference` has no variation after two-way ",
                          "centring; the additive identity is vacuous")
      } else {
        slope    <- sum(Sc * Lc) / denom
        resid    <- Sc - slope * Lc
        rsd      <- stats::sd(as.vector(resid))
        max_res  <- max(abs(resid))
        unexp    <- sum(resid * resid) / tot
        pairs    <- .deposit_identical_column_pairs(R)
        rsd_ok <- is.finite(rsd) && rsd <= 10^(-tolerance_digits)
        ok <- rsd_ok && length(pairs) == 0L
        fail_bits <- c(
          if (!rsd_ok)
            sprintf("residual SD %.3g exceeds %.0e", rsd,
                    10^(-tolerance_digits)),
          if (length(pairs))
            sprintf("%d value-identical column pair(s) in `reference`",
                    length(pairs))
        )
        a3 <- list(
          pass = ok,
          slope = unname(slope),
          residual_sd = unname(rsd),
          max_resid = unname(max_res),
          unexplained_var = unname(unexp),
          value_identical_columns = pairs,
          threshold = 10^(-tolerance_digits),
          note = if (ok) {
            "residuals are at file-rounding scale and no two reference columns are identical"
          } else {
            paste(fail_bits, collapse = "; ")
          }
        )
      }
    }
  }

  ## =========================================================================
  ## Permutation recovery
  ## =========================================================================
  ## A bijective match between row tuples is positive evidence about the row
  ## order itself. It is attempted only when `reference_intensity` is present;
  ## the other reference is a different representation, so exact row tuples
  ## need not be comparable.
  rec <- list(ok = FALSE, pi = NULL, is_identity = NA,
              source = "reference_intensity",
              reason = "no `reference_intensity` supplied")
  if (!is.null(RI)) rec <- .deposit_recover_permutation(I, RI)

  ## =========================================================================
  ## Verdict
  ## =========================================================================
  hard <- c(
    A1 = isFALSE(a1$pass),
    A2 = isFALSE(a2$pass),
    A3 = isFALSE(a3$pass),
    permutation = isTRUE(rec$ok) && !isTRUE(rec$is_identity)
  )
  marginal <- c(
    A1 = isTRUE(a1$pass) && is.finite(a1$worst_rho) &&
      a1$worst_rho < min_rho + 0.05,
    A2 = isTRUE(a2$pass) && is.finite(a2$rho) &&
      a2$rho < min_fingerprint + 0.03,
    A3 = isTRUE(a3$pass) && is.finite(a3$residual_sd) &&
      a3$residual_sd > 10^(-tolerance_digits - 1)
  )
  executed <- c(A1 = !is.na(a1$pass), A2 = !is.na(a2$pass),
                A3 = !is.na(a3$pass))
  exec_names <- names(executed)[executed]
  skip_names <- setdiff(checks, exec_names)

  verdict <- if (any(hard, na.rm = TRUE)) {
    "fail"
  } else if (any(marginal, na.rm = TRUE) || !any(executed)) {
    "suspect"
  } else {
    "ok"
  }

  reasons <- character(0)
  if (isFALSE(a1$pass))
    reasons <- c(reasons, sprintf("A1: %d/%d column(s) below |rho| = %.2f",
                                  a1$n_fail, length(a1$rho_per_column), min_rho))
  if (isFALSE(a2$pass))
    reasons <- c(reasons, if (is.finite(a2$rho)) {
      sprintf("A2: rho = %.3f below %.2f", a2$rho, min_fingerprint)
    } else {
      "A2: row means are constant, so rho is undefined"
    })
  if (isFALSE(a3$pass))
    reasons <- c(reasons, sprintf("A3: %s", a3$note))
  if (isTRUE(rec$ok) && !isTRUE(rec$is_identity))
    reasons <- c(reasons, sprintf("row order: non-identity permutation recovered (%s)",
                                  rec$reason))
  if (any(marginal, na.rm = TRUE))
    reasons <- c(reasons, sprintf("near threshold: %s",
                                  paste(names(marginal)[marginal],
                                        collapse = ", ")))
  reason_txt <- if (length(reasons))
    paste0("（", paste(reasons, collapse = "；"), "）") else ""

  ## =========================================================================
  ## Printable report
  ## =========================================================================
  status_txt <- function(p) {
    if (is.na(p)) "SKIP" else if (isTRUE(p)) "PASS" else "FAIL"
  }
  verdict_txt <- switch(verdict, ok = "OK", suspect = "SUSPECT", fail = "FAIL")

  a1_line <- if (is.na(a1$pass)) {
    sprintf("  A1 强度-检出P值一致性: SKIP（%s）", a1$note)
  } else {
    sprintf(paste0("  A1 强度-检出P值一致性: %s；最差 |rho| = %.3f",
                   "（阈值 %.2f），%d/%d 列未达标，最差列 %s，符号约定 %s"),
            status_txt(a1$pass), a1$worst_rho, min_rho,
            a1$n_fail, length(a1$rho_per_column), a1$worst_column,
            a1$sign_convention)
  }
  a2_line <- if (is.na(a2$pass)) {
    sprintf("  A2 平台指纹: SKIP（%s）", a2$reference_note)
  } else {
    sprintf("  A2 平台指纹: %s；rho = %.3f（阈值 %.2f）",
            status_txt(a2$pass), a2$rho, min_fingerprint)
  }
  a3_line <- if (is.na(a3$pass)) {
    sprintf("  A3 加性恒等式: SKIP（%s）", a3$note)
  } else {
    sprintf(paste0("  A3 加性恒等式: %s；slope = %.4f，残差 SD = %.3g",
                   "（阈值 %.0e），最大残差 = %.3g，未解释方差 = %.3g，",
                   "逐值相同列 %d 对"),
            status_txt(a3$pass), a3$slope, a3$residual_sd,
            10^(-tolerance_digits), a3$max_resid, a3$unexplained_var,
            length(a3$value_identical_columns))
  }
  perm_line <- if (isTRUE(rec$ok)) {
    sprintf("  置换恢复: 成功；%s；π %s",
            rec$reason,
            if (isTRUE(rec$is_identity)) {
              "为恒等置换（行序与参考一致）"
            } else {
              "为非恒等置换，intensity[π, ] 与 reference_intensity 同序"
            })
  } else {
    sprintf("  置换恢复: 未完成（%s）", rec$reason)
  }
  coverage_line <- if (length(skip_names)) {
    sprintf("  覆盖: 已执行 %s；已请求但无法执行 %s",
            paste(exec_names, collapse = "/"), paste(skip_names, collapse = "/"))
  } else {
    sprintf("  覆盖: 已执行全部请求的检查（%s）",
            paste(exec_names, collapse = "/"))
  }
  verdict_line <- switch(
    verdict,
    ok = "  判定: OK —— 所有已执行的检查通过。",
    suspect = sprintf("  判定: SUSPECT —— 无硬性失败，但存在接近阈值的测量或无法执行的检查。%s",
                      reason_txt),
    fail = sprintf("  判定: FAIL —— 沉积文件未能忠实代表所分析的数据。%s",
                   reason_txt)
  )

  summary_lines <- c(
    "cade 沉积矩阵核验报告",
    sprintf("  矩阵: %d 探针 x %d 样本", nrow(I), ncol(I)),
    a1_line, a2_line, a3_line, perm_line, coverage_line,
    sprintf("  结论: %s", verdict_txt),
    verdict_line,
    paste0("  说明: 本报告只陈述数值判定，不对成因作任何推断；原案例的结论",
           "措辞为 “the deposited file does not faithfully represent the ",
           "analysed data”。")
  )

  out <- list(
    a1 = a1,
    a2 = a2,
    a3 = a3,
    verdict = verdict,
    summary = paste(summary_lines, collapse = "\n"),
    n_probes = nrow(I),
    n_samples = ncol(I),
    permutation_recoverable = isTRUE(rec$ok),
    permutation = list(pi = rec$pi, is_identity = rec$is_identity,
                       source = rec$source, reason = rec$reason),
    checks_run = exec_names,
    checks_requested = checks,
    thresholds = list(min_rho = min_rho,
                      min_fingerprint = min_fingerprint,
                      tolerance_digits = tolerance_digits,
                      residual_sd_max = 10^(-tolerance_digits))
  )
  class(out) <- "cade_deposit"
  out
}


#' Print method for a deposit verification
#'
#' @param x Output of \code{\link{verify_deposit}()}.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.cade_deposit <- function(x, ...) {
  lab <- function(p) if (is.na(p)) "SKIP" else if (isTRUE(p)) "PASS" else "FAIL"

  a1_measure <- if (is.na(x$a1$pass)) {
    x$a1$note
  } else {
    sprintf(paste0("worst |rho| = %.3f (threshold %.2f); %d of %d ",
                   "column(s) below; sign %s"),
            x$a1$worst_rho, x$a1$threshold, x$a1$n_fail,
            length(x$a1$rho_per_column), x$a1$sign_convention)
  }
  a2_measure <- if (is.na(x$a2$pass)) {
    x$a2$reference_note
  } else {
    sprintf("rho = %.3f (threshold %.2f)", x$a2$rho, x$a2$threshold)
  }
  a3_measure <- if (is.na(x$a3$pass)) {
    x$a3$note
  } else {
    sprintf(paste0("slope = %.4f; residual SD = %.3g (threshold %.0e); ",
                   "max residual = %.3g; %d identical column pair(s)"),
            x$a3$slope, x$a3$residual_sd, x$a3$threshold, x$a3$max_resid,
            length(x$a3$value_identical_columns))
  }

  cat("<cade deposit verification>\n")
  cat(sprintf("  %d probes x %d samples\n", x$n_probes, x$n_samples))
  cat(sprintf("  verdict: %s\n\n", toupper(x$verdict)))
  cat(sprintf("  %-3s %-5s %s\n", "chk", "state", "measurement"))
  cat(sprintf("  %-3s %-5s %s\n", "A1", lab(x$a1$pass), a1_measure))
  cat(sprintf("  %-3s %-5s %s\n", "A2", lab(x$a2$pass), a2_measure))
  cat(sprintf("  %-3s %-5s %s\n", "A3", lab(x$a3$pass), a3_measure))
  cat("\n")
  cat(x$summary, "\n", sep = "")
  invisible(x)
}
