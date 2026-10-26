#' Anatomical plausibility gate for deconvolved fractions
#'
#' @description
#' Reject a deconvolution whose fractions contradict basic histology. A
#' deconvolution that reports 0\% endothelium in lung tissue, or 0\% adipocytes
#' in normal breast, is not a mediocre estimate: it is a failed one, and any
#' statistic computed from it is meaningless. This gate makes that failure loud
#' instead of silent.
#'
#' @param fractions Numeric matrix or data frame of samples x cell types, each
#'   row summing to about 1.
#' @param rules A named list. Each element is a list with \code{types} (a
#'   character vector of column names or a regex), \code{min} (minimum acceptable
#'   summed fraction), and optionally \code{samples} (a logical vector or
#'   character vector selecting which rows the rule applies to) and
#'   \code{label} (a human-readable name for the rule).
#'   \code{anatomical_rules()} supplies defaults for a few tissues.
#' @param tissue Optional tissue name; if given and \code{rules} is NULL,
#'   \code{anatomical_rules(tissue)} is used.
#' @param action One of \code{"error"} (default), \code{"warn"} or
#'   \code{"report"}.
#'
#' @return Invisibly, a data frame with one row per rule: the rule label, the
#'   observed minimum fraction across the selected samples, the threshold, and
#'   whether the rule passed.
#'
#' @details
#' \strong{Why this is needed.} In four attempts to validate the composition
#' term on real bulk data, the deconvolution reported 0.05\% endothelial and
#' 0.03\% lymphoid fraction in lung tissue, and 0\% adipocyte fraction in normal
#' breast. Both are anatomically impossible. The downstream statistics were
#' computed before anyone noticed, because nothing checked.
#'
#' \strong{The rules are tissue-specific and must be justified.} Lung
#' parenchyma is roughly 10-20\% endothelium by volume; normal breast is
#' 30-80\% adipose. The defaults in \code{anatomical_rules()} are order-of-
#' magnitude floors, deliberately set well below the literature values so that
#' they catch gross failure rather than flagging ordinary variation. They are
#' not a substitute for a tissue-specific argument.
#'
#' \strong{These thresholds are not a test of the method.} Passing the gate does
#' not mean the fractions are correct; it means they are not obviously
#' impossible. A deconvolution can pass every rule here and still be badly
#' biased.
#'
#' @examples
#' set.seed(1)
#' f <- matrix(runif(36, 0.2, 1), nrow = 6,
#'             dimnames = list(NULL, c("epithelial", "endothelial",
#'                                     "fibroblast", "T_cell", "myeloid", "other")))
#' f <- f / rowSums(f)
#' # passes: every compartment is present
#' anatomical_plausibility(f, rules = list(
#'   endo = list(types = "endothelial", min = 0.01),
#'   lym  = list(types = "T_cell", min = 0.01)), action = "report")
#' # fails: no endothelium at all
#' f2 <- f
#' f2[, "endothelial"] <- 0
#' f2 <- f2 / rowSums(f2)
#' try(anatomical_plausibility(f2, rules = list(
#'   endo = list(types = "endothelial", min = 0.01))))
#'
#' @seealso \code{\link{marker_consistency}}, \code{\link{anatomical_rules}}
#' @export
anatomical_plausibility <- function(fractions, rules = NULL, tissue = NULL,
                                    action = c("error", "warn", "report")) {
  action <- match.arg(action)
  if (!is.null(tissue) && is.null(rules)) rules <- anatomical_rules(tissue)
  if (is.null(rules)) {
    stop("supply rules, or a tissue with built-in rules (",
         paste(names(anatomical_rules()), collapse = ", "), ")", call. = FALSE)
  }

  f <- as.matrix(fractions)
  if (!is.numeric(f)) stop("fractions must be numeric", call. = FALSE)
  if (anyNA(f)) stop("fractions contain NA; a gate cannot be applied to missing values", call. = FALSE)
  if (any(f < 0, na.rm = TRUE)) stop("fractions must be non-negative", call. = FALSE)

  rs <- rowSums(f)
  if (any(rs <= 0)) stop("some rows of fractions sum to zero", call. = FALSE)
  if (any(abs(rs - 1) > 0.05)) {
    warning("fractions do not sum to 1 (range ", sprintf("%.3f-%.3f", min(rs), max(rs)),
            "); normalising rows before the gate", call. = FALSE)
    f <- f / rs
  }

  out <- vector("list", length(rules))
  names(out) <- names(rules)
  for (nm in names(rules)) {
    r <- rules[[nm]]
    if (is.null(r$types) || is.null(r$min)) {
      stop("rule '", nm, "' needs both `types` and `min`", call. = FALSE)
    }
    ## 先按精确名匹配；一个都没匹配上时，再按正则匹配。
    ## 这样 "endothelial" 既能匹配名叫 endothelial 的列，
    ## 也能匹配 "endothelial cell of artery" 这类列。
    cols <- intersect(as.character(r$types), colnames(f))
    if (!length(cols)) {
      cols <- unique(unlist(lapply(as.character(r$types), function(tt) {
        grep(tt, colnames(f), value = TRUE)
      })))
    }
    if (!length(cols)) {
      stop("rule '", nm, "': none of the named cell types are present. ",
           "Available: ", paste(colnames(f), collapse = ", "), call. = FALSE)
    }
    sel <- rep(TRUE, nrow(f))
    if (!is.null(r$samples)) {
      if (is.logical(r$samples)) {
        if (length(r$samples) != nrow(f)) {
          stop("rule '", nm, "': `samples` has length ", length(r$samples),
               " but fractions has ", nrow(f), " rows", call. = FALSE)
        }
        sel <- r$samples
      } else {
        rn <- rownames(f)
        if (is.null(rn)) stop("rule '", nm, "': `samples` given by name but fractions has no rownames",
                              call. = FALSE)
        sel <- rn %in% as.character(r$samples)
      }
      if (!any(sel)) stop("rule '", nm, "': `samples` selects no rows", call. = FALSE)
    }
    val <- min(rowSums(f[sel, cols, drop = FALSE]))
    out[[nm]] <- data.frame(
      rule = if (is.null(r$label)) nm else r$label,
      types = paste(cols, collapse = ", "),
      n_samples = sum(sel),
      min_fraction = val,
      threshold = r$min,
      pass = val >= r$min,
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL

  bad <- res[!res$pass, , drop = FALSE]
  if (nrow(bad)) {
    msg <- paste0(
      "anatomically implausible fractions:\n",
      paste(sprintf("  - %s: %s is %.4f in the worst sample, floor is %.4f",
                    bad$rule, bad$types, bad$min_fraction, bad$threshold),
            collapse = "\n"),
      "\n  A proportion that contradicts basic histology is not a mediocre estimate,",
      "\n  it is a failed one. Do not compute downstream statistics from it."
    )
    if (action == "error") stop(msg, call. = FALSE)
    if (action == "warn") warning(msg, call. = FALSE)
  }
  invisible(res)
}

#' Built-in anatomical floor rules
#'
#' @param tissue One of \code{"lung"}, \code{"breast"}, \code{"aorta"}, or
#'   \code{"blood"}. If missing, returns the list of available tissues.
#'
#' @details
#' The floors are order-of-magnitude lower than published tissue composition,
#' deliberately: they catch gross failure (a compartment reported as absent when
#' it is abundant) rather than flagging ordinary biological variation.
#'
#' References for the expected ranges: lung endothelium and epithelium are
#' major parenchymal compartments, with alveolar epithelium typically 20-40\% of
#' cells in parenchyma by volume; normal breast tissue is 30-80\% adipose; aortic
#' media is smooth-muscle dominated.
#'
#' @return A named list of rules suitable for \code{\link{anatomical_plausibility}}.
#' @examples
#' names(anatomical_rules())
#' str(anatomical_rules("lung"), max.level = 2)
#' @export
anatomical_rules <- function(tissue = c("lung", "breast", "aorta", "blood")) {
  tissue <- match.arg(tissue, several.ok = FALSE)
  switch(tissue,
    lung = list(
      endothelium = list(types = "endothelial", min = 0.03,
                         label = "lung must contain endothelium"),
      lymphoid    = list(types = "T[ _]cell|B[ _]cell|NK|natural killer|plasma|mast",
                         min = 0.01, label = "lung must contain lymphocytes"),
      epithelium  = list(types = "alveolar|basal|club|goblet|ciliated|epithelial",
                         min = 0.05, label = "lung must contain epithelium")
    ),
    breast = list(
      adipocyte   = list(types = "adipocyte", min = 0.05,
                         label = "normal breast must contain adipocytes"),
      epithelium  = list(types = "epithelial|malignant|basal|luminal", min = 0.05,
                         label = "breast must contain epithelium"),
      endothelium = list(types = "endothelial", min = 0.01,
                         label = "breast must contain endothelium")
    ),
    aorta = list(
      vsmc        = list(types = "VSMC|smooth[ _]muscle", min = 0.05,
                         label = "aortic media must contain smooth muscle"),
      endothelium = list(types = "endothelial", min = 0.01,
                         label = "aorta must contain endothelium")
    ),
    blood = list(
      lymphoid    = list(types = "T[ _]cell|B[ _]cell|NK|plasma", min = 0.05,
                         label = "blood must contain lymphocytes"),
      myeloid     = list(types = "monocyte|macrophage|myeloid", min = 0.01,
                         label = "blood must contain myeloid cells")
    )
  )
}
