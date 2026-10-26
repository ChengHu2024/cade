# ============================================================================
# Consistency diagnostics (added in CADE 4.0.0)
#
# Two functions:
#   marker_consistency()      — does the composition term match the observed
#                               change on cell-type marker genes?
#   cross_cohort_agreement()  — how reproducible are the composition term and
#                               the observed change across cohorts?
#
# Why these exist: the composition term depends only on the reference and the
# estimated proportions, while the observed log2 fold change depends on the
# expression quantification. In four bulk cohorts of the same disease, the
# composition term was reproducible at r = 0.83-0.93 while the observed change
# was reproducible at only r = 0.34-0.59. Any residual inherits the latter.
# ============================================================================

#' Default marker gene sets
#'
#' A named list of cell-type marker genes used by \code{\link{marker_consistency}}.
#' These are widely used, canonical markers, not data-derived.
#'
#' @return A named list; each element is a character vector of gene symbols.
#' @examples
#' str(CADE_markers, max.level = 1)
#' @export
CADE_markers <- list(
  VSMC        = c("ACTA2","MYH11","CNN1","TAGLN","MYOCD","DES","NOTCH3","PDGFRB",
                  "LMOD1","SYNPO2","MYL9","CALD1"),
  Fibroblast  = c("COL1A1","COL1A2","COL3A1","DCN","LUM","FBLN1","PDGFRA","COL6A1"),
  Endothelial = c("PECAM1","VWF","CDH5","CLDN5","EGFL7","KDR","TEK"),
  Macrophage  = c("CD68","CD163","MRC1","CSF1R","MSR1","C1QA"),
  T_cell      = c("CD3D","CD3E","CD2","CD3G","LCK","TRBC2"),
  B_cell      = c("CD79A","MS4A1","CD19","BANK1","CD79B"),
  Plasma      = c("MZB1","JCHAIN","IGHG1","IGKC","DERL3"),
  Monocyte    = c("LYZ","FCN1","S100A8","S100A9","VCAN"),
  Mast        = c("TPSAB1","TPSB2","CPA3","KIT","MS4A2"),
  NK_cell     = c("NKG7","GNLY","KLRD1","PRF1","NCAM1")
)

#' Marker consistency between the composition term and the observed change
#'
#' For each cell type, computes the median ratio of the observed log2 fold change
#' to the composition term across that cell type's marker genes. The ratio should
#' be near 1 when the composition term is correctly scaled for that cell type.
#'
#' @param obs Named numeric vector of observed log2 fold changes.
#' @param comp Named numeric vector of composition terms, same names as
#'   \code{obs}. For a gene not in \code{names(comp)} the gene is skipped.
#' @param markers Named list of marker gene sets, for example
#'   \code{\link{CADE_markers}}.
#' @param min_markers Minimum number of usable markers per cell type.
#' @param min_abs_comp Markers whose composition term is smaller in absolute
#'   value than this are excluded, because the ratio is unstable there.
#' @param pass_range Two-element numeric range considered a pass.
#'
#' @return A data frame with one row per cell type and columns
#'   \code{celltype}, \code{n_markers}, \code{ratio_median},
#'   \code{ratio_q25}, \code{ratio_q75}, \code{obs_median},
#'   \code{comp_median}, \code{markers}, \code{verdict}.
#'
#'   \code{verdict} is \code{"pass"} when \code{ratio_median} lies inside
#'   \code{pass_range}, \code{"sign_flipped"} when it is negative, and
#'   \code{"off"} otherwise.
#'
#' @details
#' The ratio is a calibration statistic, not a test of reference quality.
#' It equals the ratio of the observed effect scale to the composition scale in
#' a cohort, so it is sensitive to platform: cohorts whose log2 fold changes are
#' compressed (typically microarray) give ratios below 1, and RNA-seq cohorts
#' with larger fold changes give ratios near or above 1. Compare cohorts of the
#' same platform, or treat the ratio as a per-cohort calibration read-out.
#'
#' A single marker gene may itself be regulated, so the median over several
#' markers is used, never a single gene.
#'
#' @examples
#' set.seed(1)
#' g <- paste0("g", 1:200)
#' comp <- setNames(rnorm(200, 0, 0.5), g)
#' obs  <- setNames(comp * 1.1 + rnorm(200, 0, 0.1), g)
#' mk <- list(A = g[1:10], B = g[11:20])
#' # give the markers large composition terms so the ratio is stable
#' comp[mk$A] <- runif(10, 0.5, 1.5); comp[mk$B] <- runif(10, -1.5, -0.5)
#' obs[mk$A]  <- comp[mk$A] * 1.05;    obs[mk$B]  <- comp[mk$B] * 0.95
#' marker_consistency(obs, comp, mk)
#' @export
marker_consistency <- function(obs, comp, markers = CADE_markers,
                               min_markers = 3, min_abs_comp = 0.05,
                               pass_range = c(0.5, 2.0)) {
  if (is.null(names(obs)) || is.null(names(comp)))
    stop("obs and comp must be named vectors")
  if (!is.list(markers) || is.null(names(markers)))
    stop("markers must be a named list of gene vectors")
  if (length(pass_range) != 2 || any(!is.finite(pass_range)))
    stop("pass_range must be two finite numbers, e.g. c(0.5, 2.0)")
  if (pass_range[1] > pass_range[2])
    stop("pass_range must be increasing")

  out <- lapply(names(markers), function(ct) {
    mg <- intersect(markers[[ct]], intersect(names(obs), names(comp)))
    if (length(mg) < min_markers) return(NULL)
    o <- obs[mg]; c <- comp[mg]
    keep <- is.finite(o) & is.finite(c) & abs(c) > min_abs_comp
    if (sum(keep) < min_markers) return(NULL)
    ratio <- o[keep] / c[keep]
    rm_ <- stats::median(ratio)
    data.frame(
      celltype = ct,
      n_markers = sum(keep),
      ratio_median = rm_,
      ratio_q25 = stats::quantile(ratio, .25, names = FALSE),
      ratio_q75 = stats::quantile(ratio, .75, names = FALSE),
      obs_median = stats::median(o[keep]),
      comp_median = stats::median(c[keep]),
      markers = paste(mg[keep], collapse = ","),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out) || !nrow(out))
    stop("no cell type had enough usable markers")
  out$verdict <- ifelse(out$ratio_median >= pass_range[1] &
                          out$ratio_median <= pass_range[2], "pass",
                 ifelse(out$ratio_median < 0, "sign_flipped", "off"))
  out
}

#' Cross-cohort agreement of the composition term and the observed change
#'
#' Compares two or more cohorts gene by gene and reports, for each pair, the
#' Pearson and Spearman correlation of the observed log2 fold changes and of the
#' composition terms, together with the fraction of genes whose signs agree.
#'
#' @param obs_list Named list; each element a named numeric vector of observed
#'   log2 fold changes for one cohort.
#' @param comp_list Named list of composition terms, in the same order and with
#'   the same names as \code{obs_list}. Optional; when absent only the observed
#'   change is compared.
#' @param min_common Minimum number of shared genes to report a pair.
#'
#' @return A list with \code{obs} and (when \code{comp_list} is given)
#'   \code{comp}, each a data frame with columns \code{a}, \code{b},
#'   \code{n}, \code{pearson}, \code{spearman}, \code{sign_agree}.
#'
#' @details
#' In four bulk cohorts of one disease, the composition term agreed at
#' r = 0.83-0.93 across cohorts while the observed change agreed at only
#' r = 0.34-0.59. The practical consequence is that a decomposition residual,
#' which is the difference of the two, is limited by the reproducibility of the
#' observed change, not by that of the composition term. Cohorts measured on the
#' same platform agree better; grouping by platform before comparing residuals
#' is advisable.
#'
#' @section A per-cohort rescaling does not repair the disagreement:
#' Standardising each cohort (z-scores, or percentile ranks) before comparing
#' them does not help. Both are monotone transformations, so Spearman correlation
#' is unchanged by construction, and Pearson after z-scoring is numerically
#' identical to Pearson on the raw values because Pearson is invariant to the
#' affine map that z-scoring applies. Verified on four bulk cohorts: raw and
#' z-scored Pearson both span 0.346-0.595, and Spearman is 0.3561 under raw,
#' z-scored and rank transforms alike.
#'
#' Percentile ranks are actively harmful for sign-based statistics: every rank is
#' positive, so a sign-agreement computed on ranks is 100 percent by
#' construction and carries no information.
#'
#' What does help is stratification. Grouping cohorts by measurement platform
#' raised the observed-change correlation from 0.346-0.400 across platforms to
#' 0.595 for array versus array and 0.590 for RNA-seq versus RNA-seq.
#'
#' @examples
#' set.seed(1)
#' g <- paste0("g", 1:500)
#' c1 <- setNames(rnorm(500, 0, .4), g)
#' o1 <- setNames(c1 + rnorm(500, 0, .2), g)
#' o2 <- setNames(c1 * 1.05 + rnorm(500, 0, .25), g)
#' res <- cross_cohort_agreement(list(A = o1, B = o2), list(A = c1, B = c1))
#' res$obs
#' res$comp
#' @export
cross_cohort_agreement <- function(obs_list, comp_list = NULL,
                                   min_common = 100) {
  if (!is.list(obs_list) || length(obs_list) < 2)
    stop("obs_list must be a named list with at least two elements")
  nm <- names(obs_list)
  if (is.null(nm) || any(!nzchar(nm)))
    stop("obs_list elements must be named")
  if (anyDuplicated(nm))
    stop("obs_list element names must be unique")

  pair_df <- function(a, b) {
    g <- intersect(names(a), names(b))
    if (length(g) < min_common) return(NULL)
    x <- a[g]; y <- b[g]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < min_common) return(NULL)
    x <- x[ok]; y <- y[ok]
    data.frame(n = length(x),
               pearson  = stats::cor(x, y),
               spearman = stats::cor(x, y, method = "spearman"),
               sign_agree = 100 * mean(sign(x) == sign(y)),
               stringsAsFactors = FALSE)
  }

  do_pairs <- function(L) {
    rows <- list()
    for (i in seq_len(length(L) - 1)) for (j in seq(i + 1, length(L))) {
      d <- pair_df(L[[i]], L[[j]])
      if (is.null(d)) next
      rows[[length(rows) + 1]] <- cbind(data.frame(a = names(L)[i], b = names(L)[j],
                                                   stringsAsFactors = FALSE), d)
    }
    if (!length(rows)) return(NULL)
    do.call(rbind, rows)
  }

  out <- list(obs = do_pairs(obs_list))
  if (!is.null(comp_list)) {
    if (length(comp_list) != length(obs_list))
      stop("comp_list must have the same length as obs_list")
    if (!identical(names(comp_list), names(obs_list)))
      stop("comp_list names must match obs_list names in the same order; ",
           "otherwise cohorts are mislabelled in the output", call. = FALSE)
    out$comp <- do_pairs(comp_list)
  }
  out
}
