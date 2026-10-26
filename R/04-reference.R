# ============================================================================
# cade: reference quality and dual-reference comparison
#
# A deconvolution reference is not a neutral ingredient. If every cell in it
# comes from diseased tissue, then the per-cell means already carry the disease
# signature, and "fewer cells" cannot be separated from "each cell changed".
#
# There is no way to detect this from the bulk data alone. The only defence is
# to know what the reference was built from, and - when two references exist -
# to run both and compare the conclusions.
# ============================================================================

#' Check a reference's composition and coverage
#'
#' @param cell_type Character vector, the cell-type label of every cell in the
#'   reference.
#' @param sample Optional character vector, the sample of origin for every cell.
#' @param min_cells Minimum acceptable cells per type, default 100.
#' @param group Optional character vector of disease group per cell
#'   (`"normal"` / `"disease"`); used to report whether the reference contains
#'   any normal cells.
#' @param markers Optional named list mapping cell type -> character vector of
#'   marker genes that must be present in the reference gene list. Pass the
#'   reference gene names via `genes`.
#' @param genes Optional character vector of gene names in the reference.
#'
#' @return A list with `by_type` (a data frame), `flags` (a character vector of
#'   problems found) and `ok` (logical, no flags raised).
#'
#' @details
#' Two failure modes are common and both are silent:
#'
#' * **A disease-only reference.** Reported through `group`.
#' * **A starved cell type.** A cell type represented by a handful of cells
#'   gives an unstable per-cell mean, and its fraction is then estimated with
#'   large error. Reported through `min_cells`.
#'
#' A third, marker coverage, is reported when `markers` and `genes` are given:
#' deconvolution on a gene absent from the reference returns exactly zero, which
#' looks like "not expressed" rather than "not measured".
#'
#' @examples
#' ct <- rep(c("A", "B", "C"), c(500, 50, 5))
#' check_reference(ct)
#' @export
check_reference <- function(cell_type, sample = NULL, min_cells = 100,
                            group = NULL, markers = NULL, genes = NULL) {
  cell_type <- as.character(cell_type)
  tab <- table(cell_type)
  n_total <- length(cell_type)
  by_type <- data.frame(cell_type = names(tab),
                        n_cells = as.integer(tab),
                        pct = as.numeric(tab) / n_total,
                        stringsAsFactors = FALSE)

  if (!is.null(sample)) {
    if (length(sample) != n_total) stop("`sample` length mismatch.")
    dom <- tapply(seq_along(sample), cell_type, function(i) {
      tt <- table(sample[i]); max(tt) / sum(tt)
    })
    by_type$max_single_sample_frac <- as.numeric(dom[by_type$cell_type])
  }

  flags <- character(0)
  starved <- by_type$cell_type[by_type$n_cells < min_cells]
  if (length(starved))
    flags <- c(flags, sprintf("cell type(s) with fewer than %d cells: %s",
                              min_cells, paste(starved, collapse = ", ")))

  if (any(is.na(by_type$max_single_sample_frac)) ||
      any(by_type$max_single_sample_frac > 0.5, na.rm = TRUE)) {
    bad <- by_type$cell_type[!is.na(by_type$max_single_sample_frac) &
                             by_type$max_single_sample_frac > 0.5]
    if (length(bad))
      flags <- c(flags, sprintf("cell type(s) dominated by a single sample (>50%%): %s",
                                paste(bad, collapse = ", ")))
  }

  if (!is.null(group)) {
    if (length(group) != n_total) stop("`group` length mismatch.")
    if (!any(group %in% c("normal", "control", "healthy")))
      flags <- c(flags, paste("reference contains no normal/control cells;",
                              "composition and residual effects cannot be",
                              "separated with such a reference"))
  }

  if (!is.null(markers) && !is.null(genes)) {
    missing <- lapply(markers, setdiff, y = genes)
    missing <- missing[vapply(missing, length, integer(1)) > 0]
    if (length(missing))
      flags <- c(flags, sprintf(
        "marker gene(s) absent from the reference (deconvolution returns 0 for these): %s",
        paste(sprintf("%s: %s", names(missing),
                      vapply(missing, paste, character(1), collapse = ", ")),
              collapse = "; ")))
  }

  list(by_type = by_type, flags = flags, ok = length(flags) == 0)
}



#' Run a decomposition with two references and compare the conclusions
#'
#' @param obs Named numeric vector of observed log2 fold changes.
#' @param frac_a,frac_b Fraction matrices from deconvolving with reference A and
#'   reference B.
#' @param mu_a,mu_b Per-cell mean matrices (linear) for the two references.
#' @param is_case Logical vector over rows of the fraction matrices.
#' @param validate Run [check_inputs()] on both references first. Default
#'   `TRUE`; the check runs once per reference rather than once per call.
#'
#' @return A list with `a` and `b` (the two `decompose_composition()` outputs),
#'   `agreement` (a one-row data frame with the correlation of the residuals and
#'   the sign-agreement rate) and `discordant` (a data frame of the genes whose
#'   residual sign flips).
#'
#' @details
#' The point of the comparison is not to pick a winner. Two references that
#' differ in status (disease versus normal) give different answers, and the
#' disagreement is itself the finding: a conclusion that flips between
#' references is not reportable as a property of the disease.
#'
#' Validation is performed once per reference, so a failure names which
#' reference is at fault (A or B) rather than merely "input validation failed".
#'
#' @examples
#' mu <- matrix(c(500, 50, 300, 40), nrow = 2,
#'              dimnames = list(c("A", "B"), c("g1", "g2")))
#' fr <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4,
#'              byrow = TRUE, dimnames = list(paste0("s", 1:4), c("A", "B")))
#' ob <- c(g1 = 0.5, g2 = -0.2)
#' ic <- c(FALSE, FALSE, TRUE, TRUE)
#' compare_references(ob, fr, fr, mu, mu, ic, validate = FALSE)$agreement
#' @export
compare_references <- function(obs, frac_a, frac_b, mu_a, mu_b, is_case,
                               validate = TRUE) {

  ## Validate both references together so the error message can name which
  ## one is at fault.
  if (validate) {
    ca <- check_inputs(obs, frac_a, mu_a, is_case)
    cb <- check_inputs(obs, frac_b, mu_b, is_case)
    if (!ca$ok)
      stop(paste0("reference A failed validation:\n",
                  paste0("  * ", ca$errors, collapse = "\n")), call. = FALSE)
    if (!cb$ok)
      stop(paste0("reference B failed validation:\n",
                  paste0("  * ", cb$errors, collapse = "\n")), call. = FALSE)
    ws <- unique(c(ca$warnings, cb$warnings))
    for (w in ws) warning(w, call. = FALSE)
  }

  da <- decompose_composition(obs, frac_a, mu_a, is_case, validate = FALSE)
  db <- decompose_composition(obs, frac_b, mu_b, is_case, validate = FALSE)

  g <- intersect(da$gene, db$gene)
  da <- da[match(g, da$gene), ]
  db <- db[match(g, db$gene), ]

  ok <- is.finite(da$resid) & is.finite(db$resid)
  ra <- da$resid[ok]; rb <- db$resid[ok]; gg <- g[ok]
  flip <- sign(ra) != sign(rb) & ra != 0 & rb != 0

  agreement <- data.frame(
    n_gene         = length(gg),
    cor_resid      = if (length(gg) > 2) stats::cor(ra, rb) else NA_real_,
    sign_agree     = mean(!flip),
    n_disagree     = sum(flip),
    cor_pred       = stats::cor(da$pred[ok], db$pred[ok]),
    stringsAsFactors = FALSE
  )
  discordant <- data.frame(gene = gg[flip], resid_a = ra[flip],
                           resid_b = rb[flip], stringsAsFactors = FALSE)
  discordant <- discordant[order(-abs(discordant$resid_a)), , drop = FALSE]

  list(a = da, b = db, agreement = agreement, discordant = discordant)
}
