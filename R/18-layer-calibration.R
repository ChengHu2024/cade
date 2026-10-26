# ============================================================================
# cade: per-layer calibration
#
# WHY THIS EXISTS
#
# A calibration fitted over pooled samples from several anatomical layers
# describes no layer. In the companion study, the same 76 AAA cases and 13
# donors measured in both media and adventitia gave an all-gene slope of
# 0.269 in media and 0.037 in adventitia, under one and the same reference.
# Pooling the two would have produced a number that holds for neither.
#
# This function fits the calibration separately in each layer and reports the
# per-layer result, never a pooled one.
# ============================================================================

#' Calibration of the composition term within each anatomical layer
#'
#' @param decomp A `cade_decomp` object from [decompose_composition()].
#' @param layer A vector assigning each *gene* to a layer, or a data frame
#'   with columns `gene` and `layer`. Because calibration is gene-wise, the
#'   layer must be a property of the gene, not of the sample. If the layers
#'   differ by sample, decompose each layer separately and pass the results
#'   here as a named list instead.
#' @param expr Optional named vector of mean expression, used for the
#'   stratified calibration. Defaults to `decomp$mean_expr`.
#' @param strata Expression strata to report. `"all"` is always included.
#' @param min_genes Minimum number of genes required in a layer x stratum
#'   cell; cells below this are reported as `NA` with a note rather than
#'   silently fitted.
#'
#' @return A list of class `cade_layercal` with a `table` data frame
#'   (layer, stratum, n, slope, slope_lo, slope_hi, r2, intercept, verdict)
#'   and a `note` string.
#'
#' @details
#' The slope is from the ordinary least-squares regression of `obs` on
#' `pred`. A layer whose slope cannot be distinguished from zero is marked
#' `"not_calibrated"`; a layer with fewer than `min_genes` genes in a stratum
#' is marked `"insufficient"`. **No layer is ever dropped silently**, because
#' the interesting result is usually the layer that fails.
#'
#' @examples
#' \dontrun{
#' d <- decompose_composition(obs, frac, mu, is_case)
#' lay <- data.frame(gene = d$gene,
#'                   layer = rep(c("media", "adventitia"), length.out = nrow(d)))
#' calibrate_by_layer(d, lay)
#' }
#' @export
calibrate_by_layer <- function(decomp, layer, expr = NULL,
                               strata = c("all", "top50", "top25"),
                               min_genes = 100L) {

  if (!inherits(decomp, "cade_decomp"))
    stop("`decomp` must be a `cade_decomp` object from decompose_composition().",
         call. = FALSE)
  if (!all(c("gene", "obs", "pred") %in% names(decomp)))
    stop("`decomp` is missing gene/obs/pred columns.", call. = FALSE)

  ## ---- accept either a data frame or a named vector ---------------------
  if (is.data.frame(layer)) {
    if (!all(c("gene", "layer") %in% names(layer)))
      stop("`layer` data frame must have columns `gene` and `layer`.",
           call. = FALSE)
    lay <- stats::setNames(as.character(layer$layer), layer$gene)
  } else {
    if (is.null(names(layer)))
      stop("`layer` must be a named vector (names = gene) or a data frame ",
           "with columns gene and layer. An unnamed vector would be matched ",
           "by position and silently mislabel every gene.", call. = FALSE)
    ## as.character() on a named vector DROPS the names, which would make
    ## every gene look unassigned. Preserve them explicitly.
    lay <- stats::setNames(as.character(layer), names(layer))
  }

  miss <- setdiff(decomp$gene, names(lay))
  if (length(miss))
    stop(sprintf(paste0("%d gene(s) have no layer assignment (e.g. %s). ",
                        "Every gene must be assigned; an unassigned gene ",
                        "would be dropped silently."),
                 length(miss), paste(utils::head(miss, 5), collapse = ", ")),
         call. = FALSE)

  if (is.null(expr)) {
    ## decomp$mean_expr is an unnamed column of the data frame, in the order
    ## of decomp$gene. Indexing it by gene name would return all NA and
    ## silently empty every expression stratum, so name it here.
    expr <- stats::setNames(as.numeric(decomp$mean_expr), decomp$gene)
  } else {
    if (is.null(names(expr)))
      stop("`expr` must be a NAMED vector (names = gene); an unnamed vector ",
           "would be matched by position and silently mislabel every gene.",
           call. = FALSE)
    expr <- stats::setNames(as.numeric(expr[decomp$gene]), decomp$gene)
  }
  if (all(!is.finite(expr)))
    stop("the expression vector has no finite values; expression strata ",
         "cannot be formed.", call. = FALSE)

  g <- decomp$gene
  lay_v <- lay[g]
  layers <- sort(unique(lay_v))
  if (length(layers) < 2L)
    warning("only one layer present; the comparison this function exists for ",
            "cannot be made.", call. = FALSE)

  fit_one <- function(idx) {
    n <- length(idx)
    if (n < min_genes)
      return(list(n = n, slope = NA_real_, slope_lo = NA_real_,
                  slope_hi = NA_real_, r2 = NA_real_, intercept = NA_real_,
                  verdict = "insufficient"))
    x <- decomp$pred[idx]; y <- decomp$obs[idx]
    ok <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)
    if (n_ok < min_genes)
      return(list(n = n_ok, slope = NA_real_, slope_lo = NA_real_,
                  slope_hi = NA_real_, r2 = NA_real_, intercept = NA_real_,
                  verdict = "insufficient"))
    x <- x[ok]; y <- y[ok]
    if (stats::var(x) == 0)
      return(list(n = n_ok, slope = NA_real_, slope_lo = NA_real_,
                  slope_hi = NA_real_, r2 = NA_real_, intercept = NA_real_,
                  verdict = "degenerate"))
    m <- stats::lm(y ~ x)
    cf <- summary(m)$coefficients
    slope <- cf[2, 1]
    se    <- cf[2, 2]
    dfree <- m$df.residual
    tcrit <- stats::qt(0.975, dfree)
    r2    <- summary(m)$r.squared
    ## A slope whose 95% interval covers zero is not calibration.
    verdict <- if (is.finite(slope) && abs(slope) > 0 &&
                   (slope - tcrit * se > 0 | slope + tcrit * se < 0))
                 "calibrated" else "not_calibrated"
    list(n = n_ok, slope = slope, slope_lo = slope - tcrit * se,
         slope_hi = slope + tcrit * se, r2 = r2,
         intercept = cf[1, 1], verdict = verdict)
  }

  rows <- list()
  for (L in layers) {
    in_layer <- g[lay_v == L]
    for (s in strata) {
      if (s == "all") {
        idx <- which(g %in% in_layer)
      } else {
        th <- as.numeric(sub("top", "", s)) / 100
        ## expression threshold computed WITHIN the layer, so that a
        ## low-expressing layer is not entirely discarded by a global cut
        e_layer <- expr[in_layer]
        cut <- stats::quantile(e_layer, 1 - th, na.rm = TRUE, names = FALSE)
        idx <- which(g %in% in_layer & is.finite(expr) & expr >= cut)
      }
      r <- fit_one(idx)
      rows[[length(rows) + 1L]] <- data.frame(
        layer = L, stratum = s, n = r$n,
        slope = r$slope, slope_lo = r$slope_lo, slope_hi = r$slope_hi,
        r2 = r$r2, intercept = r$intercept, verdict = r$verdict,
        stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)

  ## A one-line comparison that makes the point of the function visible.
  all_rows <- tab[tab$stratum == "all", ]
  note <- if (nrow(all_rows) >= 2) {
    s <- all_rows$slope
    names(s) <- all_rows$layer
    paste0("All-gene slope by layer: ",
           paste(sprintf("%s = %.3f", names(s), s), collapse = "; "),
           if (any(is.finite(s)) && length(unique(round(s[is.finite(s)], 2))) > 1)
             ". The layers do not calibrate alike; report them separately."
           else ".")
  } else {
    "Only one layer present."
  }

  out <- list(table = tab, note = note, layers = layers,
              strata = strata, min_genes = min_genes,
              convention = attr(decomp, "convention"))
  class(out) <- "cade_layercal"
  out
}

#' @export
print.cade_layercal <- function(x, ...) {
  cat("<cade_layercal>\n")
  cat(sprintf("  convention = %s   layers = %s\n",
              if (is.null(x$convention)) "not recorded" else x$convention,
              paste(x$layers, collapse = ", ")))
  cat("\n")
  t <- x$table
  t$slope <- ifelse(is.na(t$slope), "   NA", sprintf("%6.3f", t$slope))
  t$r2    <- ifelse(is.na(t$r2),    "   NA", sprintf("%5.3f", t$r2))
  t$n     <- sprintf("%6d", t$n)
  print(t[, c("layer", "stratum", "n", "slope", "r2", "verdict")],
        row.names = FALSE)
  cat("\n  ", x$note, "\n", sep = "")
  invisible(x)
}
