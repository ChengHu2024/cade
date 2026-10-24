# ============================================================================
# cade: plots
#
# Deliberately minimal: ggplot2 is in Suggests, not Imports, so the package
# installs without it. Each plot function checks for ggplot2 and fails with a
# clear message rather than a missing-namespace error.
# ============================================================================

.need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("package 'ggplot2' is required for plotting; install it with install.packages('ggplot2').")
  invisible(TRUE)
}

#' Calibration scatter plot
#'
#' @param decomp Output of [decompose_composition()].
#' @param levels Expression quantiles to facet on, default the upper quarter.
#' @param point_alpha Transparency.
#'
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' d <- data.frame(gene = paste0("g", 1:300), obs = rnorm(300),
#'                 pred = rnorm(300), mean_expr = runif(300, .1, 100))
#' plot_calibration(d)
#' }
#' @export
plot_calibration <- function(decomp, levels = c(0, 0.75), point_alpha = 0.25) {
  .need_ggplot2()
  d <- decomp[is.finite(decomp$obs) & is.finite(decomp$pred) &
              is.finite(decomp$mean_expr), , drop = FALSE]
  d$facet <- NA_character_
  for (lv in levels) {
    keep <- expression_filter(d$mean_expr, lv)
    lab <- if (lv == 0) "all genes" else sprintf("top %d%% expressed",
                                                 round((1 - lv) * 100))
    d$facet[keep] <- lab
  }
  d <- d[!is.na(d$facet), , drop = FALSE]
  d$facet <- factor(d$facet, levels = unique(d$facet))

  labs <- do.call(rbind, lapply(split(d, d$facet), function(x) {
    f <- stats::lm(obs ~ pred, data = x)
    data.frame(facet = x$facet[1],
               slope = unname(stats::coef(f)[2]),
               r2 = summary(f)$r.squared,
               n = nrow(x))
  }))
  labs$label <- sprintf("slope = %.2f\nR2 = %.2f\nn = %d",
                        labs$slope, labs$r2, labs$n)

  ggplot2::ggplot(d, ggplot2::aes(x = pred, y = obs)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    ggplot2::geom_point(alpha = point_alpha, size = 0.7, colour = "#2c7fb8") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = "#d95f02", linewidth = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_text(data = labs,
                       ggplot2::aes(x = -Inf, y = Inf, label = label),
                       hjust = -0.08, vjust = 1.15, size = 3, inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ facet, scales = "free") +
    ggplot2::labs(x = "composition-implied log2 fold change",
                  y = "observed log2 fold change",
                  title = "Calibration of the composition term") +
    ggplot2::theme_bw(base_size = 10)
}


#' Observed / composition / residual decomposition plot
#'
#' @param decomp Output of [decompose_composition()].
#' @param genes Genes to show, in the order given.
#' @param facet Optional grouping variable, e.g. cohort. Recycled to `length(genes)`
#'   when a single value is supplied.
#' @param facet_label Axis label for the facet.
#'
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' d <- data.frame(gene = c("A","B"), obs = c(1,-1), pred = c(.6,-.4),
#'                 resid = c(.4,-.6), mean_expr = 1)
#' plot_decomposition(d, genes = c("A","B"))
#' }
#' @export
plot_decomposition <- function(decomp, genes = NULL, facet = NULL,
                               facet_label = NULL) {
  .need_ggplot2()
  d <- decomp
  if (!is.null(genes)) d <- d[d$gene %in% genes, , drop = FALSE]
  if (!nrow(d)) stop("no genes to plot.")
  d$gene <- factor(d$gene, levels = unique(c(genes, d$gene)))

  long <- rbind(
    data.frame(gene = d$gene, comp = "composition", value = d$pred),
    data.frame(gene = d$gene, comp = "cell-intrinsic", value = d$resid)
  )
  long$comp <- factor(long$comp, levels = c("composition", "cell-intrinsic"))
  pts <- data.frame(gene = d$gene, value = d$obs)
  if (!is.null(facet)) {
    long$facet <- rep(facet, length.out = nrow(long))
    pts$facet  <- rep(facet, length.out = nrow(pts))
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_col(data = long,
                      ggplot2::aes(x = gene, y = value, fill = comp),
                      position = "stack", width = 0.68) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = gene, y = value),
                        shape = 18, size = 2.4, colour = "black") +
    ggplot2::scale_fill_manual(values = c(composition = "#9ecae1",
                                          `cell-intrinsic` = "#fc9272")) +
    ggplot2::labs(x = NULL, y = "log2 fold change", fill = NULL,
                  title = "Observed change decomposed (diamond = observed)") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(facet))
    p <- p + ggplot2::facet_wrap(~ facet, scales = "free") +
      ggplot2::labs(x = facet_label %||% NULL)
  p
}

`%||%` <- function(a, b) if (is.null(a)) b else a
