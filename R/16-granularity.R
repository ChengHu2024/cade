# ============================================================================
# cade: reference granularity is not a well-defined complexity variable
#
# The number of cell-type compartments in a reference is not a property of the
# data. It is a property of how the analyst merges the fine-grained profiles,
# and the merge is only defined once a within-compartment weighting convention
# is chosen. This file wraps the fully crossed experiment that makes that
# dependence visible:
#
#   factor A  true composition granularity: pseudo-bulk generated at
#             A = 6 / 10 / 15 compartments
#   factor B  reference granularity: mu aggregated to B = 6 / 10 / 15
#   factor C  aggregation weight: atlas cells / pool cells / pool RNA / equal
#   repeats   5 sampling seeds
#
# For every cell of the 3 x 3 x 4 x 5 design the composition-implied log2 fold
# change is correlated with the observed log2 fold change, and the spread of
# those correlations across the 12 conditions that share a reference level is
# the quantity of interest. A reference level whose spread exceeds
# `sd_threshold` is flagged as not identifiable: the answer depends on the
# merge, not on the data.
#
# WHY THE TRUE LEVEL MATTERS. At a coarse true level the cells that make up a
# compartment are not resolved, so the fine composition inside the compartment
# is set by the atlas cell counts (the biological abundance). At the finest
# true level the supplied fine composition is used as it is. The reference at
# level B is always built from the same fine profiles, so a reference that is
# too coarse cannot represent the fine composition that actually generated the
# bulk. That mismatch is what the experiment measures.
#
# WHY THE NOISE TERM EXISTS. A pseudo-bulk built from the mean profiles alone
# would be an exact linear combination of the reference columns whenever the
# reference is at the finest level, so r would be exactly 1 and the design
# would be degenerate. Real pseudo-bulk carries per-gene sampling and
# biological noise. A lognormal term (default SD 0.1 per gene per pseudo-bulk)
# supplies that variation; without it the fine-reference r has zero variance.
# ============================================================================


`%||%` <- function(x, y) if (is.null(x)) y else x


# ---------------------------------------------------------------------------
# Built-in merge ladder
# ---------------------------------------------------------------------------
# Rules are applied in order until the requested number of compartments is
# reached. Every rule merges one compartment into another (or two into one)
# and is skipped when its members are not all present. The first five rules
# reproduce the compartment names used by the package (B_plasma, T_NK,
# Myeloid, Epithelial); the remaining rules coarsen the immune compartment and
# allow a six-compartment vascular scheme to be reached when the atlas lacks
# epithelium. The ladder never invents data: an unknown cell type is kept as
# its own compartment and reported.
#
# Rule 5 groups Fibroblast with VSMC as "Stromal". That is the coarsening that
# matters in vascular tissue, where the two mesenchymal compartments carry the
# largest and opposite composition changes; without it a ten-compartment atlas
# cannot be reduced to six compartments at all. It is recorded in the summary.
cade_merge_rules <- function() {
  list(
    ## NOTE: Endothelial is deliberately never merged. It is a canonical
    ## compartment name that passes through every level, so it is "known"
    ## without appearing as a merge target.
    list(members = c("B_cell", "Plasma"),               target = "B_plasma"),
    list(members = c("T_cell", "NK_cell"),              target = "T_NK"),
    list(members = c("Monocyte", "Macrophage"),         target = "Myeloid"),
    list(members = c("Alveolar_epi", "Airway_epi"),     target = "Epithelial"),
    list(members = c("Fibroblast", "VSMC"),             target = "Stromal"),
    list(members = c("Mast"),                           target = "Myeloid"),
    list(members = c("Dendritic"),                      target = "Myeloid"),
    list(members = c("Neutrophil"),                     target = "Myeloid"),
    list(members = c("B_plasma", "T_NK"),               target = "Lymphoid"),
    list(members = c("Lymphoid", "Myeloid"),            target = "Immune"),
    list(members = c("Pericyte"),                       target = "Stromal")
  )
}


# Cell-type names the built-in ladder understands. Endothelial, VSMC,
# epithelial and the immune labels pass through at every level; the merge
# targets are listed so that a caller can tell "unannotated" from "already
# coarsened".
cade_known_types <- function() {
  c("B_cell", "Plasma", "T_cell", "NK_cell", "Monocyte", "Macrophage",
    "Mast", "Dendritic", "Neutrophil", "Endothelial", "Fibroblast", "VSMC",
    "Pericyte", "Alveolar_epi", "Airway_epi", "B_plasma", "T_NK", "Myeloid",
    "Stromal", "Epithelial", "Lymphoid", "Immune")
}


# Resolve a level map for one target number of compartments.
#
# Returns a list with `map` (named character vector type -> group), `groups`
# (unique group labels in first-appearance order), `unknown` (types that are
# not named by any rule) and `fallback` (logical: was a generic merge needed?).
cade_level_map <- function(types, k, level_map = NULL) {
  n <- length(types)

  if (!is.null(level_map)) {
    if (is.function(level_map)) {
      m <- level_map(types, k)
    } else if (is.list(level_map) && !is.null(names(level_map))) {
      ## group label -> member types
      mm <- character(0)
      for (g in names(level_map))
        mm[level_map[[g]]] <- g
      m <- mm
    } else if (is.character(level_map)) {
      m <- level_map
    } else {
      stop("`level_map` must be NULL, a function, a named list, or a named ",
           "character vector.", call. = FALSE)
    }
    if (is.null(names(m)))
      stop("`level_map` must return a named vector (names = cell types).",
           call. = FALSE)
    miss <- setdiff(types, names(m))
    if (length(miss))
      stop(sprintf("`level_map` does not map cell type(s): %s",
                   paste(miss, collapse = ", ")), call. = FALSE)
    m <- m[types]
    if (anyNA(m) || any(!nzchar(m)))
      stop("`level_map` contains NA or empty group labels.", call. = FALSE)
    groups <- unique(m)
    if (length(groups) != k)
      warning(sprintf("`level_map` produced %d group(s), not the requested %d.",
                      length(groups), k), call. = FALSE)
    return(list(map = m, groups = groups, unknown = character(0),
                fallback = FALSE))
  }

  if (k == n)
    return(list(map = stats::setNames(types, types), groups = types,
                unknown = character(0), fallback = FALSE))
  if (k < 1L || k > n)
    stop(sprintf("cannot aggregate %d cell type(s) to %d compartment(s).",
                 n, k), call. = FALSE)

  known <- cade_known_types()
  unknown <- setdiff(types, known)

  current <- stats::setNames(types, types)   # type -> group
  fallback <- FALSE

  for (rule in cade_merge_rules()) {
    if (length(unique(current)) <= k) break
    present <- unique(current)
    if (!all(rule$members %in% present)) next
    if (length(unique(current[current %in% rule$members])) < 2L) next
    current[current %in% rule$members] <- rule$target
  }

  ## If the named ladder cannot reach k (for example a data set full of
  ## unannotated types), merge the last group into the previous one in a
  ## deterministic order. The summary records that this happened.
  while (length(unique(current)) > k) {
    gs <- unique(current)
    current[current == gs[length(gs)]] <- gs[length(gs) - 1L]
    fallback <- TRUE
  }

  groups <- unique(current)
  list(map = current, groups = groups, unknown = unknown, fallback = fallback)
}


# ---------------------------------------------------------------------------
# Weight construction
# ---------------------------------------------------------------------------
# Returns a named list of numeric weight vectors (one per method) and a
# character vector recording whether each was supplied or proxied.
cade_granularity_weights <- function(mu, frac, types) {
  G <- ncol(mu)
  out <- list()
  notes <- character(0)

  atlas <- attr(mu, "cell_counts")
  if (!is.null(atlas)) {
    atlas <- as.numeric(atlas[types])
    if (anyNA(atlas))
      stop("attr(mu, 'cell_counts') does not cover every cell type.",
           call. = FALSE)
    notes["atlas_cells"] <- "supplied (attr(mu, 'cell_counts'))"
  } else {
    atlas <- rowMeans(mu)                     # genes x types -> per-type mean
    notes["atlas_cells"] <- "proxy: rowMeans(mu)"
    if (stats::sd(atlas) == 0)
      notes["atlas_cells"] <- paste0(notes["atlas_cells"],
                                     " (degenerate: constant across types)")
  }
  names(atlas) <- types
  out$atlas_cells <- atlas

  pool <- attr(frac, "cell_counts")
  if (!is.null(pool)) {
    pool <- as.numeric(pool[types])
    if (anyNA(pool))
      stop("attr(frac, 'cell_counts') does not cover every cell type.",
           call. = FALSE)
    notes["pool_cells"] <- "supplied (attr(frac, 'cell_counts'))"
  } else {
    pool <- colMeans(frac)
    notes["pool_cells"] <- "proxy: colMeans(frac)"
  }
  names(pool) <- types
  out$pool_cells <- pool

  rna <- attr(frac, "rna")
  if (!is.null(rna)) {
    rna <- as.numeric(rna[types])
    if (anyNA(rna))
      stop("attr(frac, 'rna') does not cover every cell type.", call. = FALSE)
    notes["pool_rna"] <- "supplied (attr(frac, 'rna'))"
  } else {
    rna <- rowSums(mu)                        # total expression per type
    notes["pool_rna"] <- "proxy: rowSums(mu)"
    if (stats::sd(rna) == 0)
      notes["pool_rna"] <- paste0(notes["pool_rna"],
                                  " (degenerate: constant across types)")
  }
  names(rna) <- types
  out$pool_rna <- rna

  out$equal <- stats::setNames(rep(1, length(types)), types)
  notes["equal"] <- "constant 1"

  list(weights = out, notes = notes)
}


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------
# mu is K x G. Rows are averaged within each group with the given weights.
cade_aggregate_mu <- function(mu, map, w) {
  groups <- unique(map)
  out <- matrix(0, length(groups), ncol(mu),
                dimnames = list(groups, colnames(mu)))
  notes <- character(0)
  for (g in groups) {
    mem <- names(map)[map == g]
    ww <- w[mem]
    if (any(!is.finite(ww)) || any(ww < 0))
      stop("aggregation weights must be finite and non-negative.", call. = FALSE)
    if (sum(ww) <= 0) {
      ww <- rep(1, length(mem))
      notes <- c(notes, sprintf("group %s: all weights zero; equal weights used", g))
    }
    ww <- ww / sum(ww)
    out[g, ] <- as.numeric(crossprod(mu[mem, , drop = FALSE], ww))
  }
  attr(out, "notes") <- notes
  out
}


# Collapse a samples x fine-types fraction matrix to the merged groups by
# summing the member proportions (fractions add; they must not be averaged).
cade_collapse_frac <- function(frac, map) {
  groups <- unique(map)
  out <- vapply(groups, function(g)
    rowSums(frac[, names(map)[map == g], drop = FALSE]),
    numeric(nrow(frac)))
  colnames(out) <- groups
  out
}


# Expand a coarse composition back to the fine types. The total of each coarse
# group is kept; the fine types inside the group share it in proportion to
# `within` (atlas cell counts by default). This is how a bulk "generated at
# level A" is defined: the data only resolve A compartments, and the
# unobserved within-compartment split follows the atlas.
cade_expand_truth <- function(frac, map, within) {
  groups <- unique(map)
  F_coarse <- cade_collapse_frac(frac, map)
  out <- matrix(0, nrow(frac), length(map),
                dimnames = list(rownames(frac), names(map)))
  for (g in groups) {
    mem <- names(map)[map == g]
    ww <- within[mem]
    if (any(!is.finite(ww)) || sum(ww) <= 0) ww <- rep(1, length(mem))
    ww <- ww / sum(ww)
    out[, mem] <- F_coarse[, g] * rep(ww, each = nrow(frac))
  }
  out
}


# Sample one pseudo-bulk from a fine composition and add lognormal noise.
# `mu` is fine-types x genes.
cade_pseudo_bulk <- function(mu, f, n_cells, noise_sd) {
  p <- as.numeric(f)
  p[p < 0] <- 0
  if (sum(p) <= 0) stop("pseudo-bulk composition sums to zero.", call. = FALSE)
  p <- p / sum(p)
  p <- as.numeric(stats::rmultinom(1, n_cells, p)) / n_cells
  bulk <- as.numeric(crossprod(mu, p))
  if (noise_sd > 0)
    bulk <- bulk * exp(stats::rnorm(length(bulk), 0, noise_sd))
  bulk
}


# Observed log2 fold change between the case and control pseudo-bulks.
cade_obs_log2fc <- function(bulk_case, bulk_ctrl) {
  scale <- max(mean(bulk_case), mean(bulk_ctrl), 1)
  eps <- 1e-9 * scale
  log2((bulk_case + eps) / (bulk_ctrl + eps))
}


# One fit: generate the pseudo-bulks at true level A, aggregate the reference
# at level B with weight w, and correlate the composition-implied log2FC with
# the observed one. `f_fine` is the true fine composition; `f_B` is that
# composition collapsed to the reference level; `mu_B` is the aggregated
# reference.
cade_gran_fit <- function(mu, f_fine, f_B, is_case, mu_B, seed,
                          n_cells, noise_sd, expr_floor) {
  set.seed(seed)
  fc <- colMeans(f_fine[is_case, , drop = FALSE])
  ft <- colMeans(f_fine[!is_case, , drop = FALSE])
  bc <- cade_pseudo_bulk(mu, fc, n_cells, noise_sd)
  bt <- cade_pseudo_bulk(mu, ft, n_cells, noise_sd)

  FBc <- colMeans(f_B[is_case, , drop = FALSE])
  FBt <- colMeans(f_B[!is_case, , drop = FALSE])
  if (sum(FBc) <= 0 || sum(FBt) <= 0) return(NA_real_)
  FBc <- FBc / sum(FBc)
  FBt <- FBt / sum(FBt)

  num <- as.numeric(crossprod(mu_B, FBc))
  den <- as.numeric(crossprod(mu_B, FBt))
  eps <- 1e-12
  pred <- log2(pmax(num, eps) / pmax(den, eps))
  obs <- cade_obs_log2fc(bc, bt)

  mean_expr <- (bc + bt) / 2
  ## Expression filter (never effect-size). The floor is one count on a CP10k
  ## scale, expressed relative to the library size of the matrix so that the
  ## same rule works for CPM or raw counts. Low-expression genes have unstable
  ## log2 fold changes and would dilute the correlation.
  keep <- is.finite(pred) & is.finite(obs) & bc > 0 & bt > 0 &
    pmin(bc, bt) >= expr_floor
  if (sum(keep) < 30L)
    keep <- is.finite(pred) & is.finite(obs) & bc > 0 & bt > 0 &
      mean_expr >= stats::median(mean_expr)
  if (sum(keep) < 30L) return(NA_real_)
  if (stats::sd(pred[keep]) == 0 || stats::sd(obs[keep]) == 0)
    return(NA_real_)
  stats::cor(pred[keep], obs[keep])
}


#' Crossed granularity experiment for composition attribution
#'
#' Repeats the composition decomposition of pseudo-bulk data across a fully
#' crossed design of true composition granularity, reference granularity and
#' within-compartment aggregation weight, and reports how much the
#' reference-implied correlation moves when only the merge convention changes.
#'
#' @param mu Per-cell mean expression, **cell types in rows and genes in
#'   columns**, on a linear scale. Row names are the fine-grained cell types;
#'   column names are genes. Optional attributes:
#'   `"cell_counts"` (named vector of cells per type in the full atlas, used by
#'   the `atlas_cells` weight), `"n_cells"` (cells sampled per pseudo-bulk,
#'   default 2000) and `"noise_sd"` (lognormal per-gene noise per pseudo-bulk,
#'   default 0.1).
#' @param frac Samples x cell types matrix or data frame of cell-type
#'   proportions. Columns are matched to `mu` by name. Optional attributes
#'   `"cell_counts"` (cells per type in this pool, used by `pool_cells`) and
#'   `"rna"` (total RNA per type in this pool, used by `pool_rna`).
#' @param is_case Logical vector over the rows of `frac`; `TRUE` marks case
#'   samples.
#' @param levels Integer vector of compartment counts to test; each level is
#'   used both as the true composition granularity and as the reference
#'   granularity, so the design has `length(levels)^2` cells before weights and
#'   seeds.
#' @param weights Aggregation weights. Any subset of `"atlas_cells"`,
#'   `"pool_cells"`, `"pool_rna"` and `"equal"`. See **Weights**.
#' @param seeds Integer vector of sampling seeds. Each seed fixes the pseudo-bulk
#'   draw through [stats::set.seed()].
#' @param level_map Optional custom merge rule. Either a function
#'   `function(types, k)` returning a named character vector
#'   (names = fine types, values = group labels), a named list mapping group
#'   labels to member types, or a named character vector type -> group. When
#'   supplied it is used for every requested level and the built-in ladder is
#'   bypassed. When `NULL`, [cade_merge_rules()] is applied in order.
#' @param sd_threshold Numeric. A reference level whose population SD of `r`
#'   across true levels x weights x seeds is below this value is reported as
#'   identifiable.
#' @param n_cores Number of workers. Values above 1 use
#'   [parallel::mclapply()]; on Windows the function falls back to serial with
#'   a warning.
#'
#' @details
#' **The merge ladder is data-set specific.** The built-in ladder is written
#' for a ten-compartment vascular reference (aorta). It is *not* the ladder
#' used in the companion paper, which starts from a fifteen-compartment lung
#' atlas. In particular, the paper keeps Fibroblast and VSMC separate at every
#' level; a ten-compartment vascular atlas cannot be reduced to six
#' compartments without merging the two mesenchymal compartments, so the
#' ladder used here groups them as `Stromal`. Consequently, for a ten-type
#' input the returned `K = 6` result must **not** be compared with the `K = 6`
#' row of the paper's Table 3: it is a different partition. The exact partition
#' is returned in `merge_rules_used`, and the summary warns whenever the input
#' has a number of cell types other than fifteen.
#'
#' **What `r` is.** The function computes
#' `r = cor(composition log2FC, observed log2FC)`, where the observed log2FC
#' comes from a pseudo-bulk sampled from `frac` and generated with the
#' finest-grained `mu`, and the composition term comes from the weighted `mu`
#' at the reference granularity. The paper's experiment additionally runs NNLS
#' deconvolution on the observed side; this function does not. Absolute `r`
#' values therefore differ (about 0.92 at `K = 10` here versus about 0.75 in
#' the paper). What is validated is the direction and magnitude of the change
#' in SD with reference granularity, not the paper's absolute numbers.
#'
#' @return A list of class `"cade_granularity"` with elements:
#'   \describe{
#'     \item{`raw`}{Data frame, one row per fit, with columns `true_level`,
#'       `ref_level`, `weight`, `seed` and `r`.}
#'     \item{`matrix_r`}{`levels` x `levels` matrix of mean `r` over weights and
#'       seeds (rows = true level, columns = reference level).}
#'     \item{`sd_by_ref`}{Named numeric vector: population SD of `r` across all
#'       true levels x weights x seeds that share a reference level.}
#'     \item{`sd_sample`}{The same spread with the sample (n - 1) denominator.}
#'     \item{`identifiable`}{Named logical vector: `sd_by_ref < sd_threshold`.}
#'     \item{`level_groups`}{List of compartment labels used at each requested
#'       level (one character vector per level, in first-appearance order).}
#'     \item{`level_maps`}{List of named character vectors mapping each fine cell
#'       type to its compartment at each requested level; this is the exact
#'       merge that was used.}
#'     \item{`merge_rules_used`}{Named list, one element per requested level;
#'       each element maps a coarse compartment to the fine cell types it
#'       contains, so a caller can see what a given `K` actually means.}
#'     \item{`r_definition`}{Character scalar:
#'       `"cor(comp, obs) without deconvolution"`.}
#'     \item{`proxies`}{Named character vector recording, for each weight, whether
#'       it was supplied or proxied.}
#'     \item{`unknown_types`}{Cell types not named by the built-in ladder; they
#'       are kept as singleton compartments and listed in the summary.}
#'     \item{`settings`}{List of the settings that affect the result (`n_cells`,
#'       `noise_sd`, `sd_threshold`, `seeds`, `weights`).}
#'     \item{`summary`}{Character vector: human-readable report.}
#'   }
#'
#' @section Weights:
#' When the fine profiles are merged into a compartment, the compartment
#' profile is `sum_k (w_k / sum(w)) * mu[k, ]` over the fine types `k` that
#' belong to it. The four conventions are:
#' \describe{
#'   \item{`atlas_cells`}{cells per type in the full atlas
#'     (`attr(mu, "cell_counts")`; proxy: `rowMeans(mu)`).}
#'   \item{`pool_cells`}{cells per type in this data pool
#'     (`attr(frac, "cell_counts")`; proxy: `colMeans(frac)`).}
#'   \item{`pool_rna`}{total RNA per type in this pool (`attr(frac, "rna")`;
#'     proxy: `rowSums(mu)`).}
#'   \item{`equal`}{all weights equal to one.}
#' }
#' Every substitution is recorded in `proxies` and in the summary. On a matrix
#' that is already normalised per cell (for example CP10k, where every cell type
#' sums to the same total), `rowMeans(mu)` and `rowSums(mu)` are constant across
#' types; both proxies then collapse to `equal`. This is a property of the
#' data, not of the function, and the summary says so.
#'
#' @section What the spread means:
#' A reference level with a large spread is not a level at which the
#' decomposition is "less accurate"; it is a level at which the answer is a
#' function of the merge convention. The fine reference is expected to be
#' stable because every fine profile is available; a coarse reference is
#' expected to be unstable because the within-compartment mix can no longer be
#' recovered from the data. The `identifiable` flag encodes the pre-specified
#' rule, not a biological conclusion.
#'
#' @examples
#' set.seed(1)
#' types <- c("B_cell", "Plasma", "T_cell", "NK_cell", "Monocyte",
#'            "Macrophage", "Mast", "Endothelial", "Fibroblast", "VSMC")
#' mu <- matrix(stats::runif(length(types) * 200, 0.1, 5), nrow = length(types),
#'              dimnames = list(types, paste0("g", 1:200)))
#' frac <- matrix(stats::runif(12 * length(types)), nrow = 12,
#'                dimnames = list(paste0("s", 1:12), types))
#' frac <- frac / rowSums(frac)
#' is_case <- rep(c(FALSE, TRUE), each = 6)
#' g <- granularity_crossed(mu, frac, is_case, levels = c(5, 10),
#'                          seeds = 1:2, sd_threshold = 0.1)
#' g$sd_by_ref
#' @export
granularity_crossed <- function(mu, frac, is_case, levels = c(6, 10),
                                weights = c("atlas_cells", "pool_cells",
                                            "pool_rna", "equal"),
                                seeds = 1:5,
                                level_map = NULL,
                                sd_threshold = 0.1,
                                n_cores = 1L) {

  ## ---------------------------------------------------------------- input
  mu <- as.matrix(mu)
  frac <- as.matrix(frac)
  if (!is.numeric(mu) || !is.numeric(frac))
    stop("`mu` and `frac` must be numeric.", call. = FALSE)
  if (is.null(rownames(mu)) || any(!nzchar(rownames(mu))))
    stop("`mu` must have non-empty cell-type row names.", call. = FALSE)
  if (anyDuplicated(rownames(mu)))
    stop("`mu` has duplicated cell-type row names.", call. = FALSE)
  if (is.null(colnames(mu)) || any(!nzchar(colnames(mu))))
    stop("`mu` must have non-empty gene column names.", call. = FALSE)
  if (anyDuplicated(colnames(mu)))
    stop("`mu` has duplicated gene names.", call. = FALSE)
  if (is.null(colnames(frac)))
    stop("`frac` must have cell-type column names.", call. = FALSE)
  if (any(!is.finite(mu)) || any(!is.finite(frac)))
    stop("`mu` or `frac` contain non-finite values.", call. = FALSE)
  if (any(mu < 0) || any(frac < 0))
    stop("`mu` and `frac` must be non-negative.", call. = FALSE)

  types <- rownames(mu)
  if (!setequal(colnames(frac), types))
    stop(sprintf(paste0("cell types in `frac` and `mu` differ.\n",
                        "  only in mu: %s\n  only in frac: %s"),
                 paste(setdiff(types, colnames(frac)), collapse = ", "),
                 paste(setdiff(colnames(frac), types), collapse = ", ")),
         call. = FALSE)
  frac <- frac[, types, drop = FALSE]

  if (!is.logical(is_case))
    stop("`is_case` must be logical.", call. = FALSE)
  if (length(is_case) != nrow(frac))
    stop(sprintf("`is_case` has length %d but `frac` has %d rows.",
                 length(is_case), nrow(frac)), call. = FALSE)
  if (anyNA(is_case) || !any(is_case) || all(is_case))
    stop("`is_case` must contain both TRUE and FALSE and no NA.", call. = FALSE)

  rs <- rowSums(frac)
  if (any(rs <= 0))
    stop("some samples have all-zero fractions.", call. = FALSE)
  frac <- frac / rs

  weights <- as.character(weights)
  allowed <- c("atlas_cells", "pool_cells", "pool_rna", "equal")
  bad <- setdiff(weights, allowed)
  if (length(bad))
    stop(sprintf("unknown weight(s): %s. Allowed: %s",
                 paste(bad, collapse = ", "), paste(allowed, collapse = ", ")),
         call. = FALSE)
  weights <- unique(weights)

  levels <- sort(unique(as.integer(levels)))
  if (!length(levels) || anyNA(levels) || any(levels < 1L))
    stop("`levels` must be positive integers.", call. = FALSE)
  if (length(levels) == 1L && levels > length(types))
    stop("`levels` exceeds the number of cell types.", call. = FALSE)
  if (any(levels > length(types)))
    stop(sprintf("level(s) %s exceed the %d cell type(s) in `mu`.",
                 paste(levels[levels > length(types)], collapse = ", "),
                 length(types)), call. = FALSE)

  seeds <- as.integer(seeds)
  if (!length(seeds) || anyNA(seeds))
    stop("`seeds` must be a non-empty integer vector.", call. = FALSE)

  sd_threshold <- as.numeric(sd_threshold)
  if (length(sd_threshold) != 1L || !is.finite(sd_threshold) ||
      sd_threshold <= 0)
    stop("`sd_threshold` must be a positive number.", call. = FALSE)

  n_cores <- as.integer(n_cores %||% 1L)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L)
    stop("`n_cores` must be a positive integer.", call. = FALSE)

  n_cells <- attr(mu, "n_cells") %||% 2000L
  n_cells <- as.integer(n_cells)
  if (is.na(n_cells) || n_cells < 10L)
    stop("attr(mu, 'n_cells') must be an integer >= 10.", call. = FALSE)
  noise_sd <- attr(mu, "noise_sd") %||% 0.1
  noise_sd <- as.numeric(noise_sd)
  if (length(noise_sd) != 1L || !is.finite(noise_sd) || noise_sd < 0)
    stop("attr(mu, 'noise_sd') must be a non-negative number.", call. = FALSE)

  ## One count on a CP10k scale, expressed relative to the library size of the
  ## matrix, so the same expression floor works for CPM and raw counts.
  expr_floor <- 1e-4 * stats::median(rowSums(mu))
  if (!is.finite(expr_floor) || expr_floor <= 0) expr_floor <- 0

  ## -------------------------------------------------------- merge + weights
  maps <- lapply(levels, function(k) cade_level_map(types, k, level_map))
  names(maps) <- as.character(levels)
  level_groups <- lapply(maps, `[[`, "groups")
  level_maps <- lapply(maps, `[[`, "map")
  unknown_types <- sort(unique(unlist(lapply(maps, `[[`, "unknown"))))

  ## Unknown names must never be silently merged or dropped: the caller has to
  ## say what they are. Aggregating to the finest level is still allowed,
  ## because it does not merge anything.
  if (length(unknown_types) && any(levels < length(types)))
    stop(sprintf(paste0(
      "cell type(s) not covered by the built-in merge ladder: %s.\n",
      "  Aggregating them would silently invent a compartment.\n",
      "  Either supply `level_map` (function or named list) or run only at ",
      "levels = %d."),
      paste(unknown_types, collapse = ", "), length(types)), call. = FALSE)

  fallback_used <- any(vapply(maps, `[[`, logical(1), "fallback"))

  wobj <- cade_granularity_weights(mu, frac, types)
  wvecs <- wobj$weights[weights]
  proxies <- wobj$notes[weights]

  ## Aggregated references: one matrix per (reference level, weight).
  mu_B <- list()
  agg_notes <- character(0)
  for (b in levels) {
    for (wn in weights) {
      key <- paste(b, wn, sep = "|")
      mu_B[[key]] <- cade_aggregate_mu(mu, maps[[as.character(b)]]$map,
                                       wvecs[[wn]])
      if (length(attr(mu_B[[key]], "notes")))
        agg_notes <- c(agg_notes, paste0("K=", b, ", ", wn, ": ",
                                         attr(mu_B[[key]], "notes")))
    }
  }

  ## True fine composition per true level. At the finest level the supplied
  ## composition is used unchanged; below it, the within-compartment split
  ## follows the atlas counts.
  ## The within-compartment composition of a coarse truth is biological, not
  ## an analyst choice: use the atlas counts whenever they exist, otherwise
  ## equal shares. This is independent of the weights requested for the fit.
  atlas_w <- wobj$weights$atlas_cells
  if (is.null(atlas_w)) atlas_w <- stats::setNames(rep(1, length(types)), types)
  truth <- list()
  for (a in levels) {
    if (a == length(types)) {
      truth[[as.character(a)]] <- frac
    } else {
      truth[[as.character(a)]] <- cade_expand_truth(frac,
                                                    maps[[as.character(a)]]$map,
                                                    atlas_w)
    }
  }
  ## Collapsed truth at each reference level.
  truth_B <- list()
  for (t in names(truth)) {
    for (b in levels) {
      truth_B[[paste(t, b, sep = "|")]] <-
        cade_collapse_frac(truth[[t]], maps[[as.character(b)]]$map)
    }
  }

  ## ------------------------------------------------------------------ fits
  grid <- expand.grid(true_level = levels, ref_level = levels,
                      weight = weights, seed = seeds,
                      stringsAsFactors = FALSE)
  grid <- grid[, c("true_level", "ref_level", "weight", "seed")]

  one <- function(i) {
    tl <- grid$true_level[i]; rl <- grid$ref_level[i]
    wn <- grid$weight[i]; sd <- grid$seed[i]
    tryCatch(
      cade_gran_fit(mu = mu,
                    f_fine = truth[[as.character(tl)]],
                    f_B = truth_B[[paste(tl, rl, sep = "|")]],
                    is_case = is_case,
                    mu_B = mu_B[[paste(rl, wn, sep = "|")]],
                    seed = sd, n_cells = n_cells, noise_sd = noise_sd,
                    expr_floor = expr_floor),
      error = function(e) NA_real_)
  }

  if (n_cores > 1L) {
    if (.Platform$OS.type == "windows") {
      warning("parallel::mclapply() is not available on Windows; ",
              "running serially.", call. = FALSE)
      r <- vapply(seq_len(nrow(grid)), one, numeric(1))
    } else {
      r <- unlist(parallel::mclapply(seq_len(nrow(grid)), one,
                                     mc.cores = n_cores))
    }
  } else {
    r <- vapply(seq_len(nrow(grid)), one, numeric(1))
  }

  raw <- data.frame(true_level = grid$true_level,
                    ref_level = grid$ref_level,
                    weight = grid$weight,
                    seed = grid$seed,
                    r = as.numeric(r),
                    stringsAsFactors = FALSE)
  n_na <- sum(!is.finite(raw$r))

  ## -------------------------------------------------------------- summaries
  matrix_r <- matrix(NA_real_, length(levels), length(levels),
                     dimnames = list(as.character(levels),
                                     as.character(levels)))
  for (a in levels) for (b in levels) {
    v <- raw$r[raw$true_level == a & raw$ref_level == b]
    matrix_r[as.character(a), as.character(b)] <-
      if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  }

  sd_by_ref <- stats::setNames(rep(NA_real_, length(levels)),
                               as.character(levels))
  sd_sample <- sd_by_ref
  n_by_ref <- stats::setNames(rep(0L, length(levels)),
                              as.character(levels))
  for (b in levels) {
    v <- raw$r[raw$ref_level == b]
    v <- v[is.finite(v)]
    n_by_ref[as.character(b)] <- length(v)
    if (length(v) >= 2L) {
      sd_by_ref[as.character(b)] <- sqrt(mean((v - mean(v))^2))
      sd_sample[as.character(b)] <- stats::sd(v)
    }
  }

  identifiable <- stats::setNames(sd_by_ref < sd_threshold,
                                  names(sd_by_ref))
  identifiable[is.na(identifiable)] <- FALSE

  ## --------------------------------------------------------------- summary
  lines <- c(
    "CADE crossed granularity experiment",
    sprintf("  cell types (finest): %d", length(types)),
    sprintf("  levels: %s", paste(levels, collapse = ", ")),
    sprintf("  achieved compartments: %s",
            paste(vapply(level_groups, length, integer(1)), collapse = ", ")),
    sprintf("  weights: %s", paste(weights, collapse = ", ")),
    sprintf("  seeds: %s", paste(seeds, collapse = ", ")),
    sprintf("  pseudo-bulk: %d cells, lognormal noise SD %.3g, expression floor %.3g",
            n_cells, noise_sd, expr_floor),
    sprintf("  conditions per reference level: %d",
            length(levels) * length(weights) * length(seeds)),
    if (n_na > 0L)
      sprintf("  fits without a usable r (too few expressed genes): %d", n_na)
    else NULL,
    "  weight sources:",
    paste0("    ", names(proxies), ": ", proxies),
    if (length(unknown_types))
      sprintf("  unknown cell type(s), kept as singletons: %s",
              paste(unknown_types, collapse = ", "))
    else "  unknown cell type(s): none",
    if (fallback_used)
      "  NOTE: the built-in ladder could not reach every requested level; a deterministic fallback merge was used (see level_groups)."
    else NULL,
    if (length(agg_notes)) paste0("  ", agg_notes) else NULL,
    "  reference-level spread of r:",
    sprintf("    K=%2d  SD=%s  n=%d  %s",
            levels,
            ifelse(is.na(sd_by_ref), "NA", sprintf("%.4f", sd_by_ref)),
            n_by_ref,
            ifelse(identifiable, "identifiable", "NOT identifiable"))
  )
  if (any(!identifiable, na.rm = TRUE)) {
    lines <- c(lines,
      sprintf("  verdict: granularity is not identifiable at K = %s (SD > %.2g).",
              paste(levels[!identifiable], collapse = ", "), sd_threshold))
  } else {
    lines <- c(lines,
      sprintf("  verdict: no reference level exceeded SD = %.2g.", sd_threshold))
  }
  if (length(types) != 15L) {
    lines <- c(lines, sprintf(paste0(
      "  WARNING: merge ladder adapted for %d input cell type(s); ",
      "not the 15-type ladder of the companion paper."), length(types)))
  }

  merge_rules_used <- lapply(level_maps, function(m)
    split(names(m), factor(m, levels = unique(m))))

  out <- list(raw = raw, matrix_r = matrix_r, sd_by_ref = sd_by_ref,
              sd_sample = sd_sample, identifiable = identifiable,
              level_groups = level_groups, level_maps = level_maps,
              merge_rules_used = merge_rules_used,
              r_definition = "cor(comp, obs) without deconvolution",
              proxies = proxies,
              unknown_types = unknown_types, n_by_ref = n_by_ref,
              settings = list(n_cells = n_cells, noise_sd = noise_sd,
                              expr_floor = expr_floor,
                              sd_threshold = sd_threshold,
                              seeds = seeds, weights = weights),
              summary = lines)
  class(out) <- "cade_granularity"
  out
}


#' Print a crossed granularity experiment
#'
#' @param x Output of [granularity_crossed()].
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.cade_granularity <- function(x, ...) {
  cat(x$summary, sep = "\n")
  cat("\n")
  cat("  mean r (rows = true level, columns = reference level):\n")
  print(round(x$matrix_r, 3))
  invisible(x)
}
