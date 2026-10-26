# ============================================================================
# Zero-model test for the composition term
#
# The composition term has no intrinsic scale: a correlation of 0.28 between
# `comp` and `obs` is only informative if it exceeds what the same arithmetic
# produces when the cell-type-to-expression-profile correspondence carries no
# information. null_permutation() permutes the cell types of the reference
# spectrum, recomputes `comp` under the requested convention and reports where
# the observed correlation sits inside that null.
#
# Two implementation choices keep the loop fast enough to be usable at
# n_perm = 20000 with pure R:
#
#   * the permutation is applied to the small K x n sample matrix, never to the
#     large mu. For a permutation `pi` of the cell-type axis,
#         mu[, pi] %*% f  ==  mu %*% f[pi^-1]
#     and `pi` and its inverse are both uniform on the symmetric group, so the
#     two formulations have exactly the same distribution while only the
#     K x n matrix is copied.
#   * the correlation is computed from a pre-centred `obs` with an explicit dot
#     product instead of stats::cor(), which re-centres and re-scans on every
#     call.
# ============================================================================

## Floor applied before log2(); log2(0) would be -Inf.
.NULLPERM_EPS <- 1e-12

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

## Validate and normalise the inputs. Stops with a message that names the input
## and the mismatch; returns the four aligned pieces plus bookkeeping.
.nullperm_prepare <- function(obs, frac, mu, is_case, genes) {

  ## ---- obs ---------------------------------------------------------------
  if (!is.numeric(obs)) {
    stop("`obs` must be a numeric vector of log2 fold changes (names = genes), not ",
         paste(class(obs), collapse = "/"), ".", call. = FALSE)
  }
  if (is.null(names(obs))) {
    stop("`obs` must be named: the names are the gene identifiers that are ",
         "aligned with `colnames(mu)`.", call. = FALSE)
  }
  if (anyDuplicated(names(obs))) {
    dup <- unique(names(obs)[duplicated(names(obs))])
    stop("`obs` has duplicated gene names: ",
         paste(utils::head(dup, 5L), collapse = ", "), ".", call. = FALSE)
  }
  if (!all(is.finite(obs))) {
    stop(sprintf("`obs` contains %d non-finite value(s); a correlation needs finite values.",
                 sum(!is.finite(obs))), call. = FALSE)
  }
  if (length(obs) < 3L) {
    stop(sprintf("`obs` has %d gene(s); at least 3 are needed for a correlation.",
                 length(obs)), call. = FALSE)
  }

  ## ---- frac --------------------------------------------------------------
  if (is.data.frame(frac)) {
    bad <- names(frac)[!vapply(frac, is.numeric, logical(1))]
    if (length(bad)) {
      stop("`frac` must contain only the numeric cell-type proportion columns; ",
           "non-numeric column(s) found: ", paste(bad, collapse = ", "),
           ". Select the cell-type columns before calling (e.g. `frac[, celltypes]`).",
           call. = FALSE)
    }
  }
  frac <- as.matrix(frac)
  if (!is.numeric(frac)) {
    stop("`frac` must be numeric (cell-type proportions, samples in rows).",
         call. = FALSE)
  }
  if (nrow(frac) < 2L) {
    stop(sprintf("`frac` has %d row(s); at least 2 samples are required.", nrow(frac)),
         call. = FALSE)
  }
  if (is.null(colnames(frac))) {
    stop("`frac` must have cell-type column names matching `rownames(mu)`.",
         call. = FALSE)
  }
  if (anyDuplicated(colnames(frac))) {
    dup <- unique(colnames(frac)[duplicated(colnames(frac))])
    stop("`frac` has duplicated cell-type column names: ",
         paste(dup, collapse = ", "), ".", call. = FALSE)
  }
  if (!all(is.finite(frac))) {
    stop(sprintf("`frac` contains %d non-finite value(s).", sum(!is.finite(frac))),
         call. = FALSE)
  }
  if (any(frac < 0)) {
    stop("`frac` contains negative values; proportions must be non-negative.",
         call. = FALSE)
  }
  rs <- rowSums(frac)
  if (any(rs <= 0)) {
    stop(sprintf("%d row(s) of `frac` sum to zero; every sample needs a non-zero total.",
                 sum(rs <= 0)), call. = FALSE)
  }
  if (any(abs(rs - 1) > 1e-6)) {
    warning(sprintf("`frac` rows do not all sum to 1 (range %.6f to %.6f); proportions are used as supplied, not renormalised.",
                    min(rs), max(rs)), call. = FALSE)
  }

  ## ---- is_case -----------------------------------------------------------
  if (!is.logical(is_case)) {
    stop("`is_case` must be a logical vector: TRUE for case (AAA) samples, ",
         "FALSE for control (normal/donor) samples.", call. = FALSE)
  }
  if (length(is_case) != nrow(frac)) {
    stop(sprintf("`is_case` has length %d but `frac` has %d rows.",
                 length(is_case), nrow(frac)), call. = FALSE)
  }
  if (anyNA(is_case)) {
    stop("`is_case` contains NA; every row of `frac` must be classified.",
         call. = FALSE)
  }
  if (!any(is_case)) {
    stop("`is_case` has no TRUE values; at least one case sample is required.",
         call. = FALSE)
  }
  if (all(is_case)) {
    stop("`is_case` is all TRUE; at least one control sample is required.",
         call. = FALSE)
  }

  ## ---- mu ----------------------------------------------------------------
  if (is.data.frame(mu)) mu <- as.matrix(mu)
  if (!is.matrix(mu)) {
    stop("`mu` must be a matrix with cell types in rows and genes in columns, not ",
         paste(class(mu), collapse = "/"), ".", call. = FALSE)
  }
  if (!is.numeric(mu)) {
    stop("`mu` must be numeric (per-cell mean expression, linear scale).",
         call. = FALSE)
  }
  if (!all(is.finite(mu))) {
    stop(sprintf("`mu` contains %d non-finite value(s).", sum(!is.finite(mu))),
         call. = FALSE)
  }
  if (any(mu < 0)) {
    stop("`mu` contains negative values; per-cell mean expression must be non-negative.",
         call. = FALSE)
  }
  if (is.null(rownames(mu)) || is.null(colnames(mu))) {
    stop("`mu` must carry both dimnames: rownames = cell types, colnames = genes.",
         call. = FALSE)
  }
  if (anyDuplicated(rownames(mu))) {
    dup <- unique(rownames(mu)[duplicated(rownames(mu))])
    stop("`mu` has duplicated cell-type rownames: ",
         paste(dup, collapse = ", "), ".", call. = FALSE)
  }
  if (anyDuplicated(colnames(mu))) {
    dup <- unique(colnames(mu)[duplicated(colnames(mu))])
    stop("`mu` has duplicated gene colnames: ",
         paste(utils::head(dup, 5L), collapse = ", "), ".", call. = FALSE)
  }

  ## ---- orientation -------------------------------------------------------
  ## The package convention is cell types x genes, so colnames(frac) must be the
  ## same set as rownames(mu) and names(obs) must be the same set as
  ## colnames(mu). A genes x cell types matrix (the layout of the CSV written by
  ## the analysis pipeline) is detected and transposed with a warning.
  ct      <- colnames(frac)
  ct_rows <- setequal(ct, rownames(mu))
  ct_cols <- setequal(ct, colnames(mu))
  gene_rows <- length(intersect(rownames(mu), names(obs)))
  gene_cols <- length(intersect(colnames(mu), names(obs)))
  transposed <- FALSE

  if (ct_rows && !ct_cols) {
    ## already the package convention
  } else if (!ct_rows && ct_cols) {
    transposed <- TRUE
  } else if (ct_rows && ct_cols) {
    ## cell-type names match both margins (symmetric test); use gene overlap
    transposed <- gene_rows > gene_cols
  } else {
    only_frac <- setdiff(ct, rownames(mu))
    only_mu   <- setdiff(rownames(mu), ct)
    stop(sprintf(
      "`colnames(frac)` and `rownames(mu)` must be the same set of cell types. In `frac` but not in `rownames(mu)`: [%s]; in `rownames(mu)` but not in `frac`: [%s]. If `mu` is genes x cell types, supply its transpose.",
      paste(only_frac, collapse = ", "), paste(only_mu, collapse = ", ")),
      call. = FALSE)
  }

  if (transposed) {
    warning(sprintf(
      "`mu` is arranged as genes x cell types rather than the package convention (cell types x genes): `rownames(mu)` has %d entries of which %d are genes in `obs`, while `colnames(mu)` (%d entries) matches `colnames(frac)` (nrow(mu) = %d, ncol(mu) = %d). Transposing `mu`; pass `t(mu)` to silence this warning.",
      nrow(mu), gene_rows, length(ct), nrow(mu), ncol(mu)), call. = FALSE)
    mu <- t(mu)
  }

  if (!setequal(ct, rownames(mu))) {
    only_frac <- setdiff(ct, rownames(mu))
    only_mu   <- setdiff(rownames(mu), ct)
    stop(sprintf(
      "`colnames(frac)` and `rownames(mu)` must be the same set of cell types. In `frac` but not in `rownames(mu)`: [%s]; in `rownames(mu)` but not in `frac`: [%s].",
      paste(only_frac, collapse = ", "), paste(only_mu, collapse = ", ")),
      call. = FALSE)
  }
  ## align both margins to one order
  mu   <- mu[ct, , drop = FALSE]
  frac <- frac[, ct, drop = FALSE]

  ## ---- genes -------------------------------------------------------------
  genes_arg <- genes
  if (is.null(genes)) {
    missing <- setdiff(names(obs), colnames(mu))
    if (length(missing)) {
      warning(sprintf("%d gene(s) in `obs` are absent from `colnames(mu)` and are dropped, e.g. %s.",
                      length(missing),
                      paste(utils::head(missing, 5L), collapse = ", ")),
              call. = FALSE)
    }
    genes <- intersect(names(obs), colnames(mu))
  } else {
    if (!is.character(genes)) {
      stop("`genes` must be a character vector of gene names.", call. = FALSE)
    }
    if (anyDuplicated(genes)) {
      dup <- unique(genes[duplicated(genes)])
      stop("`genes` contains duplicates: ", paste(utils::head(dup, 5L), collapse = ", "),
           ".", call. = FALSE)
    }
    miss_obs <- setdiff(genes, names(obs))
    if (length(miss_obs)) {
      stop(sprintf("%d gene(s) in `genes` are absent from `names(obs)`: %s.",
                   length(miss_obs),
                   paste(utils::head(miss_obs, 5L), collapse = ", ")), call. = FALSE)
    }
    miss_mu <- setdiff(genes, colnames(mu))
    if (length(miss_mu)) {
      stop(sprintf("%d gene(s) in `genes` are absent from `colnames(mu)`: %s.",
                   length(miss_mu),
                   paste(utils::head(miss_mu, 5L), collapse = ", ")), call. = FALSE)
    }
  }
  if (length(genes) < 3L) {
    if (is.null(genes_arg)) {
      stop(sprintf("only %d gene(s) are shared by `obs` and `colnames(mu)`; at least 3 are needed for a correlation.",
                   length(genes)), call. = FALSE)
    } else {
      stop(sprintf("`genes` selects only %d gene(s); at least 3 are needed for a correlation.",
                   length(genes)), call. = FALSE)
    }
  }

  list(obs = obs[genes],
       mu = mu[, genes, drop = FALSE],
       frac_case = frac[is_case, , drop = FALSE],
       frac_ctrl = frac[!is_case, , drop = FALSE],
       genes = genes,
       n_case = sum(is_case),
       n_ctrl = sum(!is_case),
       transposed = transposed)
}

## Draw (or accept) the permutations as a K x n_perm integer matrix. A stored
## permutation `pi` always means "permute the cell-type axis of mu", i.e.
## comp = mu[, pi] %*% f. The loop applies the equivalent form on the fraction
## rows, which is mu %*% f[pi^-1], so supplied permutations are inverted here
## and every permutation source keeps one meaning.
.nullperm_draw <- function(n_perm, K, seed, perm) {
  if (!is.null(perm)) {
    perm <- as.matrix(perm)
    if (!is.numeric(perm)) {
      stop("`perm` must be an integer matrix with one permutation per row.",
           call. = FALSE)
    }
    if (ncol(perm) != K) {
      stop(sprintf("`perm` has %d column(s) but `mu` has %d cell types; every row must be a permutation of 1:%d.",
                   ncol(perm), K, K), call. = FALSE)
    }
    ok <- apply(perm, 1L, function(p) identical(sort(as.integer(p)), seq_len(K)))
    if (!all(ok)) {
      stop(sprintf("row(s) %s of `perm` are not permutations of 1:%d.",
                   paste(utils::head(which(!ok), 5L), collapse = ", "), K),
           call. = FALSE)
    }
    ## pi -> pi^-1, stored as K x n_perm
    perms <- matrix(0L, nrow = K, ncol = nrow(perm))
    for (b in seq_len(nrow(perm))) perms[, b] <- order(perm[b, ])
    return(perms)
  }
  ## R's RNG is not numpy's; a fixed `seed` makes a run reproducible within R,
  ## not identical to a Python reference. Supply `perm` to replay an exact
  ## randomisation.
  set.seed(seed)
  replicate(n_perm, sample.int(K, K))
}

## Correlation of one already-computed composition vector against the
## pre-centred `obs`.
.nullperm_r <- function(comp, obs_c, obs_ss) {
  x  <- comp - mean(comp)
  xx <- sum(x * x)
  if (!is.finite(xx) || xx <= 0) return(NA_real_)
  sum(x * obs_c) / sqrt(xx * obs_ss)
}

## null_permutation() for one permutation.
.nullperm_one <- function(pi, M, Fc, Fn, fc, fn, mean_of_log, spearman,
                          obs_c, obs_ss) {
  comp <- if (mean_of_log) {
    rowMeans(log2(pmax(M %*% Fc[pi, , drop = FALSE], .NULLPERM_EPS))) -
      rowMeans(log2(pmax(M %*% Fn[pi, , drop = FALSE], .NULLPERM_EPS)))
  } else {
    drop(log2(pmax(M %*% fc[pi], .NULLPERM_EPS)) -
           log2(pmax(M %*% fn[pi], .NULLPERM_EPS)))
  }
  if (spearman) comp <- rank(comp)
  .nullperm_r(comp, obs_c, obs_ss)
}

## Run the null loop, optionally over forked workers.
.nullperm_run <- function(perms, M, Fc, Fn, fc, fn, mean_of_log, spearman,
                          obs_c, obs_ss, progress, n_cores) {
  n_perm <- ncol(perms)
  one <- function(b) {
    .nullperm_one(perms[, b], M, Fc, Fn, fc, fn, mean_of_log, spearman,
                  obs_c, obs_ss)
  }

  if (n_cores > 1L) {
    n_cores <- min(n_cores, n_perm)
    parts <- parallel::splitIndices(n_perm, n_cores)
    if (progress) {
      cat(sprintf("  null_permutation: %d permutations on %d forked workers\n",
                  n_perm, n_cores))
    }
    res <- parallel::mclapply(parts, function(ix) {
      vapply(ix, one, numeric(1))
    }, mc.cores = n_cores)
    out <- unlist(res, use.names = FALSE)
    if (length(out) != n_perm) {
      stop("a forked worker returned the wrong number of permutations; ",
           "re-run with `n_cores = 1` to see the underlying error.", call. = FALSE)
    }
    return(out)
  }

  out   <- numeric(n_perm)
  every <- max(1L, n_perm %/% 20L)          # 20 updates, whatever n_perm is
  if (progress) {
    cat(sprintf("  null_permutation: %d permutations\n", n_perm))
  }
  for (b in seq_len(n_perm)) {
    out[b] <- one(b)
    if (progress && (b %% every == 0L || b == n_perm)) {
      cat(sprintf("\r  null_permutation: %d/%d (%3.0f%%)", b, n_perm,
                  100 * b / n_perm))
      utils::flush.console()
    }
  }
  if (progress) cat("\n")
  out
}

# ---------------------------------------------------------------------------
# Exported interface
# ---------------------------------------------------------------------------

#' Zero-model test for the composition term
#'
#' Permutes the cell types of the reference spectrum `mu` (leaving each cell
#' type's expression profile intact), recomputes the composition term, and
#' reports where the observed correlation between the composition term and the
#' observed log2 fold changes falls in the resulting null distribution.
#'
#' @param obs Named numeric vector of observed log2 fold changes; the names are
#'   gene identifiers and are aligned with `colnames(mu)`.
#' @param frac Data frame or matrix of cell-type proportions with samples in rows
#'   and cell types in columns. `colnames(frac)` must be the same set as
#'   `rownames(mu)`.
#' @param mu Numeric matrix of per-cell mean expression in the package
#'   convention **cell types x genes**: rownames are cell types, colnames are
#'   genes. A genes x cell types matrix is detected from the dimnames and
#'   transposed with a warning, so a matrix read from a genes-in-rows CSV works
#'   after `t()`.
#' @param is_case Logical vector over the rows of `frac`; `TRUE` marks case (AAA)
#'   samples, `FALSE` marks control (normal/donor) samples.
#' @param n_perm Number of permutations. Ignored when `perm` is supplied.
#' @param seed Random seed passed to [set.seed()]. Reproducible within R; it does
#'   not reproduce another language's random stream (see `perm`). Ignored when
#'   `perm` is supplied.
#' @param convention `"mean_of_log"` (default) takes `log2` of the per-sample
#'   composition and averages within group; `"log_of_mean"` averages the
#'   proportions within group first and takes one `log2`.
#' @param stat `"pearson"` (default) or `"spearman"`.
#' @param genes Optional character vector restricting the test to a gene subset.
#'   Defaults to the genes shared by `names(obs)` and `colnames(mu)`.
#' @param n_cores Number of forked workers, via [parallel::mclapply()]. Forking
#'   is not available on Windows, where the function warns and runs serially.
#' @param progress Print a progress line.
#' @param perm Optional integer matrix of pre-computed permutations, one row per
#'   permutation, each row a permutation of `seq_len(nrow(mu))` applied to the
#'   cell-type axis as `mu[, perm[b, ]] \%*\% f` (the same convention as permuting
#'   the columns of a genes x cell types matrix). Supply it to replay a fixed
#'   randomisation exactly, for example the permutation stream of a Python
#'   reference run; `n_perm` and `seed` are then ignored.
#'
#' @return A list of class `"cade_nullperm"` with
#'   \item{r_obs}{the observed correlation between the composition term and
#'     `obs`;}
#'   \item{null}{numeric vector of length `n_perm` with the null correlations;}
#'   \item{null_mean, null_sd}{mean and standard deviation of `null`;}
#'   \item{null_lo, null_hi}{2.5\% and 97.5\% quantiles of `null`;}
#'   \item{z}{`(r_obs - null_mean) / null_sd`;}
#'   \item{n_ge}{number of permutations with a correlation at least `r_obs`;}
#'   \item{p}{Monte Carlo p value, `(1 + n_ge) / (1 + n_perm)`;}
#'   \item{p_se}{Monte Carlo standard error, `sqrt(p * (1 - p) / n_perm)`;}
#'   \item{n_perm}{number of permutations;}
#'   \item{n}{number of genes used;}
#'   \item{K}{number of cell types;}
#'   \item{convention, stat, seed}{the settings used;}
#'   \item{elapsed}{elapsed time in seconds;}
#'   \item{mu_orientation}{always `"celltypes_x_genes"`, recording the layout
#'     `mu` was interpreted in.}
#'
#' @details
#' The permutation is applied to the small `K x n` fraction matrix instead of to
#' `mu` itself. For a permutation `pi` that reorders the cell-type axis,
#' `mu[, pi] \%*\% f` equals `mu \%*\% f[pi^-1]`; since `pi` and its inverse are both
#' uniform on the symmetric group, permuting the rows of the fraction matrix
#' draws from the same null distribution while copying only the small matrix.
#' A permutation supplied through `perm` always uses the first form, so a stored
#' stream from another implementation reproduces the same `null` values.
#'
#' `p` is a Monte Carlo estimate with standard error `p_se`; at `p` near 0.05 and
#' `n_perm` below 10000 the warning states that the resolution is not enough to
#' separate 0.05 from 0.06. If `p` sits at the resolution floor `1 / (1 + n_perm)`
#' no permutation reached `r_obs` and only a bound is reported.
#'
#' @examples
#' set.seed(1)
#' K <- 3L
#' G <- 40L
#' n <- 12L
#' mu <- matrix(stats::runif(K * G, 1, 50), K, G,
#'              dimnames = list(c("A", "B", "C"), paste0("g", seq_len(G))))
#' fr <- matrix(stats::runif(n * K), n, K,
#'              dimnames = list(paste0("s", seq_len(n)), c("A", "B", "C")))
#' fr <- fr / rowSums(fr)
#' is_case <- rep(c(FALSE, TRUE), each = n / 2)
#' ## an obs that tracks the composition term closely
#' comp <- rowMeans(log2(pmax(crossprod(mu, t(fr[is_case, , drop = FALSE])), 1e-12))) -
#'         rowMeans(log2(pmax(crossprod(mu, t(fr[!is_case, , drop = FALSE])), 1e-12)))
#' obs <- stats::setNames(comp + stats::rnorm(G, 0, 0.05), paste0("g", seq_len(G)))
#' out <- null_permutation(obs, fr, mu, is_case, n_perm = 50, progress = FALSE)
#' out$r_obs
#' out$p
#' out
#' @export
null_permutation <- function(obs, frac, mu, is_case,
                             n_perm = 2000L, seed = 1L,
                             convention = c("mean_of_log", "log_of_mean"),
                             stat = c("pearson", "spearman"),
                             genes = NULL,
                             n_cores = 1L,
                             progress = TRUE,
                             perm = NULL) {
  t_start <- proc.time()[["elapsed"]]
  convention <- match.arg(convention)
  stat       <- match.arg(stat)

  chk <- .nullperm_prepare(obs, frac, mu, is_case, genes)
  mu  <- chk$mu                                    # K x G, aligned
  K   <- nrow(mu)
  G   <- length(chk$genes)
  if (K < 2L) {
    stop("`mu` has fewer than 2 cell types; a permutation test needs at least 2.",
         call. = FALSE)
  }

  ## ---- settings ----------------------------------------------------------
  if (is.null(perm)) {
    if (length(n_perm) != 1L || !is.numeric(n_perm) || is.na(n_perm) || n_perm < 1) {
      stop("`n_perm` must be a single positive integer.", call. = FALSE)
    }
    n_perm <- as.integer(n_perm)
    if (length(seed) != 1L || !is.numeric(seed) || is.na(seed)) {
      stop("`seed` must be a single integer passed to set.seed().", call. = FALSE)
    }
  } else {
    if (!is.matrix(perm) && !is.data.frame(perm)) {
      stop("`perm` must be a matrix with one permutation per row.", call. = FALSE)
    }
    n_perm <- nrow(as.matrix(perm))
  }
  if (length(n_cores) != 1L || !is.numeric(n_cores) || is.na(n_cores) || n_cores < 1) {
    stop("`n_cores` must be a single positive integer.", call. = FALSE)
  }
  n_cores <- as.integer(n_cores)
  if (n_cores > 1L && .Platform$OS.type == "windows") {
    warning("parallel::mclapply() is serial on Windows, so `n_cores = ", n_cores,
            "` gives no speed-up; falling back to a serial loop.",
            call. = FALSE, immediate. = TRUE)
    n_cores <- 1L
  }

  ## ---- pre-computed pieces ----------------------------------------------
  M  <- t(mu)                                     # G x K
  Fc <- t(chk$frac_case)                           # K x n_case
  Fn <- t(chk$frac_ctrl)                           # K x n_ctrl
  fc <- rowMeans(Fc)
  fn <- rowMeans(Fn)
  mean_of_log <- convention == "mean_of_log"
  spearman    <- stat == "spearman"

  obs_g <- chk$obs
  if (spearman) {
    obs_c <- rank(obs_g)
    obs_c <- obs_c - mean(obs_c)
  } else {
    obs_c <- obs_g - mean(obs_g)
  }
  obs_ss <- sum(obs_c * obs_c)
  if (!is.finite(obs_ss) || obs_ss <= 0) {
    stop("`obs` has no variance; a correlation is undefined.", call. = FALSE)
  }

  ## ---- observed ----------------------------------------------------------
  comp_obs <- if (mean_of_log) {
    rowMeans(log2(pmax(M %*% Fc, .NULLPERM_EPS))) -
      rowMeans(log2(pmax(M %*% Fn, .NULLPERM_EPS)))
  } else {
    drop(log2(pmax(M %*% fc, .NULLPERM_EPS)) -
           log2(pmax(M %*% fn, .NULLPERM_EPS)))
  }
  ## Spearman correlates ranks on both sides; the null loop ranks each
  ## permuted composition, so the observed composition must be ranked too.
  if (spearman) comp_obs <- rank(comp_obs)
  r_obs <- .nullperm_r(comp_obs, obs_c, obs_ss)

  ## ---- null --------------------------------------------------------------
  perms <- .nullperm_draw(n_perm, K, seed, perm)
  null  <- .nullperm_run(perms, M, Fc, Fn, fc, fn, mean_of_log, spearman,
                         obs_c, obs_ss, progress, n_cores)

  if (!all(is.finite(null))) {
    stop(sprintf("%d of %d permutations produced an undefined correlation (a constant composition term). Check `mu`, `frac` and `is_case`.",
                 sum(!is.finite(null)), n_perm), call. = FALSE)
  }

  ## ---- summary -----------------------------------------------------------
  null_mean <- mean(null)
  null_sd   <- stats::sd(null)
  q         <- stats::quantile(null, c(0.025, 0.975), names = FALSE)
  n_ge      <- sum(null >= r_obs)
  p         <- (1 + n_ge) / (1 + n_perm)
  p_se      <- sqrt(p * (1 - p) / n_perm)
  z         <- if (is.finite(null_sd) && null_sd > 0) (r_obs - null_mean) / null_sd else NA_real_
  elapsed   <- as.numeric(proc.time()[["elapsed"]] - t_start)

  if (!is.finite(r_obs)) {
    stop("the observed composition term has no variance; a correlation is undefined.",
         call. = FALSE)
  }

  if (n_perm < 10000L && p >= 0.04 && p <= 0.09) {
    warning("p is close to 0.05 and n_perm is small; the Monte Carlo standard error is non-negligible. Consider n_perm >= 20000.",
            call. = FALSE)
  }

  out <- list(
    r_obs      = r_obs,
    null       = null,
    null_mean  = null_mean,
    null_sd    = null_sd,
    null_lo    = q[1L],
    null_hi    = q[2L],
    z          = z,
    n_ge       = n_ge,
    p          = p,
    p_se       = p_se,
    n_perm     = n_perm,
    n          = G,
    K          = K,
    convention = convention,
    stat       = stat,
    seed       = seed,
    elapsed    = elapsed,
    mu_orientation = "celltypes_x_genes"
  )
  class(out) <- "cade_nullperm"
  out
}

#' Print method for a zero-model permutation test
#'
#' @param x Output of [null_permutation()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' mu <- matrix(stats::runif(3 * 20, 1, 50), 3, 20,
#'              dimnames = list(c("A", "B", "C"), paste0("g", 1:20)))
#' fr <- matrix(stats::runif(8 * 3), 8, 3,
#'              dimnames = list(paste0("s", 1:8), c("A", "B", "C")))
#' fr <- fr / rowSums(fr)
#' obs <- stats::setNames(stats::rnorm(20), paste0("g", 1:20))
#' print(null_permutation(obs, fr, mu, rep(c(FALSE, TRUE), each = 4),
#'                        n_perm = 20, progress = FALSE))
#' @export
print.cade_nullperm <- function(x, ...) {
  cat(sprintf("<cade null permutation> %d genes, %d cell types\n", x$n, x$K))
  cat(sprintf("  convention %-12s statistic %s\n", x$convention, x$stat))
  cat(sprintf("  observed r          %+.4f\n", x$r_obs))
  cat(sprintf("  null mean (sd)      %+.4f (%.4f)\n", x$null_mean, x$null_sd))
  cat(sprintf("  null 95%% interval   [%+.4f, %+.4f]\n", x$null_lo, x$null_hi))
  cat(sprintf("  z                   %+.2f\n", x$z))
  cat(sprintf("  n_ge                %d / %d\n", x$n_ge, x$n_perm))
  cat(sprintf("  p (Monte Carlo SE)  %.4f (%.4f)\n", x$p, x$p_se))
  cat(sprintf("  seed                %s\n", format(x$seed)))
  cat(sprintf("  elapsed             %.1f s\n", x$elapsed))
  invisible(x)
}
