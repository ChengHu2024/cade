# ============================================================================
# cade: pseudo-bulk simulation with known ground truth
#
# Real data has no ground truth. The only way to know whether the
# decomposition recovers what it claims to recover is to construct data where
# the answer is known: mix single cells in known proportions, inject a known
# within-cell change, and see whether the framework separates the two.
# ============================================================================

#' Simulate pseudo-bulk samples with known composition and known within-cell change
#'
#' @param mu Per-cell mean expression, cell types x genes, LINEAR scale. Usually
#'   the same matrix used as the deconvolution reference.
#' @param frac_ctrl,frac_case Named numeric vectors of cell-type proportions for
#'   the two groups. Names must match `rownames(mu)`.
#' @param n_ctrl,n_case Number of samples per group.
#' @param intrinsic Named numeric vector of within-cell log2 fold changes to
#'   inject into the case group, e.g. `c(IL6 = 1)` doubles IL6 everywhere. Genes
#'   not named are unchanged.
#' @param noise Log-normal standard deviation per gene per sample, default 0.05.
#'   `0` gives a noiseless mixture and is only useful for testing the algebra.
#' @param seed Integer seed.
#'
#' @return A list with `bulk` (samples x genes, linear), `is_case`, `truth`
#'   (a data frame with the true composition term, the true within-cell term,
#'   and their sum).
#'
#' @details
#' The simulated bulk is constructed as
#' `sum_k f[k] * mu[k, g]` and the within-cell change is applied AFTER mixing,
#' multiplied on the case samples. This means the within-cell term is exactly
#' the injected log2 fold change, and the composition term is exactly
#' `log2(sum_k f_case mu / sum_k f_ctrl mu)` computed on the noiseless `mu`.
#' Both are therefore known to machine precision.
#'
#' The simulation is a joint test of the algebra and of the deconvolution: if
#' the recovered residual does not match `truth$intrinsic`, either the algebra
#' is wrong or the deconvolution did not recover the composition.
#'
#' @examples
#' mu <- matrix(c(100, 10, 50, 5, 20, 2), nrow = 3,
#'              dimnames = list(c("A", "B", "C"), c("g1", "g2")))
#' sim <- simulate_pseudobulk(mu,
#'          frac_ctrl = c(A = 0.6, B = 0.3, C = 0.1),
#'          frac_case = c(A = 0.2, B = 0.5, C = 0.3),
#'          n_ctrl = 5, n_case = 5, intrinsic = c(g1 = 1))
#' head(sim$truth)
#' @export
simulate_pseudobulk <- function(mu, frac_ctrl, frac_case, n_ctrl, n_case,
                                intrinsic = NULL, noise = 0.05, seed = 1) {
  mu <- as.matrix(mu)
  if (is.null(rownames(mu))) stop("`mu` must have cell-type row names.")
  ct <- rownames(mu)

  check_frac <- function(f, nm) {
    f <- as.numeric(f[ct])
    if (any(is.na(f))) stop(sprintf("`%s` is missing cell type(s).", nm))
    if (any(f < 0)) stop(sprintf("`%s` has negative values.", nm))
    if (sum(f) <= 0) stop(sprintf("`%s` sums to zero.", nm))
    f / sum(f)
  }
  fc <- check_frac(frac_ctrl, "frac_ctrl")
  fa <- check_frac(frac_case, "frac_case")
  names(fc) <- ct; names(fa) <- ct

  if (!is.null(intrinsic)) {
    nf <- setdiff(names(intrinsic), colnames(mu))
    if (length(nf))
      stop(sprintf("`intrinsic` names not in `mu`: %s", paste(nf, collapse = ", ")))
  }

  set.seed(seed)
  make <- function(f, n) {
    base <- as.numeric(crossprod(f, mu))
    m <- t(replicate(n, base * exp(stats::rnorm(length(base), 0, noise))))
    colnames(m) <- colnames(mu)
    m
  }
  bc <- make(fc, n_ctrl)
  ba <- make(fa, n_case)

  if (!is.null(intrinsic)) {
    idx <- match(names(intrinsic), colnames(ba))
    ## NOTE: must use sweep(), not `ba[, idx] * v`. A plain matrix-by-vector
    ## product recycles v in column-major order, so with more than one
    ## intrinsic gene the multipliers alternate across rows instead of
    ## applying one per column, and the injected changes cancel out. That
    ## failure is silent: the returned bulk looks plausible and the truth
    ## table is still correct, so only score_simulation() exposes it.
    ba[, idx] <- sweep(ba[, idx, drop = FALSE], 2,
                       2 ^ as.numeric(intrinsic), `*`)
  }

  bulk <- rbind(bc, ba)
  rownames(bulk) <- paste0("pb", seq_len(nrow(bulk)))
  is_case <- c(rep(FALSE, n_ctrl), rep(TRUE, n_case))

  ## ground truth from the noiseless means
  num <- as.numeric(crossprod(fa, mu))
  den <- as.numeric(crossprod(fc, mu))
  good <- is.finite(num) & is.finite(den) & num > 0 & den > 0
  comp_term <- rep(NA_real_, ncol(mu)); comp_term[good] <- log2(num[good] / den[good])
  within <- rep(0, ncol(mu))
  if (!is.null(intrinsic))
    within[match(names(intrinsic), colnames(mu))] <- as.numeric(intrinsic)

  truth <- data.frame(gene = colnames(mu),
                      composition = comp_term,
                      intrinsic = within,
                      total = comp_term + within,
                      stringsAsFactors = FALSE)

  list(bulk = bulk, is_case = is_case, truth = truth,
       f_ctrl = fc, f_case = fa, intrinsic = intrinsic, noise = noise)
}


#' Score a recovered decomposition against a simulation's ground truth
#'
#' @param decomp Output of [decompose_composition()] applied to a simulated
#'   bulk.
#' @param truth The `truth` element of [simulate_pseudobulk()].
#' @param intrinsic_genes Optional character vector of genes whose within-cell
#'   term was non-zero; used for the recovery report.
#' @param tol Tolerance for calling a residual "non-zero" when computing the
#'   false-positive rate, in log2 units. Default `0.3`.
#'
#' @return A one-row data frame with the recovery metrics.
#'
#' @details
#' The informative quantities are:
#'
#' * `cor_pred` and `median_error_pred`: does the composition term recover the
#'   truth?
#' * `cor_resid`: does the residual recover the injected within-cell term?
#' * `recovery_intrinsic`: the mean of `resid / intrinsic` over the injected
#'   genes. `1` means perfect recovery; a value near `0.3` means the framework
#'   attributes most of the injected change to composition.
#' * `false_positive_rate`: among genes with no injected change, the fraction
#'   whose residual exceeds a tolerance. These are genes wrongly called
#'   cell-intrinsic.
#'
#' @examples
#' mu <- matrix(c(100, 10, 50, 5, 20, 2), nrow = 3,
#'              dimnames = list(c("A", "B", "C"), c("g1", "g2")))
#' sim <- simulate_pseudobulk(mu, c(A = .6, B = .3, C = .1),
#'          c(A = .2, B = .5, C = .3), 5, 5, intrinsic = c(g1 = 1))
#' # a perfect decomposition, to exercise the scorer
#' d <- data.frame(gene = colnames(mu), obs = sim$truth$total,
#'                 pred = sim$truth$composition,
#'                 resid = sim$truth$total - sim$truth$composition,
#'                 pct = NA, mean_expr = 1)
#' score_simulation(d, sim$truth, "g1")
#' @export
score_simulation <- function(decomp, truth, intrinsic_genes = NULL,
                             tol = 0.3) {
  m <- merge(decomp, truth, by = "gene", suffixes = c("", "_true"))
  if (!nrow(m)) stop("no genes in common between `decomp` and `truth`.")

  ok <- is.finite(m$pred) & is.finite(m$composition)
  cor_pred <- stats::cor(m$pred[ok], m$composition[ok])
  med_err  <- stats::median(abs(m$pred[ok] - m$composition[ok]))

  okr <- is.finite(m$resid)
  ## Correlation is undefined when either side has zero variance (e.g. a
  ## decomposition that attributes everything to composition gives resid == 0
  ## for every gene). Report NA rather than a spurious number or a warning.
  cor_resid <- NA_real_
  if (sum(okr) > 2 && stats::sd(m$resid[okr]) > 0 &&
      stats::sd(m$intrinsic[okr]) > 0)
    cor_resid <- stats::cor(m$resid[okr], m$intrinsic[okr])

  recov <- NA_real_
  if (!is.null(intrinsic_genes)) {
    ii <- m$gene %in% intrinsic_genes & is.finite(m$resid) & m$intrinsic != 0
    if (any(ii))
      recov <- mean(m$resid[ii] / m$intrinsic[ii])
  }

  fp <- if (is.null(intrinsic_genes)) NA_real_ else {
    bg <- !(m$gene %in% intrinsic_genes) & is.finite(m$resid)
    if (!any(bg)) NA_real_ else mean(abs(m$resid[bg]) > tol)
  }

  data.frame(n_gene = nrow(m), cor_pred = cor_pred, median_error_pred = med_err,
             cor_resid = cor_resid, recovery_intrinsic = recov,
             false_positive_rate = fp)
}
