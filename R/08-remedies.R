# ============================================================================
# cade: actionable remedies ("what to do", not just "what is wrong")
#
# Design constraint (preserves authenticity):
#   This file performs string matching and appends advice. It does NOT read,
#   modify, or recompute any numeric quantity. It is called once while
#   assembling the return value of check_inputs(), so it cannot affect a single
#   digit of obs / pred / resid.
#
# Every string in this file is ASCII, because R packages must be portable.
# A Chinese-language cheat sheet is distributed alongside the package.
# ============================================================================

.input_remedies <- function(errors, warnings, d) {
  rule <- function(pat, fn) list(pat = pat, fn = fn)

  RULES <- list(
    ## ---------------- obs ----------------
    rule("`obs` is not numeric",
         function() "as.numeric(obs) -- but confirm it really is a log2 fold change, not a linear ratio, before taking logs."),
    rule("`obs` has no names",
         function() "names(obs) <- rownames(de_table). obs must carry gene names or it cannot be aligned to mu."),
    rule("at least .* are required",
         function() "At least 10 genes are required. If you only care about a few candidates, pass the full matrix (cade is not designed for single-gene inference)."),
    rule("contains .* non-finite",
         function() "obs <- obs[is.finite(obs)], or drop NA rows from the topTable() output."),
    rule("duplicated gene names",
         function() "Collapse duplicated probes first: limma::avereps() or rowsum() by gene symbol."),

    ## ---------------- frac ----------------
    rule("`frac` has no cell-type column names",
         function() "colnames(frac) <- the cell-type labels from your deconvolution output. V1/V2 ... will not do."),
    rule("`frac` has duplicated cell-type names",
         function() "Decide whether the duplication is real. If the labels merely clash, colnames(frac) <- make.unique(colnames(frac)). If two clusters are the same population, merge them: frac <- cbind(T_cell = frac[, c('T CD4','T CD8')])."),
    rule("row\\(s\\) of `frac` sum to zero",
         function() "Drop those samples: keep <- rowSums(frac) > 0; frac <- frac[keep, ]. Apply the same row selection to obs and is_case."),
    rule("`frac` contains negative",
         function() "Deconvolution output should not be negative. If it is NNLS/CIBERSORT, clamp and renormalise: frac[frac < 0] <- 0; frac <- frac / rowSums(frac)."),

    ## ---------------- mu ----------------
    rule("`mu` has no cell-type row names",
         function() "rownames(mu) <- cell_type_labels. The rows of mu must be cell types."),
    rule("`mu` has duplicated cell-type row names",
         function() "Merge the duplicated rows: mu <- rowsum(mu, group = rownames(mu)) / as.numeric(table(rownames(mu))[unique(rownames(mu))])."),
    rule("all-zero expression in `mu`",
         function() "Those cell types have no expression data in the reference. Either remove them from frac and mu together, or use a reference that contains them."),
    rule("looks log-scaled|log-scaled \\(or otherwise compressed\\)",
         function() paste0(
           "The reference must be a LINEAR per-cell mean. Common inversions: ",
           "if you hold log2(CP10k + 1), use mu <- 2^L - 1; ",
           "if log2(CP10k), use mu <- 2^L; ",
           "if ln(CP10k + 1), use mu <- expm1(L). ",
           "After inversion max(mu) should reach several hundred ",
           "(mitochondrial and ribosomal transcripts are always high).")),
    rule("looks like raw counts",
         function() paste0(
           "That is raw counts, not a per-cell mean. Normalise by library size first: ",
           "cp10k <- sweep(counts, 2, colSums(counts), '/') * 1e4; ",
           "mu <- rowsum(cp10k, group = cell_type) / as.numeric(table(cell_type)).")),

    ## ---------------- alignment ----------------
    rule("present in `mu` but not in `frac`",
         function() sprintf(
           "The two cell-type sets differ. Present in mu but not in frac: %s.",
           paste(setdiff(d$cell_types_mu, d$cell_types_frac), collapse = ", "))),
    rule("present in `frac` but not in `mu`",
         function() sprintf(
           "The deconvolution produced cell types that mu does not have: %s.",
           paste(setdiff(d$cell_types_frac, d$cell_types_mu), collapse = ", "))),
    rule("differ in case, space or encoding",
         function() "Naming only. Use trimws() and tolower() on both sides before aligning."),

    ## ---------------- is_case ----------------
    rule("`is_case` contains .* NA",
         function() paste0(
           "NA propagates everywhere and does NOT raise an error: every gene ",
           "would be labelled 'unchanged'. Fix the grouping first: ",
           "is_case <- !is.na(is_case) & is_case, or drop those samples.")),
    rule("`is_case` is not logical",
         function() "is_case <- as.logical(is_case); if it is 0/1 integer, use is_case == 1."),
    rule("has length .* but `frac` has",
         function() "The vectors must correspond one-to-one. Check whether some step dropped samples from frac or is_case but not the other."),
    rule("no TRUE \\(no case samples\\)",
         function() "Is the contrast reversed? Try is_case <- !is_case."),
    rule("no FALSE \\(no control samples\\)",
         function() "Composition attribution needs two groups. A case-only dataset cannot be analysed this way."),
    rule("smaller group",
         function() sprintf(
           "The smaller group has only %d samples, so its group mean is unstable. Aim for at least 5 per group; below that report direction only, not magnitude.",
           min(d$n_case, d$n_control, na.rm = TRUE))),

    ## ---------------- genes ----------------
    rule("gene\\(s\\) shared between `obs` and `mu`",
         function() sprintf(
           "Only %d genes are shared between obs and mu. Check that both sides use the same gene identifier (symbol vs Ensembl).",
           d$n_genes_shared)),
    rule("are absent from `mu`",
         function() sprintf(
           "%d genes are absent from the reference and will be returned as NA. If those genes matter, the reference gene list was filtered (common in 10x-filtered references); use a different reference or state the limitation.",
           d$n_gene_dropped))
  )

  map <- function(msgs) {
    vapply(msgs, function(m) {
      for (r in RULES) if (grepl(r$pat, m, ignore.case = TRUE)) return(r$fn())
      ""
    }, character(1), USE.NAMES = FALSE)
  }

  re <- map(errors); rw <- map(warnings)
  list(
    for_errors   = data.frame(message = errors,   remedy = re, stringsAsFactors = FALSE),
    for_warnings = data.frame(message = warnings, remedy = rw, stringsAsFactors = FALSE),
    n_unmatched  = sum(re == "") + sum(rw == "")
  )
}


#' Print an actionable diagnosis of a decomposition's inputs
#'
#' Runs [check_inputs()] and prints, for every problem found, a concrete fix.
#'
#' @param obs,frac,mu,is_case As in [check_inputs()].
#' @param ... Passed to [check_inputs()].
#'
#' @return Invisibly, the list from [check_inputs()].
#'
#' @details
#' `check_inputs()` tells you what is wrong; this tells you what to do about it.
#' It never touches the numbers: it calls `check_inputs()` and formats its
#' return value.
#'
#' A Chinese-language cheat sheet accompanies the package.
#'
#' @examples
#' mu <- matrix(c(500, 50, 300, 40), nrow = 2,
#'              dimnames = list(c("A", "B"), c("g1", "g2")))
#' fr <- matrix(c(0.8, 0.2, 0.8, 0.2, 0.2, 0.8, 0.2, 0.8), nrow = 4,
#'              byrow = TRUE, dimnames = list(paste0("s", 1:4), c("A", "B")))
#' diagnose_inputs(c(g1 = .5, g2 = -.2), fr, mu, c(FALSE, FALSE, TRUE, TRUE))
#' @export
diagnose_inputs <- function(obs, frac, mu, is_case, ...) {
  chk <- check_inputs(obs, frac, mu, is_case, ...)
  cat("\n========== cade input diagnosis ==========\n\n")
  if (chk$ok && !length(chk$warnings)) {
    cat("No problems found. Inputs look usable.\n\n")
    return(invisible(chk))
  }
  show <- function(df, header) {
    if (!nrow(df)) return(invisible(NULL))
    cat(header, "\n", sep = "")
    for (i in seq_len(nrow(df))) {
      cat(sprintf("  [%d] %s\n", i, df$message[i]))
      if (nzchar(df$remedy[i])) cat(sprintf("      -> %s\n", df$remedy[i]))
      else cat("      -> (no specific remedy in the lookup table; see ?check_inputs)\n")
      cat("\n")
    }
  }
  show(chk$remedies$for_errors,   "ERRORS (computation would stop):")
  show(chk$remedies$for_warnings, "WARNINGS (computation continues):")
  if (chk$remedies$n_unmatched)
    cat(sprintf("%d message(s) had no matching remedy.\n\n", chk$remedies$n_unmatched))
  cat("==========================================\n\n")
  invisible(chk)
}
