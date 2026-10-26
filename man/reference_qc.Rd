\name{reference_qc}
\alias{reference_qc}
\title{Diagnostic for a single-cell reference matrix}
\usage{
reference_qc(mu, mu_se, n_cells = NULL, rel_se_soft = 0.25)
}
\arguments{
\item{mu}{Numeric matrix, genes in rows, cell types in columns.}
\item{mu_se}{Numeric matrix of standard errors, same shape.}
\item{n_cells}{Named integer vector of cell counts per cell type. Optional.}
\item{rel_se_soft}{Median relative standard error above which a cell type is
flagged.}
}
\value{
A data frame with one row per cell type: \code{celltype}, \code{n_cells},
\code{median_rel_se}, \code{q75_rel_se}, \code{pct_rel_se_gt_50},
\code{flagged}.
}
\description{
Reports, for each cell type, the number of cells behind the reference mean and
the median relative standard error across genes.
}
\details{
A cell type whose median relative standard error exceeds \code{rel_se_soft} is
flagged. Such a cell type still enters the decomposition, but its contribution
to the composition term carries more uncertainty than the interval from
\code{\link{decompose_uq}} assumes unless that uncertainty is passed in as
\code{mu_se}.
}
\examples{
set.seed(1)
mu <- matrix(runif(4 * 200, 1, 10), 200, 4,
             dimnames = list(NULL, c("A", "B", "C", "D")))
se <- mu * 0.2
reference_qc(mu, se, n_cells = c(A = 5000, B = 2000, C = 800, D = 300))
}
