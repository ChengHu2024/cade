\name{marker_consistency}
\alias{marker_consistency}
\title{Marker consistency between the composition term and the observed change}
\usage{
marker_consistency(obs, comp, markers = CADE_markers, min_markers = 3,
                   min_abs_comp = 0.05, pass_range = c(0.5, 2))
}
\arguments{
\item{obs}{Named numeric vector of observed log2 fold changes.}
\item{comp}{Named numeric vector of composition terms, same names as \code{obs}.}
\item{markers}{Named list of marker gene sets, for example \code{\link{CADE_markers}}.}
\item{min_markers}{Minimum number of usable markers per cell type.}
\item{min_abs_comp}{Markers whose composition term is smaller in absolute value
than this are excluded.}
\item{pass_range}{Two-element numeric range considered a pass.}
}
\value{
A data frame with one row per cell type: \code{celltype}, \code{n_markers},
\code{ratio_median}, \code{ratio_q25}, \code{ratio_q75}, \code{obs_median},
\code{comp_median}, \code{markers}, \code{verdict}.
}
\description{
For each cell type, computes the median ratio of the observed log2 fold change
to the composition term across that cell type's marker genes.
}
\details{
The ratio is a calibration statistic, not a test of reference quality. It equals
the ratio of the observed effect scale to the composition scale in a cohort, so
it is sensitive to platform: cohorts whose log2 fold changes are compressed
(typically microarray) give ratios below 1, and RNA-seq cohorts give ratios near
or above 1.

A single marker gene may itself be regulated, so the median over several markers
is used, never a single gene.
}
\examples{
set.seed(1)
g <- paste0("g", 1:200)
comp <- stats::setNames(rnorm(200, 0, 0.5), g)
obs <- stats::setNames(comp * 1.1, g)
mk <- list(A = g[1:10], B = g[11:20])
comp[mk$A] <- runif(10, 0.5, 1.5)
comp[mk$B] <- runif(10, -1.5, -0.5)
obs[mk$A] <- comp[mk$A] * 1.05
obs[mk$B] <- comp[mk$B] * 0.95
marker_consistency(obs, comp, mk)
}
