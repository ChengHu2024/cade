\name{cross_cohort_agreement}
\alias{cross_cohort_agreement}
\title{Cross-cohort agreement of the composition term and the observed change}
\usage{
cross_cohort_agreement(obs_list, comp_list = NULL, min_common = 100)
}
\arguments{
\item{obs_list}{Named list; each element a named numeric vector of observed
log2 fold changes for one cohort.}
\item{comp_list}{Named list of composition terms, same order and names as
\code{obs_list}. Optional.}
\item{min_common}{Minimum number of shared genes to report a pair.}
}
\value{
A list with \code{obs} and, when \code{comp_list} is given, \code{comp}. Each is
a data frame with columns \code{a}, \code{b}, \code{n}, \code{pearson},
\code{spearman}, \code{sign_agree}.
}
\description{
Compares two or more cohorts gene by gene and reports, for each pair, the
correlation of the observed log2 fold changes and of the composition terms.
}
\details{
In four bulk cohorts of one disease, the composition term agreed at
\eqn{r = 0.83}--\eqn{0.93} across cohorts while the observed change agreed at
only \eqn{r = 0.34}--\eqn{0.59}. A decomposition residual is the difference of
the two, so it is limited by the reproducibility of the observed change, not by
that of the composition term.

Cohorts measured on the same platform agree better. Group by platform before
comparing residuals across cohorts.

A per-cohort rescaling does not repair the disagreement. Z-scoring and
percentile ranks are monotone transformations, so Spearman correlation is
unchanged by construction, and Pearson after z-scoring equals Pearson on the raw
values. Verified on four bulk cohorts: raw and z-scored Pearson both span
0.346--0.595, and Spearman is 0.3561 under all three transforms. Percentile
ranks also make any sign-agreement statistic trivially 100\% because every rank
is positive. Stratification does help: grouping by platform raised the
observed-change correlation from 0.346--0.400 to 0.595 (array versus array) and
0.590 (RNA-seq versus RNA-seq).
}
\examples{
set.seed(1)
g <- paste0("g", 1:500)
c1 <- stats::setNames(rnorm(500, 0, 0.4), g)
o1 <- stats::setNames(c1 + rnorm(500, 0, 0.2), g)
o2 <- stats::setNames(c1 * 1.05 + rnorm(500, 0, 0.25), g)
res <- cross_cohort_agreement(list(A = o1, B = o2), list(A = c1, B = c1))
res$obs
res$comp
}
