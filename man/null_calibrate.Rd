\name{null_calibrate}
\alias{null_calibrate}
\title{Empirical null distribution for the within-cell residual}
\usage{
null_calibrate(resid_normal, resid_observed, probs = c(0.5, 0.75, 0.9, 0.95, 0.99))
}
\arguments{
\item{resid_normal}{Matrix of residuals from normal-versus-normal decompositions,
one column per comparison, genes in rows.}
\item{resid_observed}{Named numeric vector of observed residuals.}
\item{probs}{Quantile levels to report for the null.}
}
\value{
A list with \code{null_quantiles}, \code{null_pooled}, and \code{genes}, a data
frame with \code{gene}, \code{resid}, \code{abs_resid}, \code{pct},
\code{exceeds_95} and \code{exceeds_99}.
}
\description{
Builds an empirical null for the residual by decomposing comparisons among
normal samples of the same study, then reports each gene's position in that null.
}
\details{
Residuals from normal-versus-normal comparisons are compared by \strong{pooling},
not by averaging. Averaging over complementary splits cancels the residual
exactly, because a 4-versus-4 split and its complement give opposite residuals.

The pooled null describes a single group-level comparison of the same size as
the splits used to build it. A case-control comparison with larger groups has
less noise and therefore a narrower null; match group sizes, or use
\code{\link{scale_null}}.

The null does not correct for a systematic difference between the reference and
the study cohort, which affects cell-type-specific genes in a consistent
direction. It controls donor-to-donor and batch variation only.
}
\examples{
set.seed(1)
n_g <- 500
null_mat <- matrix(rnorm(n_g * 20, 0, 0.2), n_g, 20,
                   dimnames = list(paste0("g", 1:n_g), NULL))
observed <- stats::setNames(rnorm(n_g, 0, 0.2), paste0("g", 1:n_g))
observed[1:10] <- observed[1:10] * 5
out <- null_calibrate(null_mat, observed)
head(out$genes)
table(out$genes$exceeds_99)
}
