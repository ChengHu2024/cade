\name{scale_null}
\alias{scale_null}
\title{Scale a normal-versus-normal null to a different group size}
\usage{
scale_null(n_a, n_b, n_0 = 4, probs = 0.99)
}
\arguments{
\item{n_a, n_b}{Group sizes of the case-control comparison.}
\item{n_0}{Side size of the normal-versus-normal splits used for the null.}
\item{probs}{Unused; retained for interface stability.}
}
\value{
A numeric scalar, the factor to multiply the null quantiles by.
}
\description{
Rescales a null built at a reference group size to the size of the comparison
being evaluated.
}
\details{
The factor is \code{sqrt(((1/n_a + 1/n_b) / (2/n_0)))}. It assumes independent,
identically distributed samples within each group and equal variance across
groups. It is an approximation; when the data allow, build the null at the
matched group size instead.
}
\examples{
scale_null(4, 4, 4)     # 1
scale_null(76, 13, 4)   # about 0.42
}
