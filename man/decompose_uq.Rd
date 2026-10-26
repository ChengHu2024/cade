\name{decompose_uq}
\alias{decompose_uq}
\title{Uncertainty-aware decomposition of a bulk log2 fold change}
\usage{
decompose_uq(obs, f_ctrl, f_case, mu, mu_se, f_se_ctrl = NULL,
             f_se_case = NULL, n_boot = 2000, seed = 1, ci = 0.95)
}
\arguments{
\item{obs}{Named numeric vector of observed log2 fold changes per gene.}
\item{f_ctrl, f_case}{Named numeric vectors of group-mean cell-type proportions.}
\item{mu}{Numeric matrix of reference means (genes in rows, cell types in columns), linear scale.}
\item{mu_se}{Numeric matrix of standard errors of \code{mu}, same shape.}
\item{f_se_ctrl, f_se_case}{Standard errors of the proportions.}
\item{n_boot}{Number of bootstrap replicates.}
\item{seed}{Random seed.}
\item{ci}{Credible level of the interval.}
}
\value{
A data frame with one row per gene: \code{gene}, \code{obs}, \code{comp_med},
\code{resid_med}, \code{resid_lo}, \code{resid_hi}, \code{p_same_sign},
\code{direction}, \code{call}.
}
\description{
Propagates uncertainty in the reference matrix and in the estimated cell-type
proportions through the composition term, and returns an interval for the
within-cell residual.
}
\details{
The decomposition identity is \code{resid = obs - comp = d_case - d_ctrl},
where \code{d} is a group's deviation from the reference. The residual is a
difference of two deviations, not the case group's deviation.

The interval is a parametric bootstrap. Coverage at nominal 95\% was verified at
95.0--95.2\% on three single-cell references with direction accuracy of
99.9--100\%.

\strong{Limitation.} The bootstrap captures estimation noise only. It assumes
the point estimate is unbiased, so systematic bias such as a reference whose
cell states differ from the studied tissue is not reflected in the interval.
}
\examples{
set.seed(1)
mu <- matrix(runif(4 * 100, 1, 10), 100, 4,
             dimnames = list(paste0("g", 1:100), c("A", "B", "C", "D")))
mu_se <- mu * 0.1
f1 <- c(A = .5, B = .2, C = .2, D = .1)
f2 <- c(A = .1, B = .4, C = .3, D = .2)
obs <- log2(as.vector(f2 \%*\% t(mu)) / as.vector(f1 \%*\% t(mu)))
names(obs) <- rownames(mu)
out <- decompose_uq(obs, f1, f2, mu, mu_se,
                    f_se_ctrl = rep(0.005, 4), f_se_case = rep(0.005, 4),
                    n_boot = 200)
table(out$call)
}
