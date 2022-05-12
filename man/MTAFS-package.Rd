\name{MTAFS-package}
\alias{MTAFS-package}
\docType{package}
\title{
MTAFS
}
\description{
This package provides functions for multi-trait analysis of GWAS summary statistics.
}
\details{
\tabular{ll}{
Package: \tab MTAFS\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2022-05-11\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains two main functions,
\code{eigenDecomp} and \code{MTAFS}.
\code{eigenDecomp} conducts eigendecompositions and its output is used by \code{MTAFS}. \code{MTAFS} analyzes multiple GWAS summary statistics and return p-values.
}
\author{
Qiaolan Deng, Chi Song, Shili Lin

Maintainer: Qiaolan Deng <deng.295@osu.edu>
}
\references{
Deng, Qiaolan, Chi Song, and Shili Lin. "An Adaptive and Robust Method for Multi-trait Analysis of Genome-wide Association Studies Using Summary Statistics." arXiv preprint arXiv:2204.09751 (2022).
}
\keyword{ package }
\seealso{
\code{\link{MTAFS}}, \code{\link{eigenDecomp}}.
}
\examples{
# load a simulated dataset
data("simulation")
Z <- simulation

# estimate the variance-covariance matrix
Rhat <- cov(Z)


## the default setting

# conduct eigendecomposition
eigRhat <- eigenDecomp(Rhat)

# use MTAFS to analyze the summary statistics
pMTAFS <- MTAFS(Z, eigRhat)
head(pMTAFS)


## user-defined setting

# conduct eigendecomposition
eigRhat2 <- eigenDecomp(Rhat, varExp = c(0.4,0.6,0.8))

# use MTAFS to analyze the summary statistics
pMTAFS2 <- MTAFS(Z, eigRhat2)
head(pMTAFS2)


}
