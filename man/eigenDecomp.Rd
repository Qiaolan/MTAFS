\name{eigenDecomp}
\alias{eigenDecomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Eigendecomposition for MTAFS
}
\description{
 Eigendecomposition on the estimated variance-covariance matrix R.
}
\usage{
eigenDecomp(estR, varExp=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estR}{
  The variance-covariance matrix R (T by T) estimated by linkage disequillibrium score regression (recommended) or sample covariance.
}
  \item{varExp}{
  A vector of user-defined percentage of variance explained, e.g., \code{varExp=c(0.2,0.5,0.9)}. The default is NULL, which selects five percentages: the percentage explained by first eigenvector (denoted as v0), 1, and three percentages evenly distributed between v0 and 1.
}
}
\details{
 \code{eigenDecomp} takes the estimated R as the input. It calculates the number of eigenvectors (q's) used in MTAFS and their corresponding variance explained (v's).
}
\value{
 \code{eigenDecomp} returns a list.
  \item{q}{Number of eigenvectors included in MTAFS.}
  \item{v}{Percentage of variance explained by the selected eigenvectors.}
  \item{eigenD}{A list containing all eigenvectors and eigenvalues}
}
\references{
Deng, Q., Song, C., & Lin, S. (2022). An Adaptive and Robust Method for Multi-trait Analysis of Genome-wide Association Studies Using Summary Statistics. arXiv preprint arXiv:2204.09751.
}
\author{Qiaolan Deng}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(gwasArea)
Rhat <- cov(gwasArea)
eigR <- eigenDecomp(estR = Rhat)
str(eigR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

