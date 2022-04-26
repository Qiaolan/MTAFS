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
eigenDecomp(estR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estR}{
  The variance-covariance matrix R (T by T) estimated by linkage disequillibrium score regression (recommended) or sample covariance.
}
}
\details{
 \code{eigenDecomp} takes the estimated R as the input. It calculates the number of eigenvectors (q's) used in MTAFS and their corresponding variance explained (v's).
}
\value{
 \code{eigenDecomp} returns a list.
  \item{q}{Number of eigenvectors included in MTAFS. That is, \eqn{q_1,q_2,q_3}}
  \item{v}{Percentage of variance explained by \eqn{2,q_1,q_2,q_3,T}}
  \item{eigenD}{A list containing all eigenvectors and eigenvalues}
}
\references{
Deng, Q., Song, C., & Lin, S. (2022). An Adaptive and Robust Method for Multi-trait Analysis of Genome-wide Association Studies Using Summary Statistics. arXiv preprint arXiv:2204.09751.
}
\author{Qiaolan Deng}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(simulation)
Rhat <- cov(simulation)
eigR <- eigenDecomp(estR = Rhat)
str(eigR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

