\name{MTAFS}
\alias{MTAFS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-trait analysis of GWAS summary statistics using MTAFS
}
\description{
Analyze multi-trait GWAS summary statistics and get the p-values of MTAFS
}
\usage{
MTAFS(Z, eigR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
  A matrix of z scores obtained from GWAS summary statistics, where its row and column correspond to SNP and trait, respectively.
}
  \item{eigR}{
  A list returned by \code{eigenDecomp}.
}
}
\details{
\code{MTAFS} applies the MTAFS method to GWAS summary statistics. It requires z scores from GWAS summary statistics and returns p-values of MTAFS.
}
\value{
  \code{MTAFS} returns a matrix of p-values of MTAFS, where its row and column correspond to SNP and trait, respectively.
}
\references{
Deng, Q., Song, C., & Lin, S. (2022). An Adaptive and Robust Method for Multi-trait Analysis of Genome-wide Association Studies Using Summary Statistics. arXiv preprint arXiv:2204.09751.
}
\author{
Qiaolan Deng
}



\seealso{
  \code{eigenDecomp}
}
\examples{
data("simulation")
Z <- simulation
Rhat <- cov(Z)
eigR <- eigenDecomp(estR = Rhat)
p_MTAFS <- MTAFS(Z,eigR)
head(p_MTAFS)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{method}
