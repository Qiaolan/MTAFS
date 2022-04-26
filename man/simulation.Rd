\name{simulation}
\alias{simulation}
\docType{data}
\title{
  Simulation data for MTAFS
}
\description{
 This is an simulation dataset.
}
\usage{data("simulation")}
\format{
\code{simulation} is a matrix of z scores, where its row and column correspond to SNP and trait, respectively. It includes 300 SNPs and 212 traits.
}

\references{
Deng, Q., Song, C., & Lin, S. (2022). An Adaptive and Robust Method for Multi-trait Analysis of Genome-wide Association Studies Using Summary Statistics. arXiv preprint arXiv:2204.09751.
}
\examples{
data(simulation)
dim(simulation)
}
\keyword{datasets}
