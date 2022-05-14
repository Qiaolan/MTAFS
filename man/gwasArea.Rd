\name{gwasArea}
\alias{gwasArea}
\docType{data}
\title{
  An example dataset for MTAFS
}
\description{
 This is an example dataset for MTAFS and it is a subset of the real dataset of 212 Area image-derived phenotypes.
}
\usage{data("gwasArea")}
\format{
\code{gwasArea} is a matrix of z scores, where its row and column correspond to SNP and trait, respectively. It includes 300 SNPs and 212 traits.
}

\references{
Deng, Q., Song, C., & Lin, S. (2022). An Adaptive and Robust Method for Multi-trait Analysis of Genome-wide Association Studies Using Summary Statistics. arXiv preprint arXiv:2204.09751.
}
\examples{
data(gwasArea)
dim(gwasArea)
}
\keyword{datasets}
