#' @export
# function to conduct eigendecomposition

eigenDecomp <- function(estR, varExp=NULL){

  # conduct eigen-decomposition
  eig <- eigen(estR)

  # obtain U and lambda
  U <- eig$vectors
  lambda <- eig$values

  # get variance explained by eigenvalues
  percentExp <- cumsum(lambda/sum(lambda))

  if (!is.null(varExp)) {

    varExp <- sort(varExp)

    if (min(percentExp) > min(varExp)) {
      stop(paste0("The minimum percentage of variance explained is ", round(min(percentExp),2)))
    }

    v <- varExp

  } else {
    # percentage explained by the first two eigenvalues
    minExp <- percentExp[2]

    # set v1, v2, and v3
    v <- seq(minExp,1,length.out = 5)
  }




  # find the number of eigenvector corresponding to variance explained: q1, q2, q3
  nEigen <- vector()
  percentQ <- vector() # percentage actually explained by q's
  for (per in v) {
    tmp <- percentExp <= per
    if (any(tmp)) {
      percentQ <- c(percentQ,per)
      nEigen <- c(nEigen,max(which(tmp)))
    }
  }


  return(list(q=nEigen, v=percentQ, eigenD=eig))
}
