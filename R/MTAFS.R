#' @export
MTAFS <- function(Z, eigR){

  # load values from eigendecompostion
  U <- eigR$eigenD$vectors
  lambda <- eigR$eigenD$values
  n.pcs <- ncol(U)
  n.snps <- nrow(Z)
  d <- eigR$q

  # to save p-values
  results <- matrix(NA,nrow = n.snps, ncol = length(d)+1)

  # For each level E in {2,q1,q2,q3,T}, get p_{AF(E)}

  j=1
  for (i in d) {
    lambda_inv <- diag(1/sqrt(lambda[1:i]))
    z_pcs <- Z %*% U[,1:i] %*% lambda_inv # transformed z scores
    results[,j] <- af(z_pcs, weights = NULL)
    j=j+1
  }

  # MTAFS
  results[,j] <- cct(results[,1:(j-1)])
  colnames(results) <- c(paste0(as.character(eigR$q),"(",as.character(round(eigR$v,2)*100),"%",")"),"MTAFS")

  return(results)
}










