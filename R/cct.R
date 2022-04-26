#' @title Cauchy's method.
#' @description  Cauchy's method combines p-values regardless of dependency structure. We improved Liu's code to handle matrix of p-values.
#' @details Input a matrix of p-values and return a vector of p-values.
#' @param pvals A matrix of p-values (#SNPs by #traits).
#' @param weights Weights for the Cauchy's method. Null represents uniform weights.
#' @return A numeric vector of p-values.
#' @import pracma stats
#' @importFrom pracma isempty
#' @noRd
#' @keywords internal
cct <- function(pvals,weights=NULL){

  k = ncol(pvals)

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/k,k)
  } else if ((sum(weights < 0) > 0)){
    stop("All the weights must be positive and sum to 1")
  }

  weights = matrix(weights, ncol = 1)

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (base::sum(is.small) == 0){
    cct.stat <- tan((0.5-pvals)*pi) %*% weights
  }else{
    cct.stat <- tan((0.5-pvals)*pi) %*% weights
    index_small <- which(Rfast::rowsums(is.small,parallel = TRUE)!=0)        # the row that has nonzero

    for (i in index_small) {
      cct.stat[i] <- sum((weights[is.small[i,]]/pvals[i,is.small[i,]])/pi)
      cct.stat[i] <- cct.stat[i] + sum(weights[!is.small[i,]]*tan((0.5-pvals[i,!is.small[i,]])*pi))
    }

  }

  #### check if the test statistic is very large.

  index_large <- which(cct.stat > 1e+15)
  if(isempty(index_large)){
    pval <- 1-pcauchy(cct.stat)
  }else{
    pval <- 1-pcauchy(cct.stat)
    pval[index_large] <- (1/cct.stat[index_large])/pi
  }

  return(pval)
}


