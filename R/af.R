#' @title Adaptive Method in the spirit of Adaptive Fisher.
#' @description  Adaptive Fisher considers the sum of negative log ordered p-values. The p-value of its test statistic can be approximated by Davies's method or Saddlepoint approximation
#' @details Input a matrix of p-values and return a vector of p-values.
#' @param Z A matrix of independent z scores.
#' @param weights Weights for the Cauchy's method. Null represents uniform weights.
#' @return A numeric vector of p-values.
#' @importFrom  CompQuadForm davies
#' @importFrom  Rfast rowSort colCumSums
#' @import  mvtnorm stats
#' @noRd
#' @keywords internal
af <- function(Z, weights=NULL){

  n_col = 1:ncol(Z)
  n.snp = nrow(Z)
  K = ncol(Z)

  # transform to log two-sided p-value
  logp <- log(2) + pnorm(abs(Z),lower.tail = FALSE,log.p = TRUE)

  # sorted
  logp_sort <- t(rowSort(-logp,descending = TRUE,parallel = TRUE))
  rm(logp)

  # cumsum of sorted log p-values
  logp_cumsum <- colCumSums(logp_sort)
  rm(logp_sort)

  # davies method
  davies_v <- Vectorize(davies,"q")



  p_cumsum_davies = matrix(NA,ncol=length(n_col),nrow=n.snp)

  j=1
  for (n in n_col) {
    w=ifelse(n/1:K <= 1,n/1:K,1)
    lambda=0.5*w
    suppressWarnings({tmp = davies_v(logp_cumsum[j,],lambda=lambda,h=rep(2,K),lim=1e5,acc=1e-8)})
    pval = unlist(tmp[3,])

    index.fault = which(unlist(tmp[2,])==1)
    index.wrongp = which((pval<=0 | pval>1))
    index = union(index.fault,index.wrongp)

    # if davies fail to converge, turn to saddlepoint approximation
    if(length(index)>0) pval[index] = saddle_approx(logp_cumsum[n,index],df=rep(2,K),a=lambda,lower.tail = FALSE)

    p_cumsum_davies[,j] = pval
    j=j+1
  }

  # test statistic
  p_T <- cct(p_cumsum_davies, weights = weights)


  return(p_T)
}


# saddlepoint functions, revised from survey package

saddle_approx <- function (x, df, a, lower.tail = TRUE)
{
  satterthwaite <- function(a, df) {
    if (any(df > 1)) {
      a <- rep(a, df)
    }
    tr <- mean(a)
    tr2 <- mean(a^2)/(tr^2)
    list(scale = tr * tr2, df = length(a)/tr2)
  }

  sat <- satterthwaite(a, df)
  guess <- pchisq(x/sat$scale, sat$df, lower.tail = lower.tail)


  lambda <- rep(a, df)
  sad <- sapply(x, saddle, lambda = lambda)
  if (lower.tail){
    sad <- 1 - sad
  }
  guess <- ifelse(is.na(sad), guess, sad)

  return(guess)

}


saddle <- function (x, lambda)
{
  d <- max(lambda)
  lambda <- lambda/d
  x <- x/d
  k0 <- function(zeta) -sum(log(1 - 2 * zeta * lambda))/2
  kprime0 <- function(zeta) sapply(zeta, function(zz) sum(lambda/(1 -
                                                                    2 * zz * lambda)))
  kpprime0 <- function(zeta) 2 * sum(lambda^2/(1 - 2 * zeta *
                                                 lambda)^2)
  n <- length(lambda)
  if (any(lambda < 0)) {
    lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
  }
  else if (x > sum(lambda)) {
    lmin <- -0.01
  }
  else {
    lmin <- -length(lambda)/(2 * x)
  }
  lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
  hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
                     upper = lmax, tol = 1e-08)$root
  w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
  v <- hatzeta * sqrt(kpprime0(hatzeta))
  if (abs(hatzeta) < 1e-04)
    NA
  else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}

