#' @title Logit-Beta density function
#' @aliases rlogitbeta
#' @rdname dlogitbeta
#' @description Calculates the PDF of the logit-beta distribution. The logit-beta distribution is defined as the 
#' logit transformation of a beta distributed random variable
#' @param x vector of quantiles
#' @param a first shape parameter of a beta distribution
#' @param b second shape parameter of a beta distribution
#' @param n number of observations
#' @return density or random variable generation
#' @author Devin S. Johnson
#' @export

dlogitbeta = function(x, a, b, log=FALSE){
  if(!log){
    return(dbeta(plogis(x), a, b)*plogis(x)/(1+exp(x)))
  } else{
    return(dbeta(plogis(x), a, b, log=TRUE) + plogis(x, log=TRUE) - log((1+exp(x))))
  }
}

#' @rdname dlogitbeta
#' @export
rlogitbeta = function(n, a, b){
  return(qlogis(rbeta(n,a,b)))
}