#' Pre-defined dwell time distributions for HSMM models
#' 
#' Various discrete dwell time distributions defined in Zucchini et al second edition.
#' Note: if you use unstructured without a geometric tail (gt), the number of dwell times must exceed the longest 
#' observed sequence in the state or -lnl could become undefined (Nan). 
#' 
#' @param theta parameter vector for dwell time distributions 
#' @param n number of discrete states for dwell time  
#' @return vector of probabilities for well times
#' @author Jeff Laake
#' @export geometric shifted_poisson shifted_negbinomial shifted_binomial unstructured unstructured_gt
#' @keywords utility
#' @references Zucchini, W., MacDonald, I.L. and Langrock, R. 2016. Hidden Markov Models for Time Series: 
#' An introduction using R, 2nd ed. CRC Press.
geometric=function(theta,n)
	dgeom((1:n)-1,plogis(theta))

shifted_poisson=function(theta,n)
	dpois((1:n)-1,exp(theta))

shifted_negbinomial=function(theta,n)
	dnbinom((1:n)-1, exp(theta[1]), plogis(theta[2]))

shifted_binomial=function(theta,n)
	dbinom((1:n)-1,plogis(theta))

unstructured_gt=function(theta,n)
{
	if(length(theta)!=n) stop("unstructured mismatch on parameters")
	p=c(exp(theta),1)
	p=p/sum(p)
	return(p[1:n])
}

unstructured=function(theta,n)
{
	if(length(theta)!=n) stop("unstructured mismatch on parameters")
	p=c(exp(theta),1)
	p=p/sum(p)
	return(p)
}




