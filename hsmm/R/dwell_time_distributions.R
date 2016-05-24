#' Pre-defined dwell time distributions for HSMM models
#' 
#' Various discrete dwell time distributions defined in Zucchini et al second edition.
#' Note: if you use unstructured without a geometric tail (gt), the number of dwell times must exceed the longest 
#' observed sequence in the state or -lnl could become undefined (Nan). 
#' 
#' @param theta parameter vector for dwell time distributions 
#' @param n number of discrete states for dwell time  
#' @return vector of probabilities for dwell times
#' @author Jeff Laake
#' @export geometric shifted_poisson shifted_negbinomial shifted_binomial unstructured unstructured_gt 
#' @keywords utility
#' @aliases shifted_poisson shifted_negbinomial shifted_binomial unstructured_gt unstructured
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
#' Compute dwell time distribution function values for HSMM models
#' 
#' @param theta parameter vector for dwell time distributions 
#' @param n number of discrete states for dwell time  
#' @param dtf dwell time distribution function
#' @return vector of probabilities for dwell times from approximation using geometric tail
#' @export
dt.values=function(dtf,n,theta)
{
	m=dtf$steps
	values=dtf$fct(theta,n)
	if(m>1)
		z=1-values[m]/(1-sum(values[1:(m-1)]))
	else
		z=1-values[m]
	if(n>m+1)
	   values[(m+1):n]=values[m]*z^((m+1):n-m)
	return(values[1:n])
}
#' Plot dwell time distributions for HSMM models
#' 
#' @param object fitted hsmm model object object
#' @param range single value or vector of values with length = number of states; range value used from 1:range to plot distribution
#' @param dm vector for combining parameter values for dt distribution
#' @param labels labels to be used for states in plotting
#' @param ... plot parameters
#' @return none
#' @export
plot_dt=function(object,range,dm=NULL,labels=NULL,...)
{
	if(is.null(dm))dm=lapply(1:length(object$dtf),function(x) matrix(c(1),nrow=1))
	if(is.null(labels))labels=paste("state",1:length(object$dtf))
	pars=split(object$results$par,rep(1:length(object$type),object$type))
	if(length(range)==1)range=rep(range,length(object$dtf))
	for (i in 1:length(object$dtf))
	{
		plot(1:range[i],dt.values(object$dtf[[i]],range[i],sum(pars[[i]]*dm[[i]][1,])),xlab="Time",ylab=paste("Dwell time distribution:",labels[i]),type="b",...) 
		if(nrow(dm[[i]])>1)
			for(j in 2:nrow(dm[[i]]))
				lines(1:range[i],dt.values(object$dtf[[i]],range[i],sum(pars[[i]]*dm[[i]][j,])),lty=j,pch=j,type="b",...) 
	}
	invisible()
}




