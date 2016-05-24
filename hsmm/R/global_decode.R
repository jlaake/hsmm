#' Global decoding of HMM 
#' 
#' Computes sequence of state predictions for each individual
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param ddl design data list
#' @param state.names names for states used to label output; if NULL uses strata.labels + Dead state
#' @author Jeff Laake
#' @return matrix of state predictions
#' @export global_decode
#' @keywords utility
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 82.
global_decode=function(object,ddl=NULL,state.names=NULL)
{  	
	parmlist=hsmm_likelihood(object$results$par,object$type,object$data,ddl,object$dtf,object$fct_gamma,object$fct_dmat,
			            object$fct_delta,object$pformula,object$omega,debug=FALSE,mat=TRUE)
	dmat=parmlist$dmat
	gamma=parmlist$gamma
	delta=parmlist$delta
	x=object$data
	T=ncol(object$data)
	m=ncol(dmat[1,1,,])
	first=rep(1,nrow(x))
	states=matrix(NA,nrow=nrow(x),ncol=T)
	state.names=do.call("c",sapply(1:length(object$dtf),function(x) paste(state.names[x],1:object$dtf[[x]]$steps,sep="")))
	for(i in 1:nrow(x))
	{
		psi=matrix(NA,nrow=T,ncol=m)
		# Assign psi value at first occasion
		psi[first[i],]=delta[i,]%*%diag(dmat[i,first[i],x[i,first[i]],]) 
		for(t in (first[i]+1):T)
		{
			psi[t,]=apply(psi[t-1,]*gamma[i,t-1,,],2,max)%*%diag(dmat[i,t,x[i,t],])
		}
		states[i,T]=which.max(psi[T,])
		for(t in (T-1):first[i])
			states[i,t]=which.max(psi[t,]*gamma[i,t,,states[i,t+1]])
	}
	states=t(apply(states,1,function(x) state.names[x]))
	return(states)
}


