#' Hidden Markov Model Functions
#' 
#' R implementation of HMMs described in processed report except function HMMLikelihood renamed to R_HMMLikelihood and changed to compute
#' values for all capture histories and return  lnl, alpha, phi, v, dmat, and gamma values. loglikelihood is called with a fitted hmm model
#' and then computes the gamma,dmat and delta matrices and calls R_HMMLikelihood function. These are not used by the fitting code.
#'  
#' @param x single observed sequence (capture history) 
#' @param first occasion to initiate likelihood calculation for sequence 
#' @param m number of states
#' @param T number of occasions; sequence length
#' @param dmat observation probability matrices
#' @param gamma transition matrices
#' @param delta initial distribution
#' @export R_HMMLikelihood
#' @return return log-likelihood, alpha, v and phi arrays 
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. See page 45. 
R_HMMLikelihood=function(x,first,m,T,dmat,gamma,delta)
{  
	# Arguments:
	# x: observed sequence (capture (encounter) history)
	# first: occasion to start sequence
	# m: number of states
	# T: number of occasions; sequence length
	# dmat: array of occasion specific observation probabilty matrices
	# gamma: array of occasion specific transition matrices
	# delta: initial state distribution
	# Other variables:
	# lnl: log likelihood value
	# phi: alpha/sum(alpha) sequence as defined in Zucchini/MacDonald
	# v: temp variable to hold phi calculations
	# u: sum(v)
	alpha=array(NA,dim=c(nrow(x),T,m))
	phimat=array(NA,dim=c(nrow(x),T,m))
	vmat=array(NA,dim=c(nrow(x),T,m))
	lnl=rep(NA,nrow(x))
	for(i in 1:nrow(x))
	{
		#occ=1
		# Assign prob state vector for initial observation: delta*p(x_first)
		v=delta[i,]%*%diag(dmat[i,first[i],x[i,first[i]],]) 
		# Compute log-likelihood contribution for first observation; for
		# models that condition on first observation u=1,lnl=0
		u=sum(v)
		phi=v/u
		alpha[i,first[i],]=v
		vmat[i,first[i],]=v
		phimat[i,first[i],]=phi
		lnl[i]=log(u)
		# Loop over occasions for this encounter history (x)
		for(t in (first[i]+1):T)
		{
			#occ=occ+1
			# Compute likelihood contribution for this occasion
			v=phi%*%gamma[i,t-1,,]%*%diag(dmat[i,t,x[i,t],])  
			u=sum(v)
			lnl[i]=lnl[i]+log(u)
			# Compute updated state vector
			phi=v/u
			# Compute alpha vector
			vmat[i,t,]=v
			alpha[i,t,]=alpha[i,t-1,]%*%gamma[i,t-1,,]%*%diag(dmat[i,t,x[i,t],])
			phimat[i,t,]=phi
		}
	}
	return(list(lnl=lnl,alpha=alpha,phi=phimat,v=vmat,dmat=dmat,gamma=gamma))
} 
