#' Compute transition matrix for HMM from HSMM
#' 
#' Computes expanded transition matrix for HMM from aggregated state transition
#' matrix and dwell time distributions.
#' 
#' @param omega transition matrix for aggregated states
#' @param dm list of dwell time distribution vectors; each list element is vector of dwell time probabilities in order of aggregated states  
#' @param eps epsilon value for setting small probabilities to 0; anything less than eps is set to exactly 0.  
#' @return transition matrix for hmm to approximate hsmm
#' @author Walter Zucchini
#' @export hsmm2hmm
#' @keywords utility
#' @references Zucchini, W., MacDonald, I.L. and Langrock, R. 2016. Hidden Markov Models for Time Series: 
#' An introduction using R, 2nd ed. CRC Press.
###############################################################################
hsmm2hmm<-function(omega,dm,eps=1e-10){
	mv<- sapply(dm,length)
	m <- length(mv)
	G <- matrix(0,0,sum(mv))
	for (i in 1:m){
		mi <- mv[[i]]
		F <- cumsum(c(0,dm[[i]][-mi]))
		ci <- ifelse(abs(1-F)>eps,dm[[i]]/(1-F),1)
		cim <- ifelse(1-ci>0,1-ci,0)
		Gi <- matrix(0,mi,0)
		for (j in 1:m){
			if(i==j){ if(mi==1)
				{ Gi <- cbind(Gi,c(rep(0,mv[[j]]-1),cim))} else
				{ Gi <- cbind(Gi,rbind(cbind(rep(0,mi-1),diag(cim[-mi],mi-1,mi-1)),
									c(rep(0,mi-1),cim[[mi]])))}
			} else { if(mi==1)
				{  Gi <- cbind(Gi,matrix(c(omega[[i,j]]*ci,rep(0,mv[[j]]-1)),1))}else
				{  Gi <- cbind(Gi,cbind(omega[[i,j]]*ci,matrix(0,mv[[i]],mv[[j]]-1)))}
			}
		}
		G <- rbind(G,Gi)
	}
	G
}
#' Computes dwell time distribution values and transition matrix
#' 
#' Create expanded transition matrix for HMM from aggregated state transition
#' matrix (omega) and dwell time distributions (dtf) calling hsmm2hmm.
#' 
#' @param pars list of parameter vectors; length of list is m - number of states
#' @param dtf list of length m with a function for each state
#' @param mv vector with number of steps (m_i) for each dwell time distribution; length of mv is m - the number of states
#' @param omega known m by m state transition matrix
#' @return expanded state transition matrix defining HSMM
#' 
gamma_dtd=function(pars, dtf, mv, omega)
{
	if(length(pars)!=length(dtf)) stop("Internal error: pars and dtf do not match")
	if(length(mv)!=length(dtf)) stop("Internal error: mv and dtf do not match")
	dm=vector("list",length=length(mv))
	for(i in 1:length(mv))
		dm[[i]]=dtf[[i]](pars[[i]], mv[i])
	return(hsmm2hmm(omega, dm))
}
#' Create expanded state-dependent observation matrix for HMM from HSMM
#' 
#' Creates expanded state-dependent matrix for HMM from aggregated state-dependent observation
#' matrix.
#' 
#' @param dmat state-dependent observation matrix for aggregated states
#' @param mv vector with number of steps in each dwell time distribution  
#' @return expanded state-dependent observation matrix for hmm to approximate hsmm
#' @author Jeff Laake
#' @export dmat_hsmm2hmm
#' @keywords utility
###############################################################################
dmat_hsmm2hmm<-function(dmat,mv)
{
	m <- length(mv)
	return(dmat[,rep(1:m,mv)])
}
