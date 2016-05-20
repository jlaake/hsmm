#' Hidden semi-Markov likelihood function Pre-defined dwell time distributions for HSMM models
#' 
#' Various discrete dwell time distributions defined in Zucchini et al second edition
#' 
#' @param pars parameter vector  
#' @param type vector with number of parameters in each of the mv dwell time distributions followed by number of parameters in state dependent observation matrix (dmat) and initial vector (delta)  
#' @param data matrix of observations; row is a sequence; column is an occasion; value is 1 to k where k is number of possible observations which match rows of dmat
#' @param mv vector with number of steps (m_i) for each dwell time distribution; length of mv is m - the number of states
#' @param dtf list of length m with a function for each state
#' @param fct_dmat function that specifies the state dependent observation matrix (dmat)
#' @param fct_delta function that specifies the initial probability vector
#' @param omega known m by m state transition matrix
#' @param debug if TRUE, prints the parameter values and -log likelihood
#' @param mat if TRUE, return matrices rather than -log likelihood
#' @return negative log-likelihood value or list of matrices depending on value of mat
#' @author Jeff Laake
#' @export hsmm_likelihood
#' @keywords utility
#' @useDynLib hsmm
#' @references Zucchini, W., MacDonald, I.L. and Langrock, R. 2016. Hidden Markov Models for Time Series: 
#' An introduction using R, 2nd ed. CRC Press.

hsmm_likelihood=function(pars,type,data,ddl,dtf,fct_gamma,fct_dmat,fct_delta,pformula,omega,debug=FALSE,mat=FALSE)
{
	T=ncol(data)
	pars=split(pars,rep(1:length(type),type))
	mv=sapply(dtf,function(x) x$steps)
	# create transition matrix calling list of dwell time distributions and then hsmm2hmm in fct_gamma
	gamma=fct_gamma(pars=pars[1:length(mv)],type=type[1:length(mv)],dtf=dtf,ddl=ddl,mv=mv,omega=omega)
	# create state dependent observation matrix
	dmat=fct_dmat(pars=pars[[length(mv)+1]],pformula,ddl,mv)
	# create initial state matrix
	if(length(pars)==length(mv)+2)
    	delta=fct_delta(pars=pars[length(mv)+2],mv,nrow(data))
	else
		delta=fct_delta(pars=NULL,mv,nrow(data))
	if(mat)return(list(gamma=gamma,dmat=dmat,delta=delta))
	if(debug) sapply(pars,function(x) cat("\n",x))   
	neglnl=hmm.lnl(x=data,start=rep(1,nrow(data)),m=sum(mv),T=T,dmat,gamma,delta,freq=rep(1,nrow(data)),debug=FALSE)
	#neglnl=-sum(R_HMMLikelihood(x=data,first=rep(1,nrow(data)),m=sum(mv),T=T,dmat,gamma,delta)$lnl)
	if(debug)cat("\n -lnl=",neglnl)
	return(neglnl)
}

