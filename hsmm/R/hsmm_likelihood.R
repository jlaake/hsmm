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
#' @return negative log-likelihood value
#' @author Jeff Laake
#' @export hsmm_likelihood
#' @keywords utility
#' @references Zucchini, W., MacDonald, I.L. and Langrock, R. 2016. Hidden Markov Models for Time Series: 
#' An introduction using R, 2nd ed. CRC Press.
#' @examples
#' {
#' # get female sea lion attendance data
#' data(attendance)
#' # detection can occur when sea lion is on land (state 1) but is 0 when sea lion at sea (state 2). Observations
#' # are not seen (1) and seen (2). 
#' fct_dmat=function(pars)
#'   matrix(c(1-plogis(pars),plogis(pars),1,0),ncol=2)
#' # assume intial state has equal probability across all expanded states
#' fct_delta=function(pars,mv,nr)
#'   matrix(rep(1/sum(mv),sum(mv)*nr),byrow=T,ncol=sum(mv))  
#' # call optim to fit unstructured order 1 (q=2) and geometric tail for time on land and shifted_poisson for time at sea; first 2 parameters 
#' # are for unstructured_gt, third is log of poisson rate, and fourth is logit of detection probability - probability of seeing sea lion given it is on land. 
#'  optim(c(0,0,log(3),0),type=c(2,1,1),hsmm_likelihood,data=attendance,dtf=list(unstructured_gt,shifted_poisson),mv=c(2,6),fct_dmat=fct_dmat,fct_delta=fct_delta)
#' 
#'# Four state model - pre-birth, birth, at sea, on land. Detection can occur in states 1,2 and 4 but is 0 in state 3 (at sea). Observations
#'# are not seen (1) and seen (2). 
#' library(optimx)
#'# observation matrix has 3 detection probabilities; one for pre-birth, one for birth period and one for on land (nursing).  Turns out last 2 can probably
#'# be equal.
#'fct_dmat=function(pars)
#'  matrix(c(1-plogis(pars[1]),plogis(pars[1]),1-plogis(pars[2]),plogis(pars[2]),1,0, 1-plogis(pars[3]),plogis(pars[3])),nrow=2)
#'# assume intial state is pre-birth
#'fct_delta=function(pars,mv,nr)
#'  matrix(c(1,rep(0,sum(mv)-1)),byrow=T,ncol=sum(mv),nrow=nr)
#'# 4 state model - pre-birth, birth period, at sea, on land
#'omega=matrix(c(0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0),byrow=TRUE,ncol=4)
#'# call optim to fit geometric (state 1), shifted_poisson (state 2), shifted_poisson for state 3 time at sea; and unstructured_gt (q=2) for on land 
#'# are for unstructured_gt, third is log of poisson rate, and fourth is logit of detection probability - probability of seeing sea lion given it is on land. 
#'# from initial run that ran out of iterations
#'initial=c( -1.61158059 , 2.11453537 , 1.15941438 ,-0.83525677, -0.22348341 ,-2.08088003  ,0.06749222 ,-0.51567059)
#'optimx(initial,type=c(1,1,1,2,3),hsmm_likelihood,data=attendance,
#' dtf=list(geometric,shifted_poisson,shifted_poisson,unstructured_gt),omega=omega,
#' mv=c(1,8,6,2),fct_dmat=fct_dmat,fct_delta=fct_delta)

#'  }
hsmm_likelihood=function(pars,type,data,ddl,dtf,fct_gamma,fct_dmat,fct_delta,pformula,omega,debug=FALSE)
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
	if(debug) sapply(pars,function(x) cat("\n",x))   
	neglnl=-sum(R_HMMLikelihood(x=data,first=rep(1,nrow(data)),m=sum(mv),T=T,dmat,gamma,delta)$lnl)
	if(debug)cat("\n -lnl=",neglnl)
	return(neglnl)
}

