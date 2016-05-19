#' Hidden Markov Model likelihood 
#' 
#' hmm.lnl which is a wrapper for the FORTRAN code hmm_like.f to compute negative log-likelihood for the hmm.
#'  
#' For an R version of the HMMLikelihood and related code see \code{\link{R_HMMLikelihood}}
#'  
#' @param x   matrix of observed sequences (row:id; column:occasion/time)
#' @param start the first occasion to start sequence for each x value 
#' @param m  number of states
#' @param T number of occasions; sequence length
#' @param dmat observation probability matrices
#' @param gamma transition matrices
#' @param delta initial distribution
#' @param freq vector of history frequencies or 1 
#' @param debug if TRUE, print out par values and -log-likelihood
#' @return negative log-likelihood value for each capture history.
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @seealso R_HMMLikelihood
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
hmm.lnl=function(x,start,m,T,dmat,gamma,delta,freq,debug)
{
	lnl=.Fortran("hmmlike",as.integer(x),as.integer(nrow(x)),as.integer(m),as.integer(T),
			as.integer(nrow(dmat[1,1,,])),as.integer(start),as.double(freq),as.double(dmat),
			as.double(gamma),as.double(delta),lnl=double(nrow(x)),PACKAGE="hsmm")$lnl
	if(any(is.nan(lnl) | is.infinite(lnl)))
	{
		if(debug)
		{
			cat("lnl is nan for these rows\n")
			cat(which(is.nan(lnl)))
			cat("lnl is infinite for these rows\n")
			cat(which(is.infinite(lnl)))
		}
		lnl[is.nan(lnl)]=-1e5
		lnl[is.infinite(lnl)]=-1e5
	}
	-sum(lnl)
}



