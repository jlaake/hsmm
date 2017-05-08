#' Fitting code for sea lion attendance data using hsmm 
#' 
#' Computes MLEs of parameters for dwell time distributions and sighting probabilities to data in df, based on
#' model defined by ddl (design data list), dtf (dwell time distributions), pformula (formula for sighting probabilities),
#' and omega (known state transition matrix). 
#' 
#' @param df encounter history matrix where each row is sequence for each occasion (column) 
#' @param ddl design data list of dwell time dataframe (dt) and sighting probability dataframe (p)
#' @param dtf list of length s where s is the number of states. The list elements specify the states in order
#' and contain a dwell time function (fct), a formula for the dwell time function parameter (at present only a
#' single parameter is allowed), and the number of time steps for the distribution (steps)
#' @param pformula formula to be applied to ddl$p to create the design matrix for the sighting probability model
#' @param omega an s by s matrix with known transition values between the states
#' @param initial vector of parameter initial values; must equal number of parameters if provided
#' @param method method to be used with optim; currently only a single method is supported
#' @param debug if TRUE will print out parameter values and negative log-likelihood value at current parameter values
#' @param nll if TRUE will return -log_likelihood function instead of optimizing 
#' @param mat if TRUE, return matrices rather than -log likelihood
#' for MLE infrence. Arguments 'method' and 'debug' will be ignored.
#' @author Jeff Laake
#' @import optimx
#' @return fitted list of model results
#' @export 
#' @keywords utility
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 82.
#' @examples
#' { 
#' data(attendance)
#' # use first 10 records from 1999 to keep execution time reasonable
#' attendance=attendance[attendance$year==1999,][1:10,]
#' # 4 state model - pre-birth, birth period, at sea, on land
#'omega=matrix(c(0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0),byrow=TRUE,ncol=4)
#' # strata.labels only used to create values in design data; values in encounter history are just 
#' # 0 (not seen) or 1 (seen) and not strata values;  encounter history is converted to
#' # 1 (not seen) and 2 (seen)
#' # process data
#' xp=process.data(attendance,model="ATTEND",strata.labels=c("P","B","S","L"))
#' # make design data and then delete fields that aren't used including fix 
#' # which has default settings that
#' # are not useful for this model
#' ddl=make.design.data(xp)
#' # fixed real values
#' ddl$p$fix=NULL
#' ddl$dt$fix=NULL
#' ddl$dt$cohort=NULL
#' ddl$dt$Cohort=NULL
#' ddl$dt$age=NULL
#' ddl$dt$Age=NULL
#' ddl$p$cohort=NULL
#' ddl$p$Cohort=NULL
#' ddl$p$age=NULL
#' ddl$p$Age=NULL
#' # set p=0 for "S" at sea 
#' ddl$p$fix=ifelse(ddl$p$stratum%in%c("P","L","B"),NA,0)
#' # set p=0 for occasions in which not a single animal was seen; presumably no effort
#' ddl$p$fix[ddl$p$time%in%which(colSums(xp$ehmat)==nrow(xp$ehmat))]=0
#' # create covariate for birth and on-land p
#' ddl$p$birth=ifelse(ddl$p$stratum=="P",0,1)
#' # specify shifted poisson distributions for each state with different numbers of steps
#' state1=list(fct=shifted_poisson,steps=25,formula=~1)
#' state2=list(fct=shifted_poisson,steps=12,formula=~1)
#' state3=list(fct=shifted_poisson,steps=6,formula=~1)
#' state4=list(fct=shifted_poisson,steps=3,formula=~1)
#' dtf=list(state1,state2,state3,state4)
#' system.time(model<-fit_attendance(xp$ehmat,ddl=ddl,dtf=dtf,pformula=~birth,
#'                 omega=omega,debug=TRUE,initial=c(2,2,1,0,-2,3)))
#' par(mfrow=c(2,2))
#' plot_dt(model,range=c(30,18,10,6),labels=c("pre-birth","birth","at sea","on land"))
#' }
fit_attendance=function(df,ddl,dtf,pformula,omega,initial=NULL,method="nlminb",debug=FALSE, nll=FALSE, mat=FALSE)
{
  # define functions fct_gamma,fct_dmat,fct_delta
  fct_dmat=function(pars,pformula,ddl,mv)
  {
    nobs= length(levels(ddl$p$id))
    nocc=length(levels(ddl$p$time))
    dm=model.matrix(pformula,ddl$p)
    p=plogis(dm%*%pars)
    p[!is.na(ddl$p$fix)]=ddl$p$fix[!is.na(ddl$p$fix)]
    xx=lapply(split(p,list(ddl$p$id,ddl$p$occ)),function(x) 
      dmat_hsmm2hmm(matrix(c(1-x[1],x[1],1-x[2],x[2],1-x[3],x[3],1-x[4],x[4]),nrow=2),mv=mv))
    zz=aperm(array(as.vector(do.call("cbind",xx)),dim=c(2,sum(mv),nobs,nocc)),c(3,4,1,2))
    return(zz)	
  }
  fct_gamma=function(pars,type,dtf,ddl,mv,omega)
  {
    nobs= length(levels(ddl$dt$id))
    nocc=length(levels(ddl$dt$time))
    dt=matrix(NA,nrow=nobs*nocc,ncol=length(dtf))
    for(i in 1:length(mv))
    {
      dm=model.matrix(dtf[[i]]$formula,ddl$dt[as.numeric(ddl$dt$stratum)==i,])
      dt[,i]=dm%*%pars[[i]]
    }
    xx=apply(dt,1, function(x) as.vector(hsmm:::gamma_dtd(split(x,1:length(dtf)),dtf,mv,omega)))
    xx=aperm(array(xx,dim=c(sum(mv),sum(mv),nobs,nocc)),c(3,4,1,2))
    return(xx)
  }
  
  fct_delta=function(pars,mv,nr)
    matrix(c(1,rep(0,sum(mv)-1)),byrow=T,ncol=sum(mv),nrow=nr)
  # compute number of parameters
  type=vector("numeric",length(dtf))
  cnames=NULL
  for (i in 1:length(dtf))
  {
    dm=model.matrix(dtf[[i]]$formula,ddl$dt[as.numeric(ddl$dt$stratum)==i,])
    if(!is.null(ddl$dt$fix))dm[!is.na(ddl$dt$fix[as.numeric(ddl$dt$stratum)==i]),]=0
    dm=dm[,colSums(dm)!=0,drop=FALSE]
    cnames=c(cnames,paste("state",i,":",colnames(dm),sep=""))
    type[i]=ncol(dm)
  }
  dm=model.matrix(pformula,ddl$p)
  if(!is.null(ddl$p$fix))dm[!is.na(ddl$p$fix),]=0
  dm=dm[,colSums(dm)!=0,drop=FALSE]
  cnames=c(cnames,paste("p:",colnames(dm),sep=""))
  type=c(type,ncol(dm))
  if(is.null(initial))
    initial=rep(0,sum(type))
  else
    if(length(initial)!=sum(type)) 
      stop(paste("initial parameter vector not correct length. Should be length",sum(type)))
  # call optim to fit model
  if(nll){
    foo = function(pars){
      hsmm:::hsmm_likelihood(pars,type=type,data=df,ddl=ddl,dtf=dtf,
                             fct_gamma=fct_gamma,fct_dmat=fct_dmat,fct_delta=fct_delta,
                             pformula=pformula,omega=omega,debug=debug,mat=mat)
    }
    return(foo)
  } else{
    fit=optimx(initial,type=type,hsmm_likelihood,data=df,ddl=ddl,
               dtf=dtf,omega=omega,pformula=pformula,fct_gamma=fct_gamma,fct_dmat=fct_dmat,fct_delta=fct_delta,
               method=method,debug=debug)
    par=as.matrix(fit[,1:sum(type)])[1,]
    names(par)=cnames
    vc=solve(attr(fit,"details")[,"nhatend"][[1]])
    rownames(vc)=cnames
    colnames(vc)=cnames
    results=list(data=df,ddl=ddl,dtf=dtf,pformula=pformula,type=type,omega=omega,fct_gamma=fct_gamma,fct_dmat=fct_dmat,fct_delta=fct_delta,
                 results=list(par=par,vc=vc,fit=fit))
    return(results)
  }
}


