
library(CIPinnipedAnalysis)
library(hsmm)
library(optimx)

fit_attendance=function(df,ddl,dtf,pformula,omega,initial=NULL,method="nlminb",debug=FALSE)
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
		dt=matrix(NA,nrow=nobs*nocc,ncol=sum(type[1:length(mv)]))
		for(i in 1:length(mv))
		{
			dm=model.matrix(dtf[[i]]$formula,ddl$dt[as.numeric(ddl$dt$stratum)==i,])
			dt[,i]=dm%*%pars[[i]]
		}
		xx=apply(dt,1, function(x) as.vector(hsmm:::gamma_dtd(split(x,rep(1:length(type),type)),dtf,mv,omega)))
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
  fit=optimx(initial,type=type,hsmm_likelihood,data=df,ddl=ddl,
					dtf=dtf,omega=omega,pformula=pformula,fct_gamma=fct_gamma,fct_dmat=fct_dmat,fct_delta=fct_delta,method=method,debug=debug)
  par=as.matrix(fit[,1:sum(type)])[1,]
  names(par)=cnames
  vc=solve(attr(fit,"details")[,"nhatend"][[1]])
  rownames(vc)=cnames
  colnames(vc)=cnames
  results=list(df=df,ddl=ddl,dtf=dtf,pformula=pformula,omega=omega,results=list(par=par,vc=vc,fit=fit))
  return(results)
}

plot_dt=function(object,range, labels=NULL,...)
{
	if(is.null(labels))labels=paste("state",1:length(object$dtf))
	for (i in 1:length(object$dtf))
	{
		plot(1:range,sapply(range, function(x) object$dtf[[i]]$fct(object$results$par[i],x)),xlab="Time",ylab=paste("Dwell time distribution:",labels[i]),type="b",...) 
	}
}

create_attendance=function(years,firstday,trim=TRUE)
{
	#Firstday of June to 25 July
	resights=getCalcurData(db="Zc", tbl="Alive")
	df=NULL
	for(year in years)
	{
		xx=get_attendance(resights,year=year,firstday=1)
		x=t(apply(do.call("rbind",strsplit(xx$ch,"")),1,function(x) as.numeric(x)))
		attendance=x
		#only use those seen in last 7 days
		if(trim) attendance=attendance[which(!apply(attendance[,49:55],1,paste,collapse="")=="0000000"),]
		df=rbind(df,attendance)
	}
	return(df)
}

# 4 state model - pre-birth, birth period, at sea, on land
omega=matrix(c(0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0),byrow=TRUE,ncol=4)

year=2014
df=create_attendance(year,1,trim=FALSE)
x=data.frame(ch=apply(df,1,paste,collapse=""),stringsAsFactors=F)
xp=process.data(x,model="ATTEND",strata.labels=c("P","B","S","L"))
ddl=make.design.data(xp)
# fixed real values
ddl$p$fix=NULL
ddl$dt$fix=NULL
ddl$p$fix=ifelse(ddl$p$stratum%in%c("P","L","B"),NA,0)
ddl$p$fix[ddl$p$time%in%which(colSums(xp$ehmat)==nrow(xp$ehmat))]=0
# create covariate for birth and on-land p
ddl$p$birth=ifelse(ddl$p$stratum=="P",0,1)

# geometric for pre-birth and on-land but shifted poisson for birth and at-sea
state1=list(fct=geometric,steps=1,formula=~1)
state2=list(fct=shifted_poisson,steps=8,formula=~1)
state3=list(fct=shifted_poisson,steps=6,formula=~1)
state4=list(fct=geometric,steps=1,formula=~1)
dtf=list(state1,state2,state3,state4)

model1=fit_attendance(xp$ehmat,ddl=ddl,dtf=list(state1,state2,state3,state4),pformula=~birth,omega=omega,debug=TRUE,initial=c(-2,2,2,-1,-2,2.5))
dev.new()
par(mfrow=c(2,2))
plot_dt(model1,15,labels=c("pre-birth","birth","at sea","on land"))

#std geometric model for all states
state1=list(fct=geometric,steps=1,formula=~1)
state2=list(fct=geometric,steps=1,formula=~1)
state3=list(fct=geometric,steps=1,formula=~1)
state4=list(fct=geometric,steps=1,formula=~1)
dtf=list(state1,state2,state3,state4)

model0=fit_attendance(xp$ehmat,ddl=ddl,dtf=list(state1,state2,state3,state4),pformula=~birth,omega=omega,debug=TRUE,initial=model1.2012$results$par)
par(mfrow=c(2,2))
plot_dt(model0,15,labels=c("pre-birth","birth","at sea","on land"))

model0.1999=model0
model1.1999=model1

model0.2012=model0
model1.2012=model1

model0.2014=model0
model1.2014=model1


