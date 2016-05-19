
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


plot_dt=function(object,range,dm,labels=NULL,...)
{
	if(is.null(labels))labels=paste("state",1:length(object$dtf))
	pars=split(object$results$par,rep(1:length(object$type),object$type))
	for (i in 1:length(object$dtf))
	{
		   plot(1:range,sapply(range, function(x) object$dtf[[i]]$fct(sum(pars[[i]]*dm[[i]][1,]),x)),xlab="Time",ylab=paste("Dwell time distribution:",labels[i]),type="l",...) 
		   if(nrow(dm[[i]])>1)
		   lines(1:range,sapply(range, function(x) object$dtf[[i]]$fct(sum(pars[[i]]*dm[[i]][2,]),x)),lty=2,pch=2,type="l",...) 
	}
}


# code that created data; requires database and CIPinnipedAnalysis and CalcurData packages
create.data=FALSE
if(create.data)
{
	library(CIPinnipedAnalysis)
	create_attendance=function(years,firstday=20,trim=TRUE)
	{
		#Firstday of June to 25 July
		resights=getCalcurData(db="Zc", tbl="Alive")
		df=NULL
		for(year in years)
		{
			xx=get_attendance(resights,year=year,firstday=firstday)
			x=t(apply(do.call("rbind",strsplit(xx$ch,"")),1,function(x) as.numeric(x)))
			attendance=x
			#only use those seen in last 7 days
			rownames(attendance)=xx$brand
			if(trim) attendance=attendance[which(!apply(attendance[,49:55],1,paste,collapse="")=="0000000"),]
			df=rbind(df,attendance)
		}
		return(df)
	}
	data=NULL
	firstday=25
	for(year in 1999:2002)
	{
		df=create_attendance(year,firstday,trim=FALSE)
		x=data.frame(ch=apply(df,1,paste,collapse=""),year=year,period="early",brand=rownames(df),stringsAsFactors=F)
		data=rbind(data,x)
	}
	for(year in 2012:2015)
	{
		df=create_attendance(year,firstday,trim=FALSE)
		x=data.frame(ch=apply(df,1,paste,collapse=""),year=year,period="late",brand=rownames(df),stringsAsFactors=F)
		data=rbind(data,x)
	}
	data$year=factor(data$year)
	data$period=factor(data$period)
}



library(hsmm)
library(optimx)

data(attendance)
# 4 state model - pre-birth, birth period, at sea, on land
omega=matrix(c(0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0),byrow=TRUE,ncol=4)
# strata.labels only used to create values in design data; values in encounter history are just 
# 0 (not seen) or 1 (seen) and not strata values;  encounter history is converted to 1 (not seen) and 2 (seen)
# process data
xp=process.data(attendance,model="ATTEND",groups=c("year","period"),strata.labels=c("P","B","S","L"))
# make design data and then delete fields that aren't used including fix which has default settings that
# are not useful for this model
ddl=make.design.data(xp)
# fixed real values
ddl$p$fix=NULL
ddl$dt$fix=NULL
ddl$dt$cohort=NULL
ddl$dt$Cohort=NULL
ddl$dt$age=NULL
ddl$dt$Age=NULL
ddl$p$cohort=NULL
ddl$p$Cohort=NULL
ddl$p$age=NULL
ddl$p$Age=NULL
# set p=0 for "S" at sea 
ddl$p$fix=ifelse(ddl$p$stratum%in%c("P","L","B"),NA,0)
# set p=0 for occasions in which not a single animal was seen; presumably no effort
for (year in c(1999:2002,2012:2015))
   ddl$p$fix[ddl$p$year==year][ddl$p$time[ddl$p$year==year]%in%which(colSums(xp$ehmat[xp$data$year==year,])==nrow(xp$ehmat[xp$data$year==year,]))]=0
# create covariate for birth and on-land p
ddl$p$birth=ifelse(ddl$p$stratum=="P",0,1)

# geometric for pre-birth and on-land but shifted poisson for birth and at-sea
state1=list(fct=shifted_poisson,steps=25,formula=~1)
state2=list(fct=shifted_poisson,steps=12,formula=~1)
state3=list(fct=shifted_poisson,steps=6,formula=~period)
state4=list(fct=shifted_poisson,steps=3,formula=~period)
dtf=list(state1,state2,state3,state4)

model3<-fit_attendance(xp$ehmat,ddl=ddl,dtf=list(state1,state2,state3,state4),pformula=~birth:year,omega=omega,debug=TRUE,initial=init)

# plot distributions
dev.new()
par(mfrow=c(2,2))
plot_dt(model2,15,dm=list(matrix(c(1,0),nrow=1),matrix(c(1,0),nrow=1),matrix(c(1,0,1,1),nrow=2,byrow=T),matrix(c(1,0,1,1),nrow=2,byrow=T)),labels=c("pre-birth","birth","at sea","on land"))
 
		
#std geometric model for all states
state1=list(fct=geometric,steps=1,formula=~1)
state2=list(fct=geometric,steps=1,formula=~1)
state3=list(fct=geometric,steps=1,formula=~1)
state4=list(fct=geometric,steps=1,formula=~1)
dtf=list(state1,state2,state3,state4)

