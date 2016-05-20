

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

plot_dt(model2,range=c(45,15,15,10),dm=list(matrix(c(1,0),nrow=1),matrix(c(1,0),nrow=1),matrix(c(1,0,1,1),nrow=2,byrow=T),matrix(c(1,0,1,1),nrow=2,byrow=T)),labels=c("pre-birth","birth","at sea","on land"))


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


