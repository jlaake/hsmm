# script that analyzes and displays graphs used in ISEC 2016 presentation
library(hsmm)
# get attendance data 
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
ddl$dt$cohort=NULL
ddl$dt$Cohort=NULL
ddl$dt$age=NULL
ddl$dt$Age=NULL
ddl$p$cohort=NULL
ddl$p$Cohort=NULL
ddl$p$age=NULL
ddl$p$Age=NULL
# fixed real values
ddl$p$fix=NULL
ddl$dt$fix=NULL
# set p=0 for "S" at sea 
ddl$p$fix=ifelse(ddl$p$stratum%in%c("P","L","B"),NA,0)
# set p=0 for occasions in which not a single animal was seen; presumably no effort
for (year in c(1999:2002,2012:2015))
   ddl$p$fix[ddl$p$year==year][ddl$p$time[ddl$p$year==year]%in%which(colSums(xp$ehmat[xp$data$year==year,])==nrow(xp$ehmat[xp$data$year==year,]))]=0
# create covariate for birth and on-land p
ddl$p$birth=ifelse(ddl$p$stratum=="P",0,1)
# used shifted poisson functions for each state dwell time distribution
state1=list(fct=shifted_poisson,steps=25,formula=~1)
state2=list(fct=shifted_poisson,steps=12,formula=~1)
state3=list(fct=shifted_poisson,steps=6,formula=~period)
state4=list(fct=shifted_poisson,steps=3,formula=~period)
dtf=list(state1,state2,state3,state4)
#NOTE:  It takes several hours to fit the model 
init=c(  2.61279500,  1.99485100,  0.94510640,  0.25324590,  0.09476032,  0.22813880, -2.53685900,
       3.76270700,  3.23724300,  3.28566300,  2.87230900,  2.26122200 , 2.75874700,  2.69134900, 2.9666820)
model=fit_attendance(xp$ehmat,ddl=ddl,dtf=list(state1,state2,state3,state4),
		  pformula=~birth:year,omega=omega,debug=TRUE,initial=init)
# plot dwell time distributions
dev.new()
png("plots.png",pointsize=10,width=960,height=960)
par(mfrow=c(2,2),cex=1.5,cex.axis=1.5,cex.lab=1.5)
plot_dt(model,range=c(45,15,15,10),dm=list(matrix(c(1,0),nrow=1),matrix(c(1,0),nrow=1),matrix(c(1,0,1,1),nrow=2,byrow=T),matrix(c(1,0,1,1),nrow=2,byrow=T)),labels=c("pre-birth","birth","at sea","on land"))
dev.off()

# global decode to get state predictions
status=global_decode(model,ddl,state.names=c("P","B","S","L"))
# display histogram of birth timing
dev.new()
hist(apply(status,1,function(x) which(x=="B1")),xlab="Days from 25 May",ylab="Number of females giving birth",main="",cex.lab=1.5,cex.axis=1.5)
# display barplot of proportion on land in July 
dev.new()
barplot(apply(status[,37:61],2,function(x) {
					tab=table(substr(x,1,1))
					1-tab[names(tab)=="S"]/sum(tab)
				}),names=1:25, xlab="July date",ylab="Proportion on land",cex.axis=1.25,cex.lab=1.25)
