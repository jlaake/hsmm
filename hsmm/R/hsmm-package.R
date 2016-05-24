#' Attendance data
#' 
#' Observations of 54 female sea lions with pups over 51 consecutive days. 
#' 
#' @name attendance
#' @docType data
#' @format The data are a matrix with 54 rows and 51 columns containing a 2 if the female was seen at the rookery and 1 if not seen.
#' @keywords datasets
#' @examples 
#' 
#'\dontrun{
#'# code that created data; requires database and CIPinnipedAnalysis and CalcurData packages
#'create.data=FALSE
#'if(create.data)
#'{
#'	library(CIPinnipedAnalysis)
#'	create_attendance=function(years,firstday=20,trim=TRUE)
#'	{
#'		#Firstday of June to 25 July
#'		resights=getCalcurData(db="Zc", tbl="Alive")
#'		df=NULL
#'		for(year in years)
#'		{
#'			xx=get_attendance(resights,year=year,firstday=firstday)
#'			x=t(apply(do.call("rbind",strsplit(xx$ch,"")),1,function(x) as.numeric(x)))
#'			attendance=x
#'			#only use those seen in last 7 days
#'			rownames(attendance)=xx$brand
#'			if(trim) attendance=attendance[which(!apply(attendance[,49:55],1,paste,collapse="")=="0000000"),]
#'			df=rbind(df,attendance)
#'		}
#'		return(df)
#'	}
#'	data=NULL
#'	firstday=25
#'	for(year in 1999:2002)
#'	{
#'		df=create_attendance(year,firstday,trim=FALSE)
#'		x=data.frame(ch=apply(df,1,paste,collapse=""),year=year,period="early",brand=rownames(df),stringsAsFactors=F)
#'		data=rbind(data,x)
#'	}
#'	for(year in 2012:2015)
#'	{
#'		df=create_attendance(year,firstday,trim=FALSE)
#'		x=data.frame(ch=apply(df,1,paste,collapse=""),year=year,period="late",brand=rownames(df),stringsAsFactors=F)
#'		data=rbind(data,x)
#'	}
#'	data$year=factor(data$year)
#'	data$period=factor(data$period)
#'}
#'}


NULL
