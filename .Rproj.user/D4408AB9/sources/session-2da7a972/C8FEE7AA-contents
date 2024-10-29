#' getNormValues
#'
#' Get a set of values with a specific (approximative) mean and SD.
#'
#' @param m desired mean
#' @param dsd desired standard deviation
#' @param nb number of values
#'
#' @return Returns a numeric vector of (deterministic) values with a specific mean and SD.
#'
#' @export
getNormValues <- function(m, dsd, nb){
    switch(paste("n",nb,sep=""), 
        n1=m,
        n2=c(m-(0.71*dsd),m+(0.71*dsd)),
        n3=c(m-dsd,m,m+dsd),
        n4=rep(c(m-(0.87*dsd),m+(0.87*dsd)),2),
        n5=c(rep(c(m-(1.19*dsd),m-(1*dsd),m,m+(1.19*dsd),m-(1*dsd)),1)),
        n6=rep(c(m-(1.155*dsd),m-(0.92*dsd),m+(1.155*dsd),m-(0.92*dsd)),2),
        n7=c(rep(c(m-(1.2*dsd),m-(dsd),m+(1.2*dsd),m-(dsd)),2),m),
        ifelse(nb != round(nb), c(m,getNormValues(m,dsd,nb-1)), rep(c(m-(0.9*dsd),m+(0.9*dsd)),nb/2))
        )
}

#' getGlobalVar
#'
#' Returns the global variance (removing differences between groups)
#'
#' @param y A numeric vector containing the values
#' @param x A numeric vector of same length as y, containing the groups to which the corresponding value in y belongs.
#'
#' @return Returns the group-adjusted variance.
#'
#' @export
getGlobalVar <- function(y,x){
	x <- x[which(!is.na(y))]
	y <- y[which(!is.na(y))]
	if(length(y)==0)	return(NA)
	return(sum((y-sapply(x,FUN=function(z){ mean(y[which(x==z)],na.rm=T) }))^2))
}

# imputes missing values as random values between 0 and the minimum value
imputeZero2Min <- function(x,nb=NULL){
	x <- as.numeric(x)
	if(is.null(nb))	nb <- sum(is.na(x))
	minVal <- min(x,na.rm=T)
	runif(nb,min(0,minVal),max(0,minVal))
}

# imputes missing values as draws from a normal distribution with the same global variance, 
# and with a mean equal to minValue-(SD*sd.offset)
imputeNorm <- function(y,x=NULL,sd.offset=1,trim.negative=F){
	y <- as.numeric(y)
	if(sum(is.na(y))==0)   return(y)
	if(is.null(x)){
		SD <- sd(y,na.rm=T)
	}else{
		SD <- sqrt(getGlobalVar(y,x))
	}
	if(is.na(SD))	return(y)
	y[which(is.na(y))] <- rnorm(sum(is.na(y)),min(y,na.rm=T)-(sd.offset*SD),SD)
	return(y)
}

# imputes missing values with values having exactly a mean and SD
imputeNormD <- function(y,x=NULL,sd.offset=1,trim.negative=F){
    y <- as.numeric(y)
    if(is.null(x)){
            SD <- sd(y,na.rm=T)
    }else{
            SD <- sqrt(getGlobalVar(y,x))
    }
    if(is.na(SD))	return(y)
    dmean <- min(y,na.rm=T)-(sd.offset*SD)
    if(is.null(x))	x <- rep(1,length(y))
    for(onegroup in unique(x)){
        if(any(is.na(y[x==onegroup]))){
            y[which(x==onegroup & is.na(y))] <- getNormValues(dmean, SD, sum(x==onegroup & is.na(y)))
        }
    }
    return(y)
}

#' imputeAll
#'
#' Impute all missing values in a matrix or data.frame.
#'
#' @param e The numeric matrix.
#' @param x The groups to which each column in e belongs.
#' @param method The imputation method; either: "norm", "normD", or "zero2min",
#' @param sd.offset The number of standard deviation to substract from the mean.
#' @param log.transform Logical; whether to log-transform for the purpose of calculating the normal distribution.
#'
#' @return Returns a numeric vector of (deterministic) values with a specific mean and SD.
#'
#' @export
imputeAll <- function(e, x=NULL, method="normD", sd.offset=1, log.transform=T){
    if(method=="normD1") method <- "normD"
    if(method=="normD2"){
        method <- "normD"
        sd.offset <- 2
    }
    if(is.null(nrow(e)))	return(e)
    if(nrow(e)==0 | !any(is.na(e)))	return(e)
    if(method=="zero2min"){
            e <- t(apply(e,1,FUN=function(z){ z[which(is.na(z))] <- imputeZero2Min(z); return(z) }))
    }else{
        if(log.transform)	e <- log2(e+0.01)
        if(method=="normD"){
                    e <- t(apply(e,1,FUN=function(z){ imputeNormD(z, x, sd.offset=sd.offset) }))
            }else{
                    e <- t(apply(e,1,FUN=function(z){ imputeNorm(z, x, sd.offset=sd.offset) }))
            }
        if(log.transform)	e <- 2^e - 0.01
    }
    return(e)
}

