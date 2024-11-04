#' Calculate RIA-based k_loss
#'
#' Models the rate of loss of the light isotope (k_loss) on the basis of the relative isotope abundance (RIA).
#'
#' @param x The timepoints.
#' @param y The RIA values corresponding to the timepoints.
#' @param a The starting value of the k_loss for model-fitting.
#'
#' @return A numeric vector of length 4 containing: 
#'	    the k_loss,
#'	    the standard error of the k_loss, 
#'	    the sum of residuals of the linear model, 
#'	    the number of points on which the model was built.
#'
#' @export
getRIAmod <- function(x,y,a=0.1){
  if(length(x) != length(y))	stop("x and y are of different lengths!")
  f <- function(x,a){ exp(-a*x) }
  y <- as.numeric(y)
  x2 <- x[which(!is.na(y))]
  y <- y[which(!is.na(y))]
  if(length(y)==1){
    if(x2==0)	return(c(NA,NA,NA,1))
    return( c(-log(y)/x2,NA,NA,1) )
  }
  fm <- try(suppressWarnings(nls(y~f(x2,a), start = c(a=a))), silent=T)
  if(class(fm)=="nls"){
    return(c(summary(fm)$coefficients[1:2],sum(residuals(fm)^2),length(x2)))
  }
  return(c(NA,NA,NA,length(x2)))
}

#' Calculate NLI-based k_loss
#'
#' Models the rate of loss of the light isotope (k_loss) on the basis of the normalized light intensities
#'
#' @param x The timepoints.
#' @param y The normalized intensities values corresponding to the timepoints.
#' @param y0 The intersect of the model. If NULL (default), will attempt to fit it.
#' @param a The starting value of the k_loss for model-fitting.
#'
#' @return A numeric vector of length 4 containing: 
#'	    the k_loss,
#'	    the standard error of the k_loss, 
#'	    the sum of residuals of the linear model, 
#'	    the number of points on which the model was built.
#'
#' @export
getNLImod <- function(x,y,y0=NULL,a=0.1){
  if(length(x) != length(y))	stop("x and y are of different lengths!")
  #f <- function(x,a,b){ b*(1-a)^x}
  f <- function(x,a,b){ b*exp(-a*x) }
  y <- as.numeric(y)
  x2 <- x[which(!is.na(y))]
  y <- y[which(!is.na(y))]
  if(length(y)<1) return(c(NA,NA,NA,1))
  if(length(y)<2){
    if(is.null(y0) | x2==0 | y0==y) return(c(NA,NA,NA,1))
    return( c(-log(y/y0)/x2,NA,NA,1) )
  }
  if(is.null(y0)){
    fm <- try(suppressWarnings(nls(y~f(x2,a,b), start = c(a=a,b=max(y)))), silent=T)
  }else{
    fm <- try(suppressWarnings(nls(y~f(x2,a,y0), start = c(a=a))), silent=T)
  }
  if(class(fm)=="nls"){
    return(c(summary(fm)$coefficients[1:2],sum(residuals(fm)^2),length(x2)))
  }
  return(c(NA,NA,NA,length(x2)))
}

#' Calculate H/L-based k_loss
#'
#' Models the rate of loss of the light isotope (k_loss) on the basis of the slope of the log-transformed isotope ratio.
#'
#' @param x The timepoints.
#' @param y The ln(H/L+1) values corresponding to the timepoints.
#' @param tryRobust Whether to try to fit a robust linear model (requires the 'MASS' library). Ignored with less than 10 data points.
#'
#' @return A numeric vector of length 4 containing: 
#'	    the k_loss, i.e. ln(H/L+1)slope, 
#'	    the standard error of the slope, 
#'	    the sum of residuals of the linear model, 
#'	    the number of points on which the model was built.
#'
#' @export
getHoLmod_old <- function(x,y,tryRobust=F){
  if(length(x) != length(y))	stop("x and y are of different lengths!")
  y <- as.numeric(y)
  x2 <- x[which(!is.na(y))]
  y <- y[which(!is.na(y))]
  if(length(y)==1){
    if(x2==0)	return(c(NA,NA,NA,1))
    return(c(y/x2,NA,NA,1))
  }
  if(tryRobust & length(y) >= 10){
    fm <- try(rlm(y~x2+0), silent=T)
    if(!"lm" %in% class(fm))	fm <- try(lm(y~x2+0), silent=T)
  }else{
    fm <- try(lm(y~x2+0), silent=T)
  }
  if("lm" %in% class(fm)){
    return(c(as.numeric(summary(fm)$coefficients[1:2]),sum(residuals(fm)^2),length(x2)))
  }
  return(c(NA,NA,NA,length(x2)))
}

