#' normalizeLightChannel
#' 
#' Normalizes the light channel so that the channel sum is constant. Will be used for the NLI method. 
#'
#' @param o a pSILAC object.
#' @param method normalization method, either 'linear' (default) or 'geometric'.
#' @param timepoints an optional vector indicating the timepoint for each sample. If NULL (default), falls back to the value of the `time` column of the design.
#' @param removeNAs logical; whether to remove NA values before normalizing (default TRUE)
#'
#' @return data.frame.
#'
#' @export
normalizeLightChannel <- function(o, method="linear", timepoints=NULL, removeNAs=T){
  if(class(o) != "pSILAC")	stop("o should be a pSILAC object.")
  method <- match.arg(method, c("geometric","linear"))
  if(is.null(timepoints)) timepoints <- unique(o$design$time)
  csu <- .getChannelSum(o,T)
  p <- row.names(csu)
  if(removeNAs) csu <- csu[which( apply(o$heavy[p,which(o$design$time %in% timepoints)],1,FUN=function(x){!any(is.na(x))}) & 
                                    apply(o$light[p,which(o$design$time %in% timepoints)],1,FUN=function(x){!any(is.na(x))})
  ),]
  sizeFactors <- apply(log(csu),2,na.rm=T,FUN=median)
  if(method=="geometric"){
    normFactors <- 1/(sizeFactors/mean(sizeFactors))
    normFactors <- normFactors/exp(mean(log(normFactors)))
    o$NCS <- as.data.frame(exp(t(t(log(.getChannelSum(o,F)))*normFactors)))
    o$NLI <- as.data.frame(exp(t(t(log(o$light))*normFactors)))
  }else{
    normFactors <- exp(max(sizeFactors)-sizeFactors)
    o$NCS <- as.data.frame(t(t(.getChannelSum(o,F))*normFactors))
    o$NLI <- as.data.frame(t(t(o$light)*normFactors))
  }
  return(o)
}

.getChannelSum <- function(o, inBothChannels=T){
  p <- row.names(o$light)
  if(inBothChannels) p <- p[which(p %in% row.names(o$heavy))]
  csu <- o$light[p,]
  tmp <- as.data.frame(o$heavy)[p,]
  for(i in 1:ncol(o$light)){ csu[,i] <- apply(cbind(csu[,i],tmp[,i]),1,na.rm=T,FUN=sum) }
  return(csu)
}