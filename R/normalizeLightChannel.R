#' @title normalizeLightChannel
#' 
#' @description Normalizes the light channel so that the channel sum is constant. 
#' The normalized light channel can be used for the NLI method. 
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
  
  # Ensure the object 'o' is of class 'pSILAC'
  if(class(o) != "pSILAC") stop("o should be a pSILAC object.")
  
  # Match the method argument to 'geometric' or 'linear'; use 'linear' by default if unspecified
  method <- match.arg(method, c("geometric","linear"))
  
  # Assign timepoints based on provided vector or default to unique time values from o$design
  if(is.null(timepoints)) timepoints <- unique(o$design$time)
  
  # Calculate the channel sum for the light channel (with optional NA removal)
  csu <- .getChannelSum(o, T)
  
  # Get row names (protein identifiers) of the channel sum
  p <- row.names(csu)
  
  # Filter out rows with NA values in either the heavy or light channels for the specified timepoints
  if(removeNAs) csu <- csu[which(apply(o$heavy[p, which(o$design$time %in% timepoints)], 1, 
                                       FUN = function(x) {!any(is.na(x))}) & 
                                   apply(o$light[p, which(o$design$time %in% timepoints)], 1, 
                                         FUN = function(x) {!any(is.na(x))})
  ), ]
  
  # Calculate size factors as the median of log-transformed channel sums, removing NAs as needed
  sizeFactors <- apply(log(csu), 2, na.rm = T, FUN = median)
  
  # Normalize based on the specified method
  if(method == "geometric"){
    message(paste(Sys.time(), "Normalizing the light channel using geometric method...", sep = " "))
    # Calculate normalization factors using geometric scaling
    normFactors <- 1 / (sizeFactors / mean(sizeFactors))
    normFactors <- normFactors / exp(mean(log(normFactors)))  # Scale to stabilize mean
    
    # Update the normalized channel sum (NCS) and normalized light intensity (NLI) in the pSILAC object
    o$NCS <- as.data.frame(exp(t(t(log(.getChannelSum(o, F))) * normFactors)))
    o$NLI <- as.data.frame(exp(t(t(log(o$light)) * normFactors)))
    
  } else {
    message(paste(Sys.time(), "Normalizing the light channel using linear method...", sep = " "))
    # Calculate normalization factors using linear scaling
    normFactors <- exp(max(sizeFactors) - sizeFactors)
    
    # Update the normalized channel sum (NCS) and normalized light intensity (NLI) in the pSILAC object
    o$NCS <- as.data.frame(t(t(.getChannelSum(o, F)) * normFactors))
    o$NLI <- as.data.frame(t(t(o$light) * normFactors))
  }
  
  # Return the modified pSILAC object with normalized fields
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