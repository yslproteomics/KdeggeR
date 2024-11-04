#' @title Calculate k_loss for All Peptides Based on Normalized Light Intensities
#'
#' @description This function applies the `getNLImod` function to all peptides in a pSILAC dataset to calculate the decay rate 
#' constant \( k_{loss} \) for normalized light intensities. It uses a specified starting intensity as the intersect for fitting 
#' the model, and supports parallel computation to speed up processing.
#'
#' @param o A `pSILAC` object containing the data to be analyzed.
#' @param startIntensity Character. Specifies the method for determining the starting intensity for model fitting:
#' \itemize{
#'   \item \code{"model"} - Fits the model with an intersect.
#'   \item \code{"max"} - Uses the maximum H+L intensity as the intersect (default).
#'   \item \code{"median"} - Uses the median of non-null H+L intensities as the intersect.
#' }
#' @param ncores Integer. The number of cores to use for parallel computation. If \code{NULL} (default), the function 
#' detects available cores and uses all but one.
#'
#' @return The input `pSILAC` object with an added `NLI.kloss` element, containing the calculated \( k_{loss} \), 
#' standard error, sum of squared residuals, and the number of points used in the fitting for each sample.
#'
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Uses parallel processing, if multiple cores are available, to apply `getNLImod` on each peptide across different samples.
#'   \item Calculates \( k_{loss} \), the standard error, and other fit metrics for each peptide, based on the specified starting intensity method.
#' }
#'
#' @export
calcNLIkloss <- function(o, startIntensity="max", ncores=NULL){
  
  if(class(o) != "pSILAC")	stop("o should be a pSILAC object.")
  
  # if(!all(c("NLI","NCS") %in% names(o))){
  #   message("Light channel not yet normalized; performing linear normalization...")
  #   o <- normalizeLightChannel(o, method="linear")
  #   message("Done; now fitting model on each peptide...")
  # }
  
  startIntensity <- match.arg(startIntensity, c("median","max","model"))
  
  freeIntersect <- startIntensity == "model"
  
  o$NLI.kloss <- NULL
  
  x <- unique(o$design$time)
  
  if(is.null(ncores)){
    library(parallel)
    ncores <- detectCores() - 1
  }else{
    if(ncores>1)	library(parallel)
  }
  if(ncores>1){
    library(parallel)
    cl <- makeCluster(ncores)
    clusterExport(cl, c("getNLImod",".medianNonNull"))
    clusterExport(cl, c("x","o","freeIntersect"), environment())
  }
  for(p in unique(o$design$sample)){
    if(ncores > 1){
      d <- as.data.frame(t(parApply(cl, cbind(
        suppressWarnings(apply(o$NCS[, which(o$design$sample == p)], 1, na.rm = TRUE, 
                               FUN = ifelse(startIntensity == "median", ".medianNonNull", startIntensity))),
        o$NLI[, which(o$design$sample == p)]
      ), 1, FUN = function(y) { 
        getNLImod(x, y[-1], ifelse(freeIntersect, NULL, y[1])) 
      })))
    } else {
      d <- as.data.frame(t(apply(cbind(
        suppressWarnings(apply(o$NCS[, which(o$design$sample == p)], 1, na.rm = TRUE, 
                               FUN = ifelse(startIntensity == "median", ".medianNonNull", startIntensity))),
        o$NLI[, which(o$design$sample == p)]
      ), 1, freeIntersect = freeIntersect, FUN = function(y, freeIntersect) { 
        getNLImod(x, y[-1], ifelse(freeIntersect, NULL, y[1])) 
      })))
    }
    colnames(d) <- paste(p, c("kloss","kloss.stderr","kloss.SSR","nbpoints"), sep=".")
    row.names(d) <- row.names(o$NLI)
    message(paste(" ...calculated ", sum(!is.na(d[,1])), " k_loss values for sample ", p, 
                  " (", sum(is.na(d[,1])), " missing)", sep=""))
    
    if(is.null(o$NLI.kloss)){
      o$NLI.kloss <- d
    }else{
      o$NLI.kloss <- cbind(o$NLI.kloss,d)
    }
  }
  if(ncores > 1) stopCluster(cl)
  return(o)
}
