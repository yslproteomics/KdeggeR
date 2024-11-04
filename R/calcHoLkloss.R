#' Calculate ln(H/L+1)-based k_loss for all peptides
#'
#' A wrapper to apply getHoLmod to all peptides of the dataset.
#'
#' @param e a pSILAC object.
#'
#' @return The updated pSILAC object.
#'
#' @export
calcHoLkloss <- function(e, ncores=NULL){
  if(class(e) != "pSILAC")	stop("e should be a pSILAC object.")
  e$hol.kloss <- NULL
  x <- unique(e$design$time)
  if(is.null(ncores)){
    library(parallel)
    ncores <- detectCores() - 1
  }else{
    if(ncores>1)	library(parallel)
  }
  if(ncores>1){
    library(parallel)
    cl <- makeCluster(ncores)
    clusterExport(cl, c("getHoLmod"))
    clusterExport(cl, c("x","e"), environment())
  }
  for(p in unique(e$design$sample)){
    if(ncores > 1){
      d <- as.data.frame(t(parApply(cl,e$hol[,which(e$design$sample==p)],1,FUN=function(y){ getHoLmod(x,y) })))
    }else{
      d <- as.data.frame(t(apply(e$hol[,which(e$design$sample==p)],1,FUN=function(y){ getHoLmod(x,y) })))
    }
    colnames(d) <- paste(p,c("kloss","kloss.stderr","kloss.SSR","nbpoints"),sep=".")
    row.names(d) <- row.names(e$hol)
    message(paste(" ...calculated ",sum(!is.na(d[,1]))," k_loss values for sample ",p," (",sum(is.na(d[,1]))," missing)",sep=""))
    if(is.null(e$hol.kloss)){
      e$hol.kloss <- d
    }else{
      e$hol.kloss <- cbind(e$hol.kloss,d)
    }
  }
  if(ncores > 1) stopCluster(cl)
  return(e)
}
