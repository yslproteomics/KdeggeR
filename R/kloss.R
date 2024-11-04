#' Calculate the Median of Non-Null, Non-Zero Values
#'
#' This function computes the median of a numeric vector, ignoring values 
#' that are less than or equal to zero. It returns \code{NA} if no values
#' greater than zero are present in the vector.
#'
#' @param x A numeric vector from which the median of non-zero values is to be calculated.
#' @param na.rm A logical value indicating whether \code{NA} values should be stripped before 
#'        the computation proceeds. Default is \code{TRUE}.
#' 
#' @return The median of the non-zero values in \code{x}. If no values greater than zero are
#'         found, the function returns \code{NA}.
#'
#' @details This function first removes all values from the input vector \code{x} that are 
#'          less than or equal to zero. It then calculates the median of the remaining 
#'          values. If the resulting vector is empty, \code{NA} is returned.
#'
#' @export
.medianNonNull <- function(x, na.rm = TRUE) {
  # Filter out values that are less than or equal to zero
  x <- x[which(x > 0)]
  
  # If there are remaining values, return the median
  if (length(x) > 0) {
    return(median(x, na.rm = na.rm))
  }
  
  # If no non-zero values are present, return NA
  return(NA)
}

#' Calculate all peptide and protein kloss.
#'
#' A wrapper that calls, in order, calcRIAkloss(), getHoLkloss() and calcProteinsKloss().
#'
#' @param e a pSILAC object.
#' @param method the method used to calculate protein-level rates (default 'combined'). See ?calcProteinsKloss
#' @param ag.metric the aggregation metric used for protein-level rates (default 'mean'). See ?calcProteinsKloss
#' @param ag.weights the method to calculate weights for aggregation (default 'variance'). Ignored if ag.metric != 'mean'. See ?calcProteinsKloss
#' @param in.all whether to use only peptides quantified in all samples. See ?calcProteinsKloss
#' @param ncores number of cores to use
#'
#' @return The updated pSILAC object.
#'
#' @export
calcAllRates <- function(e, method="combined", ag.metric="mean", ag.weights="both", in.all=2, ncores = 1){
  if(class(e) != "pSILAC")	stop("e should be a pSILAC object.")
  message("Calculating peptide-wise k_loss using the RIA-based method...")
  e <- calcRIAkloss(e)
  message("Calculating peptide-wise k_loss using the ln(H/L+1)-based method...")
  e <- calcHoLkloss(e)
  message("Calculating peptide-wise k_loss using the normalized light channel...")
  e <- calcNLIkloss(e)
  message(paste("Calculating protein-wise k_loss using..."))
  e <- calcProteinsKloss(e, method, ag.metric, ag.weights, in.all)
  e$info$protk.call[["in.all"]] <- in.all
  return(e)
}

#' Calculate all protein kloss.
#'
#' Aggregates peptide-level kloss for all proteins. This requires that calcRIAkloss, calcHoLkloss, and/or calcNLIkloss be run first.
#'
#' @param x a pSILAC object.
#' @param method the method used to calculate protein-level rates (default 'combined'); either 'RIA', 'hol', 'NLI', or 'combined' (uses both RIA- and NLI-based peptide kloss).
#' @param ag.metric the aggregation metric used for protein-level rates (default 'mean'). Either 'mean', 'median' or 'remodel'. If 'mean' is chosen, and if ag.weights is not null, a weighted mean will be performed.
#' @param ag.weights the method to calculate weights for weighted mean (default 'variance'). Ignored if ag.metric != 'mean'. Either 'variance' (1/variance of the model's fit), 'nbpoints' (number of datapoints used for the fit), or 'both' (nbpoints/variance).
#' @param unique.weights whether unique weights should be used for each peptides (default T), or if false a weight is calculated individually for each sample.
#' @param in.all Starting with how many peptides should peptides not quantified in all samples be ignored (default 2), or FALSE to use all peptides, even if not quantified in all samples. 
#' @param removeOutliers Starting with how many peptides should outlier peptides be removed, or FALSE to disable outlier removal (default 5). See ?removeOutlierPeptides
#' @param ncores Number of cores to use (only for remodeling). Defaults to all cores minus 1.
#' @param tryRobust Whether to try fitting a robust linear model when remodeling on the basis of H/L ratio (requires the 'MASS' library). Disabled by default.
#' @param returnKlossTableOnly Logical; whether to return only the kloss table, rather than the whole pSILAC object (default).
#' @param returnSD Logical; whether to return also the SD of the kloss (default FALSE).
#'
#' @return The updated pSILAC object, or the protein kloss table if returnKlossTableOnly is TRUE.
#'
#' @export
calcProteinsKloss <- function(x, method="combined", ag.metric="mean", ag.weights="both", unique.weights=T, in.all=2, removeOutliers=5, ncores=NULL, tryRobust=F, returnKlossTableOnly=F, returnSD=F){
  if(class(x) != "pSILAC")	stop("x should be a pSILAC object.")
  if(is.null(removeOutliers))	removeOutliers <- F
  if( (is.logical(removeOutliers) & removeOutliers) |
      (is.numeric(removeOutliers) & !(removeOutliers >= 3)) ){
	stop("removeOutliers should be either FALSE, or an integer >= 3 indicating the minimum number of peptides for outlier removal.")
  }
  method <- match.arg(method, c("combined","complement","RIA","hol","NLI"))
  if(ag.metric=="remodel" & method != "complement")	return(modelProteinsKloss(x, method=method, ag.weights=ag.weights, in.all=in.all, removeOutliers=removeOutliers, ncores=ncores, tryRobust=tryRobust,returnKlossTableOnly=returnKlossTableOnly, returnSD=returnSD))
  ag.metric <- match.arg(ag.metric, c("remodel","median","mean"))
  if(method %in% c("combined","complement")){
    if(is.null(x$RIA.kloss) | is.null(x$NLI.kloss) | (method=="combined2" & is.null(x$hol.kloss)))	stop("Peptide-level k_loss should first be calculated.")
    k1 <- calcProteinsKloss(x, method="RIA", ag.metric=ag.metric, ag.weights=ag.weights, unique.weights=unique.weights, in.all=in.all, removeOutliers=removeOutliers, ncores=ncores, tryRobust=tryRobust, returnKlossTableOnly=T, returnSD=method=="combined")
    k2 <- calcProteinsKloss(x, method="NLI", ag.metric=ag.metric, ag.weights=ag.weights, unique.weights=unique.weights, in.all=in.all, removeOutliers=removeOutliers, ncores=ncores, tryRobust=tryRobust, returnKlossTableOnly=T, returnSD=method=="combined")
    if(method=="complement"){
	message("Adding complementary NLI-based estimates to the RIA-based rates of loss")
	k1 <- k1[which(!apply(k1,1,FUN=function(k){ any(is.na(k) | k < 0) })),]
	k2 <- k2[which(!(row.names(k2) %in% row.names(k1))),]
	x$protein.kloss <- as.data.frame(rbind(k1,k2))
	x$protein.kloss$source <- rep(c("RIA","NLI"),c(nrow(k1),nrow(k2)))
    }else{
	message("Combining RIA-based and NLI-based estimates of rates of loss")
	k1 <- k1[which(!apply(k1[,2*1:(ncol(k1)/2)-1],1,FUN=function(k){ any(is.na(k) | k < 0) })),]
	k2comp <- k2[which(!(row.names(k2) %in% row.names(k1))),2*1:(ncol(k2)/2)-1]
	k2 <- k2[row.names(k1),]
	msd1 <- apply(k1[,2*1:(ncol(k1)/2)],1,FUN=function(x){ abs(mean(x[2*1:(length(x)/2)],na.rm=T)/median(x[2*1:(length(x)/2)-1],na.rm=T)) })
	msd2 <- apply(k2[,2*1:(ncol(k2)/2)],1,FUN=function(x){ abs(mean(x[2*1:(length(x)/2)],na.rm=T)/median(x[2*1:(length(x)/2)-1],na.rm=T)) })
	k1 <- k1[,2*1:(ncol(k1)/2)-1]
	k2 <- k2[,2*1:(ncol(k2)/2)-1]
	k1[which(msd1>msd2),] <- k2[which(msd1>msd2),]
	x$protein.kloss <- as.data.frame(rbind(k1,k2comp))
	x$protein.kloss$source <- rep(c("RIA","NLI"),c(nrow(k1),nrow(k2comp)))
	x$protein.kloss$source[which(msd1>msd2)] <- "NLI"
    }
    x$info$protk.method <- paste("Protein-level k_loss calculated throught the",method,"method with",ag.metric,"of peptides",ifelse(ag.metric=="mean" & !is.null(ag.weights),paste("weighted by",ifelse(ag.weights=="both","nbpoints/variance",ag.weights)),""))
    x$info$protk.call <- match.call()
    if(returnKlossTableOnly) return(x$protein.kloss)
    return(x)
  }else{
    e <- x[[paste(method,"kloss",sep=".")]]
    if(is.null(e))	stop("Peptide-level k_loss should first be calculated using the desired method.")
    proteins <- as.character(x$peptides[row.names(e),"protein"])
  }
  if(returnSD){
	x$protein.kloss <- as.data.frame(matrix(NA,nrow=length(unique(proteins)),ncol=2*length(unique(x$design$sample))), row.names=unique(proteins))
	colnames(x$protein.kloss) <- paste(rep(unique(x$design$sample),each=2),rep(c("kloss","stderr"),length(unique(x$design$sample))),sep=".")	
  }else{
	x$protein.kloss <- as.data.frame(matrix(NA,nrow=length(unique(proteins)),ncol=length(unique(x$design$sample))), row.names=unique(proteins))
	colnames(x$protein.kloss) <- unique(x$design$sample)
  }
  i <- 0;
  for(p in unique(proteins)){
    i <- i+1;
    if(i>100){
      cat(".")
      i <- 0
    }
    x$protein.kloss[p,] <- aggregateKloss(e[which(proteins==p),], unique(x$design$sample), ag.metric, ag.weights, unique.weights, in.all, removeOutliers, returnSD=returnSD)
  }
  cat("\n")
  if(returnKlossTableOnly)	return(x$protein.kloss)
  x$protein.kloss <- x$protein.kloss[which(apply(x$protein.kloss,1,FUN=function(y){ !all(is.na(y)) })),]
  x$info$protk.method <- paste("Protein-level k_loss calculated throught the",method,"method with",ag.metric,"of peptides",ifelse(ag.metric=="mean" & !is.null(ag.weights),paste("weighted by",ifelse(ag.weights=="both","nbpoints/variance",ag.weights)),""))
  x$info$protk.call <- match.call()
  return(x)
}

#' Prepare petides for protein-level kloss calculation
#'
#' @param e a matrix/data.frame of peptides values pSILAC object.
#' @param in.all Starting with how many peptides should peptides not quantified in all samples be ignored (default 2), or FALSE to use all peptides, even if not quantified in all samples. 
#' @param removeOutliers Starting with how many peptides should outlier peptides be removed, or FALSE to disable outlier removal (default 5). See ?removeOutlierPeptides
#'
#' @return A filtered version of e.
#'
#' @export
.modProtPre <- function(e, in.all, removeOutliers){
  if(is.null(e))	return(e)
  if(is.null(dim(e))) e <- as.data.frame(matrix(e,nrow=1))
  if(!is.matrix(e) & !is.data.frame(e)) return(e)
  if(in.all){
    nna <- apply(e,1,FUN=function(x){ sum(is.na(x)) })
    if(sum(nna > 0) <= (nrow(e)-in.all)){
	e <- e[which(nna==0),]
    }
  }
  if(removeOutliers)	e <- removeOutlierPeptides(e, removeOutliers)
  if(is.null(dim(e))) e <- as.data.frame(matrix(e,nrow=1))
  return(e)
}

#' Models the protein-level kloss
#'
#' This is a subroutine modeling a single protein's kloss, called by modelProteinsKloss.
#'
#' @param o a pSILAC object.
#' @param method the method used to calculate protein-level rates (default 'combined'); either 'RIA', 'hol' or 'combined' (uses both RIA- and H/L-based peptide kloss).
#' @param ag.weights the method to calculate weights for weighted mean (default 'variance'). Ignored if ag.metric != 'mean'. Either 'variance' (1/variance of the model's fit) or 'nbpoints' (number of datapoints used for the fit).
#' @param unique.weights whether unique weights should be used for each peptides (default T), or if false a weight is calculated individually for each sample.
#' @param in.all Starting with how many peptides should peptides not quantified in all samples be ignored (default 2), or FALSE to use all peptides, even if not quantified in all samples. 
#' @param removeOutliers Starting with how many peptides should outlier peptides be removed, or FALSE to disable outlier removal (default 5). See ?removeOutlierPeptides
#' @param tryRobust Whether to try fitting a robust linear model when remodeling on the basis of H/L ratio (requires the 'MASS' library). Disabled by default.
#' @param freeIntersect Logical; whether to fit also the intersect for the NLI-based 
#' @param returnSD Logical; whether to return also the SD of the kloss (default FALSE).
#'
#' @return A numeric vector with the protein's kloss across samples.
#'
#' @export
.modProt <- function(o, p, method, ag.weights, unique.weights, in.all, removeOutliers, tryRobust, freeIntersect=F, returnSD=F){
  if(!("pSILAC" %in% class(o)))	stop("'o' should be a pSILAC object... and you probably shouldn't be calling the .modProt function directly... see ?modelProteinsKloss")
  if(method %in% c("combined","RIA")){
    e <- .modProtPre(getPeptides(o, protein=p, returnValues="RIA",search=F), in.all, removeOutliers)
    if(is.null(e)){
      ria <- as.numeric(rep(NA,length(unique(o$design$sample))*3))
    }else{
      ria <- as.numeric(sapply( unique(o$design$sample), FUN=function(s){
	w <- which(o$design$sample==s)
	getRIAmod(rep(o$design$time[w],nrow(e)),as.numeric(t(e[,w])))[c(1,2,4)]
      }))
    }
    names(ria) <- paste(rep(unique(o$design$sample), each=3), rep(c("kloss","stderr","nbpoints"), length(unique(o$design$sample))), sep=".")
  }
  if(method %in% c("combined","NLI")){
    pep <- row.names(getPeptides(o, protein=p, search=F))
    psf <- as.numeric(apply(o$NCS[pep,],1,na.rm=T,FUN=median))
    if(!freeIntersect)  csu <- mean(psf,na.rm=T)*o$NCS[pep,]/psf
    e <- .modProtPre(mean(psf,na.rm=T)*o$NLI[pep,]/psf, in.all, F)
    if(is.null(e)){
      nli <- as.numeric(rep(NA,length(unique(o$design$sample))*3))
    }else{
      if(nrow(e)==0){
        nli <- as.numeric(rep(NA,length(unique(o$design$sample))*3))
        }else{
            nli <- as.numeric(sapply( unique(o$design$sample), e=e, csu=csu, freeIntersect=freeIntersect, FUN=function(s, e, csu, freeIntersect){
                w <- which(o$design$sample==s)
                getNLImod(rep(o$design$time[w],nrow(e)),as.numeric(t(e[,w])),ifelse(freeIntersect,NULL,median(as.matrix(csu[,w]),na.rm=T)))[c(1,2,4)]
            }))
        }
    }
    names(nli) <- paste(rep(unique(o$design$sample), each=3), rep(c("kloss","stderr","nbpoints"), length(unique(o$design$sample))), sep=".")
  }
  if(method == "hol"){
    e <- .modProtPre(getPeptides(o, protein=p, returnValues="hol",search=F), in.all, removeOutliers)
    if(is.null(e)){
      hol<- as.numeric(rep(NA,length(unique(o$design$sample))*3))
    }else{
      hol <- as.numeric(sapply( unique(o$design$sample), FUN=function(s){
	w <- which(o$design$sample==s)
	getHoLmod(rep(o$design$time[w],nrow(e)),as.numeric(t(e[,w])), tryRobust)[c(1,2,4)]
      }))
    }
    names(hol) <- paste(rep(unique(o$design$sample), each=3), rep(c("kloss","stderr","nbpoints"), length(unique(o$design$sample))), sep=".")
  }
  if(method == "combined"){
    return(aggregateKloss(rbind(ria,nli), unique(o$design$sample), ag.metric="mean", ag.weights=ag.weights, unique.weights=unique.weights, in.all=in.all, removeOutliers=removeOutliers))
  }
  if(returnSD){
	cols <- paste(unique(o$design$sample),"kloss",sep=".")
	cn <- unique(o$design$sample)
  }else{
	cols <- paste(rep(unique(o$design$sample),each=2),c("kloss","stderr"),sep=".")
	cn <- paste(rep(unique(o$design$sample),each=2),c("",".stderr"),sep="")
  }  
  dat <- switch(method,
    "RIA" = ria[cols],
    "hol" = hol[cols],
    "NLI" = nli[cols],
    NULL)
  names(dat) <- cn
  return(dat)
}

#' Calculate all protein kloss by fitting all the proteins' peptides on a single model
#'
#' Aggregates peptide-level kloss for all proteins. This requires that calcRIAkloss and/or calcHoLkloss be run first.
#'
#' @param o a pSILAC object.
#' @param method the method used to calculate protein-level rates (default 'combined'); either 'RIA', 'hol' or 'combined' (uses both RIA- and H/L-based peptide kloss).
#' @param ag.weights the method to calculate weights for weighted mean (default 'variance'). Ignored if ag.metric != 'mean'. Either 'variance' (1/variance of the model's fit), 'nbpoints' (number of datapoints used for the fit), or 'both' (nbpoints/variance).
#' @param unique.weights whether unique weights should be used for each peptides (default T), or if false a weight is calculated individually for each sample.
#' @param in.all Starting with how many peptides should peptides not quantified in all samples be ignored (default 2), or FALSE to use all peptides, even if not quantified in all samples. 
#' @param removeOutliers Starting with how many peptides should outlier peptides be removed, or FALSE to disable outlier removal (default 5). See ?removeOutlierPeptides
#' @param ncores Number of cores to use (only for remodeling). Defaults to all cores minus 1.
#' @param tryRobust Whether to try fitting a robust linear model when remodeling on the basis of H/L ratio (requires the 'MASS' library). Disabled by default.
#' @param returnKlossTableOnly Logical; whether to return only the kloss table, rather than the whole pSILAC object (default).
#' @param freeIntersect Logical; whether to fit also the intersect for the NLI-based method.
#' @param returnSD Logical; whether to return also the SD of the kloss (default FALSE).
#'
#' @return The updated pSILAC object, or the protein kloss table if returnKlossTableOnly is TRUE.
#'
#' @export
modelProteinsKloss <- function(o, method="combined", ag.weights="variance", unique.weights=T, in.all=2, removeOutliers=5, ncores=NULL, tryRobust=F, returnKlossTableOnly=F, freeIntersect=F, returnSD=F){
  if(class(o) != "pSILAC")	stop("o should be a pSILAC object.")
  method <- match.arg(method, c("combined","RIA","hol","NLI"))
  if(is.null(removeOutliers))	removeOutliers <- F
  if( (is.logical(removeOutliers) & removeOutliers) |
      (is.numeric(removeOutliers) & !(removeOutliers >= 3)) ){
	stop("removeOutliers should be either FALSE, or an integer >= 3 indiciating the minimum number of peptides for outlier removal.")
  }
  if(method=="NLI" & is.null(o$NLI)) stop(paste("Cannot find NLI values! first, run `normalizeLightChannel`."))
  if(is.null(ncores)){
    library(parallel)
    ncores <- detectCores() - 1
  }else{
    if(ncores>1)	library(parallel)
  }
  if(ncores>1){
    cl <- makeCluster(ncores)
    clusterExport(cl, c("o","in.all","removeOutliers","tryRobust","method","unique.weights","ag.weights"), environment())
    clusterExport(cl, c("getHoLmod","getRIAmod","getNLImod","getPeptides","removeOutlierPeptides","aggregateKloss",".mwm",".modProt",".modProtPre","medianNorm"))
    ret <- as.data.frame(t(parSapply(cl, unique(o$peptides$protein), FUN=function(p){
      .modProt(o, as.character(p), method, ag.weights, unique.weights, in.all, removeOutliers, tryRobust, freeIntersect)
    })))
    stopCluster(cl)
  }else{
    ret <- as.data.frame(t(sapply(unique(o$peptides$protein), FUN=function(p){
      .modProt(o, as.character(p), method, ag.weights, unique.weights, in.all, removeOutliers, tryRobust, freeIntersect)
    })))
  }
  row.names(ret) <- unique(o$peptides$protein)
  ret <- replace(ret, is.na(ret), NA)
  if(returnKlossTableOnly)	return(ret)
  o$protein.kloss <- ret[which(apply(ret,1,FUN=function(x){ !all(is.na(x))})),]
  o$info$protk.method <- paste("Protein-level k_loss calculated through remodeling the",ifelse(method=="combined","combined RIA and H/L",method),"values of peptides",ifelse(method=="combined" & !is.null(ag.weights),paste(", weighted by",ag.weights),""))
  o$info$protk.call <- match.call()
  return(o)
}

.mwm <- function(v,w,minVals=3,na.rm=T){
	if(!(length(v)>0))	return(NA)
	if(sum(w>0,na.rm=T) < minVals)	return(mean(v,na.rm=na.rm))
	w[which(is.na(w))] <- 0
	weighted.mean(v,w,na.rm=na.rm)
}

#' aggregateKloss
#'
#' Aggregate peptides' kloss into a protein-level kloss. Rather than calling this function directly, use calcProteinsKloss
#'
#' @param e a data.frame of peptide-level kloss, as in the RIA.kloss or hol.kloss (or NLI.kloss) slots of a pSILAC object.
#' @param samples the vector of samples' names.
#' @param ag.metric the aggregation metric used for protein-level rates (default 'mean'). Either 'mean' or 'median'. If 'mean' is chosen, and if ag.weights is not null, a weighted mean will be performed.
#' @param ag.weights the method to calculate weights for weighted mean (default 'variance'). Ignored if ag.metric != 'mean'. Either 'variance' (1/variance of the model's fit), 'nbpoints' (number of datapoints used for the fit), or 'both' (nbpoints/variance).
#' @param unique.weights whether unique weights should be used for each peptides (default T), or if false a weight is calculated individually for each sample.
#' @param in.all Starting with how many peptides should peptides not quantified in all samples be ignored (default 2), or FALSE to use all peptides, even if not quantified in all samples. 
#' @param removeOutliers Starting with how many peptides should outlier peptides be removed, or FALSE to disable outlier removal (default 5). See ?removeOutlierPeptides
#' @param returnSD Logical; whether to return also the SD of the kloss (default FALSE).
#'
#' @return A vector of length=length(samples) with protein-level kloss corresponding to each sample.
#'
#' @export
aggregateKloss <- function(e, samples, ag.metric="mean", ag.weights="variance", unique.weights=T, in.all=2, removeOutliers=5, returnSD=F){
  if(!is.na(ag.weights))	ag.weights <- match.arg(ag.weights, c("variance","nbpoints","both"))
  ag.metric <- match.arg(ag.metric, c("mean","median"))
  if(in.all){
    nna <- apply(e,1,FUN=function(x){ sum(is.na(x)) })
    if(sum(nna > 0) <= (nrow(e)-in.all)){
	e <- e[which(nna==0),]
    }
  }
  if(nrow(e)==1){
	if(returnSD){
		return(e[,paste(rep(samples,each=2),c("kloss","kloss.stderr"),sep=".")])
		
	}else{
		return(e[,paste(samples,"kloss",sep=".")])
	}
  }
  if(nrow(e)==0)	return(rep(NA,ifelse(returnSD,2,1)*length(samples)))
  if(removeOutliers)	e <- removeOutlierPeptides(e, removeOutliers)
  if(ag.metric=="mean" & !is.null(ag.weights) & !is.na(ag.weights)){
    if(ag.weights=="nbpoints"){
      weights <- as.matrix(e[,paste(samples,"nbpoints",sep=".")])
    }else{
      weights <- as.matrix(1/e[,grep("stderr",colnames(e),fixed=T)])
      weights[which(is.na(as.matrix(weights)) | is.nan(as.matrix(weights)))] <- 1
      if(ag.weights=="both")	weights <- as.matrix(weights*e[,paste(samples,"nbpoints",sep=".")])
    }
    weights[which(is.infinite(weights))] <- min(c(0,as.numeric(weights[which(!is.infinite(weights) & weights>0)])),na.rm=T)
    if(unique.weights){
      weights <- apply(weights,1,na.rm=T,FUN=median)
      k <- (sapply(samples,e=e,weights=weights,FUN=function(x,e,weights){ .mwm(e[,paste(x,"kloss",sep=".")], weights) }))
    }else{
      colnames(weights) <- samples
      k <- (sapply(samples,e=e,weights=weights,FUN=function(x,e,weights){ .mwm(e[,paste(x,"kloss",sep=".")], weights[,x]) }))
    }
  }else{
	if(ag.metric=="mean"){
		k <- sapply(samples,e=e,FUN=function(x,e){ mean(e[,paste(x,"kloss",sep=".")], na.rm=T) })
	}else{
		k <- sapply(samples,e=e,FUN=function(x,e){ median(e[,paste(x,"kloss",sep=".")], na.rm=T) })
	}
  }
  if(returnSD){
  	k <- rep(k,each=2)
  	k[2*(1:(length(k)/2))] <- sapply(samples, e=e, FUN=function(x, e){
		.errProp(e[,paste(x,".kloss.stderr",sep="")], e[,paste(x,".kloss",sep="")])
	})
  }
  return(as.numeric(replace(k, is.na(k), NA)))
}


#' compareParameters
#'
#' Runs a series of protein-level kloss calculations and compares the results. See ?compareParametersResults to plot the results.
#'
#' @param o a pSILAC object, with all peptide-level rates calculated (see ?calcAllRates)
#' @param replicates A vector of length equal to the number of samples, or a character designating the column in the design data.frame, indicating which samples are replicates of the same condition. This will be used to calculate distance/correlation between replicates.
#' @param ncores The number of cores to use. If NULL, will default to detected cores minus 1.
#' @param m A list of the comparisons to make.
#'
#' @return A list of the results, each with the following slots:
#' 'call': the function and parameters called.
#' 'kloss': the protein kloss table
#' 'params': the list of input parameters
#' 'cor': the pearson correlations between replicates.
#' 'dist': the euclidean distance between replicates.
#' 'MdAE': the median absolute error between replicates.
#'
#' @export
compareParameters <- function(o, replicates=NULL, ncores=NULL, m=NULL){
  if(class(o) != "pSILAC")	stop("o should be a pSILAC object.")
  if(is.null(o$RIA.kloss) | is.null(o$hol.kloss))	stop("Peptide-level k_loss should first be calculated.")
  if(is.null(ncores)){
    library(parallel)
    ncores <- detectCores() - 1
  }
  if(is.null(replicates))	replicates <- rep(1,length(unique(o$design$sample)))
  if(length(replicates)==1){
    if(replicates %in% colnames(o$design))	replicates <- o$design[which(!duplicated(o$design$sample)),replicates]
  }
  if(length(replicates) != length(unique(o$design$sample)))	stop("'replicates' should either indicate a column of the design data.frame, or should be a vector of length equal to the number of samples.")
  mnames <- c("method","ag.metric","ag.weights","unique.weights","in.all","removeOutliers","ncores","tryRobust")
  if(is.null(m)){
    m <- list(list("mixed","median","both",T,2,5,3,T),
	      list("mixed","mean","both",T,2,5,3,T),
	      list("RIA","remodel",NA,T,2,5,3,F),
	      list("hol","remodel",NA,T,2,5,3,T),
	      list("NLI","remodel",NA,T,2,5,3,T),
	      list("NLI","remodel",NA,T,2,5,3,F),
	      list("NLI","remodel",NA,T,2,5,3,T),
	      list("combined","remodel","variance",T,2,5,3,T),
	      list("combined","remodel","nbpoints",T,2,5,3,T),
	      list("RIA","median",NA,F,2,5,3,F),
	      list("RIA","median",NA,T,2,5,3,F),
	      list("RIA","mean","variance",T,2,5,3,F),
	      list("RIA","mean","nbpoints",T,2,5,3,F),
	      list("RIA","mean","both",T,F,5,3,F),
	      list("RIA","mean","both",T,2,5,3,F),
	      list("hol","median",NA,T,2,5,3,T),
	      list("hol","mean","variance",T,2,5,3,T),
	      list("hol","mean","nbpoints",T,2,5,3,T),
	      list("combined","median",NA,T,2,5,3,T),
	      list("combined","mean","variance",T,F,F,3,T),
	      list("combined","mean","both",T,F,F,3,T),
	      list("combined","mean","nbpoints",T,F,F,3,T),
	      list("combined","mean","variance",T,2,5,3,T),
	      list("combined","mean","both",T,2,5,3,T)
	 )
  }
  v <- cbind(1:length(replicates),1:length(replicates))
  ll <- list()
  llcor <- list()
  lldist <- list()
  llmdae <- list()
  for(i in 1:length(m)){
    params <- m[[i]]
    names(params) <- mnames
    if(params[[1]]=="mixed"){
      message(paste("calcMixedProteinsKloss(o, ",paste(paste(mnames[c(-1,-3,-8)],params[c(-1,-3,-8)],sep="="),collapse=", "),")",sep=""))
      calcMixedProteinsKloss(o, m[[i]][[2]], m[[i]][[4]], m[[i]][[5]], m[[i]][[6]], m[[i]][[7]])
    }else{
      message(paste("calcProteinsKloss(o, ",paste(paste(mnames,params,sep="="),collapse=", "),")",sep=""))
      o <- calcProteinsKloss(o, m[[i]][[1]], m[[i]][[2]], m[[i]][[3]], m[[i]][[4]], m[[i]][[5]], m[[i]][[6]], m[[i]][[7]], m[[i]][[8]])
    }
    cc <- .getReplicatesValues(cor(o$protein.kloss, use="pairwise"), replicates)
    ce <- .getReplicatesValues(as.matrix(dist(t(o$protein.kloss))), replicates)
    mdae <- .getReplicatesValues(MdAE(o$protein.kloss), replicates)
    llcor[[paste(m[[i]],collapse=".")]] <- cc
    lldist[[paste(m[[i]],collapse=".")]] <- ce
    llmdae[[paste(m[[i]],collapse=".")]] <- mdae
    ll[[i]] <- list( call=o$info$protk.call, kloss=o$protein.kloss, params=params, cor=cc, dist=ce, MdAE=mdae)
  }
  layout(matrix(1:3,nrow=3))
  boxplot(llcor,ylab="Correlation",main="Comparison across replicates",las=3)
  boxplot(lldist,ylab="Euclidean distance",las=3)
  boxplot(llmdae,ylab="Median absolute error",las=3)
  return(ll)
}

#' compareParametersResults
#'
#' Plots the results of the 'compareParameters' function, taking only the proteins which are given a value in all analyses.
#'
#' @param ll a list, as produced by the compareParameters function.
#' @param replicates A vector of length equal to the number of samples, indicating which samples are replicates of the same condition. 
#' If this is contained in the design data.frame of a pSILAC object, you may retrieve it through: o$design[which(!duplicated(o$design$sample)),"replicateColumnNames"]
#'
#' @return Nothing.
#'
#' @export
compareParametersResults <- function(ll, replicates){
  g <- row.names(ll[[1]][["kloss"]])
  llcor <- list()
  lldist <- list()
  llmdae <- list()
  for(i in 2:length(ll))	g <- g[which(g %in% row.names(ll[[i]][["kloss"]]))]
  lnames <- gsub("TRUE","T",gsub("FALSE","F",as.character(lapply(ll, FUN=function(x){ paste(x$params,collapse=".") })),fixed=T),fixed=T)
  llcor <- lapply(ll, FUN=function(x){ .getReplicatesValues(cor(x$kloss[g,], use="pairwise"), replicates) })
  lldist <- lapply(ll, FUN=function(x){ .getReplicatesValues(as.matrix(dist(t(x$kloss[g,]))), replicates) })
  llmdae <- lapply(ll, FUN=function(x){ .getReplicatesValues(MdAE(x$kloss[g,]), replicates) })
  names(llcor) <- lnames
  names(lldist) <- lnames
  names(llmdae) <- lnames
  layout(matrix(1:3,nrow=3))
  boxplot(llcor,ylab="Correlation",main="Comparison across replicates",las=3,outline=F)
  for(i in 1:length(llcor)) points(rep(i,length(llcor[[i]])), llcor[[i]], pch=20)
  boxplot(lldist,ylab="Euclidean distance",las=3,outline=F)
  for(i in 1:length(llcor)) points(rep(i,length(lldist[[i]])), lldist[[i]], pch=20)  
  boxplot(llmdae,ylab="Median absolute error",las=3,outline=F)
  for(i in 1:length(llcor)) points(rep(i,length(llmdae[[i]])), llmdae[[i]], pch=20)  
}

.getReplicatesValues <- function(x, replicates){
  x[cbind(1:ncol(x), 1:ncol(x))] <- NA
  return(as.numeric(unlist(lapply(unique(replicates),FUN=function(r){
    w <- which(replicates==r)
    x <- x[w,w]
    as.numeric(x[lower.tri(x)])
  }))))
}

#' Calculates median absolute error (MdAE)
#'
#' @param x Either a numeric vector, or a numeric matrix/data.frame with samples as columns, and features as rows.
#' @param y Should be NULL if x is a matrix/data.frame, or (if x is a vector) a numeric vector of length equal to that of x.
#'
#' @return Either the median absolute error (if x and y are vectors), or a symmetric matrix (of size=ncol(x)) of pair-wise MdAE if x is a matrix/data.frame.
#'
#' @export
MdAE <- function(x,y=NULL){
  if(!is.null(y)){
    if(is.matrix(x) | is.data.frame(x) | length(x) != length(y))	stop("If y is defined, x and y should be numeric vectors of the same length. If y is NULL, x should be a numeric matrix of data.frame.")
    return(median(abs(y-x),na.rm=T))
  }
  m <- matrix(0,ncol=ncol(x),nrow=ncol(x))
  for(i in 1:ncol(x)){
    for(j in 1:ncol(x)){
      if(i!=j)	m[i,j] <- MdAE(x[,i],x[,j])
    }
  }
  colnames(m) <- colnames(x)
  rownames(m) <- colnames(x)
  return(m)
}

#' Calculates median relatibe error (MdRE)
#'
#' @param x Either a numeric vector, or a numeric matrix/data.frame with samples as columns, and features as rows.
#' @param y Should be NULL if x is a matrix/data.frame, or (if x is a vector) a numeric vector of length equal to that of x.
#'
#' @return Either the median relative error (if x and y are vectors), or a symmetric matrix (of size=ncol(x)) of pair-wise MdRE if x is a matrix/data.frame.
#'
#' @export
MdRE <- function(x,y=NULL){
  if(!is.null(y)){
    if(is.matrix(x) | is.data.frame(x) | length(x) != length(y))	stop("If y is defined, x and y should be numeric vectors of the same length. If y is NULL, x should be a numeric matrix of data.frame.")
    return(median(2*abs(y-x)/(y+x),na.rm=T))
  }
  m <- matrix(0,ncol=ncol(x),nrow=ncol(x))
  for(i in 1:ncol(x)){
    for(j in 1:ncol(x)){
      if(i!=j)	m[i,j] <- MdRE(x[,i],x[,j])
    }
  }
  colnames(m) <- colnames(x)
  rownames(m) <- colnames(x)
  return(m)
}

getPeptideConcordance <- function(o, method="RIA", ag.func=mad){
	slot <- paste(method,"kloss",sep=".")
	wr <- intersect(row.names(o[[slot]]),row.names(o$peptides)[which(o$peptides$nbProteins==1)])
	wc <- grep("kloss$",colnames(o[[slot]]))	
	row.names(o$peptides)[which(o$peptides$nbProteins==1)]
	e <- aggregate(o[[slot]][wr,wc], by=list(protein=o$peptides[wr,"protein"]), na.rm=T, FUN=ag.func)
	row.names(e) <- e$protein
	e$protein <- NULL
	return(e)
}


# to double-check
.errProp <- function(sds, values=NULL){
	sds <- sqrt( sum(sds^2,na.rm=T)/(sum(!is.na(sds))^2) )
	if(is.null(values))	return(sds)
	mm <- suppressWarnings(max(sd(values,na.rm=T), sds, na.rm=T))
	if(is.infinite(mm)) return(NA)
	return(mm)
}

