#' Calculate All Protein k_loss by Fitting All Peptides on a Single Model
#'
#' This function aggregates peptide-level k_loss values to calculate protein-level rates by fitting all peptides 
#' associated with each protein into a single model. It requires that `calcRIAkloss()` and/or `calcHoLkloss()` 
#' have been run first to generate peptide k_loss values.
#'
#' @param o A `pSILAC` object containing peptide-level k_loss data for each protein.
#' @param method Character string specifying the method used for calculating protein-level rates. Options are `"RIA"`, `"hol"`, or `"combined"` (default), which combines both RIA- and H/L-based peptide k_loss values.
#' @param ag.weights Character string specifying the method to calculate weights for the weighted mean. Default is `"variance"`, which uses the inverse variance of the model's fit. Other options include `"nbpoints"` (based on the number of datapoints used for the fit) and `"both"` (combining both datapoint count and variance). Ignored if `ag.metric` is not set to `"mean"`.
#' @param unique.weights Logical; if `TRUE` (default), unique weights are used for each peptide. If `FALSE`, weights are calculated individually for each sample.
#' @param in.all Numeric or logical, specifying the minimum number of peptides that must be quantified in all samples for inclusion. Default is `2`. Set to `FALSE` to use all peptides, regardless of sample quantification.
#' @param removeOutliers Numeric or logical, specifying the threshold for removing outlier peptides. If a numeric value (default is `5`), outliers are removed based on this threshold. Set to `FALSE` to disable outlier removal. See `?removeOutlierPeptides` for details.
#' @param ncores Numeric, specifying the number of cores to use for parallel processing (default is all available cores minus one).
#' @param tryRobust Logical; if `TRUE`, attempts to fit a robust linear model when modeling based on the H/L ratio. Requires the `MASS` library. Default is `FALSE`.
#' @param returnKlossTableOnly Logical; if `TRUE`, returns only the k_loss table instead of the full `pSILAC` object. Default is `FALSE`.
#' @param freeIntersect Logical; if `TRUE`, also fits the intercept when using the NLI-based method.
#' @param returnSD Logical; if `TRUE`, includes the standard deviation (SD) of the k_loss in the output. Default is `FALSE`.
#'
#' @return If `returnKlossTableOnly` is `FALSE` (default), returns the updated `pSILAC` object with the added protein-level k_loss table. If `returnKlossTableOnly` is `TRUE`, returns only the protein k_loss table as a data frame.
#'
#' @details
#' This function calculates protein-level k_loss by aggregating peptide-level values, using the specified method (`"RIA"`, `"hol"`, or `"combined"`). The weighting method for aggregation can be specified, allowing flexibility in how the final protein k_loss rate is calculated. Peptides can be filtered based on the `in.all` criterion (number of samples with quantification) and by removing outliers if `removeOutliers` is enabled.
#' 
#' @export
modelProteinsKloss <- function(o, method = "combined", ag.weights = "variance", unique.weights = TRUE, in.all = 2, removeOutliers = 5, ncores = NULL, tryRobust = FALSE, returnKlossTableOnly = FALSE, freeIntersect = FALSE, returnSD = FALSE) {
  if (class(o) != "pSILAC") stop("o should be a pSILAC object.")
  method <- match.arg(method, c("combined", "RIA", "hol", "NLI"))
  if (is.null(removeOutliers)) removeOutliers <- FALSE
  if ((is.logical(removeOutliers) & removeOutliers) || (is.numeric(removeOutliers) & !(removeOutliers >= 3))) {
    stop("removeOutliers should be either FALSE, or an integer >= 3 indicating the minimum number of peptides for outlier removal.")
  }
  if (method == "NLI" & is.null(o$NLI)) stop("Cannot find NLI values! First, run `normalizeLightChannel`.")
  
  # Set up parallel processing
  if (is.null(ncores)) {
    library(parallel)
    ncores <- detectCores() - 1
  } else {
    if (ncores > 1) library(parallel)
  }
  
  # Parallel or sequential processing of proteins
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    clusterExport(cl, c("o", "in.all", "removeOutliers", "tryRobust", "method", "unique.weights", "ag.weights"), environment())
    clusterExport(cl, c("getHoLmod", "getRIAmod", "getNLImod", "getPeptides", "removeOutlierPeptides", "aggregateKloss", ".mwm", ".modProt", ".modProtPre", "medianNorm"))
    ret <- as.data.frame(t(parSapply(cl, unique(o$peptides$protein), FUN = function(p) {
      .modProt(o, as.character(p), method, ag.weights, unique.weights, in.all, removeOutliers, tryRobust, freeIntersect)
    })))
    stopCluster(cl)
  } else {
    ret <- as.data.frame(t(sapply(unique(o$peptides$protein), FUN = function(p) {
      .modProt(o, as.character(p), method, ag.weights, unique.weights, in.all, removeOutliers, tryRobust, freeIntersect)
    })))
  }
  
  # Set row names and replace NA values
  row.names(ret) <- unique(o$peptides$protein)
  ret <- replace(ret, is.na(ret), NA)
  
  # Return either the k_loss table or the full pSILAC object
  if (returnKlossTableOnly) return(ret)
  
  # Store results in pSILAC object
  o$protein.kloss <- ret[which(apply(ret, 1, FUN = function(x) { !all(is.na(x)) })), ]
  o$info$protk.method <- paste("Protein-level k_loss calculated through remodeling the", ifelse(method == "combined", "combined RIA and H/L", method), "values of peptides", ifelse(method == "combined" & !is.null(ag.weights), paste(", weighted by", ag.weights), ""))
  o$info$protk.call <- match.call()
  
  return(o)
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

#' Prepare Peptides for Protein-Level k_loss Calculation
#'
#' This function filters a matrix or data frame of peptide values in preparation for protein-level k_loss calculation. 
#' It allows for the optional exclusion of peptides not quantified in all samples and removal of outlier peptides.
#'
#' @param o A matrix or data frame of peptide values, typically from a pSILAC object.
#' @param in.all Numeric or logical, specifying the threshold for peptides not quantified in all samples.
#'   If a numeric value (default is 2), only peptides quantified in `in.all` or more samples are retained. 
#'   Set to `FALSE` to include all peptides, regardless of quantification across samples.
#' @param removeOutliers Numeric or logical, specifying the threshold for outlier removal.
#'   If a numeric value (default is 5), removes outlier peptides starting at the specified threshold.
#'   Set to `FALSE` to disable outlier removal. See `?removeOutlierPeptides` for details on outlier identification.
#'
#' @return A filtered version of `o` with peptides potentially removed based on the `in.all` and `removeOutliers` criteria.
#'
#' @details
#' The function applies two main filtering steps:
#' \enumerate{
#'   \item If `in.all` is specified as a numeric value, peptides not quantified in all samples will be ignored if their 
#'   missing values exceed the threshold.
#'   \item If `removeOutliers` is specified as a numeric value, outlier peptides are removed using `removeOutlierPeptides` 
#'   when the peptide count meets the threshold.
#' }
#' The function returns the input object `o` with the applied filters or in its original form if filtering criteria are not met.
#'
#' @export
.modProtPre <- function(o, in.all = 2, removeOutliers = 5) {
  if (is.null(o)) return(o)
  if (is.null(dim(o))) o <- as.data.frame(matrix(o, nrow = 1))
  if (!is.matrix(o) & !is.data.frame(o)) return(o)
  
  # Filter peptides not quantified in all samples based on `in.all` threshold
  if (in.all) {
    nna <- apply(o, 1, FUN = function(x) { sum(is.na(x)) })
    if (sum(nna > 0) <= (nrow(o) - in.all)) {
      o <- o[which(nna == 0), ]
    }
  }
  
  # Remove outliers based on `removeOutliers` threshold
  if (removeOutliers) o <- removeOutlierPeptides(o, removeOutliers)
  if (is.null(dim(o))) o <- as.data.frame(matrix(o, nrow = 1))
  
  return(o)
}


#' Calculate Modified Weighted Mean
#'
#' This function calculates a weighted mean of the values in `v` using weights in `w`. If the number of non-zero weights 
#' is below a specified threshold (`minVals`), the function returns the arithmetic mean of `v` instead. 
#' Missing values in `w` are treated as zero.
#'
#' @param v A numeric vector of values for which the weighted mean is to be calculated.
#' @param w A numeric vector of weights corresponding to the values in `v`. Must be the same length as `v`.
#'   Missing values in `w` are set to zero.
#' @param minVals An integer specifying the minimum number of non-zero weights required to calculate the weighted mean.
#'   If the count of non-zero weights is below this threshold, the function returns the arithmetic mean of `v`.
#'   Default is 3.
#' @param na.rm Logical, indicating whether NA values should be removed from `v` and `w` before calculating the mean.
#'   Default is `TRUE`.
#'
#' @return A numeric value representing either the weighted mean of `v` or the arithmetic mean of `v` 
#'   (if there are fewer than `minVals` non-zero weights).
#'
#' @details
#' This function provides flexibility when calculating a weighted mean, returning an unweighted mean if there are 
#' insufficient non-zero weights. This is useful in cases where an adequate number of reliable weights is required 
#' to apply weighting. Missing values in the weight vector `w` are automatically set to zero.
#'
#'
#' @export
.mwm <- function(v, w, minVals = 3, na.rm = TRUE) {
  if (!(length(v) > 0)) return(NA)
  if (sum(w > 0, na.rm = na.rm) < minVals) return(mean(v, na.rm = na.rm))
  
  # Replace NAs in weights with 0
  w[is.na(w)] <- 0
  
  # Calculate weighted mean
  weighted.mean(v, w, na.rm = na.rm)
}
