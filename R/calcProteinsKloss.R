#' Calculate Protein-Level k_loss from Peptide-Level Data
#'
#' Aggregates peptide-level k_loss values to calculate protein-level rates for all proteins in a pSILAC dataset.
#' Requires prior calculation of peptide-level k_loss values using `calcRIAkloss`, `calcHoLkloss`, and/or `calcNLIkloss`.
#'
#' @param o A `pSILAC` object containing peptide-level k_loss data for proteins.
#'
#' @param method Character string specifying the method for calculating protein-level rates. Options are:
#'   \describe{
#'     \item{`"RIA"`}{Uses only RIA-based peptide k_loss values.}
#'     \item{`"hol"`}{Uses only hol-based peptide k_loss values.}
#'     \item{`"NLI"`}{Uses only NLI-based peptide k_loss values.}
#'     \item{`"combined"`}{Combines both RIA-based and NLI-based peptide k_loss values.}
#'     \item{`"complement"`}{Adds complementary NLI-based k_loss values to the RIA-based results.}
#'   }
#'
#' @param ag.metric Character string specifying the aggregation metric for protein-level rates. Options are:
#'   \describe{
#'     \item{`"mean"` (default)}{Calculates the mean of peptide-level k_loss values; if `ag.weights` is provided, a weighted mean is calculated.}
#'     \item{`"median"`}{Calculates the median of peptide-level k_loss values.}
#'     \item{`"remodel"`}{Uses the `modelProteinsKloss` function to re-fit protein models using all peptide-level data. Does not work in the current version of KdeggeR.}
#'   }
#'
#' @param ag.weights Character string defining the weighting method for the weighted mean. Ignored if `ag.metric` is not set to `"mean"`. Options are:
#'   \describe{
#'     \item{`"variance"` (default)}{Uses the inverse variance of the model fit as weights.}
#'     \item{`"nbpoints"`}{Uses the number of datapoints for each peptide as weights.}
#'     \item{`"both"`}{Combines inverse variance and number of datapoints as weights.}
#'   }
#'
#' @param unique.weights Logical; if `TRUE` (default), unique weights are used for each peptide. If `FALSE`, weights are calculated individually for each sample.
#'
#' @param in.all Numeric or logical, specifying the minimum number of samples in which peptides must be quantified to be included. Default is `2`. Set to `FALSE` to include all peptides, regardless of quantification across samples.
#'
#' @param removeOutliers Numeric or logical, specifying the threshold for outlier removal. If a numeric value (default is `5`), removes outlier peptides based on this threshold. Set to `FALSE` to disable outlier removal. See `?removeOutlierPeptides` for details.
#'
#' @param ncores Numeric, specifying the number of cores to use for parallel processing. Defaults to all available cores minus one.
#'
#' @param tryRobust Logical; if `TRUE`, attempts to fit a robust linear model when remodeling based on the H/L ratio. Requires the `MASS` library. Ignored if `method` is not set to `"remodel"`. Default is `FALSE`.
#'
#' @param returnKlossTableOnly Logical; if `TRUE`, only the k_loss table is returned, instead of the full `pSILAC` object. Default is `FALSE`.
#'
#' @param returnSD Logical; if `TRUE`, includes the standard deviation (SD) of the k_loss values in the output. Default is `FALSE`.
#'
#' @return If `returnKlossTableOnly` is `FALSE` (default), returns the updated `pSILAC` object with the added protein-level k_loss table.
#'   If `returnKlossTableOnly` is `TRUE`, returns only the protein k_loss table as a data frame.
#'
#' @export
calcProteinsKloss <- function(o, method = "combined", ag.metric = "mean", ag.weights = "both", unique.weights = TRUE, in.all = 2, removeOutliers = 5, ncores = NULL, tryRobust = FALSE, returnKlossTableOnly = FALSE, returnSD = FALSE) {
  # Check if the input is a valid pSILAC object
  if (class(o) != "pSILAC") stop("o should be a pSILAC object.")
  
  # Validate and set default for removeOutliers
  if (is.null(removeOutliers)) removeOutliers <- FALSE
  if ((is.logical(removeOutliers) & removeOutliers) || (is.numeric(removeOutliers) & !(removeOutliers >= 3))) {
    stop("removeOutliers should be either FALSE, or an integer >= 3 indicating the minimum number of peptides for outlier removal.")
  }
  
  # Set the method for calculating protein-level k_loss and validate it
  method <- match.arg(method, c("combined", "complement", "RIA", "hol", "NLI"))
  
  # If 'remodel' method is used, call modelProteinsKloss and return results
  if (ag.metric == "remodel" & method != "complement") {
    stop("Protein k_loss remodeling not yet enabled.")
    return(modelProteinsKloss(o, method = method, ag.weights = ag.weights, in.all = in.all, removeOutliers = removeOutliers, ncores = ncores, tryRobust = tryRobust, returnKlossTableOnly = returnKlossTableOnly, returnSD = returnSD))
  }
  
  # Set and validate the aggregation metric
  ag.metric <- match.arg(ag.metric, c("remodel", "median", "mean"))
  
  # If 'combined' or 'complement' method is chosen, aggregate k_loss using both RIA and NLI methods
  if (method %in% c("combined", "complement")) {
    
    # Check if peptide-level k_loss has been calculated for RIA and NLI
    if (is.null(o$RIA.kloss) | is.null(o$NLI.kloss)) stop("Peptide-level k_loss should first be calculated.")
    
    # Calculate protein k_loss using both RIA and NLI-based methods
    k1 <- calcProteinsKloss(o, method = "RIA", ag.metric = ag.metric, ag.weights = ag.weights, unique.weights = unique.weights, in.all = in.all, removeOutliers = removeOutliers, ncores = ncores, tryRobust = tryRobust, returnKlossTableOnly = TRUE, returnSD = method == "combined")
    k2 <- calcProteinsKloss(o, method = "NLI", ag.metric = ag.metric, ag.weights = ag.weights, unique.weights = unique.weights, in.all = in.all, removeOutliers = removeOutliers, ncores = ncores, tryRobust = tryRobust, returnKlossTableOnly = TRUE, returnSD = method == "combined")
    
    # Combine or complement RIA- and NLI-based estimates
    if (method == "complement") {
      message("Adding complementary NLI-based estimates to the RIA-based rates of loss")
      # Filter out incomplete or invalid RIA results and add complementary NLI results
      k1 <- k1[which(!apply(k1, 1, FUN = function(k) { any(is.na(k) | k < 0) })), ]
      k2 <- k2[which(!(row.names(k2) %in% row.names(k1))), ]
      o$protein.kloss <- as.data.frame(rbind(k1, k2))
      o$protein.kloss$source <- rep(c("RIA", "NLI"), c(nrow(k1), nrow(k2)))
    } else {
      message("Combining RIA-based and NLI-based estimates of rates of loss")
      # Use more stable estimates between RIA and NLI methods based on mean-to-median ratio of standard deviations
      k1 <- k1[which(!apply(k1[, 2 * 1:(ncol(k1) / 2) - 1], 1, FUN = function(k) { any(is.na(k) | k < 0) })), ]
      k2comp <- k2[which(!(row.names(k2) %in% row.names(k1))), 2 * 1:(ncol(k2) / 2) - 1]
      k2 <- k2[row.names(k1), ]
      msd1 <- apply(k1[, 2 * 1:(ncol(k1) / 2)], 1, FUN = function(x) { abs(mean(x[2 * 1:(length(x) / 2)], na.rm = TRUE) / median(x[2 * 1:(length(x) / 2) - 1], na.rm = TRUE)) })
      msd2 <- apply(k2[, 2 * 1:(ncol(k2) / 2)], 1, FUN = function(x) { abs(mean(x[2 * 1:(length(x) / 2)], na.rm = TRUE) / median(x[2 * 1:(length(x) / 2) - 1], na.rm = TRUE)) })
      k1 <- k1[, 2 * 1:(ncol(k1) / 2) - 1]
      k2 <- k2[, 2 * 1:(ncol(k2) / 2) - 1]
      # Substitute less stable RIA-based estimates with NLI-based estimates when applicable
      k1[which(msd1 > msd2), ] <- k2[which(msd1 > msd2), ]
      o$protein.kloss <- as.data.frame(rbind(k1, k2comp))
      o$protein.kloss$source <- rep(c("RIA", "NLI"), c(nrow(k1), nrow(k2comp)))
      o$protein.kloss$source[which(msd1 > msd2)] <- "NLI"
    }
    
    # Record method information
    o$info$protk.method <- paste("Protein-level k_loss calculated through the", method, "method with", ag.metric, "of peptides", ifelse(ag.metric == "mean" & !is.null(ag.weights), paste("weighted by", ifelse(ag.weights == "both", "nbpoints/variance", ag.weights)), ""))
    o$info$protk.call <- match.call()
    
    # Return k_loss table only if specified
    if (returnKlossTableOnly) return(o$protein.kloss)
    return(o)
    
  } else {
    # Retrieve peptide-level k_loss data for the selected method
    e <- o[[paste(method, "kloss", sep = ".")]]
    if (is.null(e)) stop("Peptide-level k_loss should first be calculated using the desired method.")
    proteins <- as.character(o$peptides[row.names(e), "protein"])
  }
  
  # Initialize k_loss data frame with appropriate columns for k_loss and SD if requested
  if (returnSD) {
    o$protein.kloss <- as.data.frame(matrix(NA, nrow = length(unique(proteins)), ncol = 2 * length(unique(o$design$sample))), row.names = unique(proteins))
    colnames(o$protein.kloss) <- paste(rep(unique(o$design$sample), each = 2), rep(c("kloss", "stderr"), length(unique(o$design$sample))), sep = ".")
  } else {
    o$protein.kloss <- as.data.frame(matrix(NA, nrow = length(unique(proteins)), ncol = length(unique(o$design$sample))), row.names = unique(proteins))
    colnames(o$protein.kloss) <- unique(o$design$sample)
  }
  
  # Calculate aggregated k_loss for each protein
  i <- 0
  for (p in unique(proteins)) {
    i <- i + 1
    if (i > 100) {
      cat(".")  # Print progress every 100 iterations
      i <- 0
    }
    # Aggregate peptide-level k_loss for the current protein
    o$protein.kloss[p, ] <- aggregateKloss(e[which(proteins == p), ], unique(o$design$sample), ag.metric, ag.weights, unique.weights, in.all, removeOutliers, returnSD = returnSD)
  }
  cat("\n")  
  
  # Return only k_loss table if specified
  if (returnKlossTableOnly) return(o$protein.kloss)
  
  # Remove rows with all NA values
  o$protein.kloss <- o$protein.kloss[which(apply(o$protein.kloss, 1, FUN = function(y) { !all(is.na(y)) })), ]
  
  # Record method information in pSILAC object
  o$info$protk.method <- paste("Protein-level k_loss calculated through the", method, "method with", ag.metric, "of peptides", ifelse(ag.metric == "mean" & !is.null(ag.weights), paste("weighted by", ifelse(ag.weights == "both", "nbpoints/variance", ag.weights)), ""))
  o$info$protk.call <- match.call()
  
  return(o)
}

#' Aggregate Peptide k_loss Values to Protein-Level k_loss
#'
#' This function aggregates peptide-level k_loss values to calculate a protein-level k_loss for each sample.
#' While this function can be used directly, it is generally recommended to use `calcProteinsKloss` for calculating protein-level rates.
#'
#' @param o A data frame containing peptide-level k_loss values, as found in the `RIA.kloss`, `hol.kloss`, or `NLI.kloss` slots of a pSILAC object.
#' @param samples A character vector of sample names for which to aggregate peptide k_loss values.
#' @param ag.metric The aggregation metric for protein-level rates. Options are `"mean"` (default) or `"median"`. If `"mean"` is selected and `ag.weights` is not `NULL`, a weighted mean is performed.
#' @param ag.weights The weighting method for calculating the weighted mean. Options are `"variance"` (default), `"nbpoints"`, or `"both"`. Ignored if `ag.metric` is not `"mean"`. 
#'   - `"variance"`: uses the inverse of the variance of the model's fit.
#'   - `"nbpoints"`: weights by the number of datapoints used for the fit.
#'   - `"both"`: combines both the inverse variance and number of datapoints as weights.
#' @param unique.weights Logical; if `TRUE` (default), unique weights are used for each peptide. If `FALSE`, weights are calculated individually for each sample.
#' @param in.all Numeric or logical; specifies the minimum number of peptides that must be quantified across all samples for inclusion. Default is `2`. Set to `FALSE` to include all peptides, regardless of sample quantification.
#' @param removeOutliers Numeric or logical; specifies the threshold for removing outlier peptides. If a numeric value (default is `5`), removes outliers when the peptide count meets or exceeds this threshold. Set to `FALSE` to disable outlier removal. See `?removeOutlierPeptides` for details.
#' @param returnSD Logical; if `TRUE`, returns the standard deviation (SD) of the k_loss values in addition to the mean or median. Default is `FALSE`.
#'
#' @return A numeric vector of length equal to `samples`, containing the aggregated protein-level k_loss for each sample. 
#'   If `returnSD = TRUE`, the vector length is doubled, with the SD following each k_loss estimate.
#'
#' @details
#' This function aggregates k_loss values from peptides to a protein-level k_loss for each sample. Depending on the specified `ag.metric` and `ag.weights`, the function computes either a mean or median of the peptide values, with an option for weighted means if variance or datapoint count information is provided. Outliers can be removed if the peptide count meets the `removeOutliers` threshold.
#'
#' @examples
#' # Assuming 'o' contains peptide-level k_loss data and 'samples' specifies sample names:
#' aggregateKloss(o, samples = c("Sample1", "Sample2"), ag.metric = "mean", ag.weights = "variance", in.all = 3, returnSD = TRUE)
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
