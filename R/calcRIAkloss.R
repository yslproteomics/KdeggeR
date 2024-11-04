#' Calculate RIA-based k_loss for All Peptides in a Dataset
#'
#' This function serves as a wrapper to apply `getRIAmod` to each peptide in a pSILAC dataset, 
#' calculating RIA-based (Relative Isotopic Abundance) `k_loss` values. The results are stored 
#' in the pSILAC object's `RIA.kloss` field. Multi-core processing is supported for faster computation.
#'
#' @param o A `pSILAC` object containing the dataset to be processed. 
#'          This object should include RIA data and experimental design information.
#' @param ncores The number of cores to use for parallel processing. 
#'               Defaults to one less than the total number of detected cores on the system.
#'               If set to 1 or NULL, single-core processing is used.
#'
#' @return The updated `pSILAC` object with the calculated `k_loss` values added in the 
#'         `RIA.kloss` field. This field includes `k_loss`, standard error, sum of squared residuals, 
#'         and the number of data points used per sample.
#'
#' @details
#' The function iterates over each unique sample in the dataset, applying `getRIAmod` 
#' to calculate the RIA-based `k_loss` for each peptide. Results are stored in a matrix 
#' with columns representing `k_loss`, standard error (`kloss.stderr`), sum of squared residuals (`kloss.SSR`), 
#' and the number of points used (`nbpoints`). If parallel processing is enabled, `parApply` 
#' is used for faster computation.
#'
#'
#' @importFrom parallel makeCluster stopCluster detectCores parApply clusterExport
#' @export
calcRIAkloss <- function(o, ncores = NULL) {
  
  if (class(o) != "pSILAC") stop("o should be a pSILAC object.")
  o$RIA.kloss <- NULL
  x <- unique(o$design$time)
  
  if (is.null(ncores)) {
    library(parallel)
    ncores <- detectCores() - 1
  } else if (ncores > 1) {
    library(parallel)
  }
  
  # Create cluster if needed
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl, c("getRIAmod", "x", "o"), envir = environment())
  }
  
  # Loop through each unique sample
  for (p in unique(o$design$sample)) {
    sample_indices <- which(o$design$sample == p)
    if (ncores > 1) {
      d <- as.data.frame(t(parApply(cl, o$RIA[, sample_indices], 1, function(y) getRIAmod(x, y))))
    } else {
      d <- as.data.frame(t(apply(o$RIA[, sample_indices], 1, function(y) getRIAmod(x, y))))
    }
    
    # Label columns and rows in the output
    colnames(d) <- paste(p, c("kloss", "kloss.stderr", "kloss.SSR", "nbpoints"), sep = ".")
    row.names(d) <- row.names(o$RIA)
    
    # Display message on progress
    message(paste(" ...calculated ", sum(!is.na(d[, 1])), " k_loss values for sample ", p, " (", sum(is.na(d[, 1])), " missing)", sep = ""))
    
    # Combine results
    if (is.null(o$RIA.kloss)) {
      o$RIA.kloss <- d
    } else {
      o$RIA.kloss <- cbind(o$RIA.kloss, d)
    }
  }
  
  return(o)
}
