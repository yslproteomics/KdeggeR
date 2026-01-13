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
calcRIAkloss_v2 <- function(o, ncores = NULL) {
  
  if (!inherits(o, "pSILAC")) stop("o should be a pSILAC object.")
  o$RIA.kloss <- NULL
  
  # Safe wrapper around getRIAmod: enforce same length + consistent NA handling
  safe_getRIAmod <- function(y, time_vec) {
    # y: numeric vector of RIA values for one protein (columns = timepoints)
    # time_vec: numeric vector of timepoints, same length as y
    if (length(y) != length(time_vec)) {
      stop("Internal error in calcRIAkloss_v2: length(y) != length(time_vec). ",
           "Check that design and RIA columns are aligned for this sample.")
    }
    
    ok   <- is.finite(y) & is.finite(time_vec)
    n_ok <- sum(ok)
    
    if (n_ok < 2L) {
      # Return same structure as getRIAmod(): c(kloss, stderr, SSR, nbpoints)
      return(c(NA_real_, NA_real_, NA_real_, n_ok))
    }
    
    getRIAmod(time_vec[ok], y[ok])
  }
  
  # ncores handling
  if (is.null(ncores)) {
    library(parallel)
    ncores <- detectCores() - 1
  } else if (ncores > 1) {
    library(parallel)
  }
  
  # Cluster if needed
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    # Export modeling functions; time_vec is passed as an argument, not exported
    clusterExport(cl, c("getRIAmod", "safe_getRIAmod"), envir = environment())
  }
  
  # Loop through each unique sample
  for (p in unique(o$design$sample)) {
    sample_indices <- which(o$design$sample == p)
    
    # Order columns for this sample by time (sample-specific!)
    sample_times   <- o$design$time[sample_indices]
    ord            <- order(sample_times)
    sample_indices <- sample_indices[ord]
    time_vec       <- sample_times[ord]  # time vector for this sample
    
    # Subset RIA matrix for this sample, in the same order
    ria_mat <- o$RIA[, sample_indices, drop = FALSE]
    
    # Model k_loss per protein
    if (ncores > 1) {
      d <- as.data.frame(t(
        parApply(cl, ria_mat, 1, safe_getRIAmod, time_vec = time_vec)
      ))
    } else {
      d <- as.data.frame(t(
        apply(ria_mat, 1, safe_getRIAmod, time_vec = time_vec)
      ))
    }
    
    # Label columns and rows in the output
    colnames(d) <- paste(p, c("kloss", "kloss.stderr", "kloss.SSR", "nbpoints"),
                         sep = ".")
    row.names(d) <- row.names(o$RIA)
    
    # Progress message
    message(paste(" ...calculated ", sum(!is.na(d[, 1])),
                  " k_loss values for sample ", p,
                  " (", sum(is.na(d[, 1])), " missing)", sep = ""))
    
    # Combine results
    if (is.null(o$RIA.kloss)) {
      o$RIA.kloss <- d
    } else {
      o$RIA.kloss <- cbind(o$RIA.kloss, d)
    }
  }
  
  return(o)
}

