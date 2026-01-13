#' Calculate Logarithmic H/L-based k_loss for all Peptides
#'
#' This function calculates the ln(H/L + 1)-based k_loss for each peptide in a `pSILAC` dataset by applying the `getHoLmod` function. 
#' It handles data across multiple time points and samples, with optional parallel processing to speed up computations.
#'
#' @param o A `pSILAC` object, which contains the data for peptide H/L ratios and associated metadata.
#' @param ncores An integer specifying the number of cores for parallel processing. If `NULL`, it defaults to using all available cores minus one.
#' @param tryRobust Logical; if `TRUE`, attempts to fit a robust linear model (requires the `MASS` package). Ignored if fewer than 10 data points are available.
#'
#' @return The updated `pSILAC` object, with a new `hol.kloss` slot containing the calculated k_loss values and associated statistics (standard error, sum of squared residuals, and number of points) for each peptide.
#'
#' @details
#' The function iterates over each unique sample in the dataset and calculates k_loss values for each peptide, either in parallel or sequentially. 
#' For each sample, the calculated metrics include the k_loss, its standard error, sum of squared residuals, and the number of points used in the calculation.
#' A summary message is displayed for each sample, indicating the number of successfully calculated k_loss values and any missing values.
#'
#' @keywords internal
.calcHoLkloss_v0 <- function(o, ncores = NULL, tryRobust = FALSE) {
  if (class(o) != "pSILAC") stop("o should be a pSILAC object.")
  
  o$hol.kloss <- NULL
  x <- unique(o$design$time)
  
  if (is.null(ncores)) {
    library(parallel)
    ncores <- detectCores() - 1
  } else if (ncores > 1) {
    library(parallel)
  }
  
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    clusterExport(cl, c("getHoLmod"))
    clusterExport(cl, c("x", "o"), environment())
  }
  
  for (p in unique(o$design$sample)) {
    if (ncores > 1) {
      d <- as.data.frame(t(parApply(cl, o$hol[, which(o$design$sample == p)], 1, FUN = function(y) { getHoLmod(x, y) })))
    } else {
      d <- as.data.frame(t(apply(o$hol[, which(o$design$sample == p)], 1, FUN = function(y) { getHoLmod(x, y) })))
    }
    
    colnames(d) <- paste(p, c("kloss", "kloss.stderr", "kloss.SSR", "R_squared", "nbpoints"), sep = ".")
    row.names(d) <- row.names(o$hol)
    
    message(sprintf("...calculated %d k_loss values for sample %s (%d missing)", sum(!is.na(d[, 1])), p, sum(is.na(d[, 1]))))
    
    if (is.null(o$hol.kloss)) {
      o$hol.kloss <- d
    } else {
      o$hol.kloss <- cbind(o$hol.kloss, d)
    }
  }
  
  if (ncores > 1) stopCluster(cl)
  
  return(o)
}

#' Calculate Logarithmic H/L-based k_loss for all Peptides
#'
#' This function calculates the ln(H/L + 1)-based k_loss for each peptide in a `pSILAC` dataset by applying the `getHoLmod` function. 
#' It handles data across multiple time points and samples, with optional parallel processing to speed up computations.
#'
#' @param o A `pSILAC` object, which contains the data for peptide H/L ratios and associated metadata.
#' @param ncores An integer specifying the number of cores for parallel processing. If `NULL`, it defaults to using all available cores minus one.
#' @param tryRobust Logical; if `TRUE`, attempts to fit a robust linear model (requires the `MASS` package). Ignored if fewer than 10 data points are available.
#'
#' @return The updated `pSILAC` object, with a new `hol.kloss` slot containing the calculated k_loss values and associated statistics (standard error, sum of squared residuals, and number of points) for each peptide.
#'
#' @details
#' The function iterates over each unique sample in the dataset and calculates k_loss values for each peptide, either in parallel or sequentially. 
#' For each sample, the calculated metrics include the k_loss, its standard error, sum of squared residuals, and the number of points used in the calculation.
#' A summary message is displayed for each sample, indicating the number of successfully calculated k_loss values and any missing values.
#'
#' @export
calcHoLkloss <- function(o, ncores = NULL, tryRobust = FALSE) {
  if (!inherits(o, "pSILAC")) stop("o should be a pSILAC object.")
  
  o$hol.kloss <- NULL
  
  # Safe wrapper around getHoLmod: enforce same length + consistent NA handling
  safe_getHoLmod <- function(y, time_vec, tryRobust = FALSE) {
    # y: numeric vector of HoL values for one protein (cols = timepoints)
    # time_vec: numeric vector of timepoints, same length as y
    if (length(y) != length(time_vec)) {
      stop("Internal error in calcHoLkloss_v2: length(y) != length(time_vec). ",
           "Check that design and HoL columns are aligned for this sample.")
    }
    
    ok   <- is.finite(y) & is.finite(time_vec)
    n_ok <- sum(ok)
    
    if (n_ok < 2L) {
      # Match getHoLmod() output structure:
      # c(kloss, kloss.stderr, kloss.SSR, R_squared, nbpoints)
      return(c(NA_real_, NA_real_, NA_real_, NA_real_, n_ok))
    }
    
    getHoLmod(time_vec[ok], y[ok], tryRobust = tryRobust)
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
    clusterExport(cl, c("getHoLmod", "safe_getHoLmod"), envir = environment())
  }
  
  # Loop through each unique sample
  for (p in unique(o$design$sample)) {
    sample_indices <- which(o$design$sample == p)
    
    # Order columns for this sample by time (sample-specific!)
    sample_times   <- o$design$time[sample_indices]
    ord            <- order(sample_times)
    sample_indices <- sample_indices[ord]
    time_vec       <- sample_times[ord]  # time vector for this sample
    
    # Subset HoL matrix for this sample, in the same order
    hol_mat <- o$hol[, sample_indices, drop = FALSE]
    
    # Fit k_loss per protein
    if (ncores > 1) {
      d <- as.data.frame(t(
        parApply(cl, hol_mat, 1, safe_getHoLmod,
                 time_vec = time_vec, tryRobust = tryRobust)
      ))
    } else {
      d <- as.data.frame(t(
        apply(hol_mat, 1, safe_getHoLmod,
              time_vec = time_vec, tryRobust = tryRobust)
      ))
    }
    
    colnames(d) <- paste(
      p,
      c("kloss", "kloss.stderr", "kloss.SSR", "R_squared", "nbpoints"),
      sep = "."
    )
    row.names(d) <- row.names(o$hol)
    
    message(sprintf("...calculated %d k_loss values for sample %s (%d missing)",
                    sum(!is.na(d[, 1])), p, sum(is.na(d[, 1]))))
    
    if (is.null(o$hol.kloss)) {
      o$hol.kloss <- d
    } else {
      o$hol.kloss <- cbind(o$hol.kloss, d)
    }
  }
  
  return(o)
}
