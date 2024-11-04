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
