#' @title Calculate k_loss for All Peptides Based on Normalized Light Intensities
#'
#' @description This function applies `getNLImod` to each peptide in the dataset, calculating the decay rate constant \( k_{loss} \) 
#' based on normalized light intensities (NLI). Users can specify how to set the starting intensity for the model, 
#' with options to use a model-fitted intercept, maximum observed intensity, or median of non-missing intensities.
#'
#' @param o A `pSILAC` object containing the data to be analyzed.
#' @param startIntensity Character. Specifies the method for determining the starting intensity for model fitting:
#' \itemize{
#'   \item \code{"model"} - Fit the model with an intercept, allowing the initial intensity to be estimated during the fitting. This option is not functional in the current version. 
#'   \item \code{"max"} - Use the maximum of observed H+L intensities as the starting intensity (default).
#'   \item \code{"median"} - Use the median of non-null H+L intensities as the starting intensity. Consider this option if the data are not filtered using `filterMonotone`
#' }
#' @param ncores Integer. The number of cores to use for parallel computation. If \code{NULL} (default), the function detects the 
#' available cores and uses all but one. Single-core processing is used if \code{ncores = 1}.
#'
#' @return The updated `pSILAC` object with an additional `NLI.kloss` component. This component contains calculated values for each peptide, including:
#' \itemize{
#'   \item \code{kloss} - The decay rate constant.
#'   \item \code{kloss.stderr} - Standard error of the decay rate.
#'   \item \code{kloss.SSR} - Sum of squared residuals from the model fit.
#'   \item \code{nbpoints} - Number of points used in the fitting.
#' }
#'
#' @export
#'
calcNLIkloss <- function(o, startIntensity = "max", ncores = NULL) {
  
  if (class(o) != "pSILAC") stop("o should be a pSILAC object.")
  
  # Match the startIntensity argument
  startIntensity <- match.arg(startIntensity, c("median", "max", "model"))
  
  # Determine if we're in free intercept mode
  freeIntersect <- startIntensity == "model"
  
  o$NLI.kloss <- NULL
  x <- unique(o$design$time)
  
  # Set up parallel computation if specified
  if (is.null(ncores)) {
    library(parallel)
    ncores <- detectCores() - 1
  } else if (ncores > 1) {
    library(parallel)
  }
  
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    clusterExport(cl, c("getNLImod", ".medianNonNull"))
    clusterExport(cl, c("x", "o", "freeIntersect"), environment())
  }
  
  for (p in unique(o$design$sample)) {
    if (ncores > 1) {
      d <- as.data.frame(t(parApply(cl, cbind(
        suppressWarnings(apply(o$NCS[, which(o$design$sample == p)], 1, function(x) {
          if (startIntensity == "median") {
            .medianNonNull(x)
          } else if (startIntensity == "max") {
            max(x, na.rm = TRUE)
          } else if (startIntensity == "model") {
            NULL # Pass NULL to fit the intercept within getNLImod
          }
        })), o$NLI[, which(o$design$sample == p)]), 1, FUN = function(y) {
          getNLImod(x, y[-1], ifelse(freeIntersect, NULL, y[1]))
        })))
    } else {
      d <- as.data.frame(t(apply(cbind(
        suppressWarnings(apply(o$NCS[, which(o$design$sample == p)], 1, function(x) {
          if (startIntensity == "median") {
            .medianNonNull(x)
          } else if (startIntensity == "max") {
            max(x, na.rm = TRUE)
          } else if (startIntensity == "model") {
            NULL # Pass NULL to fit the intercept within getNLImod
          }
        })), o$NLI[, which(o$design$sample == p)]), 1, freeIntersect = freeIntersect, FUN = function(y, freeIntersect) {
          getNLImod(x, y[-1], ifelse(freeIntersect, NULL, y[1]))
        })))
    }
    
    colnames(d) <- paste(p, c("kloss", "kloss.stderr", "kloss.SSR", "nbpoints"), sep = ".")
    row.names(d) <- row.names(o$NLI)
    message(paste("...calculated ", sum(!is.na(d[, 1])), " k_loss values for sample ", p, 
                  " (", sum(is.na(d[, 1])), " missing)", sep = ""))
    
    if (is.null(o$NLI.kloss)) {
      o$NLI.kloss <- d
    } else {
      o$NLI.kloss <- cbind(o$NLI.kloss, d)
    }
  }
  
  if (ncores > 1) stopCluster(cl)
  
  return(o)
}

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

