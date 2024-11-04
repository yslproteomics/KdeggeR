#' Calculate the Rate of Isotope Loss (k_loss) Based on Normalized Light Intensities (NLI)
#'
#' This function models the rate of loss of a light isotope (`k_loss`) by fitting a non-linear model 
#' to normalized light intensity (NLI) values across specified timepoints. The model assumes exponential 
#' decay of NLI over time, allowing estimation of `k_loss`.
#'
#' @param x A numeric vector of timepoints associated with the NLI measurements.
#' @param y A numeric vector of normalized light intensity values at each timepoint in `x`. Must be the same length as `x`.
#' @param y0 The intercept of the model, representing the initial intensity. If `NULL` (default), it will be estimated during model fitting.
#' @param a A numeric value providing the initial guess for `k_loss` during model fitting. Default is 0.1.
#'
#' @return A numeric vector of length 4 containing the following elements:
#' \itemize{
#'   \item \code{k_loss}: The estimated rate of loss for the light isotope.
#'   \item \code{SE_k_loss}: The standard error of the estimated \code{k_loss}.
#'   \item \code{residual_sum}: The sum of residuals from the fitted model, indicating the goodness of fit.
#'   \item \code{n_points}: The number of data points used in the model, excluding any missing values in `y`.
#' }
#'
#' @details
#' This function applies a non-linear least squares method to fit an exponential decay model to the NLI values across timepoints.
#' If there is only one valid timepoint, a simplified calculation of `k_loss` is provided, using the intercept if available.
#' If model fitting fails, `NA` values are returned for `k_loss` and `SE_k_loss`, while `n_points` indicates the number of valid points.
#'
#' @export
getNLImod <- function(x, y, y0 = NULL, a = 0.1) {
  if (length(x) != length(y)) stop("x and y must be of the same length.")
  
  # Exponential decay function with parameters 'a' (k_loss) and 'b' (initial intensity)
  f <- function(x, a, b) { b * exp(-a * x) }
  
  # Remove NA values in y, adjust x accordingly
  y <- as.numeric(y)
  x2 <- x[!is.na(y)]
  y <- y[!is.na(y)]
  
  # If no or only one point is available, return NA or a simplified k_loss estimate
  if (length(y) < 1) return(c(NA, NA, NA, 1))
  if (length(y) == 1) {
    if (is.null(y0) || x2 == 0 || y0 == y) return(c(NA, NA, NA, 1))
    return(c(-log(y / y0) / x2, NA, NA, 1))
  }
  
  # Fit the non-linear model, handling cases with or without a specified intercept (y0)
  if (is.null(y0)) {
    fm <- try(suppressWarnings(nls(y ~ f(x2, a, b), start = c(a = a, b = max(y)))), silent = TRUE)
  } else {
    fm <- try(suppressWarnings(nls(y ~ f(x2, a, y0), start = c(a = a))), silent = TRUE)
  }
  
  # Return results if model fitting is successful
  if (inherits(fm, "nls")) {
    return(c(summary(fm)$coefficients[1:2], sum(residuals(fm)^2), length(x2)))
  }
  
  # If model fitting fails, return NA values for model parameters and residuals
  return(c(NA, NA, NA, length(x2)))
}
