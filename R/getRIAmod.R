#' Calculate the Rate of Isotope Loss (k_loss) Based on Relative Isotope Abundance (RIA)
#'
#' This function models the rate of loss for a light isotope (`k_loss`) using a non-linear model based on 
#' the relative isotope abundance (RIA) values across specified timepoints. The model assumes an exponential decay 
#' of RIA over time, allowing estimation of `k_loss`.
#'
#' @param x A numeric vector of timepoints corresponding to the measurements of RIA.
#' @param y A numeric vector of RIA values associated with each timepoint in `x`. Must be the same length as `x`.
#' @param a A numeric value specifying the initial guess for `k_loss` during model fitting. Default is 0.1.
#'
#' @return A numeric vector of length 4 containing the following elements:
#' \itemize{
#'   \item \code{k_loss}: The estimated rate of loss for the light isotope.
#'   \item \code{SE_k_loss}: The standard error of the estimated \code{k_loss}.
#'   \item \code{residual_sum}: The sum of residuals from the fitted model, indicating goodness of fit.
#'   \item \code{n_points}: The number of data points used in the model, excluding any missing values in `y`.
#' }
#'
#' @details
#' This function uses non-linear least squares to fit an exponential decay model to the RIA values over time.
#' If only one valid timepoint is available, a simplified calculation is used. If the model fitting fails,
#' `NA` values are returned for \code{k_loss} and \code{SE_k_loss}, while the number of valid points (`n_points`) is provided.
#'
#' @export
getRIAmod <- function(x, y, a = 0.1) {
  if (length(x) != length(y)) stop("x and y must be of the same length.")
  
  # Exponential decay function with k_loss parameter 'a'
  f <- function(x, a) { exp(-a * x) }
  
  # Remove NA values in y, adjust x accordingly
  y <- as.numeric(y)
  x2 <- x[!is.na(y)]
  y <- y[!is.na(y)]
  
  # If only one point is available, return a simplified k_loss estimate
  if (length(y) == 1) {
    if (x2 == 0) return(c(NA, NA, NA, 1))
    return(c(-log(y) / x2, NA, NA, 1))
  }
  
  # Fit the non-linear model and handle any fitting errors
  fm <- try(suppressWarnings(nls(y ~ f(x2, a), start = c(a = a))), silent = TRUE)
  
  # Return results if model fitting is successful
  if (inherits(fm, "nls")) {
    return(c(summary(fm)$coefficients[1:2], sum(residuals(fm)^2), length(x2)))
  }
  
  # If model fitting fails, return NA values for model parameters and residuals
  return(c(NA, NA, NA, length(x2)))
}
