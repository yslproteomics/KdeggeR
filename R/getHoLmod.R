#' Calculate H/L-based k_loss
#'
#' This function models the rate of loss of the light isotope (k_loss) by fitting a linear model to the log-transformed isotope ratio (ln(H/L + 1)) over specified time points.
#' Optionally, a robust linear model can be used if there are at least 10 data points.
#'
#' @param x A numeric vector of time points.
#' @param y A numeric vector of ln(H/L + 1) values corresponding to the time points in `x`.
#' @param tryRobust Logical; if `TRUE`, attempts to fit a robust linear model (requires the `MASS` package). Ignored if fewer than 10 data points are available.
#'
#' @return A numeric vector of length 5 containing:
#'   \itemize{
#'     \item The k_loss, calculated as the slope of ln(H/L + 1) over time.
#'     \item The standard error of the slope.
#'     \item The sum of squared residuals from the linear model.
#'     \item The R-squared value of the fit.
#'     \item The number of data points used to build the model.
#'   }
#'
#' @details
#' The function first removes any missing values from `y` and their corresponding time points in `x`.
#' If only one data point remains, the function returns the ratio of `y` to `x` if `x` equals zero, otherwise returning `NA` values for the slope and associated statistics.
#' For datasets with at least 10 data points and `tryRobust = TRUE`, a robust linear model is attempted, falling back to a standard linear model if the robust fit fails.
#'
#' @export
getHoLmod <- function(x, y, tryRobust = FALSE) {
  if (length(x) != length(y)) stop("x and y are of different lengths!")
  
  # Load the MASS package for the robust linear model
  if (tryRobust) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("The 'MASS' package is required for robust linear modeling. Please install it.")
    }
  }
  
  y <- as.numeric(y)
  x2 <- x[!is.na(y)]
  y <- y[!is.na(y)]
  
  if (length(y) == 1) {
    if (x2 == 0) return(c(NA, NA, NA, NA, 1))
    return(c(y / x2, NA, NA, NA, 1))
  }
  
  if (tryRobust && length(y) >= 10) {
    fm <- try(MASS::rlm(y ~ x2 + 0), silent = TRUE)
    if (!"lm" %in% class(fm)) fm <- try(lm(y ~ x2 + 0), silent = TRUE)
  } else {
    fm <- try(lm(y ~ x2 + 0), silent = TRUE)
  }
  
  if ("lm" %in% class(fm)) {
    summary_fm <- summary(fm)
    return(c(
      as.numeric(summary_fm$coefficients[1:2]),  # k_loss and standard error
      sum(residuals(fm)^2),                      # Sum of squared residuals
      summary_fm$r.squared,                      # R-squared
      length(x2)                                 # Number of data points
    ))
  }
  
  return(c(NA, NA, NA, NA, length(x2)))
}


