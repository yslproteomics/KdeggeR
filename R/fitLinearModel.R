#' @title Filter Data Based on Monotone Linear Regression
#'
#' @description This function fits a linear regression model to each row in a pSILAC dataset, 
#' optionally excluding early time points. The function filters data based on specified criteria, 
#' including a minimum number of valid values and trends in the data. It returns model summary 
#' statistics, including R-squared values, slopes, intercepts, and p-values from Grubbs' test to 
#' identify outliers.
#'
#' @param data A data frame where each row represents a pSILAC labeling series, with columns 
#' representing time points.
#' @param skip_time_point Integer specifying the number of early time points to exclude from 
#' the analysis. Must be 0 or 1; default is 1.
#' @param time_points A numeric vector indicating the time points corresponding to each column 
#' in `data`. The length must match the number of columns in `data`.
#' 
#' @details
#' This function iteratively fits a linear model to each row in `data`, focusing on the `ln(H/L + 1)` 
#' transformation of pSILAC values. Two models are fit:
#' - **Reduced model**: Excludes the first time point if `skip_time_point` is set to 1.
#' - **Full model**: Includes all time points, if the first time point has a valid value.
#'
#' After fitting, each row is evaluated based on:
#' - R-squared values from both reduced and full models (`R2` and `R2_full`).
#' - Estimated intercept and slope for both models.
#' - Grubbs' test p-value (`grubbs_pval`), which identifies outliers at the first time point 
#'   if it has the largest residual in absolute terms.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{R2}{R-squared of the reduced model, excluding the first time point if skipped.}
#'   \item{intercept}{Intercept of the reduced model.}
#'   \item{slope}{Slope of the reduced model.}
#'   \item{grubbs_pval}{P-value from Grubbs' test for the residual at the first time point, indicating 
#'                      potential outliers if it has the maximum residual.}
#'   \item{time_point1}{Logical, indicating if time point 1 was included in the full model.}
#'   \item{R2_full}{R-squared of the full model, including all time points if time point 1 has data.}
#'   \item{intercept_full}{Intercept of the full model.}
#'   \item{slope_full}{Slope of the full model.}
#' }
#' 
#' @importFrom dplyr filter select
#' @importFrom outliers grubbs.test
#'
#' @export
fitLinearModel <- function(data, skip_time_point = 1, time_points = time_points){
  
  # to start with time point number 2
  t <- skip_time_point + 1
  n <- ncol(data)
  
  # add na.omit to make the function robust against missing values
  data$valid_values <- apply(data[t:n], 1, function(x){sum(!is.na(x))})
  
  # for linear fitting, only select those with at least valid values
  data_filtered <- data %>%
    dplyr::filter(valid_values >= 3) %>%
    dplyr::select(-valid_values)
  
  # transform the data into a list of numeric vectors
  data_to_numeric <- suppressWarnings(split(data_filtered, row(data_filtered)))
  data_to_numeric <- suppressWarnings(lapply(data_to_numeric, as.numeric))
  names(data_to_numeric) <- row.names(data_filtered)
  
  n_precursors <- length(data_to_numeric)
  
  # initialize a result data.frame
  results_linear_fit <- data.frame(
    R2 = vector(mode = "numeric", length = n_precursors), 
    intercept = vector(mode = "numeric", length = n_precursors), 
    slope = vector(mode = "numeric", length = n_precursors), 
    grubbs_pval = vector(mode = "numeric", length = n_precursors),
    time_point1 = vector(mode = "logical", length = n_precursors), 
    R2_full = vector(mode = "numeric", length = n_precursors), 
    intercept_full = vector(mode = "numeric", length = n_precursors), 
    slope_full = vector(mode = "numeric", length = n_precursors)
  )
  
  residuals_list <- vector(mode = "list", length = n_precursors)
  
  row.names(results_linear_fit) <- row.names(data_filtered)
  
  for(i in 1:n_precursors){
    
    # remove the first time point
    time_points_later <- time_points[-1]
    
    # select observed values
    observed_y <- data_to_numeric[[i]]
    
    # prepare a data frame for model fitting, reduced model
    data_fit <- data.frame(time = time_points_later, 
                           lnHL = observed_y[-1])
    
    # fit linear model, reduced model
    # it can tolerate missing values, but does not report NA for the missing residuals
    model_fit <- lm(lnHL~time, data_fit)
    
    # model summary
    model_summary <- summary(model_fit)
    
    # export reduced model residuals
    residuals <- model_fit$residuals
    
    # rsqr <- model_summary$adj.r.squared
    rsqr <- model_summary$r.squared
    
    # estimated model coefficients, reduced model
    b0 <- as.numeric(model_fit$coefficients[1])
    b1 <- as.numeric(model_fit$coefficients[2])
    
    results_linear_fit$R2[i] <- rsqr
    results_linear_fit$intercept[i] <- b0
    results_linear_fit$slope[i] <- b1
    
    # data for the full model, not skipping time point 1, if there is tp1
    if(!is.na(observed_y[1])){
      
      data_fit_full <- data.frame(time = time_points, 
                                  lnHL = observed_y)
      
      # full model
      model_fit_full <- lm(lnHL~time, data_fit_full)
      model_summary_full <- summary(model_fit_full)
      
      # rsqr <- model_summary$adj.r.squared
      rsqr_full <- model_summary_full$r.squared
      
      # estimated model coefficients, reduced model
      b0_full <- as.numeric(model_fit_full$coefficients[1])
      b1_full <- as.numeric(model_fit_full$coefficients[2])
      
      results_linear_fit$time_point1[i] <- TRUE
      results_linear_fit$R2_full[i] <- rsqr_full
      results_linear_fit$intercept_full[i] <- b0_full
      results_linear_fit$slope_full[i] <- b1_full
      
    } else{
      
      results_linear_fit$time_point1[i] <- FALSE
      results_linear_fit$R2_full[i] <- rsqr
      results_linear_fit$intercept_full[i] <- b0
      results_linear_fit$slope_full[i] <- b1
      
    }
    
    # Residual for time point 1 if present
    tp1 <- time_points[1]
    observed_y_tp1 <- observed_y[1]
    
    if(!is.na(observed_y_tp1)){
      
      predicted_y_tp1 <- b0 + b1*tp1
      
      # what is the residual of the first time point data as compared to the fit
      residual_tp1 <- observed_y_tp1 - predicted_y_tp1
      
      # combine them together
      residuals_combined <- c(residual_tp1, residuals)
      
      residuals_list[[i]] <- residuals_combined
      
      # Check if time point 1 has the maximum residual
      tp1_is_max_residual <- sum(abs(residual_tp1) > abs(residuals)) == length(residuals)
      
      # Grubbs' test only if time point 1 has the maximum residual
      if(tp1_is_max_residual){
        
        outlier.test <- outliers::grubbs.test(residuals_combined, two.sided = FALSE, opposite = FALSE)
        results_linear_fit$grubbs_pval[i] <- outlier.test$p.value
        
      } else {
        
        results_linear_fit$grubbs_pval[i] <- NA
      }
      
    } else{
      
      # Handle missing time point 1 case
      residuals_combined <- c(NA, residuals)
      residuals_list[[i]] <- residuals_combined
      results_linear_fit$grubbs_pval[i] <- NA
      
    }
    
  }
  
  return(results_linear_fit)
  
}