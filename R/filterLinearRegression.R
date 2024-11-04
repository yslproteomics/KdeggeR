#' @title Filter Data Based on Linear Regression of ln(H/L + 1) Values
#'
#' @description Applies linear regression filtering to the ln(H/L + 1) data in a `pSILAC` object.
#' This function performs linear regression across each labeling series in the dataset, removes specific 
#' data points based on goodness-of-fit and outlier thresholds, and returns an updated `pSILAC` object.
#'
#' @param o A `pSILAC` object containing the experimental data, specifically the `hol` and `ria` data frames.
#' @param skip_time_point Integer. Specifies how many early time points should be excluded from the linear regression analysis.
#' @param R2_cutoff R2 cutoff to remove the outliers (0.9 by default), only curves below this cutoff will be filtered
#' @param p_cutoff Grubb's test P value cutoff to remove the outliers (0.05 by default)
#' Must be 0 or 1. Default is 1.
#' @importFrom dplyr filter select full_join
#' @importFrom purrr reduce
#' @return A modified `pSILAC` object with filtered data. The `hol` and `ria` data frames are updated by removing rows 
#' (based on specific criteria in the linear regression results), and an additional `test.outliers.hol` list contains 
#' the regression model fits for each sample.
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Fits a linear model to each series (excluding an initial time point if specified) and stores regression metrics.
#'   \item Identifies data points where the full model R-squared is below 0.9 and the p-value from Grubbs' test (outlier test) is less than 0.05.
#'   \item Sets the selected data points to `NA` in both the `hol` and `ria` data frames and recombines the series into their original format.
#' }
#' 
#' @examples
#' # Assuming `pSILAC_obj` is a valid pSILAC object with appropriate data
#' filtered_pSILAC <- filterLinearRegression(pSILAC_obj, skip_time_point = 1)
#'
#' @export
filterLinearRegression <- function(o, skip_time_point = 1, R2_cutoff = 0.9, p_cutoff = 0.05){
  
  if (class(o) != "pSILAC") stop("Input data should be a pSILAC object.")
  
  if (skip_time_point > 1) stop("Error: 'skip_time_point' cannot be greater than 1.")
  
  # Define samples and extract data
  samples <- o$design$sample
  hol_data <- o$hol
  ria_data <- o$RIA
  original_order <- colnames(hol_data)
  time_points <- unique(o$design$time)
  
  # Split hol and ria data by sample, maintaining list structure
  split_hol <- lapply(split.default(hol_data, samples), as.data.frame)
  split_ria <- lapply(split.default(ria_data, samples), as.data.frame)
  
  # Perform linear model fitting
  message(paste(Sys.time(), "Performing the linear fit of the ln(H/L) data...", sep = " "))
  split_hol_fit <- lapply(split_hol, fitLinearModel, skip_time_point = skip_time_point, time_points = time_points)
  names(split_hol_fit) <- names(split_hol)
  
  # Filter data based on valid values
  message(paste(Sys.time(), "Performing the data filtering...", sep = " "))
    
    for (i in seq_along(split_hol_fit)) {
      data_fit <- split_hol_fit[[i]]
      data_hol <- split_hol[[i]]
      data_ria <- split_ria[[i]]
      
      # Filter rows based on conditions and get IDs to remove
      remove_ids <- data_fit %>%
        filter(time_point1 == TRUE, R2_full < R2_cutoff, grubbs_pval < p_cutoff) %>%
        rownames_to_column("id") %>%
        pull(id)
      
      # Print message on removed data points
      message(sprintf("Data points removed from tp1: %d, sample: %s", length(remove_ids), names(split_hol_fit)[i]))
    

    data_hol[rownames(data_hol) %in% remove_ids, 1] <- NA
    data_ria[rownames(data_ria) %in% remove_ids, 1] <- NA
    
    split_hol[[i]] <- data_hol
    split_ria[[i]] <- data_ria
    
    }

  
  split_hol <- lapply(split_hol, function(x) { x$id <- row.names(x); return(x) })
  
  filter_hol <- purrr::reduce(split_hol, dplyr::full_join, by = "id")
  
  row.names(filter_hol) <- filter_hol$id
  
  o$hol <- filter_hol %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  split_ria <- lapply(split_ria, function(x) { x$id <- row.names(x); return(x) })
  
  filter_ria <- purrr::reduce(split_ria, dplyr::full_join, by = "id")
  
  row.names(filter_ria) <- filter_ria$id
  
  o$RIA <- filter_ria %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  o$test.outliers.hol <- split_hol_fit
  
  return(o)
}

