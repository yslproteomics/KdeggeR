#' Filter Data Based on Valid Values, Skipping Early Time Points
#'
#' This function filters rows in a data frame to retain only those with a specified minimum number of valid (non-NA) values. It skips early time points before applying the filtering criteria, based on a specified number.
#'
#' The input data frame is grouped and filtered based on the `values_cutoff` parameter, which defines the minimum number of valid values required to keep a row. An optional parameter, `skip_time_point`, allows you to exclude early time points from this filtering process, effectively starting the valid value counting from a later time point.
#'
#' @param data A data frame containing values from one pSILAC labeling series.
#' @param values_cutoff Integer specifying the minimum number of valid values required to retain a row in the data frame. Default is 2.
#' @param skip_time_point Integer specifying the number of early time points to exclude from the valid value counting. For example, if set to 1, the first column will be skipped. Default is 1.
#' @importFrom dplyr filter select
#' @return A filtered data frame, with only rows meeting the minimum valid values threshold, starting from the specified time point.
filterValidValuesSkipTimePoints <- function(data, values_cutoff = 2, skip_time_point = 1) {
  # Check that the values_cutoff is within a valid range
  if (values_cutoff > ncol(data)) {
    stop("Error: 'values_cutoff' cannot exceed the number of columns in 'data'.")
  }
  
  # Define the starting column index after skipping the specified time points
  t <- skip_time_point + 1
  n <- ncol(data)
  
  # Calculate the number of valid values per row, ignoring NA values in skipped columns
  data$valid_values <- apply(data[, t:n], 1, function(x) { sum(!is.na(x)) })
  
  # Filter rows based on the minimum number of valid values
  data_filtered <- data %>%
    dplyr::filter(valid_values >= values_cutoff) %>%
    dplyr::select(-valid_values)
  
  # Return the filtered data
  return(data_filtered)
}
