#' @title Filter data based on valid values
#'
#' @description A wrapper to apply the filter_valid_skip_early_tp function to all pSILAC labeling series in the dataset.
#'
#' @param o a pSILAC object.
#' @param values_cutoff Specifies the minimum number of valid values to keep per pSILAC labeling series (default is 2).
#' @param skip_time_point Specifies how many early time points should be excluded from valid value filtering (default is 1). 
#' @importFrom dplyr select full_join
#' @importFrom purrr reduce
#' @return The updated pSILAC object with filtered data.
#' @export
filterValidValues <- function(o, values_cutoff = 2, skip_time_point = 1) {
  
  if (class(o) != "pSILAC") stop("Input data should be a pSILAC object.")
  
  # Define individual samples based on the design table
  samples <- unique(o$design$sample)
  
  #############################################################################
  # Use the RIA data frame from the pSILAC object
  ria_data <- o$RIA
  
  # Original column order
  original_order <- colnames(o$RIA)
  
  # Split data into individual pSILAC labeling series
  split_ria <- lapply(split.default(ria_data, samples), as.data.frame)
  
  # Filter data based on valid values
  message("Filtering the RIA data based on valid values.")
  
  split_ria <- lapply(split_ria, filterValidValues_skipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point)
  
  split_ria <- lapply(split_ria, function(x) { x$id <- row.names(x); return(x) })
  
  filter_ria <- purrr::reduce(split_ria, dplyr::full_join, by = "id")
  
  row.names(filter_ria) <- filter_ria$id

  o$RIA <- filter_ria %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))

  message("RIA filtering completed.")
  
  #############################################################################
  # Use the HOL data frame from the pSILAC object
  hol_data <- o$hol
  
  # Original column order
  original_order <- colnames(o$hol)
  
  # Split data into individual pSILAC labeling series
  split_hol <- lapply(split.default(hol_data, samples), as.data.frame)
  
  # Filter data based on valid values
  message("Filtering the Ln H/L data based on valid values.")
  
  split_hol <- lapply(split_hol, filterValidValues_skipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point)
  
  split_hol <- lapply(split_hol, function(x) { x$id <- row.names(x); return(x) })
  
  filter_hol <- purrr::reduce(split_hol, dplyr::full_join, by = "id")
  
  row.names(filter_hol) <- filter_hol$id
  
  o$hol <- filter_hol %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  message("Ln H/L data filtering completed.")
  
  #############################################################################
  # Use the NLI data frame from the pSILAC object
  NLI_data <- o$NLI
  
  # Original column order
  original_order <- colnames(o$NLI)
  
  # Split data into individual pSILAC labeling series
  split_NLI <- lapply(split.default(NLI_data, samples), as.data.frame)
  
  # Filter data based on valid values
  message("Filtering the NLI data based on valid values.")
  
  split_NLI <- lapply(split_NLI, filterValidValues_skipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point)
  
  split_NLI <- lapply(split_NLI, function(x) { x$id <- row.names(x); return(x) })
  
  filter_NLI <- purrr::reduce(split_NLI, dplyr::full_join, by = "id")
  
  row.names(filter_NLI) <- filter_NLI$id
  
  o$NLI <- filter_NLI %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  message("NLI data filtering completed.")
  
  return(o)
}

#' Filter data based on valid values, skipping early time points
#'
#' The input is a data frame  frame split into a list based on grouping defined in the design object.
#'
#' @param data A data frame containing values from one pSILAC labeling series.
#' @param values_cutoff Specifies the minimum number of valid values to keep per pSILAC labeling series. Default is 2.
#' @param skip_time_point Specifies how many early time points should be excluded from valid value filtering. Default is 1.
#' @importFrom dplyr filter select
#' @return a filtered data frame.
#' @export
filterValidValues_skipTimePoints <- function(data, values_cutoff = 2, skip_time_point = 1) {
  
  # Start with the specified time point
  t <- skip_time_point + 1
  n <- ncol(data)
  
  # Calculate the number of valid values for each row, ignoring NA values
  data$valid_values <- apply(data[, t:n], 1, function(x) { sum(!is.na(x)) })
  
  # Filter based on the required minimum number of valid values 
  data_filtered <- data %>%
    dplyr::filter(valid_values >= values_cutoff) %>%
    dplyr::select(-valid_values)
  
  # Return the filtered data
  return(data_filtered)
}