#' @title Filter Data Based on Valid Values
#'
#' @description Applies the `filterValidValues_skipTimePoints` function to filter each pSILAC labeling series in the dataset. 
#' This function processes multiple data frames within the pSILAC object (RIA, HOL, NLI) and filters data based on a specified 
#' minimum number of valid values, while optionally skipping a defined number of initial time points.
#'
#' @param o A pSILAC object.
#' @param values_cutoff Integer. The minimum number of valid values required to retain data within each pSILAC labeling series.
#' Default is 2.
#' @param skip_time_point Integer. The number of early time points to exclude from filtering based on valid values. Default is 1.
#' @importFrom dplyr select full_join
#' @importFrom purrr reduce
#'
#' @details This function performs filtering on three specific data frames (`RIA`, `HOL`, and `NLI`) within a pSILAC object. 
#'
#' @return The pSILAC object `o`, with the `RIA`, `HOL`, and `NLI` data frames updated to reflect the applied filtering.
#' @export
filterValidValues <- function(o, values_cutoff = 2, skip_time_point = 1) {
  
  if (class(o) != "pSILAC") stop("Input data should be a pSILAC object.")
  if (!is.numeric(values_cutoff) || values_cutoff < 0) stop("values_cutoff should be a non-negative integer.")
  if (!is.numeric(skip_time_point) || skip_time_point < 0) stop("skip_time_point should be a non-negative integer.")
  
  # Define individual samples based on the design table
  samples <- o$design$sample
  
  #############################################################################
  # Use the RIA data frame from the pSILAC object
  ria_data <- o$RIA
  
  # Original column order
  original_order <- colnames(o$RIA)
  
  # Split data into individual pSILAC labeling series
  split_ria <- lapply(split.default(ria_data, samples), as.data.frame)
  
  message(paste(Sys.time(), "Filtering the RIA data based on valid values...", sep = " "))
  
  split_ria <- lapply(split_ria, filterValidValuesSkipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point)
  
  split_ria <- lapply(split_ria, function(x) { x$id <- row.names(x); return(x) })
  
  filter_ria <- purrr::reduce(split_ria, dplyr::full_join, by = "id")
  
  row.names(filter_ria) <- filter_ria$id

  o$RIA <- filter_ria %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  #############################################################################
  # Use the HOL data frame from the pSILAC object
  hol_data <- o$hol
  
  # Original column order
  original_order <- colnames(o$hol)
  
  # Split data into individual pSILAC labeling series
  split_hol <- lapply(split.default(hol_data, samples), as.data.frame)
  
  # Filter data based on valid values
  message(paste(Sys.time(), "Filtering the ln(H/L+1) data based on valid values...", sep = " "))
  
  split_hol <- lapply(split_hol, filterValidValuesSkipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point)
  
  split_hol <- lapply(split_hol, function(x) { x$id <- row.names(x); return(x) })
  
  filter_hol <- purrr::reduce(split_hol, dplyr::full_join, by = "id")
  
  row.names(filter_hol) <- filter_hol$id
  
  o$hol <- filter_hol %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  #############################################################################
  # Use the NLI data frame from the pSILAC object
  NLI_data <- o$NLI
  
  # Original column order
  original_order <- colnames(o$NLI)
  
  # Split data into individual pSILAC labeling series
  split_NLI <- lapply(split.default(NLI_data, samples), as.data.frame)
  
  # Filter data based on valid values
  message(paste(Sys.time(), "Filtering the NLI data based on valid values...", sep = " "))
  
  split_NLI <- lapply(split_NLI, filterValidValuesSkipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point)
  
  split_NLI <- lapply(split_NLI, function(x) { x$id <- row.names(x); return(x) })
  
  filter_NLI <- purrr::reduce(split_NLI, dplyr::full_join, by = "id")
  
  row.names(filter_NLI) <- filter_NLI$id
  
  NLI <- filter_NLI %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  o$NLI <- NLI
  
  # Filter the NCS data frame based on the filtered NLI
  NCS_data <- o$NCS
  
  filter_NCS <- NCS_data[row.names(NLI), ]
  filter_NCS[is.na(NLI)] <- NA
  
  o$NCS <- filter_NCS
  
  message(paste(Sys.time(), "Valid value filtering completed.", sep = " "))
  
  return(o)
}
