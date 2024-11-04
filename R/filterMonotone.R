#' @title Filter pSILAC Data Based on Monotonic Trend
#'
#' @description This function is a wrapper for applying the `filter_monotone_trend_tp1excluded` function 
#' to all pSILAC labeling series in the input dataset. It filters data frames within a pSILAC object (`RIA`, `hol`, and `NLI`)
#' according to a monotonic trend, excluding specified early time points.
#'
#' @param o A `pSILAC` object containing the data to be filtered.
#' @param skip_time_point Integer. Specifies the number of early time points to exclude from the valid value filtering process. 
#' Default is 1.
#'
#' @importFrom dplyr filter select
#' @importFrom purrr reduce
#' @return The modified `pSILAC` object with filtered data frames `RIA`, `hol`, and `NLI`, maintaining the original column order.
#'
#' @export
filterMonotone <- function(o, skip_time_point = 1){
  
  if (class(o) != "pSILAC") stop("Input data should be a pSILAC object.")
  
  # Define individual samples based on the design table
  samples <- o$design$sample
  
  #############################################################################
  # Use the RIA data frame from the pSILAC object
  ria_data <- o$RIA
  
  # Original column order
  original_order <- colnames(o$RIA)
  
  # Split data into individual pSILAC labeling series
  split_ria <- lapply(split.default(ria_data, samples), as.data.frame)
  
  # Filter data based on valid values
  message(paste(Sys.time(), "Filtering the RIA data based on monotone decrease...", sep = " "))
  
  split_ria <- lapply(split_ria, filterMonotoneSkipTimePoints, skip_time_point = skip_time_point, mode = "RIA")
  
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
  message(paste(Sys.time(), "Filtering the ln(H/L+1) data based on monotone increase...", sep = " "))
  
  split_hol <- lapply(split_hol, filterMonotoneSkipTimePoints, skip_time_point = skip_time_point, mode = "HOL")
  
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
  message(paste(Sys.time(), "Filtering the NLI data based on monotone decrease...", sep = " "))
  
  split_NLI <- lapply(split_NLI, filterMonotoneSkipTimePoints, skip_time_point = skip_time_point, mode = "NLI")
  
  split_NLI <- lapply(split_NLI, function(x) { x$id <- row.names(x); return(x) })
  
  filter_NLI <- purrr::reduce(split_NLI, dplyr::full_join, by = "id")
  
  row.names(filter_NLI) <- filter_NLI$id
  
  o$NLI <- filter_NLI %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  message(paste(Sys.time(), "Monotone trend filtering completed.", sep = " "))
  
  return(o)
}
