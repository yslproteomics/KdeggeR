#' @title Monotone Filtering for the First Time Point in pSILAC Labeling Series
#'
#' @description
#' This function is a wrapper that applies monotone trend filtering to the first time point in each pSILAC labeling series 
#' within a dataset. It processes three data frames from a pSILAC object (`RIA`, `HOL`, and `NLI`) and ensures consistent 
#' column ordering after filtering.
#'
#' @param o a pSILAC object
#' @importFrom dplyr filter select all_of
#' @importFrom purrr reduce
#' @importFrom stats time
#' @return A filtered data frame.
#' @export
filterMonotoneTimePoint1 <- function(o){
  
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
  message(paste(Sys.time(), "Filtering the first time point of RIA data based on monotone decrease...", sep = " "))
  
  split_ria <- lapply(split_ria, filterMonotoneTp1Only, mode = "RIA")
  
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
  message(paste(Sys.time(), "Filtering the first time point of ln H/L data based on monotone increase...", sep = " "))
  
  split_hol <- lapply(split_hol, filterMonotoneTp1Only, mode = "HOL")
  
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
  message(paste(Sys.time(), "Filtering the first time point of NLI data based on monotone decrease...", sep = " "))
  
  split_NLI <- lapply(split_NLI, filterMonotoneTp1Only,  mode = "NLI")
  
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
  
  message(paste(Sys.time(), "Filtering completed", sep = " "))
  
  return(o)
}
