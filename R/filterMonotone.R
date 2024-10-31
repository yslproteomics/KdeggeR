#' Filter data based on monotone filtering - wrapper function
#'
#' @description A wrapper to apply the filter_monotone_trend_tp1excluded function to all pSILAC labeling series in the dataset.
#'
#' @param o a pSILAC object
#' @param values_cutoff Specifies the minimum number of valid values to keep per pSILAC labeling series. Default is 2.
#' @param skip_time_point Specifies how many early time points should be excluded from valid value filtering. Default is 1.
#' @importFrom dplyr filter select
#' @return A filtered data frame.
#' @export
filterMonotone <- function(o, values_cutoff = 2, skip_time_point = 1){
  
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
  message("Filtering the RIA data based on monotone decrease.")
  
  split_ria <- lapply(split_ria, filterMonotone_skipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point, mode = "RIA")
  
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
  message("Filtering the Ln H/L data based on monotone increase.")
  
  split_hol <- lapply(split_hol, filterMonotone_skipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point, mode = "HOL")
  
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
  message("Filtering the NLI data based on monotone decrease.")
  
  split_NLI <- lapply(split_NLI, filterMonotone_skipTimePoints, values_cutoff = values_cutoff, skip_time_point = skip_time_point, mode = "NLI")
  
  split_NLI <- lapply(split_NLI, function(x) { x$id <- row.names(x); return(x) })
  
  filter_NLI <- purrr::reduce(split_NLI, dplyr::full_join, by = "id")
  
  row.names(filter_NLI) <- filter_NLI$id
  
  o$NLI <- filter_NLI %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  message("NLI data filtering completed.")
  
  return(o)
}

#' Filter data based on monotone filtering
#'
#' The input is a data frame  frame split into a list based on grouping defined in the design object. Do everything with HOL, and then filter ther RIA data frane based on it?
#'
#' @param data A data frame containing values from one pSILAC labeling series.
#' @param values_cutoff Specifies the minimum number of valid values to keep per pSILAC labeling series. Default is 2.
#' @param skip_time_point Specifies how many early time points should be excluded from valid value filtering. Default is 1.
#' @param mode Specifies whether the filtering is performed with the "RIA", "HOL", or "NLI" values. 
#' @importFrom dplyr filter select
#' @return A filtered data frame.
#' @export
filterMonotone_skipTimePoints <- function(data, values_cutoff = 2, skip_time_point = 1, mode = NULL){
  
  # to start with time point number 2
  t <- skip_time_point + 1
  n <- ncol(data)
  
  # add na.omit to make the function robust against missing values
  data$valid_values <- apply(data[t:n], 1, function(x){sum(!is.na(x))})
  data$increasing_trend <- apply(data[t:n], 1, function(x){sum(diff(na.omit(x)) > 0)})
  data$decreasing_trend <- apply(data[t:n], 1, function(x){sum(-diff(na.omit(x)) > 0)})
  
  # only keep those with at least some values
  data <- data %>%
    dplyr::filter(valid_values > 0)
  
  # generate two vectors with 
  values <- unique(data$valid_values)
  trend <- values - 1
  
  # split the data based on the valid values for a custom filtering
  split_data_valid_values <- split.data.frame(data, f = as.factor(data$valid_values))
  split_data_filter_values <- vector(mode = "list", length = length(split_data_valid_values))
  
  if(mode == "RIA"){
    
    for(i in 1:length(split_data_valid_values)){
      
      split_data_filter_values[[i]] <- split_data_valid_values[[i]] %>%
        filter(valid_values == values[i] & decreasing_trend == trend[i])
      
    }
    
    data_filtered <- purrr::reduce(split_data_filter_values, rbind) %>%
      select(-valid_values, -increasing_trend, -decreasing_trend)
    
  } else if(mode == "NLI"){
    
    for(i in 1:length(split_data_valid_values)){
      
      split_data_filter_values[[i]] <- split_data_valid_values[[i]] %>%
        filter(valid_values == values[i] & decreasing_trend == trend[i])
      
    }
    
    data_filtered <- purrr::reduce(split_data_filter_values, rbind) %>%
      select(-valid_values, -increasing_trend, -decreasing_trend)
    
  } else if(mode == "HOL"){
    
    for(i in 1:length(split_data_valid_values)){
      
      split_data_filter_values[[i]] <- split_data_valid_values[[i]] %>%
        filter(valid_values == values[i] & increasing_trend == trend[i])
      
    }
    
    data_filtered <- purrr::reduce(split_data_filter_values, rbind) %>%
      select(-valid_values, -increasing_trend, -decreasing_trend)
    
  }
  
  return(data_filtered)
  
}


filter_monotone_trend_tp1only <- function(data, mode = NULL){
  
  if(is.null(mode)){
    stop("Please specify the datatype (ria, hol) to use the correct assumption")
  }
  
  # report the first value from the vector generated by diff function per row
  # in those with non-missing values in tp1, this values corresponds to the difference from the next value
  data$compare_timepoints <- apply(data, 1, function(x){diff(na.omit(x))[1]})
  data$isna_tp1 <- is.na(data[,1])
  
  data_tp1_missing <- data %>%
    filter(isna_tp1 == TRUE)
  
  data_tp1_decide <- data %>%
    filter(isna_tp1 == FALSE)
  
  if(mode == "RIA"){
    
    data_tp1_decide[,1] <- ifelse(data_tp1_decide$compare_timepoints < 0, data_tp1_decide[,1], NA)
    
    data_filtered <- data_tp1_decide %>%
      bind_rows(data_tp1_missing) %>%
      select(-compare_timepoints, -isna_tp1)
    
  } else if(mode == "HOL"){
    
    data_tp1_decide[,1] <- ifelse(data_tp1_decide$compare_timepoints > 0, data_tp1_decide[,1], NA)
    
    data_filtered <- data_tp1_decide %>%
      bind_rows(data_tp1_missing) %>%
      select(-compare_timepoints, -isna_tp1)
    
  }
  
  return(data_filtered)
  
}

