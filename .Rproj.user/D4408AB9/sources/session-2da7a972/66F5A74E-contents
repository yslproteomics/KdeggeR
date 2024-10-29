#' Filter data based on valid values, wrapper
#'
#' A wrapper to apply the filter_valid_skip_early_tp function to all pSILAC labeling series in the dataset.
#'
#' @param input_data A pSILAC object.
#' @param values_cutoff Specifies the minimum number of valid values to keep per pSILAC labeling series. Default is 2.
#' @param skip_time_point Specifies how many early time points should be excluded from valid value filtering. Default is 1.
#' @importFrom dplyr select full_join
#' @importFrom purrr reduce
#' @return The updated pSILAC object with filtered data.
#' @export
filter_valid <- function(input_data, values_cutoff = 2, skip_time_point = 1) {
  
  if (class(input_data) != "pSILAC") stop("Input data should be a pSILAC object.")
  
  # Define individual samples based on the design table
  samples <- unique(input_data$design$sample)
  
  # Use the RIA data frame from the pSILAC object
  ria_data <- input_data$RIA
  
  # Split data into individual pSILAC labeling series
  split_ria <- lapply(split.default(ria_data, samples), as.data.frame)
  
  # Filter data based on valid values
  message("Filtering the data based on valid values.")
  split_ria <- lapply(split_ria, filter_valid_skip_early_tp, values_cutoff = values_cutoff, skip_time_point = skip_time_point)
  
  split_ria <- lapply(split_ria, function(x) { x$id <- row.names(x); return(x) })
  
  filter_ria <- purrr::reduce(split_ria, dplyr::full_join, by = "id")
  
  row.names(filter_ria) <- filter_ria$id
  
  input_data$RIA <- filter_ria %>% dplyr::select(-id)
  
  message("Filtering completed.")
  
  return(input_data)
}

#' Filter data based on valid values, skipping early time points
#'
#' The input is a data frame  frame split into a list based on grouping defined in the design object.
#'
#' @param data A data frame containing values from one pSILAC labeling series.
#' @param values_cutoff Specifies the minimum number of valid values to keep per pSILAC labeling series. Default is 2.
#' @param skip_time_point Specifies how many early time points should be excluded from valid value filtering. Default is 1.
#' @importFrom dplyr filter select
#' @return A filtered data frame.
#' @export
filter_valid_skip_early_tp <- function(data, values_cutoff = 2, skip_time_point = 1) {
  
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

#' Filter data based on monotone filtering - combine into one step
#'
#' The input is a data frame  frame split into a list based on grouping defined in the design object. Do everything with HOL, and then filter ther RIA data frane based on it?
#'
#' @param data A data frame containing values from one pSILAC labeling series.
#' @param values_cutoff Specifies the minimum number of valid values to keep per pSILAC labeling series. Default is 2.
#' @param skip_time_point Specifies how many early time points should be excluded from valid value filtering. Default is 1.
#' @param mode Specifies whether the filtering is performed with the RIA or HOL values. 
#' @importFrom dplyr filter select
#' @return A filtered data frame.
#' @export
filter_monotone_trend_tp1excluded <- function(data, values_cutoff = 2, skip_time_point = 1, mode = NULL){
  
  if(is.null(mode)){
    stop("Please specify the datatype (ria, hol) to use the correct assumption")
  }
  
  # to start with time point number 2
  t <- skip_time_point + 1
  n <- ncol(data)
  
  # add na.omit to make the function robust against missing values
  data$valid_values <- apply(data[t:n], 1, function(x){sum(!is.na(x))})
  data$increasing_trend <- apply(data[t:n], 1, function(x){sum(diff(na.omit(x)) > 0)})
  data$decreasing_trend <- apply(data[t:n], 1, function(x){sum(-diff(na.omit(x)) > 0)})
  
  # generate two vectors with 
  values <- seq(values_cutoff,n-1,1)
  trend <- values - 1
  
  # split the data based on the valid values for a custiom filtering
  split_data_valid_values <- split.data.frame(data, f = as.factor(data$valid_values))
  split_data_filter_values <- vector(mode = "list", length = length(split_data_valid_values))
  
  if(mode == "ria"){
    
    for(i in 1:length(split_data_valid_values)){
      
      split_data_filter_values[[i]] <- split_data_valid_values[[i]] %>%
        filter(valid_values == values[i] & decreasing_trend == trend[i])
      
    }
    
    data_filtered <- purrr::reduce(split_data_filter_values, rbind) %>%
      select(-valid_values, -increasing_trend, -decreasing_trend)
    
  } else if(mode == "hol"){
    
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
  
  if(mode == "ria"){
    
    data_tp1_decide[,1] <- ifelse(data_tp1_decide$compare_timepoints < 0, data_tp1_decide[,1], NA)
    
    data_filtered <- data_tp1_decide %>%
      bind_rows(data_tp1_missing) %>%
      select(-compare_timepoints, -isna_tp1)
    
  } else if(mode == "hol"){
    
    data_tp1_decide[,1] <- ifelse(data_tp1_decide$compare_timepoints > 0, data_tp1_decide[,1], NA)
    
    data_filtered <- data_tp1_decide %>%
      bind_rows(data_tp1_missing) %>%
      select(-compare_timepoints, -isna_tp1)
    
  }
  
  return(data_filtered)
  
}

