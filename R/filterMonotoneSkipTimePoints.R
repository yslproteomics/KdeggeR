#' @title Filter Data Based on Monotonic Trend with Skip Time Points
#'
#' @description This function filters a single pSILAC labeling series by applying a monotonic trend 
#' filter to data columns, with an option to exclude a specified number of early time points from the filtering process.
#' The input data frame is split into subsets based on the number of valid values in each row.
#'
#' @param data A data frame containing values for one pSILAC labeling series. Each row should represent an individual observation, and each column represents a time point.
#' @param skip_time_point Integer. Specifies the number of early time points to exclude from the filtering process. Default is 1.
#' @param mode Character. Specifies the type of monotonic trend filtering to apply:
#' \itemize{
#'   \item \code{"RIA"} - filters for a monotonic decrease.
#'   \item \code{"NLI"} - filters for a monotonic decrease.
#'   \item \code{"HOL"} - filters for a monotonic increase.
#' }
#' @importFrom dplyr filter select
#' @importFrom purrr reduce
#' @return A filtered data frame with rows that meet the specified monotonic trend criteria. Columns `valid_values`, `increasing_trend`, and `decreasing_trend` are removed from the output.
filterMonotoneSkipTimePoints <- function(data, skip_time_point = 1, mode = NULL){
  
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
  
  # generate two values vectors
  values <- sort(unique(data$valid_values))
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
