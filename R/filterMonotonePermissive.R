#' @title Filter pSILAC data based on a permissive monotone trend
#'
#' @description
#' Applies a more permissive monotone trend filter to all pSILAC labeling
#' series in a `pSILAC` object. The filter starts from the last time point
#' (assumed most accurate) and walks backwards in time, replacing violations
#' of the expected monotone trend with `NA`. After filtering, series with
#' fewer than `min_points` valid values are fully discarded (set to `NA`).
#'
#' @param o A `pSILAC` object containing the data to be filtered.
#' @param skip_time_point Integer. Number of earliest time points to **exclude**
#'   from monotone filtering. Default = 0 (all time points are checked).
#' @param min_points Integer. Minimal number of valid (non-`NA`) values
#'   required per series after filtering; series with fewer points are set to
#'   `NA`. Default is 3.
#'
#' @return The modified `pSILAC` object with filtered data frames `RIA`, `hol`,
#'   `NLI`, and `NCS`, maintaining the original column order.
#'
#' @keywords internal
.filterMonotonePermissive <- function(o,
                                     skip_time_point = 0L,
                                     min_points = 3L) {
  
  if (!inherits(o, "pSILAC")) stop("Input data should be a pSILAC object.")
  
  samples <- o$design$sample
  
  #############################################################################
  # RIA
  ria_data       <- o$RIA
  original_order <- colnames(o$RIA)
  split_ria      <- lapply(split.default(ria_data, samples), as.data.frame)
  
  message(paste(Sys.time(),
                "Filtering the RIA data based on permissive monotone decrease...",
                sep = " "))
  
  split_ria <- lapply(
    split_ria,
    filterMonotonePermissiveSkipTimePoints,
    skip_time_point = skip_time_point,
    mode            = "RIA",
    min_points      = min_points
  )
  
  split_ria <- lapply(split_ria, function(x) { x$id <- row.names(x); x })
  filter_ria <- purrr::reduce(split_ria, dplyr::full_join, by = "id")
  row.names(filter_ria) <- filter_ria$id
  
  o$RIA <- filter_ria %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  #############################################################################
  # HOL
  hol_data       <- o$hol
  original_order <- colnames(o$hol)
  split_hol      <- lapply(split.default(hol_data, samples), as.data.frame)
  
  message(paste(Sys.time(),
                "Filtering the ln(H/L+1) data based on permissive monotone increase...",
                sep = " "))
  
  split_hol <- lapply(
    split_hol,
    filterMonotonePermissiveSkipTimePoints,
    skip_time_point = skip_time_point,
    mode            = "HOL",
    min_points      = min_points
  )
  
  split_hol <- lapply(split_hol, function(x) { x$id <- row.names(x); x })
  filter_hol <- purrr::reduce(split_hol, dplyr::full_join, by = "id")
  row.names(filter_hol) <- filter_hol$id
  
  o$hol <- filter_hol %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  #############################################################################
  # NLI
  NLI_data       <- o$NLI
  original_order <- colnames(o$NLI)
  split_NLI      <- lapply(split.default(NLI_data, samples), as.data.frame)
  
  message(paste(Sys.time(),
                "Filtering the NLI data based on permissive monotone decrease...",
                sep = " "))
  
  split_NLI <- lapply(
    split_NLI,
    filterMonotonePermissiveSkipTimePoints,
    skip_time_point = skip_time_point,
    mode            = "NLI",
    min_points      = min_points
  )
  
  split_NLI <- lapply(split_NLI, function(x) { x$id <- row.names(x); x })
  filter_NLI <- purrr::reduce(split_NLI, dplyr::full_join, by = "id")
  row.names(filter_NLI) <- filter_NLI$id
  
  NLI <- filter_NLI %>%
    dplyr::select(-id) %>%
    dplyr::select(all_of(original_order))
  
  o$NLI <- NLI
  
  #############################################################################
  # NCS filtered to match NLI structure (as in your original function)
  NCS_data  <- o$NCS
  filter_NCS <- NCS_data[row.names(NLI), ]
  filter_NCS[is.na(NLI)] <- NA
  o$NCS <- filter_NCS
  
  message(paste(Sys.time(), "Permissive monotone trend filtering completed.", sep = " "))
  
  return(o)
}


# Helper: monotone from the last time point
#' @keywords internal
.filter_monotone_from_last <- function(x,
                                       direction = c("decreasing", "increasing"),
                                       skip_time_point = 0L,
                                       min_points = 3L) {
  direction <- match.arg(direction)
  n <- length(x)
  
  if (n == 0L) return(x)
  
  # Indices we actually apply the monotone rule to:
  # 1..skip_time_point are left unchanged (more permissive, early noisy tps).
  start_idx <- skip_time_point + 1L
  if (start_idx > n) {
    # Nothing to filter, just check min_points on the original vector
    if (sum(is.finite(x)) < min_points) x[] <- NA
    return(x)
  }
  
  y <- x
  
  # Find last non-NA value to anchor on
  last_idx <- max(which(is.finite(y)), na.rm = TRUE)
  if (!is.finite(last_idx)) {
    # all NA
    return(y)
  }
  
  ref_val <- y[last_idx]
  
  # Walk backwards from (last_idx - 1) down to start_idx
  for (i in seq(from = last_idx - 1L, to = start_idx, by = -1L)) {
    val <- y[i]
    if (!is.finite(val)) next
    
    if (direction == "decreasing") {
      # Expected: earlier >= later; violation if val < ref_val
      if (val < ref_val) {
        y[i] <- NA
      } else {
        ref_val <- val
      }
    } else { # direction == "increasing"
      # Expected: earlier <= later; violation if val > ref_val
      if (val > ref_val) {
        y[i] <- NA
      } else {
        ref_val <- val
      }
    }
  }
  
  # Check minimal number of valid points AFTER filtering
  if (sum(is.finite(y)) < min_points) {
    y[] <- NA
  }
  
  return(y)
}


# Per sample helper
#' @keywords internal
filterMonotonePermissiveSkipTimePoints <- function(df,
                                                   skip_time_point = 0L,
                                                   mode = c("RIA", "HOL", "NLI"),
                                                   min_points = 3L) {
  mode <- match.arg(mode)
  direction <- switch(mode,
                      "RIA" = "decreasing",
                      "NLI" = "decreasing",
                      "HOL" = "increasing")
  
  if (!is.data.frame(df)) df <- as.data.frame(df)
  
  res <- t(apply(df, 1, .filter_monotone_from_last,
                 direction = direction,
                 skip_time_point = skip_time_point,
                 min_points = min_points))
  
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  colnames(res) <- colnames(df)
  rownames(res) <- rownames(df)
  res
}
