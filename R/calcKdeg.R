#' @title Calculate Protein Degradation Rates
#' 
#' @description This function calculates protein degradation rates based on the protein-level `kloss` values from a pSILAC object. 
#' It supports two types of calculations: estimating theoretical degradation rates using a specified quantile (`kperc`) 
#' or using custom rates provided in a data frame (`kcd`).
#' 
#' @param o A pSILAC object containing protein degradation data and a `protein.kloss` matrix. 
#' The object must include `protein.kloss` values, which are required for degradation rate calculation.
#' 
#' @param rate_df An optional data frame containing custom degradation rate values (`kcd`) for each sample. 
#' This data frame should contain columns:
#' - `sample`: Sample names that must match with those in the pSILAC object.
#' - `kcd`: The custom degradation rate values for each sample.
#' Required only when `type = "kcd"`.
#' 
#' @param type A character string specifying the rate calculation type. 
#' Options are:
#' - `kperc`: Estimate theoretical kcd values based on the specified quantile (`perc_neg`).
#' - `kcd`: Use custom kcd values provided in `rate_df`.
#' Default is `kperc`.
#' 
#' @param perc_neg A numeric value specifying the quantile to use when estimating theoretical kcd values if `type = kperc`. 
#' Default is `0.01`, representing the 1st percentile.
#' 
#' @return A modified pSILAC object with an added `protein.kdeg` matrix, representing the calculated degradation rates for each protein.
#' The object may also include `kperc` or `kcd` data, depending on the calculation type used.
#' 
#' @importFrom dplyr select ends_with filter pull slice
#' @export
calcKdeg <- function(o, rate_df = NULL, type = "kperc", perc_neg = 0.01) {
  
  if (!inherits(o, "pSILAC")) 
    stop("o should be a pSILAC object.")
  
  if (is.null(o$protein.kloss)) 
    stop("Protein kloss must be calculated first.")
  
  prot_kloss <- o$protein.kloss %>%
    dplyr::select(-dplyr::ends_with(".stderr"), -dplyr::any_of("source"))
  
  names_kloss <- prot_kloss %>%
    colnames() %>%
    gsub(".kloss", "", ., fixed = TRUE)
  
  if (is.null(rate_df) & type == "kperc"){  
    message(paste(Sys.time(), "Estimating theoretical kcd values based on the selected assumption..."))
    
    kperc <- apply(prot_kloss, 2, quantile, probs = perc_neg, na.rm = TRUE)
    
    o$kperc <- data.frame(sample = names_kloss, kperc = kperc)
    
    message(paste(Sys.time(), "Applying theoretical kcd values to estimate kdeg..."))
    
    prot_kdeg <- sweep(prot_kloss, 2, kperc, FUN = "-")
    
  } else if(type == "kcd"){  
    
    if(is.null(rate_df))
      stop("Please provide a data frame with kcd values.")
    
    if(nrow(rate_df) != length(unique(o$design$sample)))  
      stop("Rates must be provided for every sample.")
    
    if(!identical(sort(rate_df$sample), sort(names_kloss)))
      stop("Sample names in rate_df do not match the expected samples.")
    
    kcd <- rate_df %>%
      dplyr::filter(sample %in% names_kloss) %>% 
      dplyr::slice(match(names_kloss, sample)) %>% 
      dplyr::pull(kcd)
    
    message(paste(Sys.time(), "Using the provided kcd values to estimate kdeg..."))
    
    prot_kdeg <- sweep(prot_kloss, 2, kcd, FUN = "-")
    
    o$kcd <- rate_df
  }
  
  prot_kdeg <- prot_kdeg %>%
    dplyr::rename_with(., cols = everything(), ~gsub(".kloss", ".kdeg", .))
  
  o$protein.kdeg <- prot_kdeg
  
  return(o)
}
