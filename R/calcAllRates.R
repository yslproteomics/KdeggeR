#' Calculate Peptide and Protein k_loss Rates
#'
#' This function serves as a wrapper to sequentially call the functions `calcRIAkloss()`, `getHoLkloss()`, 
#' `calcNLIkloss()`, and `calcProteinsKloss()`, estimating k_loss rates for peptides and proteins within a pSILAC dataset.
#'
#' @param o A `pSILAC` object containing the data on which k_loss calculations will be performed.
#' @param method A character string specifying the method used for calculating protein-level k_loss rates. 
#'   Default is `combined`. For more options and details, see the documentation for `calcProteinsKloss`.
#' @param ag.metric A character string specifying the aggregation metric used for calculating protein-level rates.
#'   Default is `mean`. Additional options can be found in `calcProteinsKloss`.
#' @param ag.weights A character string defining the method for calculating weights in aggregation (default is `variance`). 
#'   This parameter is ignored if `ag.metric` is not set to `mean`. Refer to `calcProteinsKloss` for details.
#' @param in.all Logical or numeric indicating whether to use only peptides quantified in all samples.
#'   For more information, refer to `calcProteinsKloss`.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param tryRobust Logical, indicating if a robust fitting method should be used when applicable. Default is `FALSE`.
#' @param startIntensity Character string, specifying the initial intensity value used in calculations. Default is `max`.
#'
#' @return The function returns the updated `pSILAC` object with added k_loss estimates for both peptides and proteins.
#'   The updated object includes protein-level k_loss rates calculated according to the specified method and aggregation metric.
#'
#' @details
#' This function performs four main steps to estimate k_loss rates:
#' \enumerate{
#'   \item Peptide-wise k_loss rates are calculated using the RIA-based method via `calcRIAkloss`.
#'   \item Peptide-wise k_loss rates are calculated using the ln(H/L+1)-based method via `getHoLkloss`.
#'   \item Peptide-wise k_loss rates are calculated using the normalized light intensity (NLI) channel via `calcNLIkloss`.
#'   \item Protein-level k_loss rates are calculated using the specified aggregation settings via `calcProteinsKloss`.
#' }
#'
#' @export
calcAllRates <- function(o, method = "combined", ag.metric = "mean", ag.weights = "variance", in.all = 2, ncores = 1, tryRobust = FALSE, startIntensity = "max") {
  if (class(o) != "pSILAC") stop("o should be a pSILAC object.")
  
  message("Calculating peptide-wise k_loss using the RIA-based method...")
  o <- calcRIAkloss(o)
  
  message("Calculating peptide-wise k_loss using the ln(H/L+1)-based method...")
  o <- calcHoLkloss(o)
  
  message("Calculating peptide-wise k_loss using the normalized light channel...")
  o <- calcNLIkloss(o)
  
  message("Calculating protein-wise k_loss using specified aggregation method...")
  o <- calcProteinsKloss(o, method, ag.metric, ag.weights, in.all)
  
  # Store settings used in the calculation
  o$info$protk.call[["in.all"]] <- in.all
  
  return(o)
}
