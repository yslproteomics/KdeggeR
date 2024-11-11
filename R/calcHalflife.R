#' @title Calculate Protein Halflives
#' 
#' @description This function calculates the halflives of proteins based on their degradation rates (`kdeg` values) 
#' from a pSILAC object. The calculated halflives are stored in the `prot.halflife` attribute of the object.
#' 
#' @param o A pSILAC object containing protein degradation data. The object must include `protein.kdeg` values, 
#' which are required for halflife calculation. If these values are not available, run `calcKdeg()` first.
#' 
#' @return The modified pSILAC object with an added `protein.halflife` matrix, representing the calculated protein halflives.
#' 
#' @export
calcHalflife <- function(o) {
  if (!inherits(o, "pSILAC")) 
    stop("o should be a pSILAC object.")
  
  if (is.null(o$protein.kdeg)) 
    stop("Protein kdeg must be calculated first. To do so, run `calcKdeg()`.")
  
  prot_kdeg <- o$protein.kdeg
  prot_halflife <- log(2) / prot_kdeg
  
  o$protein.halflife <- prot_halflife
  return(o)
}
