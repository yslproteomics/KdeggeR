#' @title Calculate Protein Halflives
#'
#' @description
#' This function calculates the halflives of proteins based on their degradation
#' rates (`kdeg` values) from a pSILAC object. The calculated halflives are stored
#' in the `protein.halflife` attribute of the object.
#'
#' @param o
#'   A pSILAC object containing protein degradation data. The object must include
#'   `protein.kdeg` values, which are required for halflife calculation. If these
#'   values are not available, run `calcKdeg()` first.
#'
#' @return
#'   The modified pSILAC object with an added `protein.halflife` matrix,
#'   representing the calculated protein halflives. Any kdeg values <= 0 are set
#'   to NA prior to computing half-life.
#'
#' @export
calcHalflife <- function(o) {
  if (!inherits(o, "pSILAC"))
    stop("o should be a pSILAC object.")
  
  if (is.null(o$protein.kdeg))
    stop("Protein kdeg must be calculated first. To do so, run `calcKdeg()`.")
  
  prot_kdeg <- o$protein.kdeg
  
  mat <- as.matrix(prot_kdeg)
  n_bad <- sum(mat <= 0, na.rm = TRUE)
  if (n_bad > 0L) {
    message(paste(Sys.time(), "Replacing", n_bad, "kdeg values <= 0 with NA prior to half-life calculation."))
    mat[mat <= 0] <- NA_real_
  }
  
  prot_halflife <- log(2) / mat
  prot_halflife <- as.data.frame(prot_halflife, stringsAsFactors = FALSE)
  rownames(prot_halflife) <- rownames(prot_kdeg)
  colnames(prot_halflife) <- colnames(prot_kdeg)
  
  prot_halflife <- prot_halflife %>%
    dplyr::rename_with(cols = dplyr::everything(), ~ gsub(".kdeg", ".t1_2", ., fixed = TRUE))
  
  o$protein.halflife <- prot_halflife
  return(o)
}
