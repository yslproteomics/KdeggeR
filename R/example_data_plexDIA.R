#' @title Example pSILAC dataset processed with DIA-NN 1.9
#'
#' @description
#' An example pSILAC dataset containing the first 5,000 unique precursors
#' from A2780 parental and A2780Cis cell lines. The experiment includes
#' three biological (dish) replicates measured at four labeling time points
#' (1, 4, 12, and 24 hours).
#'
#' The data were analyzed and exported using DIA-NN version 1.9.
#'
#' @details
#' Mandatory columns in the data frame:
#' \itemize{
#'   \item \code{Precursor.Id}: Unique precursor identifier (used as row names)
#'   \item \code{Protein.Group}: Protein group identifiers
#'   \item \code{.*\\.L$}: Light SILAC channel intensities
#'   \item \code{.*\\.H$}: Heavy SILAC channel intensities
#' }
#'
#' @format
#' A data frame of precursor-level intensities suitable for pSILAC analysis.
"example_diann"


#' @title Example design table for a pSILAC dataset analyzed with DIA-NN 1.9
#'
#' @description
#' An example design table used for pSILAC data analysis without replicate
#' aggregation. Each sample is treated independently.
#'
#' @format
#' A data frame defining sample-to-condition and labeling time-point mappings
#' for DIA-NN–derived pSILAC data.
"example_diann_design"
