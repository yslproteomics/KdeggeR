#' @title Example pSILAC dataset analyzed in Spectronaut (v19 and newer)
#'
#' @description
#' An example pSILAC dataset containing the first 5,000 unique precursors
#' from A2780 parental and A2780Cis cell lines. The experiment includes
#' three biological (dish) replicates measured at four labeling time points
#' (1, 4, 12, and 24 hours).
#'
#' The data were analyzed and exported using Spectronaut 19.
#'
#' @details
#' Mandatory columns in the data frame:
#' \itemize{
#'   \item \code{EG.PrecursorId}: Unique precursor identifier (used as row names)
#'   \item \code{PG.ProteinGroups}: Protein group identifiers
#'   \item \code{.*EG.Channel1Quantity$}: Light SILAC channel intensities
#'   \item \code{.*EG.Channel2Quantity$}: Heavy SILAC channel intensities
#' }
#'
#' @format
#' A data frame of precursor-level intensities suitable for pSILAC analysis.
"example_spectronaut"


#' @title Example design table for a pSILAC dataset analyzed with Spectronaut
#'
#' @description
#' An example design table used for pSILAC data analysis without replicate
#' aggregation. Each sample is treated independently.
#'
#' @format
#' A data frame defining sample-to-condition and labeling time-point mappings
#' for Spectronaut-derived pSILAC data.
"example_spectronaut_design"


#' @title Example design table for a pSILAC dataset analyzed with Spectronaut (replicate design)
#'
#' @description
#' An example design table used for pSILAC data analysis with replicate
#' aggregation. Replicate numbers must be explicitly specified, and replicate
#' grouping must be encoded consistently in the sample identifiers.
#'
#' @details
#' When using replicate aggregation, the replicate structure must be reflected
#' in the \code{sample} column to ensure correct grouping during analysis.
#'
#' @format
#' A data frame defining sample, replicate, condition, and labeling time-point
#' mappings for Spectronaut-derived pSILAC data.
"example_spectronaut_design_replicates"

#' @title Example pSILAC dataset analyzed in Spectronaut using the inverted spike-in workflow (pre–Spectronaut 19)
#'
#' @description
#' An example pSILAC dataset containing the first 5,000 unique precursors
#' from A2780 parental and A2780Cis cell lines. The experiment includes
#' three biological (dish) replicates measured at four labeling time points
#' (1, 4, 12, and 24 hours).
#'
#' The data were analyzed using an inverted spike-in workflow in older
#' versions of Spectronaut (version 18 and earlier).
#'
#' @details
#' In Spectronaut versions 18 and earlier, SILAC channel intensities are
#' reported using different column naming conventions:
#' \itemize{
#'   \item \code{EG.PrecursorId}: Unique precursor identifier (used as row names)
#'   \item \code{PG.ProteinGroups}: Protein group identifiers
#'   \item \code{.*EG.ReferenceQuantity\\.\\.Settings\\.$}: Light SILAC channel intensities
#'   \item \code{.*EG.TargetQuantity\\.\\.Settings\\.$}: Heavy SILAC channel intensities
#' }
#'
#' For compatibility with downstream pSILAC analysis functions, these
#' columns should be renamed to the Spectronaut 19-style conventions
#' (e.g., \code{EG.Channel1Quantity} and \code{EG.Channel2Quantity}).
#'
#' An example renaming workflow is shown below.
#' 
#' @format
#' A data frame of precursor-level intensities suitable for pSILAC analysis.
#'
#' @examples
#' example_spectronaut_isw %>%
#'   dplyr::rename_with(
#'     ~ gsub("EG.TargetQuantity\\.\\.Settings\\.", "EG.Channel2Quantity", .),
#'     .cols = dplyr::ends_with("EG.TargetQuantity..Settings.")
#'   ) %>%
#'   dplyr::rename_with(
#'     ~ gsub("EG.ReferenceQuantity\\.\\.Settings\\.", "EG.Channel1Quantity", .),
#'     .cols = dplyr::ends_with("EG.ReferenceQuantity..Settings.")
#'   )
"example_spectronaut_isw"

