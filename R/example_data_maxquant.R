#' @title Example pSILAC dataset processed with MaxQuant
#'
#' @description
#' An example dataset containing the first 5,000 unique peptides extracted
#' from a publicly available pSILAC study (PXD057850), as reported by
#' Frankenfield et al., *Molecular & Cellular Proteomics* (2025),
#' "Benchmarking SILAC Proteomics Workflows and Data Analysis Platforms"
#' (PMID: 40315959).
#'
#' The data originate from the MaxQuant results folder
#' \code{MaxQuant_dSILAC_CurveFitting.zip}, specifically the
#' \code{peptides.txt} file.
#'
#' The dataset contains a single experimental condition measured at
#' four time points, each with four biological replicates.
#'
#' @details
#' Mandatory columns in the data frame:
#' \itemize{
#'   \item \code{Sequence}: Unique peptide identifier (used as row names)
#'   \item \code{Proteins}: Protein identifiers
#'   \item \code{^Intensity\\.L\\..*}: Light SILAC channel intensities
#'   \item \code{^Intensity\\.H\\..*}: Heavy SILAC channel intensities
#' }
#'
#' @format
#' A data frame with 5,000 rows (peptides) and MaxQuant-derived intensity
#' columns suitable for pSILAC analysis.
#'
#' @source
#' ProteomeXchange accession PXD057850.
"example_maxquant"

#' @title Example design table for pSILAC analysis of MaxQuant data
#'
#' @description
#' An example design table used for pSILAC data analysis with replicate
#' aggregation. The table corresponds to a publicly available dataset
#' (ProteomeXchange accession PXD057850) reported by Frankenfield et al.,
#' *Molecular & Cellular Proteomics* (2025),
#' "Benchmarking SILAC Proteomics Workflows and Data Analysis Platforms"
#' (PMID: 40315959).
#'
#' The design table is compatible with MaxQuant output generated from
#' the \code{MaxQuant_dSILAC_CurveFitting.zip} results folder
#' (derived from the \code{peptides.txt} file).
#'
#' The experiment consists of a single condition measured at four time
#' points, each with four biological replicates.
#'
#' @format
#' A data frame defining sample-to-condition and time-point mappings
#' suitable for pSILAC analysis with replicate aggregation.
#'
#' @source
#' ProteomeXchange accession PXD057850.
"example_maxquant_design"
