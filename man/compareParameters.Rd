% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kloss.R
\name{compareParameters}
\alias{compareParameters}
\title{compareParameters}
\usage{
compareParameters(o, replicates = NULL, ncores = NULL, m = NULL)
}
\arguments{
\item{o}{a pSILAC object, with all peptide-level rates calculated (see ?calcAllRates)}

\item{replicates}{A vector of length equal to the number of samples, or a character designating the column in the design data.frame, indicating which samples are replicates of the same condition. This will be used to calculate distance/correlation between replicates.}

\item{ncores}{The number of cores to use. If NULL, will default to detected cores minus 1.}

\item{m}{A list of the comparisons to make.}
}
\value{
A list of the results, each with the following slots:
'call': the function and parameters called.
'kloss': the protein kloss table
'params': the list of input parameters
'cor': the pearson correlations between replicates.
'dist': the euclidean distance between replicates.
'MdAE': the median absolute error between replicates.
}
\description{
Runs a series of protein-level kloss calculations and compares the results. See ?compareParametersResults to plot the results.
}
