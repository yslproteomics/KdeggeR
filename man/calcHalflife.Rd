% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcHalflife.R
\name{calcHalflife}
\alias{calcHalflife}
\title{Calculate Protein Halflives}
\usage{
calcHalflife(o)
}
\arguments{
\item{o}{A pSILAC object containing protein degradation data. The object must include `protein.kdeg` values, 
which are required for halflife calculation. If these values are not available, run `calcKdeg()` first.}
}
\value{
The modified pSILAC object with an added `protein.halflife` matrix, representing the calculated protein halflives.
}
\description{
This function calculates the halflives of proteins based on their degradation rates (`kdeg` values) 
from a pSILAC object. The calculated halflives are stored in the `prot.halflife` attribute of the object.
}
