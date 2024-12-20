% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcNLIkloss.R
\name{calcNLIkloss}
\alias{calcNLIkloss}
\title{Calculate k_loss for All Peptides Based on Normalized Light Intensities}
\usage{
calcNLIkloss(o, startIntensity = "max", ncores = NULL)
}
\arguments{
\item{o}{A `pSILAC` object containing the data to be analyzed.}

\item{startIntensity}{Character. Specifies the method for determining the starting intensity for model fitting:
\itemize{
  \item \code{"model"} - Fit the model with an intercept, allowing the initial intensity to be estimated during the fitting. This option is not functional in the current version. 
  \item \code{"max"} - Use the maximum of observed H+L intensities as the starting intensity (default).
  \item \code{"median"} - Use the median of non-null H+L intensities as the starting intensity. Consider this option if the data are not filtered using `filterMonotone`
}}

\item{ncores}{Integer. The number of cores to use for parallel computation. If \code{NULL} (default), the function detects the 
available cores and uses all but one. Single-core processing is used if \code{ncores = 1}.}
}
\value{
The updated `pSILAC` object with an additional `NLI.kloss` component. This component contains calculated values for each peptide, including:
\itemize{
  \item \code{kloss} - The decay rate constant.
  \item \code{kloss.stderr} - Standard error of the decay rate.
  \item \code{kloss.SSR} - Sum of squared residuals from the model fit.
  \item \code{nbpoints} - Number of points used in the fitting.
}
}
\description{
This function applies `getNLImod` to each peptide in the dataset, calculating the decay rate constant \( k_{loss} \) 
based on normalized light intensities (NLI). Users can specify how to set the starting intensity for the model, 
with options to use a model-fitted intercept, maximum observed intensity, or median of non-missing intensities.
}
