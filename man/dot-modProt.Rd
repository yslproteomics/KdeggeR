% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelProteinsKloss.R
\name{.modProt}
\alias{.modProt}
\title{Models the protein-level kloss}
\usage{
.modProt(
  o,
  p,
  method,
  ag.weights,
  unique.weights,
  in.all,
  removeOutliers,
  tryRobust,
  freeIntersect = F,
  returnSD = F
)
}
\arguments{
\item{o}{a pSILAC object.}

\item{method}{the method used to calculate protein-level rates (default 'combined'); either 'RIA', 'hol' or 'combined' (uses both RIA- and H/L-based peptide kloss).}

\item{ag.weights}{the method to calculate weights for weighted mean (default 'variance'). Ignored if ag.metric != 'mean'. Either 'variance' (1/variance of the model's fit) or 'nbpoints' (number of datapoints used for the fit).}

\item{unique.weights}{whether unique weights should be used for each peptides (default T), or if false a weight is calculated individually for each sample.}

\item{in.all}{Starting with how many peptides should peptides not quantified in all samples be ignored (default 2), or FALSE to use all peptides, even if not quantified in all samples.}

\item{removeOutliers}{Starting with how many peptides should outlier peptides be removed, or FALSE to disable outlier removal (default 5). See ?removeOutlierPeptides}

\item{tryRobust}{Whether to try fitting a robust linear model when remodeling on the basis of H/L ratio (requires the 'MASS' library). Disabled by default.}

\item{freeIntersect}{Logical; whether to fit also the intersect for the NLI-based}

\item{returnSD}{Logical; whether to return also the SD of the kloss (default FALSE).}
}
\value{
A numeric vector with the protein's kloss across samples.
}
\description{
This is a subroutine modeling a single protein's kloss, called by modelProteinsKloss.
}
