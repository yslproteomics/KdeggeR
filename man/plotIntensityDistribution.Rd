% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotIntensityDistribution}
\alias{plotIntensityDistribution}
\title{Plot Intensity Distribution of pSILAC Channels}
\usage{
plotIntensityDistribution(o, channel = "light", ...)
}
\arguments{
\item{o}{A `pSILAC` object containing intensity data.}

\item{channel}{A character string specifying the channel to plot. Options are `"heavy"`, `"light"`, `"sum"`, or `"NLI"`.
Default is `"light"`.}

\item{...}{Additional arguments passed to `gdensity` for customizing the plot.}
}
\description{
This function plots the log-transformed intensity distribution for a specified channel in a pSILAC object.
It allows visualization of intensity distributions across different conditions, such as time points.
}
\details{
The function uses the specified channel's intensity data from the pSILAC object `o`. If `channel = "sum"`, the sum of the
`heavy` and `light` channels is plotted. If `channel = "NLI"` and `NLI` data is not already in `o`, it will be calculated 
using `normalizeLightChannel`. The plot colors are determined by the time points in `o$design$time`.
}
\examples{
plotIntensityDistribution(o, channel = "light")

}
