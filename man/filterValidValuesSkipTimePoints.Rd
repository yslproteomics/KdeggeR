% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterValidValuesSkipTimePoints.R
\name{filterValidValuesSkipTimePoints}
\alias{filterValidValuesSkipTimePoints}
\title{Filter Data Based on Valid Values, Skipping Early Time Points}
\usage{
filterValidValuesSkipTimePoints(data, values_cutoff = 2, skip_time_point = 1)
}
\arguments{
\item{data}{A data frame containing values from one pSILAC labeling series.}

\item{values_cutoff}{Integer specifying the minimum number of valid values required to retain a row in the data frame. Default is 2.}

\item{skip_time_point}{Integer specifying the number of early time points to exclude from the valid value counting. For example, if set to 1, the first column will be skipped. Default is 1.}
}
\value{
A filtered data frame, with only rows meeting the minimum valid values threshold, starting from the specified time point.
}
\description{
This function filters rows in a data frame to retain only those with a specified minimum number of valid (non-NA) values. It skips early time points before applying the filtering criteria, based on a specified number.
}
\details{
The input data frame is grouped and filtered based on the `values_cutoff` parameter, which defines the minimum number of valid values required to keep a row. An optional parameter, `skip_time_point`, allows you to exclude early time points from this filtering process, effectively starting the valid value counting from a later time point.
}
