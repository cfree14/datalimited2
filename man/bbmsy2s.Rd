% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s2bbmsy.R
\name{bbmsy2s}
\alias{bbmsy2s}
\title{Convert B/BMSY to saturation}
\usage{
bbmsy2s(bbmsy)
}
\arguments{
\item{bbmsy}{A vector of B/BMSY values}
}
\value{
A vector of saturation values
}
\description{
Converts B/BMSY to saturation (S). Note: S = 1 - depletion = B / K = 0.5 * B/BMSY.
}
\examples{
# Convert B/BMSY of 0.5 to saturation
bbmsy2s(0.5)
}
