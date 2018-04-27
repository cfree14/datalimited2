#' Irish Sea Common sole time series
#'
#' A dataset containing the catch and biomass time series for Irish Sea Common sole
#' (\emph{Solea solea}) from 1970-2014. This stock was used as an example in the CMSY/BSM
#' user manual and is used to validate the \pkg{datalimited2} package's implementation of both CMSY and BSM.
#'
#' @format A data frame with 45 rows (years) and 4 variables:
#' \describe{
#'   \item{Stock}{Stock id}
#'   \item{yr}{Year}
#'   \item{ct}{Catch, in metric tons}
#'   \item{bt}{Biomass, in metric tons}
#' }
#' @source Froese R, Demirel N, Coro G, Kleisner KM, Winker H (2017)
#' Estimating fisheries reference points from catch and resilience. \emph{Fish and Fisheries} 18(3): 506-526.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/faf.12190/abstract}
"SOLIRIS"

#' SE Australia Tiger flathead time series
#'
#' A dataset containing the catch time series for Southeast Australia Tiger flathead
#' (\emph{Neoplatycephalus richardsoni}) from 1915-2012. This stock was used as an example in the OCOM
#' paper and is used to validate the \pkg{datalimited2} package's implementation of zBRT and OCOM.
#'
#' @format A data frame with 98 rows (years) and 3 variables:
#' \describe{
#'   \item{stock}{Stock id}
#'   \item{yr}{Year}
#'   \item{catch}{Catch, in metric tons}
#' }
#' @source Zhou S, Punt AE, Smith ADM, Ye Y, Haddon M, Dichmont CM, Smith DC
#' (2017) An optimised catch-only assessment method for data poor fisheries.
#' \emph{ICES Journal of Marine Science}: doi:10.1093/icesjms/fsx226.
#' \url{https://doi.org/10.1093/icesjms/fsx226}
"TIGERFLAT"

#' USA SNE/MA Yellowtail flounder time series
#'
#' A dataset containing the catch and biomass time series for USA Southern
#' New England/Mid-Atlantic (SNE/MA) Yellowtail flounder (\emph{Pleuronectes ferruginea}) from 1973-2014.
#' This dataset is included for additional testing of the \pkg{datalimited2} package.
#'
#' @format A data frame with 42 rows (years) and 7 variables:
#' \describe{
#'   \item{stockid}{Stock id}
#'   \item{year}{Year}
#'   \item{catch}{Catch, in metric tons}
#'   \item{biomass}{Spawning stock biomass (SSB), in metric tons}
#'   \item{f}{Fishing mortality rate}
#'   \item{bbmsy}{B/BMSY}
#'   \item{ffmsy}{F/FMSY}
#' }
"YELLSNEMATL"

#' RAM Legacy Database catch-only status predictions
#'
#' A dataset containing the catch and biomass time series for USA Southern
#' New England/Mid-Atlantic (SNE/MA) Yellowtail flounder (\emph{Pleuronectes ferruginea}) from 1973-2014.
#' This dataset is included for additional testing of the \pkg{datalimited2} package.
#'
#' @format A data frame with 161 rows (stocks) and 57 variables including B/BMSY estimates
#' and status estimates.
"preds"
