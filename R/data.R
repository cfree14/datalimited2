#' Irish Sea Common sole time series
#'
#' A dataset containing the catch and biomass time series for Irish Sea Common sole
#' (Solea solea) from 1970-2014. This stock was used as an example in the cMSY/BSM
#' user manual and is used to validate this package's implementation of cMSY and BSM.
#'
#' @format A data frame with 45 rows (years) and 4 variables:
#' \describe{
#'   \item{Stock}{stock id}
#'   \item{yr}{year}
#'   \item{ct}{catch, in metric tons}
#'   \item{bt}{biomass, in metric tons}
#' }
#' @source Froese R, Demirel N, Coro G, Kleisner KM, Winker H (2017)
#' Estimating fisheries reference points from catch and resilience. Fish and Fisheries 18(3): 506-526.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/faf.12190/abstract}
"SOLIRIS"

#' SE Australia Tiger flathead time series
#'
#' A dataset containing the catch time series for SE Australia Tiger flathead
#' (Neoplatycephalus richardsoni) from 1915-2012. This stock was used as an example in the OCOM
#' paper and is used to validate this package's implementation of zBRT and OCOM.
#'
#' @format A data frame with 98 rows (years) and 3 variables:
#' \describe{
#'   \item{stock}{stock id}
#'   \item{yr}{year}
#'   \item{catch}{catch, in metric tons}
#' }
#' @source Zhou S, Punt AE, Smith ADM, Ye Y, Haddon M, Dichmont CM, Smith DC
#' (2017) An optimised catch-only assessment method for data poor fisheries.
#' ICES Journal of Marine Science: doi:10.1093/icesjms/fsx226.
#' \url{https://doi.org/10.1093/icesjms/fsx226}
"TIGERFLAT"

#' USA SNE/MA Yellowtail flounder time series
#'
#' A dataset containing the catch time series for USA SNE/MA Yellowtail flounder
#' (Pleuronectes ferruginea) from 1973-2014.
#'
#' @format A data frame with 42 rows (years) and XX variables:
#' \describe{
#'   \item{Stock}{stock id}
#'   \item{yr}{year}
#'   \item{ct}{catch, in metric tons}
#'   ...
#' }
"YELLSNEMATL"

