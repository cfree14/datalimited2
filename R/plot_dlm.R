
#' Plot data-limited stock assessment output
#'
#' Plots the results of data-limited stock assessment models.
#'
#' @param output Output from a datalimited2 model
#' @details Produces different plots for each data-limited stock assessment model.
#'
#' 1. zBRT - Produces plots showing: (A) catch time series; (B) saturation time series; and
#' (C) B/BMSY time series.
#'
#' 2. OCOM - Produces plots showing: (A) catch time series; (B) viable r/k pairs;
#' (C) saturation time series; and (D) B/BMSY time series.
#'
#' 3. cMSY - Produces plots showing: (A) catch time series; (B) viable r/k pairs;
#' (C) B/BMSY time series, (D) F/FMSY time series.
#'
#' 4. BSM - Produces plots showing: (A) catch time series; (B) viable r/k pairs;
#' (C) B/BMSY time series, (D) F/FMSY time series.
#'
#' @examples
#' output <- ocom(year=YELLSNEMATL$year, catch=YELLSNEMATL$tc, m=0.2)
#' plot_dlm(output)
#' @export
plot_dlm <- function(output){

  # Get t
  method <- output$method
  if(method=="zBRT"){plot_zbrt(output)}
  if(method=="OCOM"){plot_ocom(output)}
  if(method=="cMSY"){plot_cmsy2(output)}
  if(method=="BSM"){plot_bsm(output)}

}
