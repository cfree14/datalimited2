
#' Plot data-limited stock assessment output
#'
#' Plots the results of data-limited stock assessment models implemented in the
#' datalimited2 package.
#'
#' @param output Output from a datalimited2 model
#' @details Produces different plots for each model:
#'
#' 1. zBRT - Plots show: (A) catch time series; (B) saturation time series; and
#' (C) B/BMSY time series.
#'
#' 2. OCOM - Plots show: (A) catch time series; (B) viable r/k pairs;
#' (C) saturation time series; and (D) B/BMSY time series.
#'
#' 3. cMSY - Plots show: (A) catch time series; (B) viable r/k pairs;
#' (C) B/BMSY time series; (D) F/FMSY time series; and (E) Kobe plot.
#'
#' 4. BSM - Plots show: (A) catch time series; (B) viable r/k pairs;
#' (C) B/BMSY time series; (D) F/FMSY time series; and (E) Kobe plot.
#'
#' @examples
#' output <- ocom(year=YELLSNEMATL$year, catch=YELLSNEMATL$tc, m=0.2)
#' plot_dlm(output)
#' @export
plot_dlm <- function(output){

  # Get method
  method <- output$method
  if(method=="zBRT"){plot_zbrt(output)}
  if(method=="OCOM"){plot_ocom(output)}
  if(method%in%c("cMSY", "BSM")){plot_cmsy2(output)}

}
