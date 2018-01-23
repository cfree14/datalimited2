
#' Plot data-limited stock assessment model output
#'
#' Plots the results of data-limited stock assessment models implemented in the
#' \pkg{datalimited2} package.
#'
#' @param output Output from a \pkg{datalimited2} model
#' @details Produces different plots for each model:
#'
#' 1. \strong{zBRT} - Plots show: (A) catch time series; (B) saturation time series; and
#' (C) B/BMSY time series.
#'
#' 2. \strong{OCOM} - Plots show: (A) catch time series; (B) viable r/K pairs;
#' (C) saturation time series; and (D) B/BMSY time series.
#'
#' 3. \strong{cMSY / BSM} - Plots show: (A) catch time series; (B) viable r/K pairs;
#' (C) B/BMSY time series; (D) F/FMSY time series; and (E) Kobe plot.
#'
#' In all plots, dashed lines show the reference point target (i.e., B/BMSY = 1,
#' F/FMSY = 1, or saturation = 0.5) and dotted lines show the overfishing limit
#' (i.e., B/BMSY = 0.5 or saturation = 0.25). If MSY is estimated, the median value
#' and 95\% confidence intervals are shown in the catch time series plot as a horizontal dashed line and
#' grey rectangle, respectively.
#'
#' @examples
#' # Fit OCOM to catch time series and plot output
#' output <- ocom(year=YELLSNEMATL$year, catch=YELLSNEMATL$catch, m=0.2)
#' plot_dlm(output)
#' @export
plot_dlm <- function(output){

  # Get method
  method <- output$method
  if(method=="zBRT"){plot_zbrt(output)}
  if(method=="OCOM"){plot_ocom(output)}
  if(method%in%c("cMSY", "BSM")){plot_cmsy2(output)}

}
