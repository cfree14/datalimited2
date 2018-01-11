
#' Plot cMSY and BSM model results
#'
#' Plots cMSY and BSM model results following the example of Froese et al. (2016).
#' Produces the following six plots:
#' \itemize{
#'   \item{A - Catch time series}
#'   \item{B - Finding viable r-k pairs}
#'   \item{C - Viable r-k pairs}
#'   \item{D - Saturation (B/k) time series}
#'   \item{E - Exploitation rate (F / (r/2)) time series}
#'   \item{F - Surplus production curve}
#' }
#'
#' @param output Output from the cMSY or BSM stock assessment models (see ?cmsy2 or ?bsm)
#' @return Six plots: (1) catch time series; (2) r-k pair search; (3) viable r-k pairs;
#' (4) saturation time series; (5) exploitation rate time series; (6) surplus production curve
#' @references Froese R, Demirel N, Coro G, Kleisner KM, Winker H (2016)
#' A Simple User Guide for CMSY and BSM (version “q”). 27 October 2016.
#' \url{http://oceanrep.geomar.de/33076/}
#' @examples
#' # Fit cMSY and plot results
#' output <- cmsy2(year=SOLIRIS$yr, catch=SOLIRIS$ct, r.low=0.18, r.hi=1.02)
#' plot_cmsy2(output)
#' plot_cmsy2_mgmt(output)
#'
#' # Fit BSM and plot results
#' output <- bsm(year=SOLIRIS$yr, catch=SOLIRIS$ct, biomass=SOLIRIS$bt, btype="CPUE", r.low=0.18, r.hi=1.02)
#' plot_cmsy2(output)
#' plot_cmsy2_mgmt(output)
#' @export
plot_cmsy2 <- function(output){

  # cMSY or BSMY?
  model <- ifelse(length(output)==6, "cMSY", "BSM")
  col.use <- ifelse(model=="cMSY", "blue", "red")

  # Unpack output
  ref_pts <- output[["ref_pts"]]
  ref_ts <- output[["ref_ts"]]
  priors <- out[["priors"]]
  if(model=="cMSY"){rv <- output[["rv.all"]]}else{rv <- output[["r_out"]]}
  if(model=="cMSY"){kv <- output[["kv.all"]]}else{kv <- output[["k_out"]]}

  # A. Catch
  ##############################################################

  # Plot catch
  ymax <- ceiling(max(ref_ts$catch))
  plot(catch ~ year, ref_ts,
       type ="l", bty="l", lwd=2, ylim=c(0, ymax),
       main="A: Catch", xlab="Year", ylab="Catch in 1000 t")
  lines(x=ref_ts$year, y=ref_ts$catch_ma, col=col.use, lwd=1)
  # points(x=yr[max.yr.i], y=max.ct, col="red", lwd=2)
  # points(x=yr[min.yr.i], y=min.ct, col="red", lwd=2)

  # B. Finding r/k plot
  ##############################################################

  # Extract r/k info
  start.r <- unlist(priors[priors$param=="r", c("lo", "hi")])
  start.k <- unlist(priors[priors$param=="k", c("lo", "hi")])
  r <- ref_pts$est[ref_pts$param=="r"]
  r_lo <- ref_pts$lo[ref_pts$param=="r"]
  r_hi <- ref_pts$hi[ref_pts$param=="r"]
  k <- ref_pts$est[ref_pts$param=="k"]
  k_lo <- ref_pts$lo[ref_pts$param=="k"]
  k_hi <- ref_pts$hi[ref_pts$param=="k"]

  # Plot all r-k pairs
  plot(x=1:10, y=1:10, type="n", xlim = start.r, ylim = start.k, log="xy", xlab="r", ylab="k",
       main="B: Finding viable r-k", pch=".", cex=3, bty="l", col="gray95")
  rect(xleft=start.r[1], ybottom=start.k[1], xright=start.r[2], ytop=start.k[2], col="grey95", border=F)
  # Add viable r/k pairs
  points(x=rv, y=kv, pch=".", cex=4, col="gray")

  # Add cMSY r/k pair, with 95% CL lines
  points(x=r, y=k, pch=19, col=col.use)
  lines(x=c(r_lo, r_hi), y=c(k, k), col=col.use)
  lines(x=c(r, r), y=c(k_lo, k_hi), col=col.use)

  # C. Viable r/k plot
  ##############################################################

  # Plot viable r/k pairs
  plot(x=rv, y=kv,
       pch=16, col="gray",log="xy", bty="l",
       xlab="r", ylab="k", main="C: Analysis of viable r-k")

  # Add cMSY r/k pair, with 95% CL lines
  points(x=r, y=k, pch=19, col=col.use)
  lines(x=c(r_lo, r_hi), y=c(k, k), col=col.use)
  lines(x=c(r, r), y=c(k_lo, k_hi), col=col.use)

  # D. Predicted biomass plot
  ##############################################################

  # Extact year and prior info
  start_yr <- ref_ts$year[1]
  end_yr <- ref_ts$year[nrow(ref_ts)]
  startbio <- unlist(priors[priors$param=="startbio", c("lo", "hi")])
  intbio <- unlist(priors[priors$param=="intbio", c("lo", "hi")])
  int.yr <- priors$year[priors$param=="intbio"]
  endbio <- unlist(priors[priors$param=="endbio", c("lo", "hi")])

  # Plot relative biomass (B/k or saturation)
  ymax <- ceiling(max(ref_ts$s_hi) / 0.1) * 0.1
  plot(s ~ year, ref_ts, lwd=1.5, xlab="Year", ylab="Relative biomass B/k", type="l",
       ylim=c(0,ymax), bty="l", main="D: Biomass",col=col.use)
  lines(x=ref_ts$year, y=ref_ts$s_lo, type="l", lty="dotted", col=col.use)
  lines(x=ref_ts$year, y=ref_ts$s_hi, type="l", lty="dotted", col=col.use)
  # Add lines for 0.5 and 0.25 biomass
  abline(h=0.5, lty="dashed")
  abline(h=0.25, lty="dotted")
  # plot biomass windows
  lines(x=c(start_yr, start_yr), y=startbio, col=col.use)
  lines(x=c(int.yr, int.yr), y=intbio, col=col.use)
  lines(x=c(end_yr, end_yr), y=endbio, col=col.use)

  # E. Exploitation rate plot
  ##############################################################

  # Plot exploitation rate
  ymax <- ceiling(max(ref_ts$er) / 0.5) * 0.5
  plot(er ~ year, ref_ts, type="l", bty="l", lwd=1.5, ylim=c(0, ymax),
       xlab="Year", ylab="F / (r/2)", main="E: Exploitation rate", col=col.use)
  abline(h=1, lty="dashed")

  # F. Parabola plot
  ##############################################################

  # Catch relative to MSY
  msy <- ref_pts$est[ref_pts$param=="msy"]
  cdivmsy <- ref_ts$catch / msy

  # Y-axis limit
  ymax <- ceiling(max(cdivmsy) / 0.25) * 0.25

  # Plot parabola
  x=seq(from=0,to=2,by=0.001)
  y.c  <- ifelse(x>0.25,1,ifelse(x>0.125,4*x,exp(-10*(0.125-x))*4*x)) # correction for low recruitment below half and below quarter of Bmsy
  y=(4*x-(2*x)^2)*y.c
  plot(x=x, y=y, xlim=c(1,0), ylim=c(0,ymax), type="l", bty="l",xlab="Relative biomass B/k",
       ylab="Catch / MSY", main="F: Equilibrium curve")
  # plot catch against CMSY estimates of relative biomass
  points(x=ref_ts$s, y=cdivmsy, pch=16, col=col.use)

}


