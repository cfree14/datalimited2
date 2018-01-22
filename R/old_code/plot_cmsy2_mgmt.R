
#' Plot cMSY and BSM model results for management
#'
#' Plots cMSY and BSM model results for management following the example of Froese et al. (2016).
#' Produces the following four plots:
#' \itemize{
#'   \item{A - Catch time series}
#'   \item{B - B/BMSY time series}
#'   \item{C - F/FMSY time series}
#'   \item{D - Kobe plot}
#' }
#'
#' @param output Output from the cMSY or BSM stock assessment model (see ?cmsy2 or ?bsm)
#' @return Four plots: (1) catch time series; (2) B/BMSY time series; (3) F/FMSY time serie; and (4) Kobe plot
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
plot_cmsy2_mgmt <- function(output){

  # cMSY or BSMY?
  model <- ifelse(length(output)==6, "cMSY", "BSM")

  # Unpack output
  ref_pts <- output[["ref_pts"]]
  ref_ts <- output[["ref_ts"]]
  priors <- output[["priors"]]

  # Extract year stuff
  yr <- ref_ts$year
  nyr <- length(yr)
  yr1 <- yr[1]
  yr2 <- yr[length(yr)]
  int.yr <- priors$year[priors$param=="intbio"]

  # A. Catch
  ##############################################

  # Extract MSY stuff
  msy <- ref_pts$est[ref_pts$param=="msy"]
  msy_lo <- ref_pts$lo[ref_pts$param=="msy"]
  msy_hi <- ref_pts$hi[ref_pts$param=="msy"]

  # Plot catch time series
  ymax <- max(ref_ts$catch)
  plot(yr, 1:nyr , type="n", bty="l",
       main="Catch", xlab="", ylab="Catch in 1000 t", ylim=c(0, ymax))
  rect(yr1, msy_lo, yr2, msy_hi, col="lightgray", border=NA)
  lines(x=c(yr1, yr2), y=rep(msy,2), lty="dashed", col="black", lwd=1.5)
  lines(x=yr, y=ref_ts$catch, lwd=2)
  text("MSY", x=yr2-1.5, y=msy+msy*0.1)

  # B. B/BMSY
  ##############################################

  # Plot B/BMSY time series
  plot(yr, rep(0,nyr), type="n", ylim=c(0,max(c(2, max(ref_ts$bbmsy_hi)))),
       ylab="B / Bmsy", xlab="Year", main="Biomass", bty="l")
  polygon(c(yr,rev(yr)), c(ref_ts$bbmsy_lo, rev(ref_ts$bbmsy_hi)), col="lightgray", border=NA)
  lines(yr, ref_ts$bbmsy, lwd=2)
  lines(x=c(yr1, yr2), y=c(1,1), lty="dashed", lwd=1.5)
  lines(x=c(yr1, yr2), y=c(0.5,0.5), lty="dotted", lwd=1.5)

  # C. F/FMSY
  ##############################################

  # Plot F/FMSY over time
  ymax <- ceiling(max(ref_ts$ffmsy_hi))
  plot(yr, rep(0,nyr), type="n", ylim=c(0, ymax),
       ylab="F / Fmsy", xlab="Year", main="Exploitation", bty="l")
  polygon(c(yr,rev(yr)), c(ref_ts$ffmsy_lo, rev(ref_ts$ffmsy_hi)), col="lightgray", border=NA)
  lines(x=yr, y=ref_ts$ffmsy, lwd=2)
  lines(x=c(yr[1], yr[nyr]), y=c(1,1), lty="dashed", lwd=1.5)

  # D. Stock status plot
  ##############################################

  # Extract ref values
  bbmsy <- ref_ts$bbmsy
  last_bbmsy <- ref_ts$bbmsy[nrow(ref_ts)]
  last_bbmsy_lo <- ref_ts$bbmsy_lo[nrow(ref_ts)]
  last_bbmsy_hi <- ref_ts$bbmsy_hi[nrow(ref_ts)]
  ffmsy <- ref_ts$ffmsy
  last_ffmsy <- ref_ts$ffmsy[nrow(ref_ts)]
  last_ffmsy_lo <- ref_ts$ffmsy_lo[nrow(ref_ts)]
  last_ffmsy_hi <- ref_ts$ffmsy_hi[nrow(ref_ts)]

  # Do stuff
  log.sd.B.Bmsy = (log(last_bbmsy_hi+0.0011)-log(last_bbmsy_lo+0.001))/(2*1.96)
  log.sd.F.Fmsy = (log(last_ffmsy_hi+0.005)-log(last_ffmsy_lo+0.001))/(2*1.96)
  x.F_Fmsy= rlnorm(20000,log(last_ffmsy+0.001),log.sd.F.Fmsy)
  y.b_bmsy =rlnorm(20000,log(last_bbmsy+0.001),log.sd.B.Bmsy)

  kernelF <- ci2d(x.F_Fmsy,y.b_bmsy,nbins=201,factor=2.2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none")
  c1 <- c(-1,100)
  c2 <- c(1,1)

  max.x1   <- max(c(2, max(kernelF$contours$"0.95"$x, last_ffmsy_hi), na.rm =T))
  max.x    <- ifelse(max.x1 > 5, min(max(5,ffmsy*2),8),max.x1)
  max.y    <- max(max(2,quantile(y.b_bmsy,0.96)))

  # Create plot
  xmax <- ceiling(max(ref_ts$ffmsy, max.x))
  plot(1000,1000, type="b", xlim=c(0,xmax), ylim=c(0,max.y),lty=3,xlab="",ylab="", bty="l")
  mtext("F / Fmsy",side=1, line=2)
  mtext("B / Bmsy",side=2, line=2)

  # extract interval information from ci2d object
  # and fill areas using the polygon function
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")

  ## Add points and trajectory lines
  lines(c1, c2, lty=3, lwd=0.7)
  lines(c2, c1, lty=3, lwd=0.7)
  lines(ffmsy, bbmsy, lty=1, lwd=1.)

  points(ffmsy, bbmsy, cex=0.8, pch=4)
  points(ffmsy[1], bbmsy[1], col=1, pch=22, bg="white", cex=1.9)
  points(ffmsy[which(yr==int.yr)], bbmsy[which(yr==int.yr)], col=1, pch=21, bg="white", cex=1.9)
  points(ffmsy[nyr], bbmsy[nyr], col=1, pch=24, bg="white", cex=1.9)

  ## Add legend
  legend('topright', c(paste(yr1),paste(int.yr),paste(yr2),"50% C.I.","80% C.I.","95% C.I."),
         lty=c(1,1,1,-1,-1,-1),pch=c(22,21,24,22,22,22),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4"),
         col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),1.7,1.7,1.7),bty="n",y.intersp = 0.9)

}



