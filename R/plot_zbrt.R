
#' Plot zBRT results
plot_zbrt <- function(output){

  # Subset data
  output <- na.omit(output[["ts"]])

  # Plot settings
  par(mfrow=c(2,2), xpd=NA)

  # Catch
  #################################

  # Plot catch
  xmin <- floor(min(output$year) / 10) * 10
  xmax <- ceiling(max(output$year) / 10) * 10
  plot(catch ~ year, output, bty="n", las=1, type="l",
       xlim=c(xmin, xmax), xlab="", ylab="Catch", main="A. Catch")

  # Saturation
  #################################

  # Plot B/BMSY time series
  yr_final <- max(output$year)
  s_final <- output$s[nrow(output)]
  xmin <- floor(min(output$year) / 10) * 10
  xmax <- ceiling(max(output$year) / 10) * 10
  plot(s ~ year, output, bty="n", type="n", las=1,
       xlim=c(xmin, xmax), ylim=c(0,1), xlab="", ylab="Saturation",
       main="B. Saturation (B/K)")
  # Add CI shading
  polygon(x=c(output$year, rev(output$year)),
          y=c(output$s_lo, rev(output$s_hi)), col="grey80", border=F)
  # Add median and 95% CI trajectories
  lines(x=output$year, y=output$s, lwd=1.1)
  # lines(x=output$year, y=output$s_lo, lwd=1.1, lty=2)
  # lines(x=output$year, y=output$s_hi, lwd=1.1, lty=2)
  # Add overfished line (B/BMSY=0.25)
  lines(x=c(xmin, xmax), y=c(0.25, 0.25), lty=3)
  # Label end year saturation
  text(x=yr_final, y=s_final, label=round(s_final,2), pos=4)

  # B/BMSY
  #################################

  # Plot B/BMSY time series
  yr_final <- max(output$year)
  bbmsy_final <- output$bbmsy[nrow(output)]
  xmin <- floor(min(output$year) / 10) * 10
  xmax <- ceiling(max(output$year) / 10) * 10
  plot(bbmsy ~ year, output, bty="n", type="n", las=1,
       xlim=c(xmin, xmax), ylim=c(0,2), xlab="", ylab=expression("B / B"["MSY"]),
       main=expression(bold("C. B/B"["MSY"])))
  # Add CI shading
  polygon(x=c(output$year, rev(output$year)),
          y=c(output$bbmsy_lo, rev(output$bbmsy_hi)), col="grey80", border=F)
  # Add median and 95% CI trajectories
  lines(x=output$year, y=output$bbmsy, lwd=1.1)
  # lines(x=output$year, y=output$bbmsy_lo, lwd=1.1, lty=2)
  # lines(x=output$year, y=output$bbmsy_hi, lwd=1.1, lty=2)
  # Add overfished line (B/BMSY=0.5)
  lines(x=c(xmin, xmax), y=c(0.5, 0.5), lty=3)
  # Label end year B/BMSY
  text(x=yr_final, y=bbmsy_final, label=round(bbmsy_final,2), pos=4)

}








