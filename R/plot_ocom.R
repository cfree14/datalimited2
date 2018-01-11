


#' Plot OCOM results
#'
#' Plot OCOM results
#'
#' @param output Output from the OCOM model (see ?ocom)
#' @return Plots
#' @references Zhou S, Punt AE, Smith ADM, Ye Y, Haddon M, Dichmont CM, Smith DC
#' (2017) An optimised catch-only assessment method for data poor fisheries.
#' ICES Journal of Marine Science: doi:10.1093/icesjms/fsx226.
#' \url{https://doi.org/10.1093/icesjms/fsx226}
#' @examples
#' output <- ocom(year=YELLSNEMATL$year, catch=YELLSNEMATL$tc, m=0.2)
#' plot_ocom(output)
#' @export
plot_ocom <- function(output){

  # Unpack output
  b_ts <- output[["b_ts"]]
  b_trajs <- output[["b_trajs"]]
  s_ts <- output[["s_ts"]]
  s_trajs <- output[["s_trajs"]]
  bbmsy_ts <- output[["bbmsy_ts"]]
  bbmsy_trajs <- output[["bbmsy_trajs"]]
  krms <- output[["krms"]]
  krms_draws <- output[["krms_draws"]]
  yr_final <- max(bbmsy_ts$year)
  bbmsy_final <- bbmsy_ts$q0.5[nrow(bbmsy_ts)]
  s_final <- s_ts$q0.5[nrow(s_ts)]

  # Viable and best r/k pairs
  ##############################################################################

  # Plot viable r-k pairs
  plot(krms_draws$k, krms_draws$r, log="xy", bty="l", col="grey60", las=1,
       xlab="K", ylab='r', main="OCOM: Viable r-k pairs")

  # Add best r-k pair
  r <- krms$q0.5[krms$param=="r"]
  r_hi <- krms$q0.975[krms$param=="r"]
  r_lo <- krms$q0.025[krms$param=="r"]
  k <- krms$q0.5[krms$param=="k"]
  k_hi <- krms$q0.975[krms$param=="k"]
  k_lo <- krms$q0.025[krms$param=="k"]
  points(k, r, pch=16, cex=1.3)
  lines(x=c(k_lo, k_hi), y=c(r,r))
  lines(x=c(k, k), y=c(r_lo,r_hi))

  # Parameter estimates
  ##############################################################################

  # Parameter estimate histograms
  params <- c("r", "k", "msy", "s")
  param_names <- c("r", "K", "MSY", "Saturation (B/K)")
  for(i in 1:length(params)){
    vals <- krms_draws[,params[i]]
    median(vals)
    hist(vals, las=1, xlab=param_names[i], main=paste0("OCOM: ", param_names[i], " posterior"), col="grey70", border=F)
    abline(v=median(vals), lwd=1.5, lty=2)
  }

  # Biomass time series
  ##############################################################################

  # Setup empty plot
  b_trajs1 <- b_trajs / 1000
  xmin <- floor(min(b_ts$year) / 10) * 10
  xmax <- ceiling(max(b_ts$year) / 10) * 10
  ymax <- ceiling(max(b_trajs1) / 10) * 10
  plot(q0.5 ~ year, b_ts, type="n", bty="n", las=1,
       xlim=c(xmin, xmax), ylim=c(0, ymax), xlab="", ylab="Biomass (in 1000s)",
       main="OCOM: Biomass time series")
  # Add randomly selected trajectories
  for(i in 1:ncol(b_trajs1)){lines(x=b_ts$year, y=b_trajs1[,i], col="grey90")}
  # Add median and 95% CI trajectories
  lines(x=b_ts$year, y=b_ts$q0.5/1000, lwd=1.1)
  lines(x=b_ts$year, y=b_ts$q0.025/1000, lwd=1.1, lty=2)
  lines(x=b_ts$year, y=b_ts$q0.975/1000, lwd=1.1, lty=2)

  # Saturation time series
  ##############################################################################

  # Setup empty plot
  xmin <- floor(min(s_ts$year) / 10) * 10
  xmax <- ceiling(max(s_ts$year) / 10) * 10
  plot(q0.5 ~ year, s_ts, type="n", bty="n", las=1,
       xlim=c(xmin, xmax), ylim=c(0,1), xlab="", ylab="Saturation",
       main="OCOM: Saturation time series")
  # Add randomly selected trajectories
  for(i in 1:ncol(s_trajs)){lines(x=s_ts$year, y=s_trajs[,i], col="grey90")}
  # Add median and 95% CI trajectories
  lines(x=s_ts$year, y=s_ts$q0.5, lwd=1.1)
  lines(x=s_ts$year, y=s_ts$q0.025, lwd=1.1, lty=2)
  lines(x=s_ts$year, y=s_ts$q0.975, lwd=1.1, lty=2)
  # Add overfished line (saturation=0.25)
  lines(x=c(xmin, xmax), y=c(0.25, 0.25), lty=3)
  # Label end year saturation
  text(x=yr_final, y=s_final, label=round(s_final,2), pos=4)

  # B/BMSY time series
  ##############################################################################

  # Setup empty plot
  xmin <- floor(min(bbmsy_ts$year) / 10) * 10
  xmax <- ceiling(max(bbmsy_ts$year) / 10) * 10
  plot(q0.5 ~ year, bbmsy_ts, type="n", bty="n", las=1,
       xlim=c(xmin, xmax), ylim=c(0,2), xlab="", ylab=expression("B / B"["MSY"]),
       main="OCOM: B/BMSY time series")
  # Add randomly selected trajectories
  for(i in 1:ncol(bbmsy_trajs)){lines(x=bbmsy_ts$year, y=bbmsy_trajs[,i], col="grey90")}
  # Add median and 95% CI trajectories
  lines(x=bbmsy_ts$year, y=bbmsy_ts$q0.5, lwd=1.1)
  lines(x=bbmsy_ts$year, y=bbmsy_ts$q0.025, lwd=1.1, lty=2)
  lines(x=bbmsy_ts$year, y=bbmsy_ts$q0.975, lwd=1.1, lty=2)
  # Add overfished line (B/BMSY=0.5)
  lines(x=c(xmin, xmax), y=c(0.5, 0.5), lty=3)
  # Label end year B/BMSY
  text(x=yr_final, y=bbmsy_final, label=round(bbmsy_final,2), pos=4)

}







