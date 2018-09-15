
# Plot MS-cMSY results
# Plot MS-cMSY results
plot_ms_cmsy <- function(output){

  # Loop through species
  spp <- species
  nspp <- length(species)
  par(mfcol=c(3,nspp), mar=c(3,4,2,1), mgp=c(2.5,0.7,0), xpd=NA)
  for(i in 1:length(species)){

    # Subset data
    yrs <- out$yrs
    ids <- out$id_try
    rs <- out$r_try[,i]
    ks <- out$k_try[,i]
    id_rk_viable <- out$id_rk_v[[i]]
    b_viable <- out$b_v[[i]]
    bbmsy_viable <- out$bbmsy_v[[i]]
    er_viable <- out$er_v[[i]]
    s1_priors <- out$s1_priors
    s2_priors <-out$s2_priors
    r_priors <- out$r_priors
    k_priors <- out$k_priors
    er_vv <- out$er_vv[[i]]
    bbmsy_vv <- out$bbmsy_vv[[i]]
    bbmsy_vv_median <- out$bbmsy_vv_median[[i]]
    top_corr <- out$top_corr

    # Plot r/k pairs
    #########################################

    # Plot r/k pairs
    plot(ks ~ rs, log="xy", type="n", bty="n", las=1, pch=15, col="gray80", xpd=NA,
         xlim=r_priors[i,], ylim=k_priors[i,], xlab="r", ylab="k", main=spp[i])

    # Add viable r/k pairs
    # Potentially reduce this to unique r/k pairs
    # There could be redundancy when evaluating multiple IDs
    points(id_rk_viable$r, id_rk_viable$k, pch=15, col="grey70")

    # Add most highly correlated pairs
    id_rk_v_ind <- unlist(top_corr[,paste0("index", i)])
    rk_corr <- subset(id_rk_viable, index %in% id_rk_v_ind)
    points(x=rk_corr$r, y=rk_corr$k, pch=15, col=freeR::tcolor("darkorange", 0.6))

    # # Add most common highly correlated pair
    # rk_mode <- mode(rk_v_ind)
    # rk_corr <- subset(rk_viable, index==rk_mode)
    # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="green")

    # # Add repeatedly highly correlated r/k pairs
    # rk_corr <- subset(rk_viable, ncorr>=5)
    # points(x=rk_corr$r, y=rk_corr$k, pch=15, col="darkorange")

    # Add legend
    if(i==1){
      legend("bottomright", bty="n", pch=15, pt.cex=1.3, cex=0.9,
             col=c("grey70", "darkorange", "red"),
             legend=c("Viable", "Effort correlated"))
    }

    # Plot BBMSY trajectories
    #########################################

    # Plot BBMSY trajectories
    ymax <- freeR::ceiling1(max(bbmsy_viable, na.rm=T), 0.5)
    plot(bbmsy_viable[,1] ~ yrs, type="n", bty="n", las=2,
         ylim=c(0, ymax), xlab="", ylab=expression("B/B"["MSY"]))
    for(k in 1:ncol(bbmsy_viable)){lines(x=yrs, y=bbmsy_viable[,k], col="grey70")}
    for(k in 1:ncol(bbmsy_vv)){lines(x=yrs, y=bbmsy_vv[,k], col=freeR::tcolor("darkorange", 0.6))}
    lines(x=yrs, y=bbmsy_vv_median, lwd=1.5, col="black")
    lines(x=c(0, max(yrs)), y=c(0.5, 0.5), lty=3)
    lines(x=yrs, y=apply(bbmsy_viable,1,median), lwd=1.5, lty=3, col="green")

    # Add legend
    if(i==1){
      legend("bottomright", bty="n", lty=c(1,1), lwd=1.5, cex=0.9,
             col=c("grey70", "darkorange"),
             legend=c("Viable", "Effort correlated"))
    }

    # Plot biomass trajectories
    #########################################

    # # Plot biomass trajectories
    # plot(b_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, max(b_viable)), xlab="Year", ylab="Biomass")
    # for(k in 1:ncol(b_viable)){lines(x=yrs, y=b_viable[,k], col="grey80")}
    # lines(x=yrs, y=true$b_ts[,i+1], lwd=1.5, col="red")

    # Plot exploitation trajectories
    #########################################

    # Plot exploitation trajectories
    plot(er_viable[,1] ~ yrs, type="n", bty="n", las=2,
         ylim=c(0, 1), xlab="", ylab="Exploitation rate")
    for(k in 1:ncol(er_viable)){lines(x=yrs, y=er_viable[,k], col="grey80")}
    for(k in 1:ncol(er_vv)){lines(x=yrs, y=er_vv[,k], col=freeR::tcolor("darkorange", 0.6))}

  }


}
