
# Calculate proportional error
calc_pe <- function(obs, preds){
  (preds-obs)/abs(obs)
}

# Calculate bias (median proportional error, MPE)
calc_bias <- function(obs, preds){
  median(calc_pe(obs, preds))
}

# Calculate inaccuracy (median absolute proportional error, MAPE)
calc_inaccuracy <- function(obs, preds){
  median(abs(calc_pe(obs, preds)))
}

# Calculate rank correltation
calc_rank_corr <- function(obs, preds){
  cor(obs, preds, method="spearman", use="pairwise.complete.obs")
}

# Calculate accurcacy (categorical)
# obs <- data_eval$bbmsy_status; preds <- data_eval$super2_status
calc_accuracy <- function(obs, preds){
  obs_f <- factor(obs, levels=c("under", "fully", "over"))
  preds_f <- factor(preds, levels=c("under", "fully", "over"))
  stats <- caret::postResample(pred=preds_f, obs=obs_f)
  accuracy <- stats[1]
  return(accuracy)
}

# Calculate Cohen's kappa (categorical)
# obs <- data_eval$bbmsy_status; preds <- data_eval$super2_status
calc_kappa <- function(obs, preds){
  obs_f <- factor(obs, levels=c("under", "fully", "over"))
  preds_f <- factor(preds, levels=c("under", "fully", "over"))
  stats <- caret::postResample(pred=preds_f, obs=obs_f)
  accuracy <- stats[1]
  kappa <- stats[2]
  if(accuracy==1 & is.na(kappa)){kappa <- 1}
  return(kappa)
}

#' Compare performance of catch-only methods
#'
#' Compares the continuous and/or categorical performance of catch-only stock assessment
#' methods. Creates both a plot and table of performance metrics.
#'
#' @param status A matrix or dataframe of continuous (e.g., B/BMSY) or
#' categorical (e.g., under, fully, overexploited) status estimates in which the
#' first column is the true (or data-rich) status estimate and additional columns
#' are status estimates from catch-only methods.
#' @param methods A character vector of the names of the catch-only models being compared,
#' listed in the order they are presented in the status matrix or dataframe.
#' @return A plot and table of continuous or categorical performance metrics.
#' @details The continuous performance of each COM is evaluated by measuring
#' each method’s bias, accuracy, and ability to correctly rank or correlate
#' across populations. We measured bias as the median proportional error (MPE)
#' and accuracy as the median absolute proportional error (MAPE). Proportional
#' error is calculated as (θest-θ)/|θ|, where θest and θ represent predicted
#' and “true” (or data-rich stock assessment) B/BMSY values. The ability to
#' correctly rank populations is measured as Spearman’s rank-order correlation
#' between predicted and “true” values.
#'
#' The categorical performance of each COM is evaluated using both percentage
#' agreement (accuracy) and Cohen’s kappa. Cohen’s kappa measures inter-rate
#' agreement between categorical items and is more robust than simple percentage
#' agreement because it takes into account the probability of agreement occurring
#' by chance alone (Cohen 1968). Although there are no definitive rules for
#' interpreting Cohen’s kappa, general guidelines suggest that values >0.70
#' are ‘excellent’, 0.4-0.7 are ‘good’, 0.2-0.4 are ‘fair’, and <0.2 are
#' ‘poor’ (Landis & Koch 1977; Fleiss 1973).
#' @examples
#' bbmsy <- select(preds, bbmsy,  mprm, comsir, sscom, cmsy13, cmsy17, zbrt, ocom, super1)
#' status <- bbmsy2catg(bbmsy)
#' methods <- c("mPRM", "COMSIR", "SSCOM", "cMSY-13", "cMSY-17", "zBRT", "OCOM", "Super")
#' performance(bbmsy, methods)
#' performance(status, methods)
#' @export
performance <- function(status, methods){

  # Remove NA values
  status_clean <- na.omit(status)
  n_removed <- nrow(status) - nrow(status_clean)
  if(n_removed!=0){cat("Note:", n_removed, "stocks removed due to missing data.\n")}

  # Determine type
  type <- ifelse(is.numeric(status_clean[,1]), "cont", "catg")

  # Continuous performance
  #####################################

  if(type=="cont"){

    # Calculate continuous performance
    perf <- data.frame(method=methods, n=NA, inaccuracy=NA, rank=NA, bias=NA)
    obs <- status_clean[,1]
    for(i in 2:ncol(status_clean)){
      preds <-  status_clean[,i]
      perf$n[i-1] <- length(obs)
      perf$inaccuracy[i-1] <- calc_inaccuracy(obs, preds)
      perf$rank[i-1] <- calc_rank_corr(obs, preds)
      perf$bias[i-1] <- calc_bias(obs, preds)
    }

    # Add bias color
    summary(perf$bias)
    bmin <- floor(min(perf$bias) / 0.1) * 0.1 # lower val
    bmax <- ceiling(max(perf$bias) / 0.1) * 0.1 # upper val
    blim <- c(bmin, bmax)
    babs <- blim[which.max(abs(blim))] # which value is largest?
    blim1 <- c(babs*-1, babs) # create - and + version of largest value so color pal is even
    b_breaks <- seq(blim1[1], blim1[2], length.out=11)
    b_bins <- cut(perf$bias, breaks=b_breaks)
    b_cols <- RColorBrewer::brewer.pal(11, "RdBu")[b_bins]

    # Plot continuous performance
    xmin <- floor(min(perf$inaccuracy) / 0.2) * 0.2
    xmax <- ceiling(max(perf$inaccuracy) / 0.2) * 0.2
    ymin <- floor(min(perf$rank) / 0.2) * 0.2
    ymax <- ceiling(max(perf$rank) / 0.2) * 0.2
    plot(rank ~ inaccuracy, perf, bty="n", las=1,
         pch=21, col="black", bg=b_cols, cex=1.4,
         xlim=c(xmin, xmax), ylim=c(ymin, ymax),
         xlab="Inaccuracy (MAPE)", ylab="Rank-order correlation")
    text(x=perf$inaccuracy, y=perf$rank, labels=perf$method, pos=2, xpd=NA)
    text(x=xmax, y=ymin, pos=2, labels=paste(unique(perf$n), "stocks"))

  }

  # Categorical performance
  #####################################

  if(type=="catg"){

    # Categorical performance
    perf <- data.frame(method=methods, n=NA, accuracy=NA, kappa=NA)
    obs <- status_clean[,1]
    for(i in 2:ncol(status_clean)){
      preds <-  status_clean[,i]
      perf$n[i-1] <- length(obs)
      perf$accuracy[i-1] <- calc_accuracy(obs, preds)
      perf$kappa[i-1] <- calc_kappa(obs, preds)
    }

    # Plot categorical performance
    xmin <- floor(min(perf$accuracy) / 0.2) * 0.2
    xmax <- ceiling(max(perf$accuracy) / 0.2) * 0.2
    ymin <- floor(min(perf$kappa) / 0.2) * 0.2
    ymax <- pmax(0.4, ceiling(max(perf$kappa) / 0.2) * 0.2)
    plot(kappa ~ accuracy, perf, bty="n", las=1,
         pch=21, col="black", bg="grey80", cex=1.4,
         xlim=c(xmin, xmax), ylim=c(ymin, ymax),
         xlab="Accuracy", ylab="Cohen's kappa")
    text(x=perf$accuracy, y=perf$kappa, labels=perf$method, pos=2, xpd=NA)
    text(x=xmax, y=ymin, pos=2, labels=paste(unique(perf$n), "stocks"))

    # Add Cohen's kappa lines
    lines(x=c(xmin, xmax), y=c(0.2,0.2), lty=3, lwd=1.4, col="grey70")
    lines(x=c(xmin, xmax), y=c(0.4,0.4), lty=3, lwd=1.4, col="grey70")
    text(x=xmax, y=0.18, pos=2, labels="poor", cex=1, col="grey70", xpd=NA)
    text(x=xmax, y=0.22, pos=2, labels="fair", cex=1, col="grey70", xpd=NA)
    text(x=xmax, y=0.42, pos=2, labels="good", cex=1, col="grey70", xpd=NA)

  }

  # Return performance metrics
  return(perf)

}




