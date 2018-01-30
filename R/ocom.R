
# For testing
# load("data/TIGERFLAT.rda")
# year <- TIGERFLAT$yr;  catch <- TIGERFLAT$catch; m=0.27

# OCOM helper functions
################################################################################

# Derive satuation (B/K) prior between [0,1]
Sdistrib = function(n, s_mean) {
  nv = 0 ; n.redo = 0
  while(nv < n) {
    n.redo = n.redo+1
    if(s_mean<=0.5) {
      si1 = fGarch::rsnorm(n*n.redo, mean=max(s_mean,0)-0.072, sd=0.189, xi=0.763)
      si = si1[si1>0 & si1<1]
    } else if(s_mean>0.5) {
      si1 = fGarch::rsnorm(n*n.redo, mean=max(s_mean,0)+0.179, sd=0.223, xi=0.904)
      si = si1[si1>0 & si1<1] }
    if(length(si)>n) si = sample(si,n);
    nv = length(si) }
  return (si)
}

# Biomass dynamics model
BDM = function(K, r, S, b, C) {
  nyr = length(C)
  B = vector()
  B[1] = K*b
  for (i in 1:nyr) {
    B[i+1] = max(min(B[i]+r*B[i]*(1-B[i]/K)-C[i], K), 0)
  }
  if (all(B[-nyr]>C) & all(B<=K)) abs(B[nyr]/K-S) else max(K)*10^4
}

# Optimized catch-only model (OCOM)
################################################################################

#' Optimized catch-only model
#'
#' Estimates biomass, fishing mortality, and stock status (B/BMSY, F/FMSY) time series and
#' biological/management quantities (i.e., r, k, MSY, BMSY, FMSY) from a time
#' series of catch and a natural mortality (M) estimate using the optimized
#' catch-only model (OCOM) from Zhou et al. 2017.
#'
#' @param year A time series of years
#' @param catch A time series of catch
#' @param m Natural mortality (1/yr)
#' @return A list with the following elements:
#' \item{ref_pts}{A data frame with estimates of r, k, MSY, BMSY, FMSY and final year saturation and confidence intervals}
#' \item{ref_ts}{A data frame with time series of biomass, saturation, fishing mortality, B/BMSY, and F/FMSY estimates and 95\% confidence intervals}
#' \item{b_ts}{A data frame with time series of biomass estimates and expanded confidence intervals}
#' \item{s_ts}{A data frame with time series of saturation estimates and expanded confidence intervals}
#' \item{f_ts}{A data frame with time series of fishing mortality (F) estimates and expanded confidence intervals}
#' \item{bbmsy_ts}{A data frame with time series of B/BMSY estimates and expanded confidence intervals}
#' \item{ffmsy_ts}{A data frame with time series of F/FMSY estimates and expanded confidence intervals}
#' \item{b_trajs}{A data frame with 1000 randomly selected biomass trajectories}
#' \item{s_trajs}{A data frame with 1000 corresponding saturation trajectories}
#' \item{f_trajs}{A data frame with 1000 corresponding fishing mortality (F) trajectories}
#' \item{bbmsy_trajs}{A data frame with 1000 corresponding B/BMSY trajectories}
#' \item{ffmsy_trajs}{A data frame with 1000 corresponding F/FMSY trajectories}
#' \item{krms_draws}{A data frame with 10,000 draws underpinning the estimakes of r, k, MSY, and saturation}
#' \item{method}{Name of the method}
#' @details The "optimized catch-only model" (OCOM) developed by Zhou et al. 2017
#' employs a stock reduction analysis (SRA) using priors for r and stock depletion
#' derived from natural mortality and saturation estimated using the Zhou-BRT method, respectively.
#' The SRA employs a Schaefer biomass dynamics model and an algorithm for identifying
#' feasible parameter combinations to estimate biomass, fishing mortality, and stock status (B/BMSY, F/FMSY) time series and
#' biological/management quantities (i.e., r, K, MSY, BMSY, FMSY).
#' @references Zhou S, Punt AE, Smith ADM, Ye Y, Haddon M, Dichmont CM, Smith DC
#' (2017) An optimised catch-only assessment method for data poor fisheries.
#' \emph{ICES Journal of Marine Science}: doi:10.1093/icesjms/fsx226.
#' \url{https://doi.org/10.1093/icesjms/fsx226}
#' @examples
#' # Fit OCOM to catch time series and plot output
#' set.seed(1) # stochastic fitting
#' output <- ocom(year=TIGERFLAT$yr, catch=TIGERFLAT$catch, m=0.27)
#' plot_dlm(output)
#'
#' # Extract reference points and time series from output
#' ref_pts <- output[["ref_pts"]]
#' ref_ts <- output[["ref_ts"]]
#' @export
ocom <- function(year, catch, m){

  # Perform a few error checks
  if(sum(is.na(catch))>0){stop("Error: NA in catch time series. Fill or interpolate.")}

  # Parameters
  nsim = 10000 # number of simulations
  summ = matrix(NA, 5, 4) # matrix to hold output

  # Rename data (to be consistent w/ original OCOM code)
  C = catch
  yr = year
  M = m
  nyr = length(yr)

  # Derive r prior with lognorm dist
  # For telosts - as designated by comment in original OCOM code
  r_median = 2*0.866*M
  r_sig2 = (2*0.866)^2*(0.0012+0.23)
  # For chondrichthyans - commented out in original OCOM code
  # r_median = 2*0.41*M
  # r_sig2 = (2*0.41)^2*(0.0012+0.2
  r_sig = sqrt(r_sig2)
  r_mu = log(r_median)
  ri = stats::rlnorm(nsim, r_mu, r_sig)

  # Estimate final year saturation using zBRT
  s_out <- fit_zbrt(year=year, catch=catch)

  # Derive S prior btw 0 and 1
  s_mean = s_out$s
  si = Sdistrib(nsim, s_mean)

  # Fit biomass dynamics model
  rs = cbind(r=ri, s=si)
  k.low = max(C); k.up=max(C)*200
  opt = apply(rs, 1, function(x){stats::optimize(BDM, c(k.low, k.up), r=x[["r"]], S=x[["s"]], b=1, C=C)})

  # Combine parameter estimates
  ki = sapply(opt, '[[', 1)
  msy = ki*ri/4
  obji = sapply(opt,'[[',2)
  kr = data.frame(k=ki, r=ri, msy, s=si, obj=obji)

  # Clean parameter estimates
  kr2 = kr
  kr2[kr2$k<1.01*k.low | kr2$k>0.99*k.up | kr2$obj>0.01,] = NA # eliminate bordering effect
  kr2 = stats::na.omit(kr2)

  # Summarize parameter estimates
  summ <- apply(kr2, 2, function(x) stats::quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))[,1:4]

  # Create biomass time series
  ########################################

  # Empty vectors for B (one biomass traj) and Bmed (the median biomass traj)
  B = Bmed = vector()

  # Number of trjactories to simulate
  ntrajs = 1000

  # Clean and sample biological quantity draws
  oc0 <- kr
  oc1 = oc0[oc0$obj<0.01 & oc0$k>1.01*min(oc0$k) & oc0$k<0.99*max(oc0$k),] #  eliminate bordering effect
  oc2 = oc1[oc1$k>stats::quantile(oc1$k, 0.25) & oc1$k<stats::quantile(oc1$k, 0.75),]
  smp = sample(1:nrow(oc2), ntrajs, replace=T)

  # Calculate median r and k
  k = oc2[smp,1]
  r = oc2[smp,2]
  kmed = stats::median(k)
  rmed = stats::median(r)

  # Simulate biomass trajectory for each of the selected r-k pairs
  # Trajectories begin at unfished biomass or carrying capacity (k)
  B <- matrix(NA, nrow=nyr, ncol=ntrajs, dimnames=list(yr,1:ntrajs))
  for(j in 1:ntrajs){
    B[1, j] = k[j]
    for(t in 1:(nyr-1)){
      B[t+1, j] = B[t,j] + r[j]*B[t, j]*(1-B[t, j]/k[j]) - C[t]
    }
  }

  # Quantiles for calculations below
  quants <- c(0.025, 0.25, 0.5, 0.75, 0.975)

  # Biomass time series
  b_trajs <- B
  b_quants <- t(apply(b_trajs, 1, function(x) stats::quantile(x, quants)))
  colnames(b_quants) <- paste0("q", quants)
  b_avgs <- apply(b_trajs, 1, mean)
  b_sds <- apply(b_trajs, 1, stats::sd)
  b_ts <- data.frame(year=yr, catch=C, b_quants, avg=b_avgs, sd=b_sds)

  # Fishing mortality time series (added by CMF)
  f_trajs <- apply(b_trajs, 2, function(x) catch / x)
  f_quants <- t(apply(f_trajs, 1, function(x) stats::quantile(x, quants)))
  colnames(f_quants) <- paste0("q", quants)
  f_avgs <- apply(f_trajs, 1, mean)
  f_sds <- apply(f_trajs, 1, stats::sd)
  f_ts <- data.frame(year=yr, catch=C, f_quants, avg=f_avgs, sd=f_sds)

  # Ref points
  ########################################

  # Biological quantities and reference points (added by CMF)
  krms_draws <- kr2
  krms <- data.frame(t(summ))
  r <- krms["r",]
  k <- krms["k",]
  msy <- krms["msy",]
  bmsy <- k / 2
  fmsy <- r / 2
  s_final <- krms["s",]
  ref_pts <- data.frame(rbind(r, k, msy, bmsy, fmsy, s_final))
  ref_pts$param <- c("r", "k", "msy", "bmsy", "fmsy", "s_final")
  rownames(ref_pts) <- NULL
  colnames(ref_pts) <- c(paste0("q", quants), "param")
  ref_pts <- subset(ref_pts, select=c(param, q0.025, q0.25, q0.5, q0.75, q0.975))

  # Create B/BMSY and F/FMSY time series
  ########################################

  # Saturation time series (added by CMF)
  s_trajs <- apply(b_trajs, 2, function(x) x/x[1])
  s_quants <- t(apply(s_trajs, 1, function(x) stats::quantile(x, quants)))
  colnames(s_quants) <- paste0("q", quants)
  s_avgs <- apply(s_trajs, 1, mean)
  s_sds <- apply(s_trajs, 1, stats::sd)
  s_ts <- data.frame(year=yr, catch=C, s_quants, avg=s_avgs, sd=s_sds)

  # Build B/BMSY dataframe (added by CMF)
  bbmsy_trajs <- s2bbmsy(s_trajs)
  bbmsy_quants <- t(apply(bbmsy_trajs, 1, function(x) stats::quantile(x, quants)))
  colnames(bbmsy_quants) <- paste0("q", quants)
  bbmsy_avgs <- apply(bbmsy_trajs, 1, mean)
  bbmsy_sds <- apply(bbmsy_trajs, 1, stats::sd)
  bbmsy_ts <- data.frame(year=yr, catch=C, bbmsy_quants, avg=bbmsy_avgs, sd=bbmsy_sds)

  # Build F/FMSY dataframe (added by CMF)
  ffmsy_trajs <- f_trajs / ref_pts$q0.5[ref_pts$param=="fmsy"]
  ffmsy_quants <- t(apply(ffmsy_trajs, 1, function(x) stats::quantile(x, quants)))
  colnames(ffmsy_quants) <- paste0("q", quants)
  ffmsy_avgs <- apply(ffmsy_trajs, 1, mean)
  ffmsy_sds <- apply(ffmsy_trajs, 1, stats::sd)
  ffmsy_ts <- data.frame(year=yr, catch=C, ffmsy_quants, avg=ffmsy_avgs, sd=ffmsy_sds)

  # Create simple time series data frame
  ref_ts <- cbind(subset(b_ts, select=c(year, catch, q0.5, q0.025, q0.975)),
                  subset(s_ts, select=c(q0.5, q0.025, q0.975)),
                  subset(f_ts, select=c(q0.5, q0.025, q0.975)),
                  subset(bbmsy_ts, select=c(q0.5, q0.025, q0.975)),
                  subset(ffmsy_ts, select=c(q0.5, q0.025, q0.975)))
  colnames(ref_ts) <- c("year", "catch",
                        "b", "b_lo", "b_hi",
                        "s", "s_lo", "s_hi",
                        "f", "f_lo", "f_hi",
                        "bbmsy", "bbmsy_lo", "bbmsy_hi",
                        "ffmsy", "ffmsy_lo", "ffmsy_hi")

  # Create model output
  ########################################

  # Prep for export
  output <- list(ref_pts=ref_pts, ref_ts=ref_ts,
                 b_ts=b_ts, s_ts=s_ts, f_ts=f_ts, bbmsy_ts=bbmsy_ts, ffmsy_ts=ffmsy_ts,
                 b_trajs=b_trajs, s_trajs=s_trajs, f_trajs=f_trajs, bbmsy_trajs=bbmsy_trajs, ffmsy_trajs=ffmsy_trajs,
                 krms_draws=krms_draws, method="OCOM")
  return(output)

}







