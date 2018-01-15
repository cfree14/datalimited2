
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
#' Estimates saturation (B/K) and stock status (B/BMSY) time series and
#' other biological quantities (e.g., r, k, MSY, final year saturation) from a time
#' series of catch and a natural mortality (M) estimate using the optimized
#' catch-only model (OCOM) from Zhou et al. (2017).
#'
#' @param year A time series of years
#' @param catch A time series of catch
#' @param m Natural mortality (1/yr)
#' @return A list with the following elements: (1) time series of B/BMSY estimates;
#' (2) 1000 randomly selected biomass trajectories; (3) 1000 corresponding B/BMSY
#' trajectories; (4) estimates of biological quanties r, k, MSY, S; and (5) the
#' 10,000 draws underpinning these values.
#' @references Zhou S, Punt AE, Smith ADM, Ye Y, Haddon M, Dichmont CM, Smith DC
#' (2017) An optimised catch-only assessment method for data poor fisheries.
#' ICES Journal of Marine Science: doi:10.1093/icesjms/fsx226.
#' \url{https://doi.org/10.1093/icesjms/fsx226}
#' @examples
#' output <- ocom(year=TIGERFLAT$yr, catch=TIGERFLAT$catch, m=0.27)
#' plot_ocom(output)
#' @export
ocom <- function(year, catch, m){

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

  # Create B/BMSY time series
  ########################################

  # Convert to B/BMSY
  # S = B / K
  # B/BMSY =
  s_trajs <- apply(b_trajs, 2, function(x) x/x[1])
  bbmsy_trajs <- s2bbmsy(s_trajs)

  # Build B/BMSY dataframe
  s_quants <- t(apply(s_trajs, 1, function(x) stats::quantile(x, quants)))
  colnames(s_quants) <- paste0("q", quants)
  s_avgs <- apply(s_trajs, 1, mean)
  s_sds <- apply(s_trajs, 1, stats::sd)
  s_ts <- data.frame(year=yr, catch=C, s_quants, avg=s_avgs, sd=s_sds)

  # Build B/BMSY dataframe
  bbmsy_quants <- t(apply(bbmsy_trajs, 1, function(x) stats::quantile(x, quants)))
  colnames(bbmsy_quants) <- paste0("q", quants)
  bbmsy_avgs <- apply(bbmsy_trajs, 1, mean)
  bbmsy_sds <- apply(bbmsy_trajs, 1, stats::sd)
  bbmsy_ts <- data.frame(year=yr, catch=C, bbmsy_quants, avg=bbmsy_avgs, sd=bbmsy_sds)

  # Create model output
  ########################################

  # Format parameter estimates
  krms <- data.frame(t(summ))
  colnames(krms) <- paste0("q", quants)
  krms$param <- c("k", "r", "msy", "s")
  krms <- subset(krms, select=c(param, q0.025, q0.25, q0.5, q0.75, q0.975))

  # Prep for export
  krms_draws <- kr2
  output <- list(bbmsy_ts=bbmsy_ts, b_ts=b_ts,
                 s_trajs=s_trajs, s_ts=s_ts,
                 bbmsy_trajs=bbmsy_trajs, b_trajs=b_trajs,
                 krms=krms, krms_draws=krms_draws)
  return(output)

}







