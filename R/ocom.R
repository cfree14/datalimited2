
# Derive S = B/K prior bwt [0,1]
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

# biomass dynamics model: one parameter r, using optimize
BDM = function(K, r, S, b, C) {
  nyr = length(C)
  B = vector()
  B[1] = K*b
  for (i in 1:nyr) {
    B[i+1] = max(min(B[i]+r*B[i]*(1-B[i]/K)-C[i], K), 0)
  }
  if (all(B[-nyr]>C) & all(B<=K)) abs(B[nyr]/K-S)  else max(K)*10^4
}

#' Optimized catch-only model
#'
#' Estimates biological quantities (r, k, MSY, final year saturation) and
#' stock status (B/BMSY) from a time series of catch and natural mortality (M)
#' estimate using the optimized catch-only model from Zhou et al. (2017).
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
#' ocom(year=YELLSNEMATL$year, catch=YELLSNEMATL$tc, m=0.2)
#' @export
ocom <- function(year, catch, m){

  # Parameters
  nsim = 10000 # number of simulations
  summ = matrix(NA, 5, 4) # matrix to hold output

  # Rename data (to be consistent w/ original OCOM code)
  C = catch;
  yr = year
  M = m
  nyr = length(yr)

  # derive r prior with lognorm dist
  r_median = 2*0.866*M;              # for teleosts
  r_sig2 = (2*0.866)^2*(0.0012+0.23)
  #   r_median = 2*0.41*M;              # for chondrichthyans
  #      r_sig2 = (2*0.41)^2*(0.0012+0.2
  r_sig = sqrt(r_sig2)
  r_mu = log(r_median)
  ri = rlnorm(nsim, r_mu, r_sig)

  # Estimate final year saturation
  s_out <- fit_zbrt(year=year, catch=catch)

  # Derive S prior btw 0 and 1
  s_mean = s_out$s
  si = Sdistrib(nsim, s_mean)

  # Fit biomass dynamics model
  rs = cbind(r=ri, s=si)
  k.low = max(C); k.up=max(C)*200
  opt = apply(rs, 1, function(x){optimize(BDM, c(k.low, k.up), r=x[["r"]], S=x[["s"]], b=1, C=C)})

  # Combine model results
  ki = sapply(opt, '[[', 1)
  msy = ki*ri/4
  obji = sapply(opt,'[[',2)
  kr = data.frame(k=ki, r=ri, msy, s=si, obj=obji)

  # Clean model results
  kr2 = kr
  kr2[kr2$k<1.01*k.low | kr2$k>0.99*k.up | kr2$obj>0.01,] = NA # eliminate bordering effect
  kr2 = na.omit(kr2)
  # plot(log(kr2$k), log(kr2$r), xlab="K", ylab='r')
  # abline(h=mean(log(kr$r)), v=mean(log(kr$k)), lty=2, col=2)

  # Summarize model results
  summ <- apply(kr2, 2, function(x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)))[,1:4]

  # Calculate B/BMSY time series
  ########################################

  # Empty vectors for B (one biomass traj) and Bmed (the median biomass traj)
  B = Bmed = vector()

  # Number of trjactories to simulate
  n.sim = 1000

  # Clean and sample biological quantity draws
  oc0 <- kr
  oc1 = oc0[oc0$obj<0.01 & oc0$k>1.01*min(oc0$k) & oc0$k<0.99*max(oc0$k),] #  eliminate bordering effect
  oc2 = oc1[oc1$k>quantile(oc1$k, 0.25) & oc1$k<quantile(oc1$k, 0.75),]
  smp = sample(1:nrow(oc2), n.sim, replace=T)

  # Calculate median r and k
  k = oc2[smp,1]
  r = oc2[smp,2]
  kmed = median(k)
  rmed = median(r)

  # Simulate biomass trajectory for each of the selected r-k pairs
  B <- matrix(NA, nrow=nyr, ncol=n.sim, dimnames=list(yr,1:n.sim))
  for(j in 1:n.sim){
    B[1, j] = k[j]
    for(t in 1:(nyr-1)){
      B[t+1, j] = B[t,j] + r[j]*B[t, j]*(1-B[t, j]/k[j]) - C[t]
    }
  }

  # Convert to B/BMSY
  bbmsy_trajs <- apply(B, 2, function(x) x/x[1])
  # plot(rep(0, nyr)~yr, type='n',  ylim=c(0, 1.5),xlab='', ylab='', las=1, yaxs='i')
  # for(i in 1:n.sim){
  #   lines(yr, bmsy_trajs[,i])
  # }

  # Calculate B/BMSY data frame
  quants <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  bbmsy_quants <- t(apply(bbmsy_trajs, 1, function(x) quantile(x, quants)))
  colnames(bbmsy_quants) <- paste0("q", quants)
  bbmsy_avgs <- apply(bbmsy_trajs, 1, mean)
  bbmsy_sds <- apply(bbmsy_trajs, 1, sd)
  bbmsy <- data.frame(year=yr, catch=C, bbmsy_quants, avg=bbmsy_avgs, sd=bbmsy_sds)

  # Prep for export
  krms_draws <- kr
  krms_confint <- summ
  b_trajs <- B
  out <- list(bbmsy, b_trajs, bbmsy_trajs, krms_confint, krms_draws)
  return(out)

}







