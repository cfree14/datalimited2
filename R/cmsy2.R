
# Notes
################################################################################

# The original cMSY/BSM code was integrated in the same script.
# I seperate these two approaches and simplify the number of required parameters.

# The original cMSY model had the following parameters:
#   REQUIRED
#     yr, ct, resilience (for cMSY)
#     bt, btype, (for BSM)
#     force.cmsy (removed b/c I seperate the cMSY/BSM models)
#   OPTIONAL (but useful)
#     r.low, r.hi, stb.low, stb.hi, int.yr,
#     intb.low, intb.hi, endb.low, endb.hi, q.start, q.end
#   OPTIONAL & REMOVED (for comparison only, not used in analysis)
#     Flim, Fpa, Blim, Bpa, Bmsy, FMSY, MSYBtrigger, B40, M, Fofl, SSB
#   REMOVED - NOW DERIVED FROM CATCH DATA
#     MinOfYear, MaxOfYear, StartYear, EndYear
#   REMOVED - NOT NECESSARY FOR ANALYSIS
#     Region, Subregion, Stock, Name,
#     EnglishName, ScientificName, SpecCode, Group, Source, Comment

# For testing
################################################################################

# # Packages
# library(R2jags)
# library(coda)
# library(parallel)
# library(foreach)
# library(doParallel)
# library(gplots)
#
# # For testing
# # Required parameters
# year <- SOLIRIS$yr
# catch <- SOLIRIS$ct
# biomass <- SOLIRIS$bt
# resilience <- "Medium"
# verbose <- T
# # Optional parameters
# r.low=0.18; r.hi=1.02
# stb.low=NA; stb.hi=NA; int.yr=NA;
# intb.low=NA; intb.hi=NA; endb.low=NA; endb.hi=NA; q.start=NA; q.end=NA

# Load example data
system.file("data", "SOLIRIS.Rdata", package = "datalimited2")

# SchaeferParallelSearch
################################################################################

# Monte Carlo filtering with Schaefer Function
SchaeferParallelSearch <- function(ni, nyr, sigR, duncert, ct, int.yr, intbio, startbt, ki, i, ri, int.yr.i, nstartbt, yr, end.yr, endbio, npoints, pt){
  ptm <- proc.time()
  # create vectors for viable r, k and bt
  inmemorytable <- vector()
  # parallelised for the points in the r-k space
  inmemorytable <- foreach(i = 1 : npoints, .combine='rbind', .packages='foreach', .inorder=TRUE) %dopar%{
    nsbt = length(startbt)
    VP   <- FALSE
    for(nj in 1:nsbt) {
      # create empty vector for annual biomasses
      bt <- vector()
      j<-startbt[nj]
      # set initial biomass, including 0.1 process error to stay within bounds
      bt[1]=j*ki[i]*exp(rnorm(1,0, 0.1*sigR))  ## set biomass in first year
      # repeat test of r-k-startbt combination to allow for different random error
      for(re in 1:ni)   {
        #loop through years in catch time series
        for (t in 1:nyr)  {  # for all years in the time series
          xt=rnorm(1,0, sigR) # set new process error for every year
          zlog.sd = sqrt(log(1+(duncert)^2))
          zt=rlnorm(1,meanlog = 0, sdlog = zlog.sd) # model the catch error as a log normal distribution.
          # calculate biomass as function of previous year's biomass plus surplus production minus catch
          bt[t+1] <- ifelse(bt[t]/ki[i] >= 0.25,
                            bt[t]+ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt,
                            bt[t]+(4*bt[t]/ki[i])*ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt) # assuming reduced r at B/k < 0.25
          # if biomass < 0.01 k, discard r-k-startbt combination
          if(bt[t+1] < 0.01*ki[i]){
            break
          } # stop looping through years, go to next upper level
          # intermediate year check
          if((t+1)==int.yr.i && (bt[t+1]>(intbio[2]*ki[i]) || bt[t+1]<(intbio[1]*ki[i]))){
            break
          }
        } # end of loop of years
        # if loop was broken or last biomass falls outside of expected ranges
        # do not store results, go directly to next startbt
        if(t < nyr || bt[yr==end.yr] > (endbio[2]*ki[i]) || bt[yr==end.yr] < (endbio[1]*ki[i])){
          next
        }else{
          #each vector will be finally appended to the others found by the threads - this is done by the .combine='rbind' option
          inmemorytablerow<-c(i,j,ri[i],ki[i],bt[1:(nyr+1)]/ki[i])
          if(length(inmemorytablerow)==(4+nyr+1)){
            if(VP==FALSE){
              inmemorytable <- inmemorytablerow
            }else{
              inmemorytable <- rbind(inmemorytable,inmemorytablerow)
            }
            VP <- TRUE
          }
        }
      } # end of repetition for random error
    } # end of j-loop of initial biomasses
    # instruction necessary to make the foreach loop see the variable:
    if(length(inmemorytable)==0){
      inmemorytable <- vector(length=4+nyr+1)*NA
    }else{
      inmemorytable
    }
  }#end loop on points

  #create the output matrix
  mdat        <- matrix(data=NA, nrow = npoints*nstartbt, ncol = 2+nyr+1)
  npointsmem = dim(inmemorytable)[1]
  npointscols = dim(inmemorytable)[2]
  #reconstruction of the processing matrix after the parallel search
  if (npointsmem>0 && npointscols>0){
    for (idxr in 1:npointsmem){
      i = inmemorytable[idxr,1]
      if (!is.na(i)){
        j = inmemorytable[idxr,2]
        mdatindex<-((i-1)*nstartbt)+which(startbt==j)
        mdat[mdatindex,1]           <- inmemorytable[idxr,3]
        mdat[mdatindex,2]           <- inmemorytable[idxr,4]
        mdat[mdatindex,3:(2+nyr+1)] <- inmemorytable[idxr,5:(4+nyr+1)]
        # if(pt==T) points(x=ri[i], y=ki[i], pch=".", cex=4, col="gray")
      }
    }
  }
  ptm<-proc.time()-ptm
  mdat <- na.omit(mdat)
  return(mdat)
}

# SchaeferMC
################################################################################

SchaeferMC <- function(ri, ki, startbio, int.yr, intbio, endbio, sigR, pt, duncert, startbins, ni, yr, nyr, ct, end.yr, verbose){
  # create vector for initial biomasses
  startbt <- seq(from =startbio[1], to=startbio[2], by = (startbio[2]-startbio[1])/startbins)
  nstartbt <- length(startbt)
  npoints <- length(ri)
  # get index of intermediate year
  int.yr.i <- which(yr==int.yr)
  #loop through r-k pairs with parallel search
  mdat<-SchaeferParallelSearch(ni, nyr,sigR,duncert,ct,int.yr,intbio, startbt, ki, i, ri, int.yr.i, nstartbt, yr, end.yr, endbio, npoints,pt)
  if(verbose==T){cat("\n")}
  return(list(mdat))
}

# Moving average function
################################################################################

# Calculate moving average
ma <- function(x){
  x.1 <- stats::filter(x,rep(1/3,3),sides=1)
  x.1[1] <- x[1]
  x.1[2] <- (x[1]+x[2])/2
  return(x.1)
}

# Functions for setting priors
################################################################################

# Set R prior
r_prior <- function(r.low, r.hi){
  # initial range of r from input file
  if(is.na(r.low)==F & is.na(r.hi)==F) {
    start.r <- c(r.low,r.hi)
  }else{
    # initial range of r based on resilience
    if(res == "High") {
      start.r <- c(0.6,1.5)} else if(res == "Medium") {
        start.r <- c(0.2,0.8)}    else if(res == "Low") {
          start.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
            start.r <- c(0.015,0.1)}
  }
  return(start.r)
}

# Set start saturation prior
startbio_prior <- function(stb.low, stb.hi, start.yr){
  # use initial biomass range from input file if stated
  if(is.na(stb.low)==F & is.na(stb.hi)==F){
    startbio <- c(stb.low,stb.hi)
  }else{
    # if start year < 1960 assume high biomass
    if(start.yr < 1960){
      startbio <- c(0.5,0.9)
    }else{
      # else use medium prior biomass range
      startbio <- c(0.2,0.6)
    }
  }
  return(startbio)
}

# Set intermediate saturation prior
intbio_prior <- function(intb.low, intb.hi, int.yr, start.yr, end.yr, startbio, yr, ct){
  # get index of years with lowest and highest catch between start+3 and end-3 years
  min.yr.i <- which.min(ct[4:(length(ct)-3)])+3
  max.yr.i <- which.max(ct[4:(length(ct)-3)])+3
  min.ct <- ct[min.yr.i]
  max.ct <- ct[max.yr.i]
  # use year and biomass range for intermediate biomass from input file
  if(is.na(intb.low)==F & is.na(intb.hi)==F){
    int.yr   <- int.yr
    intbio   <- c(intb.low,intb.hi)
    # if contrast in catch is low, use initial range again in mid-year
  }else if(min(ct)/max(ct) > 0.6) {
    int.yr    <- as.integer(mean(c(start.yr, end.yr)))
    intbio    <- startbio
    # else if year of minimum catch is after max catch then use min catch
  }else if(min.yr.i > max.yr.i) {
    int.yr    <- yr[min.yr.i-1]
    if(startbio[1]>=0.5 &  (int.yr-start.yr) < (end.yr-int.yr) &
       (min.ct/max.ct) > 0.3) intbio <- c(0.2,0.6) else intbio <- c(0.01,0.4)
       # else use max catch
  } else {
    # assume that biomass range in year before maximum catch was high or medium
    int.yr    <- yr[max.yr.i-1]
    intbio    <- if((startbio[1]>=0.5 & (int.yr-start.yr) < (end.yr-int.yr))| # if initial biomass is high, assume same for intermediate
                    # ((min.ct/max.ct < 0.3 & (max.yr.i - min.yr.i) < 25))) c(0.5,0.9) else c(0.2,0.6) }
                    (((max.ct-min.ct)/max.ct)/(max.yr.i-min.yr.i) > 0.04)) c(0.5,0.9) else c(0.2,0.6) } # if incease is steep, assume high, else medium
  out <- list(intbio, int.yr)
  return(out)
}

# Set end saturation prior
endbio_prior <- function(endb.low, endb.hi, nyr, ct.raw, ct){
  # final biomass range from input file
  if(is.na(endb.low)==F & is.na(endb.hi)==F){
    endbio   <- c(endb.low,endb.hi)
  }else{
    # else use mean final catch/max catch to estimate final biomass
    rawct.ratio=ct.raw[nyr]/max(ct)
    endbio  <- if(ct[nyr]/max(ct) > 0.8) {c(0.4,0.8)} else if(rawct.ratio < 0.5) {c(0.01,0.4)} else {c(0.2,0.6)}

    # if default endbio is low (0.01-0.4), check whether the upper bound should be lower than 0.4 for depleted stocks
    if(endbio[2]==0.4){
      if(rawct.ratio< 0.05) {endbio[2] <- 0.1} else
        if(rawct.ratio< 0.15) {endbio[2] <- 0.2} else
          if(rawct.ratio< 0.35) {endbio[2] <- 0.3} else {endbio[2] <- 0.4}
    }
  }
  return(endbio)
}

# Set K prior
k_prior <- function(endbio, start.r, ct){
  # initial prior range of k values, assuming min k will be larger than max catch / prior for r
  if(mean(endbio) <= 0.5){
    start.k <- c(max(ct)/start.r[2],4*max(ct)/start.r[1])
  }else{
    start.k <- c(2*max(ct)/start.r[2],12*max(ct)/start.r[1])
  }
  return(start.k)
}

# Fit cMSY-17 model
################################################################################

#' cMSY catch-only stock assessment model
#'
#' Estimates B/BMSY time series and other biological quantities using only
#' a time series of catch and a resilience estimate using cMSY from Froese et al. (2017).
#'
#' @param year A time series of years
#' @param catch A time series of catch
#' @param resilience Resilience of the stock: "High", "Medium", "Low", "Very low"
#' @param r.low/r.hi A user-specified prior on the species intrinsic growth rate, r (optional)
#' @param stb.low/stb.hi A user-specified prior on biomass relative to unfished biomass at the beginning of the catch time series (optional)
#' @param int.yr A user-specified year of intermediate biomass (optional)
#' @param intb.low/intb.hi A user-specified prior on biomass relative to unfished biomass in the intermediate year (optional)
#' @param endb.low/endb.hi A user-specified prior on biomass relative to unfished biomass at the end of the catch time series (optional)
#' @param q.start/q.end A user-specified start and end year for estimating the catchability coefficient (optional; default is last 5 years)
#' @param verbose Set to FALSE to suppress printed updates on CMSY/BSM progress (default=TRUE)
#' @return A list of length six with the following elements
#' \itemize{
#'   \item{A dataframe with biological quantity / reference point estimates with 95% confidence intervals}
#'   \item{A dataframe with B/BMSY and reference point time series with 95% confidence intervals}
#'   \item{A dataframe with the priors used in the cMSY analysis}
#'   \item{A vector with the viable r values}
#'   \item{A vector with the viable k values}
#'   \item{A vector with the viable saturation values}
#' }
#' @references Froese R, Demirel N, Coro G, Kleisner KM, Winker H (2017)
#' Estimating fisheries reference points from catch and resilience. Fish and Fisheries 18(3): 506-526.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/faf.12190/abstract}
#' @examples
#' output <- cmsy2(year=SOLIRIS$yr, catch=SOLIRIS$ct, r.low=0.18, r.hi=1.02)
#' plot_cmsy2(output)
#' plot_cmsy2_mgmt(output)
#' @export
cmsy2 <- function(year, catch, resilience=NA,
                  r.low=NA, r.hi=NA, stb.low=NA, stb.hi=NA, int.yr=NA,
                  intb.low=NA, intb.hi=NA, endb.low=NA, endb.hi=NA, q.start=NA, q.end=NA, verbose=T){


  # Set model parameters
  #############################################################

  # # Setup parallel processing
  # # Use 3 chains in JAGS if more than 2 cores are available
  # n.cores <- detectCores()
  # n.chains <- ifelse(n.cores > 2,3,2)
  # cl <- makeCluster(n.cores)
  # registerDoParallel(cl, cores = n.cores)

  # Set model parameters
  FullSchaefer <- F # will automatically change to TRUE if enough abundance data available
  dataUncert <- 0.1  # set observation error as uncertainty in catch - default is SD=0.1
  sigmaR <- 0.1 # overall process error for CMSY; SD=0.1 is the default
  n <- 10000 # initial number of r-k pairs
  n.new <- n # initialize n.new
  ni <- 3 # iterations for r-k-startbiomass combinations, to test different variability patterns; no improvement seen above 3
  nab <- 5 # default=5; minimum number of years with abundance data to run BSM
  duncert <- dataUncert # global defaults for uncertainty
  sigR <- sigmaR # global defaults for uncertainty

  # Setup data
  #############################################################

  # Build catch data (using original cMSY variable naming convention)
  catchData <- data.frame(yr=year, ct=catch)

  # Transform catch data
  # 1. Convert to 1000s tons (or other units)
  # 2. Calculate 3-yr moving average (average of past 3 years)
  ct.raw <- catchData$ct / 1000
  if(is.na(mean(ct.raw))){
    cat("ERROR: Missing value in Catch data; fill or interpolate\n")
  }
  ct <- ma(ct.raw)

  # Identify number of years and start/end years
  yr <- catchData$yr # functions use this quantity
  nyr <- length(yr)
  start.yr <- min(yr)
  end.yr <- max(yr)

  # Determine initial ranges for parameters and biomass
  #############################################################

  # Set priors
  res <- resilience # rename resilience
  start.r <- r_prior(r.low, r.hi)
  startbio <- startbio_prior(stb.low, stb.hi, start.yr)
  int_params <- intbio_prior(intb.low, intb.hi, int.yr, start.yr, end.yr, startbio, yr, ct)
  intbio <- int_params[[1]]
  int.yr <- int_params[[2]]
  endbio <- endbio_prior(endb.low, endb.hi, nyr, ct.raw, ct)
  start.k <- k_prior(endbio, start.r, ct)

  # Record priors into dataframe
  priors <- data.frame(cbind(c("r", "k", "startbio", "intbio", "endbio"), source="default",
                             rbind(start.r, start.k, startbio, intbio, endbio)), year=NA, stringsAsFactors=F)
  colnames(priors) <- c("param", "source", "lo", "hi", "year")
  rownames(priors) <- NULL
  priors$year[priors$param=="intbio"] <- int.yr
  priors$lo <- as.numeric(priors$lo)
  priors$hi <- as.numeric(priors$hi)
  if(!is.na(r.low)){priors$source[priors$param=="r"] <- "expert"}
  if(!is.na(stb.low)){priors$source[priors$param=="startbio"] <- "expert"}
  if(!is.na(intb.low)){priors$source[priors$param=="intbio"] <- "expert"}
  if(!is.na(endb.low)){priors$source[priors$param=="endbio"] <- "expert"}

  # Print priors (if desired)
  if(verbose==T){
    cat("startbio=",startbio,ifelse(is.na(stb.low)==T,"default","expert"),
        ", intbio=",int.yr,intbio,ifelse(is.na(intb.low)==T,"default","expert"),
        ", endbio=",endbio,ifelse(is.na(endb.low)==T,"default","expert"),"\n")
  }

  # Monte Carlo procedure
  #############################################################

  # Initialize other vectors anew for each stock
  current.attempts <- NA

  # Initialize vectors for viable r, k, bt, and all in a matrix
  mdat.all <- matrix(data=vector(), ncol=2+nyr+1)

  # Get random set of r and k from log space distribution
  ri1 = exp(runif(n, log(start.r[1]), log(start.r[2])))
  ki1 = exp(runif(n, log(start.k[1]), log(start.k[2])))

  # 1 - Call CMSY-SchaeferMC function to preliminary explore the r-k space
  if(verbose==T){cat("First Monte Carlo filtering of r-k space with ",n," points...\n")}
  MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR,
                     pt=T, duncert=dataUncert, startbins=10, ni=ni, yr=yr, nyr=nyr, ct=ct, end.yr=end.yr, verbose=verbose)
  mdat.all <- rbind(mdat.all,MCA[[1]])
  rv.all   <- mdat.all[,1]
  kv.all   <- mdat.all[,2]
  btv.all  <- mdat.all[,3:(2+nyr+1)]
  # count viable trajectories and r-k pairs
  n.viable.b   <- length(mdat.all[,1])
  n.viable.pt <- length(unique(mdat.all[,1]))
  if(verbose==T){cat("Found ",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")}

  # 2 - if the lower bound of k is too high, reduce it by half and rerun
  if(length(kv.all[kv.all < 1.1*start.k[1] & rv.all < mean(start.r)]) > 10) {
    if(verbose==T){cat("Reducing lower bound of k, resampling area with",n,"additional points...\n")}
    start.k <- c(0.5*start.k[1],start.k[2])
    ri1 = exp(runif(n, log(start.r[1]), log(start.r[2])))
    ki1 = exp(runif(n, log(start.k[1]), log(start.k[2])))
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR,
                       pt=T, duncert=dataUncert, startbins=10, ni=ni, yr=yr, nyr=nyr, ct=ct, end.yr=end.yr, verbose=verbose)
    mdat.all <- rbind(mdat.all,MCA[[1]])
    rv.all   <- mdat.all[,1]
    kv.all   <- mdat.all[,2]
    btv.all  <- mdat.all[,3:(2+nyr+1)]
    n.viable.b   <- length(mdat.all[,1])
    n.viable.pt <- length(unique(mdat.all[,1]))
    if(verbose==T){cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")}
  }

  # 3 - if few points were found then resample and shrink the log k space
  if (n.viable.b <= 1000){
    log.start.k.new  <- log(start.k)
    max.attempts     <- 3
    current.attempts <- 1
    startbins        <- 10
    while (n.viable.b <= 1000 && current.attempts <= max.attempts){
      if(n.viable.pt > 0) {
        log.start.k.new[1] <- mean(c(log(start.k[1]), min(log(kv.all))))
        log.start.k.new[2] <- mean(c(log.start.k.new[2], max(log(kv.all)))) }
      n.new <- n*current.attempts #add more points
      ri1 = exp(runif(n.new, log(start.r[1]), log(start.r[2])))
      ki1 = exp(runif(n.new, log.start.k.new[1], log.start.k.new[2]))
      if(verbose==T){cat("Shrinking k space: repeating Monte Carlo in the interval [",exp(log.start.k.new[1]),",",exp(log.start.k.new[2]),"]\n")}
      if(verbose==T){cat("Attempt ",current.attempts," of ",max.attempts," with ",n.new," additional points...","\n")}
      if(current.attempts==2 & n.viable.b < 50){
        duncert   <- 2*dataUncert
        sigR      <- 2*sigmaR
        startbins <- 20
        if(verbose==T){cat("Doubling startbins, catch and process error, and number of variability patterns \n")}
      }
      MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR,
                         pt=T, duncert=duncert, startbins=startbins, ni=2*ni, yr=yr, nyr=nyr, ct=ct, end.yr=end.yr, verbose=verbose)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      if(verbose==T){cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")}
      current.attempts=current.attempts+1 #increment the number of attempts
    }
    if(n.viable.b < 5) {
      if(verbose==T){cat("Only",n.viable.pt,"viable r-k pairs found, check data and settings \n")}
      next
    }
  }

  # 4 - if tip of viable r-k pairs is 'thin', do extra sampling there
  if(length(rv.all[rv.all > 0.9*start.r[2]]) < 5) {
    l.sample.r        <- quantile(rv.all,0.6)
    add.points        <- ifelse(is.na(current.attempts)==T,n,ifelse(current.attempts==2,2*n,ifelse(length(rv.all)>500,3*n,6*n)))
    if(verbose==T){cat("Final sampling in the tip area above r =",l.sample.r,"with",add.points,"additional points...\n")}
    log.start.k.new <- c(log(0.8*min(kv.all)),log(max(kv.all[rv.all > l.sample.r])))

    ri1 = exp(runif(add.points, log(l.sample.r), log(start.r[2])))
    ki1 = exp(runif(add.points, log.start.k.new[1], log.start.k.new[2]))
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR,
                       pt=T, duncert=duncert, startbins=10, ni=ni, yr=yr, nyr=nyr, ct=ct, end.yr=end.yr, verbose=verbose)
    mdat.all <- rbind(mdat.all,MCA[[1]])
    rv.all   <- mdat.all[,1]
    kv.all   <- mdat.all[,2]
    btv.all  <- mdat.all[,3:(2+nyr+1)]
    n.viable.b   <- length(mdat.all[,1])
    n.viable.pt <- length(unique(mdat.all[,1]))
    if(verbose==T){cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")}
  }

  # Extract model results
  #############################################################

  # get estimate of most probable r as 75th percentile of mid log.r-classes
  # get unique combinations of r-k
  unique.rk         <- unique(mdat.all[,1:2])
  # get remaining viable log.r and log.k
  log.rs           <- log(unique.rk[,1])
  log.ks           <- log(unique.rk[,2])
  # get vectors with numbers of r and mid values in classes
  # determine number of classes as a function of r-width
  r.width         <- (max(unique.rk[,1])-start.r[1])/(start.r[2]-start.r[1])
  classes         <- ifelse(r.width>0.8,100,ifelse(r.width>0.5,50,ifelse(r.width>0.3,25,12)))
  hist.log.r      <- hist(x=log.rs, breaks=classes, plot=F)
  log.r.counts    <- hist.log.r$counts
  log.r.mids      <- hist.log.r$mids
  # get most probable log.r as 75th percentile of mids with counts > 0
  log.r.est       <- as.numeric(quantile(log.r.mids[which(log.r.counts > 0)],0.75))
  median.log.r    <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.50))
  lcl.log.r       <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.5125))
  ucl.log.r       <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.9875))
  sd.log.r.est    <- (ucl.log.r - log.r.est) / 1.96
  r.est           <- exp(log.r.est)
  lcl.r.est       <- exp(log.r.est-1.96*sd.log.r.est)
  ucl.r.est       <- exp(log.r.est+1.96*sd.log.r.est)

  # get r-k pairs above median of mids
  rem            <- which(unique.rk[,1] > exp(median.log.r))
  rem.log.r      <- log(unique.rk[,1][rem])
  rem.log.k      <- log(unique.rk[,2][rem])
  # do linear regression of log k ~ log r with slope fixed to -1 (from Schaefer)
  reg            <- lm(rem.log.k ~ 1 + offset(-1*rem.log.r))
  int.reg        <- as.numeric(reg[1])
  sd.reg      <- sd(resid(reg))
  # get estimate of log(k) from y where x = log.r.est
  log.k.est      <- int.reg + (-1) * log.r.est
  # get estimates of ucl of log.k.est from y + SD where x = ucl.log.r
  ucl.log.k     <- int.reg + (-1) * lcl.log.r + sd.reg
  # get estimates of sd.log.k.est from upper confidence limit of log.k.est
  sd.log.k.est   <- (ucl.log.k - log.k.est) / 1.96
  lcl.log.k      <- log.k.est - 1.96*sd.log.k.est
  ucl.log.k      <- log.k.est + 1.96*sd.log.k.est
  k.est       <- exp(log.k.est)
  lcl.k.est   <- exp(lcl.log.k)
  ucl.k.est   <- exp(ucl.log.k)

  # get MSY from remaining log r-k pairs
  log.MSY.est     <- mean(rem.log.r + rem.log.k - log(4))
  sd.log.MSY.est  <- sd(rem.log.r + rem.log.k - log(4))
  lcl.log.MSY.est <- log.MSY.est - 1.96*sd.log.MSY.est
  ucl.log.MSY.est <- log.MSY.est + 1.96*sd.log.MSY.est
  MSY.est         <- exp(log.MSY.est)
  lcl.MSY.est     <- exp(lcl.log.MSY.est)
  ucl.MSY.est     <- exp(ucl.log.MSY.est)

  # get predicted biomass vectors as median and quantiles
  # only use biomass trajectories from r-k pairs within the confidence limits
  rem.btv.all <- mdat.all[which(mdat.all[,1] > lcl.r.est & mdat.all[,1] < ucl.r.est
                                & mdat.all[,2] > lcl.k.est & mdat.all[,2] < ucl.k.est),3:(2+nyr+1)]
  median.btv <- apply(rem.btv.all,2, median)
  median.btv.lastyr  <- median.btv[length(median.btv)-1]
  nextyr.bt  <- median.btv[length(median.btv)]
  lcl.btv    <- apply(rem.btv.all,2, quantile, probs=0.025)
  q.btv      <- apply(rem.btv.all,2, quantile, probs=0.25)
  ucl.btv    <- apply(rem.btv.all,2, quantile, probs=0.975)
  lcl.median.btv.lastyr <- lcl.btv[length(lcl.btv)-1]
  ucl.median.btv.lastyr <- ucl.btv[length(lcl.btv)-1]
  lcl.nextyr.bt <- lcl.btv[length(lcl.btv)]
  ucl.nextyr.bt <- ucl.btv[length(lcl.btv)]

  # get F derived from predicted CMSY biomass
  F.CMSY      <- ct.raw/(median.btv[1:nyr]*k.est)
  Fmsy.CMSY  <- r.est/2 # Fmsy from CMSY


  # Extract management quantities
  ##################################################

  MSY   <-MSY.est; lcl.MSY<-lcl.MSY.est; ucl.MSY<-ucl.MSY.est
  Bmsy  <-k.est/2; lcl.Bmsy<-lcl.k.est/2; ucl.Bmsy<-ucl.k.est/2
  Fmsy  <-r.est/2; lcl.Fmsy<-lcl.r.est/2; ucl.Fmsy<-ucl.r.est/2
  B.Bmsy<-2*median.btv[1:nyr]; lcl.B.Bmsy<-2*lcl.btv[1:nyr]; ucl.B.Bmsy<-2*ucl.btv[1:nyr]
  B          <-B.Bmsy*Bmsy; lcl.B<-lcl.B.Bmsy*Bmsy; ucl.B<-ucl.B.Bmsy*Bmsy
  B.last     <-B[nyr]; lcl.B.last<-lcl.B[nyr]; ucl.B.last<-ucl.B[nyr]
  B.Bmsy.last<-B.Bmsy[nyr]; lcl.B.Bmsy.last<-lcl.B.Bmsy[nyr]; ucl.B.Bmsy.last<-ucl.B.Bmsy[nyr]
  Fm           <- ct.raw/B; lcl.F<-ct.raw/ucl.B; ucl.F<-ct.raw/lcl.B
  Fmsy.vec     <- ifelse(B.Bmsy>0.5,Fmsy,Fmsy*2*B.Bmsy)
  lcl.Fmsy.vec <- ifelse(B.Bmsy>0.5,lcl.Fmsy,lcl.Fmsy*2*B.Bmsy)
  ucl.Fmsy.vec <- ifelse(B.Bmsy>0.5,ucl.Fmsy,ucl.Fmsy*2*B.Bmsy)
  F.Fmsy       <- Fm/Fmsy.vec; lcl.F.Fmsy<-lcl.F/Fmsy.vec; ucl.F.Fmsy<-ucl.F/Fmsy.vec
  F.last     <-Fm[nyr];lcl.F.last<-lcl.F[nyr];ucl.F.last<-ucl.F[nyr]
  Fmsy.last  <-Fmsy.vec[nyr];lcl.Fmsy.last<-lcl.Fmsy.vec[nyr];ucl.Fmsy.last<-ucl.Fmsy.vec[nyr]
  F.Fmsy.last<-F.Fmsy[nyr];lcl.F.Fmsy.last<-lcl.F.Fmsy[nyr];ucl.F.Fmsy.last<-ucl.F.Fmsy[nyr]

  # Print results (if desired)
  ##################################################

  # Print results
  if(verbose==T){

    # Priors
    cat("---------------------------------------\n")
    cat("Catch data used from years", min(yr),"-", max(yr), "\n")
    cat("Prior initial relative biomass =", startbio[1], "-", startbio[2],ifelse(is.na(stb.low)==T,"default","expert"), "\n")
    cat("Prior intermediate rel. biomass=", intbio[1], "-", intbio[2], "in year", int.yr,ifelse(is.na(intb.low)==T,"default","expert"), "\n")
    cat("Prior final relative biomass   =", endbio[1], "-", endbio[2],ifelse(is.na(endb.low)==T,"default","expert"), "\n")
    cat("Prior range for r =", format(start.r[1],digits=2), "-", format(start.r[2],digits=2),ifelse(is.na(r.low)==T,"default","expert,"),
        ", prior range for k =", start.k[1], "-", start.k[2],"\n")

    # cMSY results
    cat("\nResults of CMSY analysis \n")
    cat("-------------------------\n")
    cat("Altogether", n.viable.b, "viable trajectories for", n.viable.pt," r-k pairs were found \n")
    cat("r   =", r.est,", 95% CL =", lcl.r.est, "-", ucl.r.est,", k =", k.est,", 95% CL =", lcl.k.est, "-", ucl.k.est,"\n")
    cat("MSY =", MSY.est,", 95% CL =", lcl.MSY.est, "-", ucl.MSY.est,"\n")
    cat("Relative biomass in last year =", median.btv.lastyr, "k, 2.5th perc =", lcl.median.btv.lastyr,
        ", 97.5th perc =", ucl.median.btv.lastyr,"\n")
    cat("Exploitation F/(r/2) in last year =", (F.CMSY/Fmsy.CMSY)[nyr],"\n\n")

    # Mgmt results
    cat("-------------------------------------------------------------\n")
    cat("Fmsy =",Fmsy,", 95% CL =",lcl.Fmsy,"-",ucl.Fmsy,"(if B > 1/2 Bmsy then Fmsy = 0.5 r)\n")
    cat("Fmsy =",Fmsy.last,", 95% CL =",lcl.Fmsy.last,"-",ucl.Fmsy.last,"(r and Fmsy are linearly reduced if B < 1/2 Bmsy)\n")
    cat("MSY  =",MSY,", 95% CL =",lcl.MSY,"-",ucl.MSY,"\n")
    cat("Bmsy =",Bmsy,", 95% CL =",lcl.Bmsy,"-",ucl.Bmsy,"\n")
    cat("Biomass in last year =",B.last,", 2.5th perc =", lcl.B.last, ", 97.5 perc =",ucl.B.last,"\n")
    cat("B/Bmsy in last year  =",B.Bmsy.last,", 2.5th perc =", lcl.B.Bmsy.last, ", 97.5 perc =",ucl.B.Bmsy.last,"\n")
    cat("Fishing mortality in last year =",F.last,", 2.5th perc =", lcl.F.last, ", 97.5 perc =",ucl.F.last,"\n")
    cat("Exploitation F/Fmsy  =",F.Fmsy.last,", 2.5th perc =", lcl.F.Fmsy.last, ", 97.5 perc =",ucl.F.Fmsy.last,"\n")

  }

  # Combine results for export
  ##################################################

  # Create output object
  ref_pts <- data.frame(rbind(c(r.est, lcl.r.est, ucl.r.est),
                              c(k.est, lcl.k.est, ucl.k.est),
                              c(MSY.est, lcl.MSY.est, ucl.MSY.est),
                              c(Fmsy, lcl.Fmsy, ucl.Fmsy),
                              c(Bmsy, lcl.Bmsy, ucl.Bmsy)))
  colnames(ref_pts) <- c("est", "lo", "hi")
  ref_pts$param <- c("r", "k", "msy", "fmsy", "bmsy")
  ref_pts <- subset(ref_pts, select=c(param, est, lo, hi))

  # Build reference point time series
  ref_ts <- data.frame(year=year, catch=ct.raw, catch_ma=ct,
                       b=B, b_lo=lcl.B, b_hi=ucl.B,
                       bbmsy=B.Bmsy, bbmsy_lo=lcl.B.Bmsy, bbmsy_hi=ucl.B.Bmsy,
                       s=B.Bmsy/2, s_lo=lcl.B.Bmsy/2, s_hi=ucl.B.Bmsy/2,
                       f=Fm, f_lo=lcl.F, f_hi=ucl.F,
                       fmsy=Fmsy.vec, fmsy_lo=lcl.Fmsy.vec, fmsy_hi=ucl.Fmsy.vec,
                       ffmsy=F.Fmsy, ffmsy_lo=lcl.F.Fmsy, ffmsy_hi=ucl.F.Fmsy,
                       er=F.CMSY/Fmsy.CMSY)

  # #stop parallel processing clusters
  # stopCluster(cl)
  # stopImplicitCluster()

  # Return
  out <- list(ref_pts=ref_pts, ref_ts=ref_ts, priors=priors,
              rv.all=rv.all, kv.all=kv.all, btv.all=btv.all)
  return(out)

}









