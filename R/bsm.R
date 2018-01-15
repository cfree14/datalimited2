
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
# load("data/SOLIRIS.Rda")
# year <- SOLIRIS$yr
# catch <- SOLIRIS$ct
# biomass <- SOLIRIS$bt
# btype <- "CPUE"
# resilience <- "Medium"
# verbose <- T
# # Optional parameters
# r.low=0.18; r.hi=1.02
# stb.low=NA; stb.hi=NA; int.yr=NA;
# intb.low=NA; intb.hi=NA; endb.low=NA; endb.hi=NA; q.start=NA; q.end=NA

# BSM function
################################################################################

#' Bayesian state-space surplus production model
#'
#' Estimates B/BMSY time series and other biological quantities using only
#' a time series of catch and a resilience estimate using the Bayesian surplus
#' produciton model from Froese et al. (2017).
#'
#' @param year A time series of years
#' @param catch A time series of catch
#' @param biomass A time series of biomass or CPUE (type is designated in btype)
#' @param btype Biomass time series type: "None", "biomass", or "CPUE"
#' @param resilience Resilience of the stock: "High", "Medium", "Low", "Very low"
#' @param r.low,r.hi A user-specified prior on the species intrinsic growth rate, r (optional)
#' @param stb.low,stb.hi A user-specified prior on biomass relative to unfished biomass at the beginning of the catch time series (optional)
#' @param int.yr A user-specified year of intermediate biomass (optional)
#' @param intb.low,intb.hi A user-specified prior on biomass relative to unfished biomass in the intermediate year (optional)
#' @param endb.low,endb.hi A user-specified prior on biomass relative to unfished biomass at the end of the catch time series (optional)
#' @param q.start,q.end A user-specified start and end year for estimating the catchability coefficient (optional; default is last 5 years)
#' @param verbose Set to FALSE to suppress printed updates on CMSY/BSM progress (default=TRUE)
#' @return A time series of B/BSMY estimates and other stuff
#' @references Froese R, Demirel N, Coro G, Kleisner KM, Winker H (2017)
#' Estimating fisheries reference points from catch and resilience. Fish and Fisheries 18(3): 506-526.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/faf.12190/abstract}
#' @examples
#' output <- bsm(year=SOLIRIS$yr, catch=SOLIRIS$ct, biomass=SOLIRIS$bt, btype="CPUE", r.low=0.18, r.hi=1.02)
#' plot_cmsy2(output)
#' plot_cmsy2_mgmt(output)
#' @export
bsm <- function(year, catch, biomass, btype, resilience=NA,
                r.low=NA, r.hi=NA, stb.low=NA, stb.hi=NA, int.yr=NA,
                intb.low=NA, intb.hi=NA, endb.low=NA, endb.hi=NA, q.start=NA, q.end=NA, verbose=T){

  # Set model parameters
  ##############################################################################

  # Setup parallel processing
  # Use 3 chains in JAGS if more than 2 cores are available
  n.cores <- parallel::detectCores()
  n.chains <- ifelse(n.cores > 2,3,2)
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl, cores = n.cores)

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
  ##############################################################################

  # Build catch data (using original cMSY variable naming convention)
  catchData <- data.frame(yr=year, ct=catch, bt=biomass)

  # Transform catch data
  # 1. Convert to 1000s tons (or other units)
  # 2. Calculate 3-yr moving average (average of past 3 years)
  ct.raw <- catchData$ct / 1000
  if(is.na(mean(ct.raw))){
    cat("ERROR: Missing value in Catch data; fill or interpolate\n")
  }
  ct <- ma(ct.raw)

  # Transform biomass data
  if(btype=="biomass" | btype=="CPUE"){
    bt <- catchData$bt / 1000  ## assumes that biomass is in tonnes, transforms to '000 tonnes
  }else{
    bt <- NA
  }

  # Identify number of years and start/end years
  yr <- catchData$yr # functions use this quantity
  nyr <- length(yr)
  start.yr <- min(yr)
  end.yr <- max(yr)

  # Determine initial ranges for parameters and biomass
  ##############################################################################

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


  # Bayesian analysis of catch & biomass (or CPUE) with Schaefer model
  ##############################################################################

  # Fit BSM if enough biomass data is available
  if(btype != "None" & length(bt[is.na(bt)==F])>=nab){

    # Indicate that BSM is being fit
    FullSchaefer <- T

    # Set inits for r-k in lower right corner of log r-k space
    init.r      <- start.r[1]+0.8*(start.r[2]-start.r[1])
    init.k      <- start.k[1]+0.1*(start.k[2]-start.k[1])

    # Vector with no penalty (=0) if predicted biomass is within viable range, else a penalty of 10 is set
    pen.bk = pen.F = rep(0,length(ct))

    # Add biomass priors
    b.yrs = c(1,length(start.yr:int.yr),length(start.yr:end.yr))
    b.prior = rbind(matrix(c(startbio[1],startbio[2],intbio[1],intbio[2],endbio[1],endbio[2]),2,3),rep(0,3)) # last row includes the 0 pen

    # Biomass-based BSM
    ##########################################################

    if(verbose==T){cat("Running MCMC analysis....\n")}
    if(btype == "biomass"){
      # Data to be passed on to JAGS
      jags.data        <- c('ct','bt','nyr', 'start.r','startbio','start.k',
                            'init.r','init.k', 'pen.bk','pen.F','b.yrs','b.prior')
      # Parameters to be returned by JAGS
      jags.save.params <- c('r','k','P') #

      # JAGS model
      Model = "model{
        # to avoid crash due to 0 values
        eps<- 0.01
        penm[1] <- 0 # no penalty for first biomass
        Pmean[1] <- log(alpha)
        P[1] ~ dlnorm(Pmean[1],itau2)

        for (t in 2:nyr) {
        Pmean[t] <- ifelse(P[t-1] > 0.25,
        log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps)),  # Process equation
        log(max(P[t-1] + 4*P[t-1]*r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps))) # assuming reduced r at B/k < 0.25
        P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
        penm[t]  <- ifelse(P[t]<(eps+0.001),log(k*P[t])-log(k*(eps+0.001)),ifelse(P[t]>1,log(k*P[t])-log(k*(0.99)),0)) # penalty if Pmean is outside viable biomass

        }

        # ><> Biomass priors/penalties are enforced as follows
        for (i in 1:3) {
        penb[i]  <- ifelse(P[b.yrs[i]]<b.prior[1,i],log(k*P[b.yrs[i]])-log(k*b.prior[1,i]),ifelse(P[b.yrs[i]]>b.prior[2,i],log(k*P[b.yrs[i]])-log(k*b.prior[2,i]),0))
        b.prior[3,i] ~ dnorm(penb[i],100)
        }

        for (t in 1:nyr){
        Fpen[t]   <- ifelse(ct[t]>(0.9*k*P[t]),ct[t]-(0.9*k*P[t]),0) #><> Penalty term on F > 1, i.e. ct>B
        pen.F[t]  ~ dnorm(Fpen[t],1000)
        pen.bk[t] ~ dnorm(penm[t],10000)
        Bm[t] <- log(P[t]*k);
        bt[t]    ~ dlnorm(Bm[t],isigma2);
        }

        # priors
        # search in the alpha space from the center of the range. Allow high variability
        log.alpha               <- log((startbio[1]+startbio[2])/2)
        sd.log.alpha            <- (log.alpha-log(startbio[1]))/5
        tau.log.alpha           <- pow(sd.log.alpha,-2)
        alpha                   ~  dlnorm(log.alpha,tau.log.alpha)

        # search in the k space from 20% of the range
        log.km              <- log(start.k[1]+0.2*(start.k[2]-start.k[1]))
        sd.log.k            <- (log.km-log(start.k[1]))/4
        tau.log.k           <- pow(sd.log.k,-2)
        k                   ~  dlnorm(log.km,tau.log.k)

        # define process (tau) and observation (sigma) variances as inversegamma priors
        itau2 ~ dgamma(2,0.01)
        tau2  <- 1/itau2
        tau   <- pow(tau2,0.5)

        isigma2 ~ dgamma(2,0.01)
        sigma2 <- 1/isigma2
        sigma <- pow(sigma2,0.5)

        log.rm              <- mean(log(start.r))
        sigma.log.r         <- abs(log.rm - log(start.r[1]))/2
        tau.log.r           <- pow(sigma.log.r,-2)
        r                   ~  dlnorm(log.rm,tau.log.r)
      }"

    # Run CPUE-based BSM
    ##########################################################

    # Run CPUE-based BSM
    }else{

      # Catchability stuff
      ########################################

      # Expert-specified catchability (q)
      if(is.na(q.start)==F & is.na(q.end)==F){
        mean.last.ct      <-mean(ct[yr >= q.start & yr <= q.end], na.rm=T) # get mean catch of indicated years
        mean.last.cpue    <-mean(bt[yr >= q.start & yr <= q.end], na.rm=T) # get mean of CPUE of indicated years
      # Default catchability (q)
      # get prior range for q from mean catch and mean CPUE in recent years
      }else{
        lyr               <- ifelse(mean(start.r)>=0.5,5,10)  # determine number of last years to use, 5 for normal and 10 for slow growing fish
        mean.last.ct      <-mean(ct[(nyr-lyr):nyr],na.rm=T) # get mean catch of last years
        mean.last.cpue    <-mean(bt[(nyr-lyr):nyr],na.rm=T) # get mean of CPUE of last years
      }
      gm.start.r      <- exp(mean(log(start.r))) # get geometric mean of prior r range
      if(mean(endbio) >= 0.5) {  # if biomass is high
        q.1           <- mean.last.cpue*0.25*gm.start.r/mean.last.ct
        q.2           <- mean.last.cpue*0.5*start.r[2]/mean.last.ct
      } else {
        q.1           <- mean.last.cpue*0.5*gm.start.r/mean.last.ct
        q.2           <- mean.last.cpue*start.r[2]/mean.last.ct
      }
      q.prior         <- c(q.1,q.2)
      init.q          <- mean(q.prior)

      # Setup JAGS model
      ########################################

      # Data to be passed on to JAGS
      jags.data        <- c('ct','bt','nyr', 'start.r', 'start.k', 'startbio', 'q.prior',
                            'init.q','init.r','init.k','pen.bk','pen.F','b.yrs','b.prior')
      # Parameters to be returned by JAGS
      jags.save.params <- c('r','k','q', 'P')

      # JAGS model
      Model = "model{
        # to reduce chance of non-convergence, Pmean[t] values are forced >= eps
        eps<-0.01
        penm[1] <- 0 # no penalty for first biomass
        Pmean[1] <- log(alpha)
        P[1] ~ dlnorm(Pmean[1],itau2)

        for (t in 2:nyr) {
        Pmean[t] <- ifelse(P[t-1] > 0.25,
        log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps)),  # Process equation
        log(max(P[t-1] + 4*P[t-1]*r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps))) # assuming reduced r at B/k < 0.25
        P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
        penm[t]  <- ifelse(P[t]<(eps+0.001),log(q*k*P[t])-log(q*k*(eps+0.001)),ifelse(P[t]>1,log(q*k*P[t])-log(q*k*(0.99)),0)) # penalty if Pmean is outside viable biomass
        }

        # ><> Biomass priors/penalties are enforced as follows
        for (i in 1:3) {
        penb[i]  <- ifelse(P[b.yrs[i]]<b.prior[1,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[1,i]),ifelse(P[b.yrs[i]]>b.prior[2,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[2,i]),0))
        b.prior[3,i] ~ dnorm(penb[i],100)
        }

        for (t in 1:nyr){
        Fpen[t]   <- ifelse(ct[t]>(0.9*k*P[t]),ct[t]-(0.9*k*P[t]),0) #><> Penalty term on F > 1, i.e. ct>B
        pen.F[t]  ~ dnorm(Fpen[t],1000)
        pen.bk[t] ~ dnorm(penm[t],10000)
        cpuem[t]  <- log(q*P[t]*k);
        bt[t]     ~ dlnorm(cpuem[t],isigma2);
        }

        # priors
        log.alpha               <- log((startbio[1]+startbio[2])/2) # needed for fit of first biomass
        sd.log.alpha            <- (log.alpha-log(startbio[1]))/4
        tau.log.alpha           <- pow(sd.log.alpha,-2)
        alpha                   ~  dlnorm(log.alpha,tau.log.alpha)

        # search in the k space starting from 20% of the range
        log.km              <- log(start.k[1]+0.2*(start.k[2]-start.k[1]))
        sd.log.k            <- (log.km-log(start.k[1]))/4
        tau.log.k           <- pow(sd.log.k,-2)
        k                   ~  dlnorm(log.km,tau.log.k)

        # set realistic prior for q
        log.qm              <- mean(log(q.prior))
        sd.log.q            <- (log.qm-log(q.prior[1]))/4
        tau.log.q           <- pow(sd.log.q,-2)
        q                   ~  dlnorm(log.qm,tau.log.q)

        # define process (tau) and observation (sigma) variances as inversegamma prios
        itau2 ~ dgamma(4,0.01)
        tau2  <- 1/itau2
        tau   <- pow(tau2,0.5)

        isigma2 ~ dgamma(2,0.01)
        sigma2 <- 1/isigma2
        sigma <- pow(sigma2,0.5)

        log.rm              <- mean(log(start.r))
        sigma.log.r         <- abs(log.rm - log(start.r[1]))/2
        tau.log.r           <- pow(sigma.log.r,-2)
        r                   ~  dlnorm(log.rm,tau.log.r)
      }" # end of JAGS model for CPUE

    } # end of else loop for Schaefer with CPUE

    # Run JAGS model
    ##########################################################

    # Write JAGS model to file
    cat(Model, file="R/r2jags.bug")

    # Initialize JAGS model?
    if(btype=="biomass") {
      j.inits     <- function(){list("r"=stats::rnorm(1,mean=init.r,sd=0.2*init.r),
                                     "k"=stats::rnorm(1,mean=init.k,sd=0.1*init.k),
                                     "itau2"=1000,
                                     "isigma2"=1000)}} else {
                                       j.inits <- function(){list("r"=stats::rnorm(1,mean=init.r,sd=0.2*init.r),
                                                                  "k"=stats::rnorm(1,mean=init.k,sd=0.1*init.k),
                                                                  "q"=stats::rnorm(1,mean=init.q,sd=0.2*init.q),
                                                                  "itau2"=1000,
                                                                  "isigma2"=1000)}}
    # Run JAGS model
    jags_outputs <- R2jags::jags.parallel(data=jags.data,
                                  working.directory=NULL, inits=j.inits,
                                  parameters.to.save=jags.save.params,
                                  model.file="R/r2jags.bug", n.chains = n.chains,
                                  n.burnin = 30000, n.thin = 10,
                                  n.iter = 60000)


    # Extract JAGS model results
    ##########################################################

    r_raw            <- as.numeric(coda::mcmc(jags_outputs$BUGSoutput$sims.list$r))
    k_raw            <- as.numeric(coda::mcmc(jags_outputs$BUGSoutput$sims.list$k))
    # Importance sampling: only accept r-k pairs where r is near the prior range
    r_out            <- r_raw[r_raw > 0.5*start.r[1] & r_raw < 1.5 * start.r[2]]
    k_out            <- k_raw[r_raw > 0.5*start.r[1] & r_raw < 1.5 * start.r[2]]

    mean.log.r.jags  <- mean(log(r_out))
    sd.log.r.jags    <- stats::sd(log(r_out))
    r.jags           <- exp(mean.log.r.jags)
    lcl.r.jags       <- exp(mean.log.r.jags - 1.96*sd.log.r.jags)
    ucl.r.jags       <- exp(mean.log.r.jags + 1.96*sd.log.r.jags)
    mean.log.k.jags  <- mean(log(k_out))
    sd.log.k.jags    <- stats::sd(log(k_out))
    k.jags           <- exp(mean.log.k.jags)
    lcl.k.jags       <- exp(mean.log.k.jags - 1.96*sd.log.k.jags)
    ucl.k.jags       <- exp(mean.log.k.jags + 1.96*sd.log.k.jags)
    MSY.posterior     <- r_out*k_out/4 # simpler
    mean.log.MSY.jags <- mean(log(MSY.posterior))
    sd.log.MSY.jags   <- stats::sd(log(MSY.posterior))
    MSY.jags          <- exp(mean.log.MSY.jags)
    lcl.MSY.jags      <- exp(mean.log.MSY.jags - 1.96*sd.log.MSY.jags)
    ucl.MSY.jags      <- exp(mean.log.MSY.jags + 1.96*sd.log.MSY.jags)

    # CPUE-based computations
    if(btype=="CPUE") {
      q_out           <- as.numeric(coda::mcmc(jags_outputs$BUGSoutput$sims.list$q))
      mean.log.q      <- mean(log(q_out))
      sd.log.q        <- stats::sd(log(q_out))
      mean.q          <- exp(mean.log.q)
      lcl.q           <- exp(mean.log.q-1.96*sd.log.q)
      ucl.q           <- exp(mean.log.q+1.96*sd.log.q)
      F.bt.cpue       <- mean.q*ct.raw/bt
      Fmsy.cpue       <- r.jags/2
    }

    # get F from observed biomass
    if(btype == "biomass"){
      F.bt       <- ct.raw/bt
      Fmsy.bt    <- r.jags/2
    }

    # get relative biomass P=B/k as predicted by BSM, including predictions for years with NA abundance
    all.P    <- jags_outputs$BUGSoutput$sims.list$P # matrix with P distribution by year
    quant.P  <- apply(all.P,2,stats::quantile,c(0.025,0.5,0.975),na.rm=T)

    # get k, r posterior ><>
    all.k  <- jags_outputs$BUGSoutput$sims.list$k # matrix with P distribution by year
    all.r  <- jags_outputs$BUGSoutput$sims.list$r # matrix with P distribution by year

    # get B/Bmys posterior
    all.b_bmsy=NULL
    for(t in 1:ncol(all.P)){
      all.b_bmsy  <- cbind(all.b_bmsy,all.P[,t]*2)}

    # get F/Fmys posterior ><>
    all.F_Fmsy=NULL
    for(t in 1:ncol(all.P)){
      all.F_Fmsy<- cbind(all.F_Fmsy,(ct.raw[t]/(all.P[,t]*all.k))/ifelse(all.P[,t]>0.25,all.r/2,all.r/2*4*all.P[,t]))}

  } # end of MCMC Schaefer loop

  # EXTRACT RESULTS
  ##############################################################################

  # Get management results
  MSY   <-MSY.jags; lcl.MSY<-lcl.MSY.jags; ucl.MSY<-ucl.MSY.jags
  Bmsy  <-k.jags/2; lcl.Bmsy<-lcl.k.jags/2; ucl.Bmsy<-ucl.k.jags/2
  Fmsy  <-r.jags/2; lcl.Fmsy<-lcl.r.jags/2; ucl.Fmsy<-ucl.r.jags/2
  B.Bmsy<-2*quant.P[2,];lcl.B.Bmsy<-2*quant.P[1,];ucl.B.Bmsy<-2*quant.P[3,]
  B          <-B.Bmsy*Bmsy;lcl.B<-lcl.B.Bmsy*Bmsy;ucl.B<-ucl.B.Bmsy*Bmsy
  Fm           <- ct.raw/B;lcl.F<-ct.raw/ucl.B;ucl.F<-ct.raw/lcl.B
  Fmsy.vec     <- ifelse(B.Bmsy>0.5,Fmsy,Fmsy*2*B.Bmsy)
  lcl.Fmsy.vec <- ifelse(B.Bmsy>0.5,lcl.Fmsy,lcl.Fmsy*2*B.Bmsy)
  ucl.Fmsy.vec <- ifelse(B.Bmsy>0.5,ucl.Fmsy,ucl.Fmsy*2*B.Bmsy)
  F.Fmsy       <- Fm/Fmsy.vec; lcl.F.Fmsy<-lcl.F/Fmsy.vec; ucl.F.Fmsy<-ucl.F/Fmsy.vec

  # Print results (if desired)
  if(verbose==T){
    # Print priors
    cat("---------------------------------------\n")
    cat("Catch data used from years", min(yr),"-", max(yr),", abundance =", btype, "\n")
    cat("Prior initial relative biomass =", startbio[1], "-", startbio[2],ifelse(is.na(stb.low)==T,"default","expert"), "\n")
    cat("Prior intermediate rel. biomass=", intbio[1], "-", intbio[2], "in year", int.yr,ifelse(is.na(intb.low)==T,"default","expert"), "\n")
    cat("Prior final relative biomass   =", endbio[1], "-", endbio[2],ifelse(is.na(endb.low)==T,"default","expert"), "\n")
    cat("Prior range for r =", format(start.r[1],digits=2), "-", format(start.r[2],digits=2),ifelse(is.na(r.low)==T,"default","expert,"),
        ", prior range for k =", start.k[1], "-", start.k[2],"\n")
    # Print BSM results
    cat("Results from Bayesian Schaefer model (BSM) using catch &",btype,"\n")
    cat("------------------------------------------------------------\n")
    if(btype == "CPUE") cat("q   =", mean.q,", lcl =", lcl.q, ", ucl =", ucl.q,"\n")
    cat("r   =", r.jags,", 95% CL =", lcl.r.jags, "-", ucl.r.jags,", k =", k.jags,", 95% CL =", lcl.k.jags, "-", ucl.k.jags,"\n")
    cat("MSY =", MSY.jags,", 95% CL =", lcl.MSY.jags, "-", ucl.MSY.jags,"\n")
    cat("Relative biomass in last year =", quant.P[2,][nyr], "k, 2.5th perc =",quant.P[1,][nyr],
        ", 97.5th perc =", quant.P[3,][nyr],"\n")
    cat("Exploitation F/(r/2) in last year =", (ct.raw[nyr]/(quant.P[2,][nyr]*k.jags))/(r.jags/2) ,"\n\n")
    # Print catchability prior (if CPUE-based BSM)
    if(btype=="CPUE") {
      cat("Prior range of q =",q.prior[1],"-",q.prior[2],"\n")
    }
  }

  # COLLECT RESULTS FOR OUTPUT
  ##############################################################################

  # Refence points dataframe
  ref_pts <- data.frame(rbind(c(est=r.jags, lo=lcl.r.jags, hi=ucl.r.jags),
                              c(k.jags, lcl.k.jags, ucl.k.jags),
                              c(MSY.jags, lcl.MSY.jags, ucl.MSY.jags),
                              c(Bmsy, lcl.Bmsy, ucl.Bmsy),
                              c(Fmsy, lcl.Fmsy, ucl.Fmsy)))
  ref_pts$param <- c("r", "k", "msy", "bmsy", "fmsy")
  ref_pts <- subset(ref_pts, select=c(param, est, lo, hi))

  # # Format F time series data
  # if(btype=="CPUE"){f <- F.bt.cpue}else{f <- F.bt}
  # if(btype=="CPUE"){ffmsy <- F.bt.cpue/Fmsy.cpue}else{ffmsy <- F.bt/Fmsy.bt}
  #
  # # Format saturation data
  # if(btype=="CPUE"){s <-bt/(mean.q*k.jags)}else{s <-bt/k.jags}
  # if(btype=="CPUE"){s_lo <- bt/(mean.q*ucl.k.jags)}else{s_lo <- bt/ucl.k.jags}
  # if(btype=="CPUE"){s_hi <- bt/(mean.q*lcl.k.jags)}else{s_hi <- bt/lcl.k.jags}

  # Reference points time series
  ref_ts <- data.frame(year=yr, catch=ct.raw, catch_ma=ct,
                       b=B, b_lo=lcl.B, b_hi=ucl.B,
                       bbmsy=B.Bmsy, bbmsy_lo=lcl.B.Bmsy, bbmsy_hi=ucl.B.Bmsy,
                       s=B.Bmsy/2, s_lo=lcl.B.Bmsy/2, s_hi=ucl.B.Bmsy/2,
                       f=Fm, f_lo=lcl.F, f_hi=ucl.F,
                       fmsy=Fmsy.vec, fmsy_lo=lcl.Fmsy.vec, fmsy_hi=ucl.Fmsy.vec,
                       ffmsy=F.Fmsy, ffmsy_lo=lcl.F.Fmsy, ffmsy_hi=ucl.F.Fmsy,
                       er=Fm/Fmsy)

  #stop parallel processing clusters
  parallel::stopCluster(cl)
  doParallel::stopImplicitCluster()

  # Assemble output
  output <- list(ref_pts=ref_pts, ref_ts=ref_ts, priors=priors,
                 r_out=r_out, k_out=k_out)
  return(output)

}
