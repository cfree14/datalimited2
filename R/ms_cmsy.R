
#' Multispecies cMSY (beta version)
#'
#' Estimates stock status (B/BMSY) time series for multispecies fisheries from
#' time series of catch and estimates of resilience uing the multispecies cMSY
#' (MS-cMSY) method of Free et al. (in prep). Note: this model is still under development.
#'
#' @param catch A dataframe with a 'year' column and additional columns containing
#' the catch for each species in the multispecies fishery.
#' @param stocks A character vector with the names of the species in the multispecies fishery.
#' Must be reported in the same order as the catch columns.
#' @param res A character vector with the resilience estimates for species in the
#' multispecies fishery. Must be reported in the same order as the catch columns.
#' Can be retrieved using the datalimited2::resilience() function.
#' @param id_fixed A boolean (TRUE or FALSE) indicating whether to (TRUE) fix initial saturation at
#' one (no depletion) or (FALSE) whether the initial saturation should be estimated.
#' @param npairs The number of r-K pairs that should be evaulated. The default is 5,000.
#' @return A list with the following elements:
#' \item{stocks}{Names of the analyzed stocks/species}
#' \item{yrs}{Years represented in the time series}
#' \item{r_priors}{R priors for each species}
#' \item{k_priors}{K priors for each species}
#' \item{s1_priors}{Initial saturation priors for each species}
#' \item{s2_priors}{Final saturation priors for each species}
#' \item{id_try}{Initial saturation values evaluated}
#' \item{r_try}{R values evaluated}
#' \item{k_try}{K values evaluated}
#' \item{id_rk_v}{Viable r-K-initial saturation combos}
#' \item{b_v}{Viable biomass trajectories}
#' \item{bbmsy_v}{Viable B/BMSY trajectories}
#' \item{bbmsy_v_median}{Median of viable B/BMSY trajectories (cMSY estimate of status)}
#' \item{er_v}{Viable exploitation rate trajectories}
#' \item{er_vv}{Most probable exploitation rate trajectories based on effort constraint}
#' \item{top_corr}{Most highly correlated r-k-initial saturation combos from each row of correlation matrix}
#' \item{bbmsy_vv}{Most probable B/BMSY trajectories based on effort constraint}
#' \item{bbmsy_vv_median}{Median of most probable B/BMSY trajectories (MS-cMSY estimate of status)}
#' \item{method}{Name of the method}
#' @references Free CM, Rudd MB, Kleisner KM, Thorson JT, Longo C, Minto C,
#' Jensen OP (in prep) Multispecies catch-only models for assessing data-limited fisheries.
#' @export
ms_cmsy <- function(catch, stocks, res, id_fixed, npairs=5000){

  # Extract yrs/catch
  yrs <- catch$year
  C_mat <- as.matrix(select(catch, -year))

  # Time series info
  nyrs <- length(yrs)
  nstocks <- length(stocks)

  # Calculate r and k priors
  # R prior based on resilience; K prior based on max catch and r prior
  r_priors <- r_priors(res)
  k_priors <- k_priors(C_mat)
  r_priors_ln <- log(r_priors)
  k_priors_ln <- log(k_priors)

  # Calculate saturation priors
  s1_priors <- sat1_priors(yrs, nstocks)
  s2_priors <- sat2_priors(C_mat)

  # Randomly sample r-k pairs in log-space
  npairs <- npairs
  ri <- sapply(1:nstocks, function(x) exp(runif(npairs, r_priors_ln[x,1], r_priors_ln[x,2])))
  ki <- sapply(1:nstocks, function(x) exp(runif(npairs, k_priors_ln[x,1], k_priors_ln[x,2])))

  # Lists to hold viable r/k pairs and associated biomass/exploitation trajectories
  id_rk_try <- list()
  id_rk_v <- list()
  b_mats_v <- list()
  bbmsy_mats_v <- list()
  er_mats_v <- list()

  # Loop through stocks to identify viable r/k pairs and trajectoris
  for(i in 1:nstocks){

    # Get info
    stock <- stocks[i]
    c_vec <- C_mat[,i]
    # print(stock)

    # Initial depletions to evaluate
    if(id_fixed==T){
      ids <- 1
    }else{
      s1_prior <- s1_priors[i,]
      ids <- seq(s1_prior[1], s1_prior[2], 0.1)
    }

    # Get r/k pairs to evaluate
    rs <- ri[,i]
    ks <- ki[,i]
    rk_pairs <- cbind(r=rs, k=ks, viable=rep(NA, npairs))

    # Build ID/r/k combos to evaluate
    id_rk_combos <- as.data.frame(do.call("rbind",
                                          lapply(ids, function(x) cbind(id=rep(x, nrow(rk_pairs)), rk_pairs))))

    # Loop through r/k pairs to see if viable
    p <- 0.2
    sigmaP <- 0.1
    b_mat <- matrix(data=NA, nrow=nyrs, ncol=nrow(id_rk_combos))
    for(j in 1:nrow(id_rk_combos)){
      id <- id_rk_combos$id[j]
      r <- id_rk_combos$r[j]
      k <- id_rk_combos$k[j]
      b_mat[1,j] <- k * id
      for(yr in 2:nyrs){
        # b_mat[yr,j] <- b_mat[yr-1,j] +  r*b_mat[yr-1,j]/p*(1-(b_mat[yr-1,j]/k)^p) - c_vec[yr-1]
        b_mat[yr,j] <- b_mat[yr-1,j] +  r*b_mat[yr-1,j]/p*(1-(b_mat[yr-1,j]/k)^p)*exp(rnorm(1,0,sigmaP)) - c_vec[yr-1]
      }
    }

    # Create saturation matrix
    s_mat <- t(t(b_mat)/ks)

    # Reduce to viable r/k pairs and trajectories
    check0 <- apply(b_mat, 2, function(x) sum(x<0, na.rm=T)==0) # check B doesn't go below 0
    s2_lo <- s2_priors[i,1]
    s2_hi <- s2_priors[i,2]
    s2_vec <- s_mat[nrow(s_mat),]
    checkS <- s2_vec > s2_lo & s2_vec < s2_hi & !is.na(s2_vec) # check final yr saturation inside prior
    viable <- check0 & checkS # merge positive biomass and fina saturation checks
    id_rk_combos[,"viable"] <- viable
    nviable <- sum(viable)
    b_mat_viable <- b_mat[,viable]
    id_rk_viable <- id_rk_combos[viable,]

    # Derive B/BMSY
    bmsy <- id_rk_viable[,"k"] * (1 / (p+1))^(1/p)
    bbmsy_mat_viable <- t(t(b_mat_viable) / bmsy)

    # # Plot viable r/k pairs
    # par(mfrow=c(1,1))
    # plot(k ~ r, rk_pairs, log="xy", bty="n", las=1, pch=15,
    #      xlim=r_priors[i,], ylim=k_priors[i,], xlab="r", ylab="k", col="gray80")
    # points(rk_viable[,1], rk_viable[,2], pch=15, col="black")
    #
    # # Plot viable biomass trajectories
    # plot(b_mat_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, max(b_mat_viable, na.rm=T)), xlab="Year", ylab="Biomass")
    # for(k in 1:ncol(b_mat_viable)){lines(x=yrs, y=b_mat_viable[,k], col="grey80")}
    #
    # # Plot viable B/BMSY trajectories
    # plot(bbmsy_mat_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, max(bbmsy_mat_viable, na.rm=T)), xlab="Year", ylab="B/BMSY")
    # for(k in 1:ncol(bbmsy_mat_viable)){lines(x=yrs, y=bbmsy_mat_viable[,k], col="grey80")}

    # Calculate exploitation rate
    # I name the rows and columns so that I can validate covariance matrix below
    er_mat_viable <- c_vec / b_mat_viable
    rownames(er_mat_viable) <- yrs
    colnames(er_mat_viable) <- paste0(LETTERS[i], 1:ncol(er_mat_viable))

    # Plot viable exploitation trajectories
    # plot(er_mat_viable[,1] ~ yrs, type="n", bty="n", las=1,
    #      ylim=c(0, 1), xlab="Year", ylab="Exploitation rate")
    # for(k in 1:ncol(er_mat_viable)){lines(x=yrs, y=er_mat_viable[,k], col="grey80")}

    # Record results
    id_rk_v[[i]] <- id_rk_viable
    b_mats_v[[i]] <- b_mat_viable
    bbmsy_mats_v[[i]] <- bbmsy_mat_viable
    er_mats_v[[i]] <- er_mat_viable

  }

  # Measure correlation
  ###################################################

  # Build index (start and end points) for each stock
  # This is used to index the correlation matrix below
  id_rk_v_n_all <- sapply(id_rk_v, nrow) # number of viable r/k pairs for each stock
  finishes <- cumsum(id_rk_v_n_all)
  starts <- c(0,finishes[1:(length(finishes)-1)])+1
  indices <- cbind(starts, finishes)

  # Merge viable effort time series for all stocks
  er_v_all <- t(do.call("cbind", er_mats_v)) # transposed so that G=F*Ft works

  # Calculate covariance matrix then correlation matrix
  cov_mat <- er_v_all %*% t(er_v_all)
  corr_mat <- cov2cor(cov_mat)

  # # Overwrite self-comparisons
  # dim(corr_mat)
  # corr_mat1 <- corr_mat
  # for(j in 1:nrow(indices)){
  #   pt1 <- indices[j,1]
  #   pt2 <- indices[j,2]
  #   corr_mat1[pt1:pt2, pt1:pt2] <- NA
  # }
  # # Plot correlation matrix
  # corr_mat2 <- t(apply(corr_mat1, 2, rev))
  # image(x=1:ncol(corr_mat2), y=1:ncol(corr_mat2), z=corr_mat2,
  #       xaxt="n", yaxt="n", xlab="", ylab="")
  # abline(v=finishes, lwd=1.5)
  # abline(h=ncol(corr_mat2)-finishes, lwd=1.5)

  # Setup containter to hold highest correlation coefficients per row
  # Columns: row, index1, index2, index3, corr12, corr13, corr23, corr_sum
  spp <- 1:nstocks # each species gets a number
  spp_vec <- rep(1:nstocks, id_rk_v_n_all) # indexes species identity (1, 2, or 3 etc)
  spp_ind <- unlist(sapply(id_rk_v_n_all, function(x) 1:x)) # indexes index of r/k pair
  hi_corr_mat <- matrix(nrow=nrow(corr_mat), ncol=1+nstocks+nstocks, data=NA,
                        dimnames=list(NULL, c("row",
                                              paste0("index", spp),
                                              paste0("corr", apply(combn(1:3, 2),2, function(x) paste(x, collapse=""))))))
  hi_corr_mat[,"row"] <- 1:nrow(hi_corr_mat)

  # Loop through rows of correlation matix and identify highest correlation in each row-block
  for(j in 1:nrow(corr_mat)){

    # Which species am I working on?
    spp_curr <- spp_vec[j]
    spp_others <- spp[!(spp %in% spp_curr)]

    # Record index of current species
    hi_corr_mat[j, paste0("index", spp_curr)] <- spp_ind[j]

    # Loop through other species blocks
    for(k in 1:length(spp_others)){
      spp_other <- spp_others[k]
      spp_other_ind <- indices[spp_other, 1]:indices[spp_other,2]
      corrs <- corr_mat[j,spp_other_ind]
      hi_corr_mat[j,paste0("index", spp_other)] <- which.max(corrs)
      corr_col <- paste0("corr", paste(sort(c(spp_curr, spp_other)), collapse=""))
      hi_corr_mat[j,corr_col] <- max(corrs)
    }

  }

  # Average correlations
  corr_cols <- colnames(hi_corr_mat)[grepl("corr", colnames(hi_corr_mat))]
  corr_avgs <- apply(hi_corr_mat[,corr_cols ], 1, mean, na.rm=T)
  hi_corr_mat <- cbind(hi_corr_mat, corr_avg=corr_avgs)

  # Mark viable r/k pairs that show high correlation
  # n_each_index <- melt(hi_corr_mat[,paste0("index", 1:nstocks)],
  #                      variable.name="stock", value.name="index") %>%
  #   select(-Var1) %>%
  #   rename(stock=Var2) %>%
  #   mutate(stock=gsub("index", "", stock)) %>%
  #   group_by(stock, index) %>%
  #   summarize(ncorr=n()) %>%
  #   ungroup()

  # Add correlation counts to viable r/k pair table
  for(i in 1:length(id_rk_v)){
    id_rk_v_do <- as.data.frame(id_rk_v[[i]])
    id_rk_v_do1 <- id_rk_v_do %>%
      mutate(index=1:nrow(id_rk_v_do)) %>%
      select(index, id, r, k)
    id_rk_v[[i]] <- id_rk_v_do1
  }

  # Identify top 10% most highly correlated effort time series
  top_p <- 0.05
  top_n <- ceiling(nrow(hi_corr_mat) * top_p)
  top_corr <- as.data.frame(hi_corr_mat) %>%
    arrange(desc(corr_avg)) %>%
    slice(1:top_n)

  # Identify combos producing > 0.4 mean correlation
  # corr_thresh <- 0.8
  # top_corr <- as.data.frame(hi_corr_mat) %>%
  #   arrange(desc(corr_avg)) %>%
  #   filter(corr_avg>=corr_thresh)
  # if(nrow(top_corr)==0){print("No highly correlated effort time series found.")}

  # Get biomass trajectories of top 10%
  # (also sneak in calculation of cMSY prediction)
  bbmsy_v_meds <- NULL
  er_mats_vv <- NULL
  bbmsy_mats_vv <- NULL
  bbmsy_vv_meds <- NULL
  for(i in 1:length(b_mats_v)){
    vv_index <- unlist(top_corr[,paste0("index", i)])
    er_mat_v <- er_mats_v[[i]]
    er_mat_vv <- er_mat_v[,vv_index]
    er_mats_vv[[i]] <- er_mat_vv
    bbmsy_mat_v <- bbmsy_mats_v[[i]]
    bbmsy_mat_vv <- bbmsy_mat_v[,vv_index]
    bbmsy_mats_vv[[i]] <- bbmsy_mat_vv
    bbmsy_v_med <- apply(bbmsy_mat_v, 1, median)
    bbmsy_vv_med <- apply(bbmsy_mat_vv, 1, median)
    bbmsy_v_meds[[i]] <- bbmsy_v_med
    bbmsy_vv_meds[[i]] <- bbmsy_vv_med
  }

  # Things to return
  out <- list(stocks=stocks,
              yrs=yrs,
              r_priors=r_priors,
              k_priors=k_priors,
              s1_priors=s1_priors,
              s2_priors=s2_priors,
              id_try=ids,
              r_try=ri,
              k_try=ki,
              id_rk_v=id_rk_v,
              b_v=b_mats_v,
              bbmsy_v=bbmsy_mats_v,
              bbmsy_v_median=bbmsy_v_meds,
              er_v=er_mats_v,
              er_vv=er_mats_vv,
              top_corr=top_corr,
              bbmsy_vv=bbmsy_mats_vv,
              bbmsy_vv_median=bbmsy_vv_meds,
              method="MS-cMSY")
  return(out)

}


################################################################################
# MS-cMSY priors (present study)
################################################################################

# R priors
r_priors <- function(res){
  for(i in 1:length(res)){
    if(res[i]=="High"){r_prior <- c(0.5, 1.25)}
    if(res[i]=="Medium"){r_prior <- c(0.01, 1.0)}
    if(res[i]=="Low"){r_prior <- c(0.01, 0.5)}
    if(res[i]=="Very low"){r_prior <- c(0.01, 0.15)}
    if(i==1){r_priors <- r_prior}else{r_priors <- rbind(r_priors, r_prior)}
  }
  colnames(r_priors) <- c("r_lo", "r_hi")
  rownames(r_priors) <- NULL
  return(r_priors)
}

# K priors
k_priors <- function(C_mat){
  cmax <- apply(C_mat, 2, max)
  k_lo <- cmax * 2
  k_hi <- cmax * 50
  k_priors <- cbind(k_lo=k_lo, k_hi=k_hi)
  return(k_priors)
}

# Initial saturation priors
# Test: plot(x=1900:2020, y=sapply(1900:2020, function(x) sat1_priors(yrs=min(x):2020, nstocks=1))[1,], ylim=c(0,2), ylab="Saturation", xlab="")
sat1_priors <- function(yrs, nstocks){
  yr1 <- min(yrs)
  s1_hi <- 1
  if(yr1<=1945){s1_lo <- 0.8}
  if(yr1>1945 & yr1<1980){s1_lo <- 0.8+(0.1-0.8)/(1980-1945)*(yr1-1945)}
  if(yr1>=1980){s1_lo <- 0.1}
  s1_lo1 <- rep(s1_lo, nstocks)
  s1_hi1 <- rep(s1_hi, nstocks)
  s1_priors <- cbind(s1_lo=s1_lo1, s1_hi=s1_hi1)
  return(s1_priors)
}

# Final saturation priors
# plot(x=seq(0,1,0.1), y=0+0.4*seq(0,1,0.1), ylim=c(0,2), type="l"); lines(x=seq(0,1,0.1), y=0.8+0.2*seq(0,1,0.1))
sat2_priors <- function(C_mat){
  c_ends <- C_mat[nrow(C_mat),]
  c_maxs <- apply(C_mat, 2, max)
  c_ratios <- c_ends / c_maxs
  s2_lo <- 0 + 0.4*c_ratios
  s2_hi <- 0.5 + 0.4*c_ratios
  s2_priors <- cbind(s2_lo=s2_lo, s2_hi=s2_hi)
  return(s2_priors)
}
