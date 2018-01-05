
# For testing
# year <- YELLSNEMATL$year
# catch <- YELLSNEMATL$tc
# fit_zbrt(year=year, catch=catch)
# s <- zbrt(year=year, catch=catch)
# sum(s$s, na.rm=T)

# Read data
system.file("data", "YELLSNEMATL.Rdata", package = "datalimited2")
system.file("data", "BRTmodelP8.RData", package = "datalimited2")
system.file("data", "BRTmodelP38.RData", package = "datalimited2")

#' Compute predictors for the Zhou-BRT catch-only stock assessment model
catchParam <- function(catchData) {
  sid = unique(as.character(catchData$stock))
  n.stock = length(sid)
  para = matrix(NA, n.stock, 57)
  n.ab = 7  # number of years at the beginning and ending time series
  a=b=a0=b0=matrix(0, n.stock, n.ab)
  for(i in 1:n.stock) {
    stk = sid[i]
    dat = subset(catchData, stock==stk)
    yr = dat$yr-min(dat$yr) +1
    midYr = mean(yr)
    C = dat$catch/max(dat$catch)   # regressions are based on scaled catch
    C05max = sum(C[-length(C)]>0.5)
    Cmean = mean(C)
    nyr = length(yr)
    nyrToCmax = min(yr[C==max(C)]-yr[1]+1)
    nyrToCmaxR = nyrToCmax/nyr
    nyrAfterCmax = nyr-nyrToCmax
    yr.cent = yr-midYr           # regressions are based on centred year
    line0 = lm(C~yr.cent)
    yr1.cent = yr[1:nyrToCmax]-mean(yr[1:nyrToCmax])
    line1 = lm(C[1:nyrToCmax]~yr1.cent)
    yr2.cent = yr[nyrToCmax:nyr]-mean(yr[nyrToCmax:nyr])
    line2 = lm(C[nyrToCmax:nyr]~yr2.cent)
    aa0 = summary(line0)$coeff[1]; bb0 = summary(line0)$coeff[2]    #all yr
    aa1 = summary(line1)$coeff[1]; bb1 = summary(line1)$coeff[2]    #before Cmax
    aa2 = summary(line2)$coeff[1]; bb2 = summary(line2)$coeff[2]    #after Cmax
    for (j in 1:n.ab) { # periodical regressions
      yrLast.cent = yr[(nyr-j):nyr]-mean(yr[(nyr-j):nyr])
      l.last = lm(C[(nyr-j):nyr]~yrLast.cent)                   # last several years
      a[i,j] = summary(l.last)$coeff[1];
      b[i,j]=summary(l.last)$coeff[2]
      yrBegin.cent = yr[1:(j+1)]-mean(yr[1:(j+1)])
      l.begin = lm(C[1:(j+1)] ~ yrBegin.cent)                  # beginning several years
      a0[i,j] = summary(l.begin)$coeff[1]
      b0[i,j] = summary(l.begin)$coeff[2] }
    # segmented regression: use yr and breakpoint is between 0-1
    f = tryCatch(segmented(line0, seg.Z=~yr, psi=list(yr=median(yr))) , error=function(err) {})
    if(is.null(f)) {
      a.spline = NA; b1.spline = NA; b2.spline = NA; breakPoint = NA
    } else {
      a.spline = summary(f)$coef[1]
      b1.spline = summary(f)$coef[2]
      slp= slope(f)
      b2.spline = slp$yr[2]
      breakPoint = (round(f$psi.history[[5]],0)-yr[1] +1)/nyr
    }
    para[i,] = c(aa0, aa1, aa2, bb0, bb1, bb2, a[i,], b[i,], a0[i,], b0[i,], a.spline, b1.spline,b2.spline, breakPoint, C[1:n.ab], C[(nyr-n.ab): nyr], Cmean, nyrToCmaxR, nyr, C05max)
  }   # end params
  colnames(para) = c('a.AllYr', 'a.BfMax', 'a.AfMax', 'b.AllYr', 'b.BfMax', 'b.AfMax', paste('a.LsY',1:n.ab, sep=""), paste('b.LsY',1:n.ab, sep=""), paste("a.BgY", 1:n.ab, sep = ""), paste('b.BgY', 1:n.ab, sep=""), 'a.seg', 'b1.seg', 'b2.seg','breakPoint',  paste('c.BgY',1:n.ab, sep=""), paste('c.LsY', n.ab:0, sep=""), 'Cmean','nyrToCmaxR','nyr', 'C05max')
  para = data.frame(stock=sid, para)
  return(para)
}

# Fit the Zhou-BRT catch-only stock assessment model
fit_zbrt <- function(year, catch){

  # Calculate predictors
  catchData <- data.frame(stock="Dummy", yr=year, catch=catch)
  sessPar <- catchParam(catchData)

  # Center the first 37 predictors
  sessParCent <- scale(sessPar[,2:38], center=BRTmodelP8$parMean, scale=F)

  # Assemble the prediction data
  stockName <- unique(as.character(catchData$stock))
  nstk <- length(stockName)
  predDat <- data.frame(stock=stockName, sessParCent, sessPar[,39:57])

  # Estimate saturation (S) using 8-predictor and 37-predictor models
  s8 <- gbm::predict.gbm(BRTmodelP8$model, predDat, n.trees=BRTmodelP8$model$gbm.call$best.trees, type="response")
  s38 <- gbm::predict.gbm(BRTmodelP38$model, predDat, n.trees=BRTmodelP38$model$gbm.call$best.trees, type="response")

  # Bias correct 8-predictor model
  slr8.a <- BRTmodelP8$slr[[1]]; slr8.b = BRTmodelP8$slr[[2]];
  sBC8 <- (s8-slr8.a)/slr8.b
  sBC8[sBC8<=0] <- 0.01

  # Bias correct 38-predictor model
  slr38.a <- BRTmodelP38$slr[[1]]; slr38.b = BRTmodelP38$slr[[2]];
  sBC38 <- (s38-slr38.a)/slr38.b
  sBC38[sBC38<=0] <- 0.01

  # Average saturation predictions
  s <- data.frame(s8=sBC8, s38=sBC38, s=mean(sBC8, sBC38))
  return(s)
}

#' Zhou-BRT catch-only assessment model
#'
#' Estimates saturation (S; S = B/BMSY / 2) from a time series of catch using the
#' boosted regression tree model from Zhou et al. (2017).
#'
#' @param year A time series of years
#' @param catch A time series of catch
#' @return A time series of saturation and B/BMSY estimates. S8 and S38
#' correspond to the saturation estimates from the 8- and 38-predictor models.
#' S, the best estimate of saturation, is the mean of these two predictions.
#' B/BMSY is this estimate doubled (B/BMSY = S * 2).
#' @references Zhou S, Punt AE, Yimin Y, Ellis N, Dichmont CM, Haddon M, Smith DC, Smith ADM
#' (2017) Estimating stock depletion level from patterns of catch history. \emph{Fish and Fisheries}.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/faf.12201/abstract}
#' @examples
#' zbrt(year=YELLSNEMATL$year, catch=YELLSNEMATL$tc)
#' @export
zbrt <- function(year, catch){
  # Create dataframe with S estimates
  d <- data.frame(year=year, catch=catch, s8=NA, s38=NA, s=NA, bbmsy=NA)
  # Loop through years and estimate S
  # zBRT requires 7 years of data
  for(i in nrow(d):8){
    yr_use <- d$year[1:i]
    catch_use <- d$catch[1:i]
    s_out <- fit_zbrt(year=yr_use, catch=catch_use)
    d[i, c("s8", "s38", "s")] <- s_out
  }
  # Calculate B/BMSY from S
  d$bbmsy <- d$s * 2
  return(d)
}



