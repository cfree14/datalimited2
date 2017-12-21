
system.file("extdata", "BRTmodelP8.RData", package = "datalimited2")
system.file("extdata", "BRTmodelP38.RData", package = "datalimited2")

#' Compute predictors for Zhou-BRT catch-only stock assessment model
#'
#' This function calculates the catch statistics used as predictors in the Zhou-BRT
#' catch-only stock assessment model.
#'
#' @param yr A time series of years
#' @param catch A time series of catch
#' @return A time series of B/BMSY estimates
#' @export
catchParam = function(catchData) {
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
    f = tryCatch(segmented::segmented(line0, seg.Z=~yr, psi=list(yr=median(yr))) , error=function(err) {})
    if(is.null(f)) {
      a.spline = NA; b1.spline = NA; b2.spline = NA; breakPoint = NA
    } else {
      a.spline = summary(f)$coef[1]
      b1.spline = summary(f)$coef[2]
      slp= segmented::slope(f)
      b2.spline = slp$yr[2]
      breakPoint = (round(f$psi.history[[5]],0)-yr[1] +1)/nyr
    }
    para[i,] = c(aa0, aa1, aa2, bb0, bb1, bb2, a[i,], b[i,], a0[i,], b0[i,], a.spline, b1.spline,b2.spline, breakPoint, C[1:n.ab], C[(nyr-n.ab): nyr], Cmean, nyrToCmaxR, nyr, C05max)
  }   # end params
  colnames(para) = c('a.AllYr', 'a.BfMax', 'a.AfMax', 'b.AllYr', 'b.BfMax', 'b.AfMax', paste('a.LsY',1:n.ab, sep=""), paste('b.LsY',1:n.ab, sep=""), paste("a.BgY", 1:n.ab, sep = ""), paste('b.BgY', 1:n.ab, sep=""), 'a.seg', 'b1.seg', 'b2.seg','breakPoint',  paste('c.BgY',1:n.ab, sep=""), paste('c.LsY', n.ab:0, sep=""), 'Cmean','nyrToCmaxR','nyr', 'C05max')
  para = data.frame(stock=sid, para)
  return(para)
}

#' Fit the Zhou-BRT catch-only stock assessment model
#'
#' This function fits the Zhou-BRT catch-only stock assessment model. It
#' only requires a time series of catch.
#'
#' @param yr A time series of years
#' @param catch A time series of catch
#' @return A time series of B/BMSY estimates
#' @export
fit_zbrt <- function(yr, catch){

  # Calculate predictors
  catchData <- data.frame(stock="Dummy", yr=yr, catch=catch)
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
  s <- mean(sBC38, sBC8)
  return(s)

}
