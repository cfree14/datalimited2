
#' Convert B/BMSY to categorical status
#'
#' Converts B/BMSY estimates to categorical statuses (e.g., underexploited,
#' fully exploited, or overexploited).
#'
#' @param bbmsy A vector, matrix, or dataframe of B/BMSY estimates.
#' @param breaks A vector of B/BMSY thresholds used to delineate status categories,
#' listed in order of ascending value (e.g., 0.2, 0.5, 1.5). Default is 0.5 and 1.5 where
#' B/BMSY < 0.5 = overexploited, 0.5 < B/BMSY < 1.5 = fully exploited, and B/BMSY > 1.5 = underexploited.
#' @param catgs A vector of names for the status categories, listed
#' in order of descending depletion (e.g, default is over, fully, under).
#' @return A vector, matrix, or dataframe of status categories.
#' @examples
#' bbmsy <- select(preds, bbmsy,  mprm, comsir, sscom, cmsy13, cmsy17, zbrt, ocom, super1)
#' status <- bbmsy2catg(bbmsy)
#' @export
bbmsy2catg <- function(bbmsy, breaks=c(0.5, 1.5), catgs=c("over", "fully", "under")){
  obj_type <- class(bbmsy)
  if(obj_type=="numeric"){
    status <- as.character(cut(bbmsy, breaks=c(0,breaks,9999), labels=catgs))
  }else{
    status <- apply(bbmsy, 2, function(x) as.character(cut(x, breaks=c(0,breaks,9999), labels=catgs)))
  }
  return(status)
}
