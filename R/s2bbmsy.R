#' Convert B/BMSY to saturation
#'
#' Converts B/BMSY to saturation (S = B/K).
#'
#' S = 1 - depletion = B / K = B/BMSY / 2
#' B/BMSY = S * 2
#'
#' @param s Saturation (B/K) value(s)
#' @return B/BMSY values(s)
#' @export
s2bbmsy <- function(s){
  bbmsy <- s * 2
  return(bbmsy)
}

#' Convert saturation to B/BMSY
#'
#' Converts saturation (S = B/K) to B/BMSY
#'
#' S = 1 - depletion = B / K = B/BMSY / 2
#' B/BMSY = S * 2
#'
#' @param bbmsy B/BMSY values(s)
#' @return Saturation (B/K) value(s)
#' @export
bbmsy2s <- function(bbmsy){
  s <- bbmsy / 2
  return(s)
}
