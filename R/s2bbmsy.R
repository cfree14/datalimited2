#' Convert saturation to B/BMSY
#'
#' Converts saturation  to B/BMSY. Note: S = 1 - depletion = B / K = 0.5 * B/BMSY.
#'
#' @param s A vector of saturation values
#' @return A vector of B/BMSY values
#' @examples
#' # Convert saturation of 0.75 to B/BMSY
#' s2bbmsy(0.75)
#' @export
s2bbmsy <- function(s){
  bbmsy <- s * 2
  return(bbmsy)
}

#' Convert B/BMSY to saturation
#'
#' Converts B/BMSY to saturation (S). Note: S = 1 - depletion = B / K = 0.5 * B/BMSY.
#'
#' @param bbmsy A vector of B/BMSY values
#' @return A vector of saturation values
#' @examples
#' # Convert B/BMSY of 0.5 to saturation
#' bbmsy2s(0.5)
#' @export
bbmsy2s <- function(bbmsy){
  s <- bbmsy / 2
  return(s)
}
