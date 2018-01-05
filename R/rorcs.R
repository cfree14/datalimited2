
# For testing
# load("data/refined_orcs_approach_bct_model.Rdata")
# scores <- c(2, # TOA 1 - Status of assessed stocks in fishery
#             NA, # TOA 3 - Behavior affecting capture (2 or 3 only)
#             1, # TOA 5 - Discard rate
#             2, # TOA 6 - Targeting intensity
#             NA, # TOA 7 - M compared to dominant species
#             2, # TOA 8 - Occurence in catch
#             1.51, # TOA 9 - Value (US$/lb) - continuous value
#             3, # TOA 10 - Recent trends in catch
#             1, # TOA 11 - Habitat loss
#             2, # TOA 12 - Recent trend in effort
#             2, # TOA 13 - Recent trend in abundance index
#             1) # TOA 14 - Proportion of population protected
# rorcs(scores)

# Read data
system.file("data", "refined_orcs_approach_bct_model.Rdata", package = "datalimited2")

#' Predict stock status using the refined ORCS approach
#'
#' This function predicts stock status (i.e., under, fully, or overexploited) using
#' the refined ORCS approach (Free et al. 2017).
#'
#' @param scores A numeric vector of length twelve containing scores for the
#' following "Table of Attributes" questions:
#' \itemize{
#'   \item{TOA 1 - Status of assessed stocks in fishery}
#'   \item{TOA 3 - Behavior affecting capture (2 or 3 only)}
#'   \item{TOA 5 - Discard rate}
#'   \item{TOA 6 - Targeting intensity}
#'   \item{TOA 7 - M compared to dominant species}
#'   \item{TOA 8 - Occurence in catch}
#'   \item{TOA 9 - Value (US$/lb) - continuous value}
#'   \item{TOA 10 - Recent trends in catch}
#'   \item{TOA 11 - Habitat loss}
#'   \item{TOA 12 - Recent trend in effort}
#'   \item{TOA 13 - Recent trend in abundance index}
#'   \item{TOA 14 - Proportion of population protected}
#' }
#' @return Stock status (i.e., under, fully, or overexploited)
#' @references Free CM, Jensen OP, Wiedenmann J, Deroba JJ (2017) The
#' refined ORCS approach: a catch-based method for estimating stock status
#' and catch limits for data-poor fish stocks. Fisheries Research 193: 60-70.
#' \url{https://doi.org/10.1016/j.fishres.2017.03.017}
#' @export
rorcs <- function(scores){

  # Check for errors in TOA scores
  error.flag <- F
  if(!scores[2]%in%c(NA,2,3)){
    cat("Error: The answer to TOA #3 must be 2, 3, or NA. It is currently ", scores[2], ".", sep="", fill=T)
    error.flag <- T
  }
  pos.to.check <- c(1,3:6,8:12)
  toa.to.check <- c(1, 5:8, 10:14)
  for(i in 1:length(pos.to.check)){
    index <- pos.to.check[i]
    if(!scores[index]%in%c(NA,1,2,3)){
      cat("Error: The answer to TOA #", toa.to.check[i],  " must be 1, 2, 3, or NA. It is currently ", scores[index], ".", sep="", fill=T)
      error.flag <- T
    }
  }

  # If there is an error, create empty output
  if(error.flag==T){
    # Create output
    out <- data.frame(status="error", p.under=NA, p.fully=NA, p.over=NA)
  # If there are no errors, predict stock status
  }else{
    # Build TOA dataframe
    scores.df <- as.data.frame(matrix(scores, nrow=1))
    colnames(scores.df) <- c("toa1", "toa3", "toa5", "toa6",
                             "toa7", "toa8", "toa9b", "toa10a",
                             "toa11", "toa12", "toa13", "toa14")

    # Predict status probabilities
    probs <- predict(rorcs_model, newdata=scores.df, type="prob", na.action=NULL)
    colnames(probs) <- c("p.under", "p.fully", "p.over")

    # Identify most-likely status
    status <- gsub("p.", "", colnames(probs)[which.max(probs)])

    # Create output
    out <- cbind(status, probs)
  }

  # Return status info
  return(out)

}

