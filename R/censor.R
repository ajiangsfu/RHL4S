#' A function rather aimed at developers
#' @author Alex Xu
#' @noRd
#' 
## a helper function which censor the data
censor <- function(data, percentile=0.99, cap=NA){
  data <- as.matrix(data)
  if(is.na(cap)){
    highlim <- quantile(data, percentile)
  } else {
    highlim <- cap
  }
  lowlim <- quantile(data, 1-percentile)
  data[data<lowlim] <- lowlim
  data[data>highlim] <- highlim
  return(as.data.table(data))
}
