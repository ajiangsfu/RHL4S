#' A function to get RHL4S calls
#' @description This is the function to calculate RHL4S LPS (Linear Predictor Score) scores and empiracal Bayes' probs,
#' and get RHL4S calls for a new IHC data set. The required prior information is:
#' 4 selected RHL4S spatial scores and their weights, LPS score means and sds for RHL4S POS and NEG groups
#' @details Make sure the new IHC data is comparable to the BCCA IHC data, if not, calibration is required before calling this function
#' @param newdat A new IHC data frame, samples are in columns, and spatial scores are in rows.
#'  Notice that data are already pre-processed
#' @param varsIn The model variable list
#' @return A data frame with LPS score, Empirical Bayesian probabilities for two groups and classification
#' @keywords RHL4S, LPS, IHC
#' @author Aixiang Jiang
#' @export

predict_RHL4S = function(newdat, varsIn = c("CXCR5_HRS_spatial_score", "PD1_CD4_spatial_score",
                                            "Mac_spatial_score", "CXCR5_B_spatial_score")){

  ## should load the weights and LPS parameter within the package
  load(system.file("extdata", "rhl4s.rda", package = "RHL4S"))

  ## the 1st challenge is to make variable comparable
  wts = rhl4s$IMCtopMean_weights

  ## I do need the varsIn should contain the following info
  ## CXCR5_HRS   PD1_CD4       Mac   CXCR5_B
  ## and should make sure that the terms are in a new data set
  ## when they are some inconsisent, deal it outside of this function
  ## actually, since I need to put all together, I should make sure the names are matched

  ords = sapply(names(wts), function(xx) {
    grep(xx, varsIn)
  })

  varsIn = varsIn[ords]

  lps = data.matrix(newdat[,varsIn]) %*% wts

  ## I should calibrate lps score
  lps = lps - rhl4s[[2]]

  ## the following line is changed accordingly
  LPS_prob = getProb(lps, groupMeans = (rhl4s[[3]])[,1], groupSds = (rhl4s[[3]])[,2])

  ## actually the probcut is in the model
  probcut = rhl4s[[4]]
  LPS_prob_class = ifelse(LPS_prob >= probcut, "High", "Low")

  ## lps, LPS_prob, LPS_prob_class
  outs = data.frame(cbind(lps, LPS_prob))
  colnames(outs) = c("LPS_score", "LPS_prob")
  outs$LPS_prob_class = LPS_prob_class[,1]
  return(outs)
}
