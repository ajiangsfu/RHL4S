
#' @title Get RHL4S risk group classification for a new MC-IF dataset
#' @description This function calculates RHL4S LPS (Linear Predictor Score) scores and empirical Bayes' probabilities, and generates RHL4S risk group classification for a new MC-IF dataset.
#' Prior information required includes the four selected RHL4S spatial scores and their weights, as well as the LPS score means and standard deviations for the RHL4S POS and NEG groups.
#' This information is embedded in the package and does not need to be provided by the user.
#' @details It is important to ensure that the new MC-IF data is comparable to the BCCA MC-IF data. If not, calibration may be required before calling this function.
#' @param newdat A pre-processed MC-IF data frame, where samples are in columns and spatial scores are in rows.
#' @param varsIn A list of model variables.
#' @return A data frame with LPS scores, empirical Bayesian probabilities for two groups, and RHL4S risk groop classification.
#' @keywords RHL4S, LPS, MC-IF
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

