#' A wrap-up function to get RHL4S calls from IHC raw data
#' @description This is the function to integrate IHC data for the given markers, calculate spatial score, filer cells, summarize to patient level data,
#'   and calculate RHL4S LPS (Linear Predictor Score) scores and empiracal Bayes' probabilities, and classification based on probability cutoff
#' @details Make sure the new IHC data is comparable to the BCCA IHC data, if not, calibration is required
#' @param datapath This data folder should contain 6 subfolder for each marker separately: CD30, CD4, CD68, CD20,PD1,CXCR5.
#' @param patientFile The xlsx or xls file with path should contain the patient information, related slide ID and low quality cores.
#' @param patientSheet The sheet name in patientFile to indicate the core used to calculate patient label data
#' @return A data frame with patient level spatial score, LPS score, Empirical Bayesian probabilities for two groups and classification
#' @keywords RHL4S, LPS, IHC
#' @author Aixiang Jiang
#' @export

RHL4Scalls = function(datapath = NULL, patientFile = NULL, patientSheet = "A"){
  dat = merge_markerFolders(path = datapath)
  scores = spatial_score(data = dat)
  patientDat = convert_patient(data = scores,  patient = patientFile, core = patientSheet)
  ## save patientDat for output
  outs = patientDat

  ## change the names
  colnames(patientDat) = gsub("pos", "", colnames(patientDat))
  colnames(patientDat) = gsub("PD1_T", "PD1_CD4", colnames(patientDat))
  colnames(patientDat) = gsub("CD68_Mac", "Mac", colnames(patientDat))

  res = predict_RHL4S(newdat = patientDat)
  colnames(res) = gsub("LPS", "RHL4S", colnames(res))

  outs = cbind(outs, res)

  return(outs)
}
