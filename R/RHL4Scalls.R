#' @title A wrapper function to get RHL4S risk group classifications from MC-IF raw data
#' @description This function integrates MC-IF data for the given markers, calculates spatial scores, filters cells, summarizes to patient-level data,
#' calculates RHL4S LPS (Linear Predictor Score) scores and empirical Bayes' probabilities, and performs RHL4S risk group classification based on a probability cutoff.
#' Prior to using this function, make sure that the new MC-IF data is comparable to the BCCA MC-IF data. Calibration may be required if the data is not comparable.
#' @param datapath A path to the directory that contains 6 subfolders for each marker separately: CD4, CD20, CD30, CD68, CXCR5, and PD1.
#' @param patientFile A path to the xlsx or xls file that contains patient information, related slide ID, and low quality cores.
#' @param patientSheet The name of the sheet in patientFile that indicates the core used to calculate patient-level data.
#' @return A data frame with patient-level spatial scores, LPS scores, empirical Bayesian probabilities for two groups, and RHL4S risk group classification.
#' @keywords RHL4S, LPS, MC-IF
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
