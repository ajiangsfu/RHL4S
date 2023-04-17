#' Empirical Bayes' probability calculation
#' @description This is a function to calculate Empirical Bayes' probability,
#'   which is called by other functions but it can also be called directly.
#' @details
#' The coding is written based on Wright 2003 for LPS (linear prediction score) approach,
#' but it can be used for other approaches such as PS and PRPS as well.
#' For LPS, empirical Bayesian probability is required since LPS itself does not have a natural cutoff for classification,
#' alternatively, if the testing and training data are comparable, we can also make classification based on LPS cutoffs
#' but this is not commonly used. For PS and PRPS, 0 is a natural cutoff for two group classification, however,
#' this probability calculation step is still useful if we allow UNCLASS group in the final classification besides the two types
#' of classes from the training.
#'
#' When a Empirical Bayes' probability is calculated, by default, the 1st group in the input mean and sd vectors is treated as the
#' test group. When we calculate the probabilities, we first calcualte probability that a sample belongs to either group,
#' and then use the following formula to get Empirical Bayes' probability:
#' \eqn{prob(x) = d_test(x)/(d_test(x) + d_ref(x))}
#' Here prob(x) is the Empirical Bayes' probability of a given sample, d_test(x) is the density value assuming that a given sample
#' belongs to the test group, d_ref(x) is the density value assuming that a given sample belongs to the reference group.
#'
#' @param inscore a classification score, which can be any types of scores
#' @param groupMeans a numeric vector of two items: two classification score means for two training groups/classes
#' @param groupSds a numeric vector of two items: two classification score standard deviations for two training groups/classes
#' @return A probability for a sample belong to a group
#' @keywords Empirical Bayes' probability
#' @author Aixiang Jiang
#' @references
#' Wright G, Tan B, Rosenwald A, Hurt EH, Wiestner A, Staudt LM. A gene expression-based method
#' to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S
#' A. 2003 Aug 19;100(17):9991-6.
#' @export
getProb = function(inscore, groupMeans, groupSds){
  ### assume groupMeans contain 2 values for 2 group, and the 1st one is for positive group
  ### assume groupSds contain 2 values for 2 group, and the 1st one is for positive group
  #### d1 and d0 are density value for the given observed score under normal assumption
  d1 = dnorm(inscore,mean = groupMeans[1], sd= groupSds[1])
  d0 = dnorm(inscore,mean = groupMeans[2], sd= groupSds[2])
  return(d1/(d1+d0))
}


