#' This function computes the switch error rate 
#' 
#' This function computes the switch error rate in the vector case only
#' First, we find a match 0/1/-1 encoding between the prediction and truth genotypes. 0 = match, 1 | -1 = mismatch
#' We then find all locations that are mismatch. 
#' For switch error rate, however, we only want to find standalone mismatches as well as the first mismatch in a stretch of mismatches.
#' Therefore, we use the diff function to find the difference between adjacent elements in the match 0/1/-1 encoding
#' We sum those for which the difference is not 1, suggesting that it is a new standalone, or start of a stretch, 
#' We also add one to this sum for the very first mismatch
#' We finally divide by num_snps to find the switch error rate
#' 
#' @param truth a vector of genotypes in 0/1 encoding, the truth
#' @param predicted a vector of genotyptes in 0/1 encoding, predicted from rhapsodi
#' @param num_snps an integer, the number of genotypes or the length of the haplotype
#' 
#' @return to_return a numeric, the switch error rate
#' 
sim_ser <- function(truth, predicted, num_snps){
  match01_encoding <- predicted - truth #match is a 0, mismatch is a 1 | -1
  which_mismatch <- which(match01_encoding != 0)
  switch_errors <- sum(diff(which_mismatch) != 1) + 1 #adding one for the first mismatch
  to_return <- switch_errors / num_snps
  return (to_return)
}