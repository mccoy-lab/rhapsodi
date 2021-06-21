#' This function computes the switch error rate 
#' 
#' This function computes the switch error rate in the vector case only
#' First, we find a match 0/1/-1 encoding between the prediction and truth genotypes. 0 = match, 1 | -1 = mismatch
#' We then find all locations that are mismatch. 
#' For switch error rate, however, we only want to find standalone mismatches as well as the first mismatch in a stretch of mismatches.
#' Therefore, we use the diff function to find the difference between adjacent elements in the match 0/1/-1 encoding
#' Then from this diff result, we pick only the locations before the mismatch locations and save this as comp_with_before. 
#' Knowing that mismatches are 1 or -1, we expect continuations of mismatches in comp_with_before to be | -1 - -1 | = 0, | -1 - 1 | 2, | 1 - 1 | = 0, | 1 - -1 | = 2. 
#' And since matches are 0, a new mismatch would be |0 - 1| = 1, or |0 - - 1| = 1.
#' The values in this vector will be 0, 1, or 2. If 0 | 2, then the value is a continuation of a mismatch stretch. If the value is 1, it's a new switch error 
#' Finally, we find the number that are equal to 1 in comp_with_before and this is the number of switch errors
#' We divide by num_snps to find the switch error rate
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
  comp_with_before <- abs(diff(match01_encoding))[which_mismatch - 1]
  switch_errors <- sum(comp_with_before == 1, na.rm=TRUE) #if comp_with_before values are 0 | 2, then the value is a continuation following another error; if 1, it's a new switch error
  to_return <- switch_errors / num_snps
  return (to_return)
}