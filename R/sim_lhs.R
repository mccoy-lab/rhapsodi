#' This function computes the largest haplotype segment 
#' 
#' This function computes the largest haplotype segment in the vector case only
#' First, we find a match 0/1/-1 encoding between the prediction and truth genotypes. 0 = match, 1 | -1 = mismatch
#' Then we use the run length encoding to find the lengths and values of runs of equal values in the match 0/1/-1 encoding,
#' specifically we're interested in stretches of 0's. 
#' So we look for values in the run length encoding that are 0
#' Then from those we extract the corresponding lengths and find the maximum length, or our largest hapltoype segment
#' 
#' @param truth a vector of genotypes in 0/1 encoding, the truth
#' @param predicted a vector of genotyptes in 0/1 encoding, predicted from rhapsodi
#' 
#' @return to_return an integer, the largest haplotype segment
#' 
sim_lhs <- function(truth, predicted){
  match01_encoding <- predicted - truth #match is a 0, mismatch is a 1 | -1
  rleres <- rle(match01_encoding) #any match location will be 0
  to_return <- max(rleres$lengths[which(rleres$values == 0)]) #find longest length of matches
  return (to_return)
}