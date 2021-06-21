#' This function computes the hamming distance (or the opposite of accuracy), ignoring NAs rather than penalizing for their presence
#' 
#' This function computes the hamming distance by finding mismatches between the truth and prediction, but ignoring NAs
#' We divide the number of mismatches by the number of SNPs and return this ratio multiplied by 100
#' Accuracy can be found by subtracting the returned value(s) from 100
#' Note that this function works for vectors of haplotypes or for matrices with
#' multiple haplotypes. In the first case, we return a single value. 
#' In the latter case, we return a vector of values.
#'
#' @param truth either a vector of true genotypes or a matrix with multiple vectors of true genotypes, where genotypes are encoded as 0/1
#' @param predicted either a vector of rhapsodi precited genotyptes or a matrix with multiple vectors of rhapsodi predicted genotypes, dimension matching that of `truth`, and genotypes are encoded as 0/1
#' @param num_snps an integer, the number of SNPs, or the number of genotypes/length of the haplotype(s)
#' 
#' @return to_return numeric(s), in percentage format, if `truth` & `predicted` are vectors, to_return is a single value; else if the inputs are matrices, to_return is a vector where the value of to_return is the hamming distance
#'
sim_hamming_distance_ignoreNA <- function(truth, predicted, num_snps){
  if (is.null(dim(truth))){
    num_mismatch <- sum((predicted - truth) != 0, na.rm=TRUE)
    to_return <- num_mismatch / num_snps * 100
  } else {
    num_mismatch_byCol <- colSums((predicted - truth) != 0, na.rm=TRUE)
    to_return <- num_mismatch_byCol / num_snps * 100
  }
  return (to_return)
}