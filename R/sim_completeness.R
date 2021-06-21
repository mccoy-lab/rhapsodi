#' This function computes the completeness of a haplotype
#' 
#' This function computes the completeness of a haplotype by finding the number of NAs
#' and then dividing the number of NAs by the total number of genotypes or SNPs in the haplotype
#' Finally, we subtract this ratio from 1 to find the ratio of complete or non-NA genotypes
#' in the hapltoype. Note that this function works for a single vector/haplotype or for a matrix with
#' multiple haplotypes. In the first case, we return a single value. 
#' In the latter case, we return a vector of values.
#' 
#' @param predicted the rhapsodi predicted haplotype vector or haplotypes matrix/data frame
#' @param num_snps an integer, the number of SNPs or the length of the haplotype
#' 
#' @return to_return the completeness of the input. If `predicted` is a vector, `to_return` is a numeric. else if `predicted` is 2-dimensional, `to_return` is a vector of numerics 
#'
sim_completeness <- function(predicted, num_snps){
  if (is.null(dim(predicted))){
    num_nas <- sum(is.na(predicted))
    to_return <- 1 - (num_nas/num_snps)
  } else {
    num_nas_byCol <- colSums(is.na(predicted))
    to_return <- 1 - (num_nas_byCol/num_snps)
  }
  return (to_return)
}