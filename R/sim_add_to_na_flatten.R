#' This function overwrites NAs from random indices with "real" data from the input, in effect sparsifying the input data
#' 
#' This function takes as input a matrix with full "real" data
#' It then creates a vector of all NAs equal in length to the number of cells in the input matrix
#' Then `num_nas` number of locations are randomly chosen in that vector and overwritten by the corresponding indices from the flattened input data 
#' Finally, the vector is reshaped to a matrix matching the size and dimension of the input and returned
#' This function is efficient for higher numbers of NAs (i.e. lower coverage samples)
#' 
#' @param to_add_from a matrix, the original data, with nrow of `num_snps` and ncol of `num_gametes`
#' @param num_nas an integer, the number of positions that should be randomly chosen to switch to NA
#' @param num_gametes an integer; the number of gametes or the number of columns for `to_return`
#' @param num_snps an integer; the number of SNPs or the number of rows for `to_return`
#'
#' @return to_return a matrix with nrow of `num_snps` and ncol of `num_gametes`. The data is a combination of locations carried over from the input `to_add_From` and NAs
#'
sim_add_to_na_flatten <- function(to_add_from, num_nas, num_gametes, num_snps){
  to_return <- rep(NA, (num_snps*num_gametes))
  coords_to_keep_genotype <- sample(1:(num_snps*num_gametes), size=((num_snps*num_gametes)-num_nas), replace= FALSE)
  to_return[coords_to_keep_genotype] <- as.vector(to_add_from)[coords_to_keep_genotype]
  to_return <- matrix(to_return, ncol=num_gametes, nrow=num_snps)
  return (to_return)
}