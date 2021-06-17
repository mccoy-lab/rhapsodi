#' This function changes random indices to NA, in effect removing data
#' 
#' This function takes as input a matrix and carries over its data to a vector
#' Then `num_nas` number of locations are randomly chosen in that vector and overwritten by NA
#' Finally, the vector is reshaped to a matrix matching the size and dimension of the input and returned
#' This function is efficient for lower numbers of NAs (i.e. higher coverage samples)
#' 
#' @param to_change a matrix, the original data, with nrow of `num_snps` and ncol of `num_gametes`
#' @param num_nas an integer, the number of positions that should be randomly chosen to switch to NA
#' @param num_gametes an integer; the number of gametes or the number of columns for `to_return`
#' @param num_snps an integer; the number of SNPs or the number of rows for `to_return`
#'
#' @return to_return a matrix with nrow of `num_snps` and ncol of `num_gametes`. The data is a combination of locations carried over from the input `to_change` and locations that were switched to NAs
#'
sim_add_na_flatten <- function(to_change, num_nas, num_gametes, num_snps){
  to_change <- as.vector(to_change)
  coords_to_change <- sample(1:(num_snps*num_gametes), size=num_nas, replace = FALSE)
  to_change[coords_to_change] <- NA
  to_return <- matrix(to_change, ncol=num_gametes, nrow=num_snps)
  return (to_return)
}