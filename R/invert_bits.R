#' A function to invert the values in a data frame or matrix
#' 
#' This function replaces 0s with 1s and 1s with 0s in a dataframe or matrix
#' 
#' @param input_data gamete alleles
#' 
#' @return input_data inverted from the actual input
#' 
invert_bits <- function(input_data) {
  return(abs(input_data-1))
}