#' This function computes the number of NAs that should be used to replace genotypes given the size of the genotype matrix and the missing genotype rate
#' 
#' This function computes the number of NAs that are needed to achieve a certain missing genotype rate
#' given the size of the input genotype matrix which is num_snps rows and num_gametes columns. 
#' 
#' @param num_gametes an integer, the number of gametes or the number of columns of the genotype matrix
#' @param num_snps an integer, the number of snps or the number of rows of the genotype matrix 
#' @param missing_genotype_rate a numeric, < 1, the missing genotype rate
#'
#' @return num_nas an integer, the product of the missing genotype rate and the number of genotypes (`num_snps` * `num_gametes`)
#'
sim_find_num_nas <- function(num_gametes, num_snps, missing_genotype_rate){
  num_genotypes <- num_gametes * num_snps
  num_nas <- floor(num_genotypes * missing_genotype_rate)
  return(num_nas)
}