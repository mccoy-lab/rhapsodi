#'
#'
#'
sim_find_num_nas <- function(num_gametes, num_snps, missing_genotype_rate){
  num_genotypes <- num_gametes * num_snps
  num_nas <- floor(num_genotypes * missing_genotype_rate)
  return(num_nas)
}