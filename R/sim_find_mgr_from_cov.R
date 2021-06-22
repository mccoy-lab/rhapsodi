#' This function computes the missing genotype rate given the coverage
#' 
#' This function computes the missing genotype rate of the simulation given the desired coverage
#' by using the poisson density function parameterized by x = 0 and lambda equal to the coverage metric
#' to return the probability of the density function, or the missing genotype rate
#' Note that the function checks to assure that the coverage is greater than 0
#' 
#' @param coverage a numeric, must be greater than 0. The desired coverage of sequencing (e.g. 0.01x, 0.1x, 1.204x, etc)
#' 
#' @return missing_genotype_rate a numeric, the corresponding missing genotype rate for the given coverage, < 1
#' 
#' @importFrom stats dpois
#' 
sim_find_mgr_from_cov <- function(coverage){
  stopifnot(coverage > 0)
  message(paste0("The coverage of this simulation is: ", coverage))
  missing_genotype_rate <- dpois(0, coverage)
  message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
  return (missing_genotype_rate)
}