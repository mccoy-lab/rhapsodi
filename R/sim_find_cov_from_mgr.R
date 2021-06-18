#' This function computes the coverage given the desired missing genotype rate
#' 
#' This function computes the coverage of the simulation given the desired missing genotype rate
#' by reversing the poisson density function and taking the negative natural log of the input missing genotype rate
#' 
#' @param missing_genotype_rate a numeric, if greater than 1, assumed to be a percentage and is divided by 100. 
#' 
#' @return coverage a numeric, the coverage of the simulation corresponding to the desired missing_genotype_rate
#'
sim_find_cov_from_mgr <- function(missing_genotype_rate){
  if (missing_genotype_rate > 1){ #if entered as a percentage
    missing_genotype_rate <- missing_genotype_rate / 100
  }
  stopifnot(missing_genotype_rate < 1)
  message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
  coverage <- -log(missing_genotype_rate)
  message(paste0("The coverage of this simulation is: ", coverage))
  return(coverage)
}