#'
#'
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