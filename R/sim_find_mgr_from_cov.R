#'
#'
#'
sim_find_mgr_from_cov <- function(coverage){
  stopifnot(coverage > 0)
  message(paste0("The coverage of this simulation is: ", coverage))
  missing_genotype_rate <- dpois(0, coverage)
  message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
  return (missing_genotype_rate)
}