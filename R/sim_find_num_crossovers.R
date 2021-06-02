#'
#'
#'
sim_find_num_crossovers <- function(num_gametes, recomb_lambda){
  stopifnot(recomb_lambda > 0)
  num_recomb_sites <- rpois(num_gametes, recomb_lambda)
  message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))
  return(num_recomb_sites)
}