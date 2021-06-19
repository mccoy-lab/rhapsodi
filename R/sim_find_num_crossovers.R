#' This function finds the number of crossovers for each gamete
#' 
#' This function uses a poisson distribution with mean equal to `recomb_lambda` to randomly generate the number of crossovers
#' for each gamete for `num_gametes` gametes
#' 
#' @param num_gametes an integer, the number of gametes for which to generate a number of crossovers
#' @param recomb_lambda a numeric, the average recombination rate, used as the mean for the poisson distribution
#' 
#' @return num_recomb_sites a vector of length `num_gametes` containing integers for the number of crossovers for each gamete
#'
sim_find_num_crossovers <- function(num_gametes, recomb_lambda){
  stopifnot(recomb_lambda > 0)
  num_recomb_sites <- rpois(num_gametes, recomb_lambda)
  message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))
  return(num_recomb_sites)
}