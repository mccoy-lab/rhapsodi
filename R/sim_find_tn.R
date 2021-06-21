#' This function finds the number of true negatives when comparing two sets of recombination breakpoints
#' 
#' This function find sthe number of true negatives by finding the gametes for which both truth and prediction say there are no 
#' recombination breakpoints (NA in start and stop columns)
#' 
#' @param truth_dt_na a data table for the truth or generated data in which the `start` and `stop` columns have NA, meaning no recombination breakpoints
#' @param pred_dt_na a data table for the predicted data from rhapsodi in which the `start` and `stop` columns have NA, meaning no predicted recombination breakpoints
#' 
#' @return tn an integer, the number of true negatives
#'
sim_find_tn <- function(truth_dt_na, pred_dt_na){
  tn <- nrow(merge(truth_dt_na, pred_dt_na, by="gam")) #both truth and predicted return NA for a given gamete meaning that they both agree that gamete has no recombination events
  return (tn)
}