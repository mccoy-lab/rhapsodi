sim_find_tn <- function(truth_dt_na, pred_dt_na){
  tn <- nrow(merge(truth_dt_na, pred_dt_na, by="gam")) #both truth and predicted return NA for a given gamete meaning that they both agree that gamete has no recombination events
  return (tn)
}