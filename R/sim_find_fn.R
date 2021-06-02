sim_find_fn <- function(truth_intersect_dt, truth_dt_nona, cons=FALSE, no_truths=FALSE, no_preds=FALSE){
  if (!no_truths & !no_preds){
    fn <- nrow(truth_intersect_dt[is.na(truth_intersect_dt$Predicted_Start),]) #number of true recombinantion spots that don't interesect any predictions
    if (cons){
      fn <- fn + sum(duplicated(paste0(truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$gam, "_", truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$Predicted_Start, "_", truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$Predicted_End))) #count the overflows when multiple truths intersect a single predicted
    }
  } else if (no_truths & !no_preds){
    fn <- 0
  } else if (!no_truths & no_preds){
    fn <- nrow(truth_dt_nona) ##no predictions, but some truths, i.e. no non-na predictions but there are non-na truths, all are false negatives
  } else {fn <- NA}
  return (fn)
}