find_fp <- function(pred_intersect_dt, pred_dt_nona, no_truths=FALSE, no_preds=FALSE){
  if (!no_truths & !no_preds){
    fp <- nrow(pred_intersect_dt[is.na(pred_intersect_dt$True_Start),]) #number predicted as recombination spots, but don't intersect with truth at all
  } else if(no_truths & !no_preds){
    fp <- nrow(pred_dt_nona) #no truths, but some predictions; i.e. no non-na truths, but there are non-na predictions, all are false positives
  } else if(!no_truths & no_preds){
    fp <- 0
  } else { fp <- NA}
  return (fp)
}