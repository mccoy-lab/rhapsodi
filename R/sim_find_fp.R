#' This function finds the number of false positives when comparing two sets of recombination breakpoints
#' 
#' This function finds the number of true negatives by finding the predictions which don't intersect with true recombination spots at all
#' If there are no non-NA truths but some non-NA predictions, the number of false positives are equal to the number of non-NA predictions
#' because predictions are present for breakpoints that aren't there
#' If there are some non-NA truths but no non-NA predictions, we have 0 false positives
#' 
#' 
#' @param pred_intersect_dt a data table from foverlaps intersecting the predictions from rhapsodi and the truth data with columns of `gam`, `True_Start`, `True_End`, `Predicted_Start`, `Predicted_End`
#' @param pred_dt_nona a data table of just the predicted data from rhapsodi with no NAs in `start` or `stop` columns
#' @param no_truths bool, default is FALSE; if TRUE, that means there were only NAs in the truth information and no real recombination break points
#' @param no_preds bool, default is FALSE; if TRUE, that means there were only NAs in the prediction information and no predicted recombination break points 
#' 
#' @return fp an integer, the number of false negatives
#'
sim_find_fp <- function(pred_intersect_dt, pred_dt_nona, no_truths=FALSE, no_preds=FALSE){
  if (!no_truths & !no_preds){
    fp <- nrow(pred_intersect_dt[is.na(pred_intersect_dt$True_Start),]) #number predicted as recombination spots, but don't intersect with truth at all
  } else if(no_truths & !no_preds){
    fp <- nrow(pred_dt_nona) #no truths, but some predictions; i.e. no non-na truths, but there are non-na predictions, all are false positives
  } else if(!no_truths & no_preds){
    fp <- 0
  } else { fp <- NA}
  return (fp)
}