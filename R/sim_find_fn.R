#' This function finds the number of false negatives when comparing two sets of recombinantion breakpoints
#'
#' This function finds the number of false negatives when comparing by finding the locations for which nnothing 
#' intersects true recombination locations.  
#' When multiple predictions intersect a single truth, in the conservative approach, we consider all but one of the
#' the intersections to be false negatives because we in essence aren't predicting a true recomb spot, or we're falsely 
#' saying there's not a recombination spot there.  
#' 
#' @param truth_intersect_dt a data table from foverlaps intersecting the truth and the predicted data with columns of `gam`, `Predicted_Start`, `Predicted_End`, `True_Start`, `True_End`
#' @param truth_dt_nona a data table of just the truth data with no NAs in `start` or `stop` columns.
#' @param cons a bool, default is FALSE; only applicable when multiple predictions intersect a single truth. if `cons` is TRUE, we take a conservative approach and consider all but one of these multiple intersections as false negatives   
#' @param no_truths a bool, default is FALSE, if true, that means there were only NAs in the truth information and no real recombination break points
#' @param no_preds a bool, default is FALSE, if true, that means there wer only NAs in the prediction information and no predicted recombination break points
#' 
#' @return fn an integer, the number of false negatives
#' 
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