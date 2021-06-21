#' This function finds the number of true positives when comparing two sets of recombination breakpoints
#' 
#' This function finds the number of true positives by finding the number of true recombinantion breakpoints
#' which are intersected by a rhapsodi predicted recombination breakpoint. In liberal mode (default), every single
#' intersection is counted as a true positive, even if a single predicted breakpoint intersects multiple true breakpoints.
#' In the conservative mode, we only count one of these multiple intersections as a true positive.
#' If there are no non-NA truths or no non-NA predictions, there are 0 true positives
#'
#' @param truth_intersect_dt a data table from foverlaps intersecting the truth and the predicted data with columns `gam`, `Predicted_Start`, `Predicted_End`, `True_Start`, `True_End`
#' @param cons a bool, default is FALSE; only applicable when multiple predictions intersect a single truth. If `cons` is TRUE, we take a conservative approach and for each multiple intersection, we consider only one of these breakpoints as a true positive
#' @param no_truths a bool, default is FALSE. If TRUE, that means there were only NAs in the truth information and no real recombination break points
#' @param no_preds a bool, default is FALSE. If TRUE, that means there were only NAs in the prediction information from rhapsodi and no predicted recombination break points
#'
sim_find_tp <- function(truth_intersect_dt, cons=FALSE, no_truths=FALSE, no_preds=FALSE){
  if (!no_truths & !no_preds){
    tp <- nrow(truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]) #want to count +1 for every truth that intersects a prediction, even if multiple truths intersect a single prediction; accomplish this by counting number of truths that intersect any prediction from the non NA intersection
    if (cons){
      tp <- tp - sum(duplicated(paste0(truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$gam, "_", truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$Predicted_Start, "_", truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$Predicted_End))) #want to count only one truth when multiple truths intersect a single prediction; accomplish this by subtracting the sum of the boolean reporting which items are duplicated in the prediction identities from the non NA intersection; note duplicated() only returns TRUE for the 2nd, 3rd, 4th, etc occurrences, not the 1st
    }
  } else if(no_truths & !no_preds){
    tp <- 0
  } else if (!no_truths & no_preds){
    tp <- 0
  } else { tp <- NA }
  return (tp)
}