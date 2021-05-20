find_tp <- function(truth_intersect_dt, cons=FALSE, no_truths=FALSE, no_preds=FALSE){
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