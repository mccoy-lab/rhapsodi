#' This function computes and returns metrics given the number of true positives, false positives, true negatives, and false negatives that were observed in a prediction task
#' 
#' This function computes metrics for precision, recall, accuracy, specificity, false discovery rate, false positive rate, and the f1 score
#' returning the metrics in a named list and computing them from inputs of the number of true positives `tp`,
#' the number of false positives `fp`, the number of true negatives `tn`, the number of false negatives `fn`. 
#' For any metric, if the denominator is equal to 0, the returned value for that metric is NA.
#'
#' @param tp an integer for the number of true positives
#' @param fp an integer for the number of false positives
#' @param tn an integer for the number of true negatives
#' @param fn an integer for the number of false negatives
#' 
#' @return metric_list a named list returning precision, recall, accuracy, specificity, fpr, fdr, and f1 (the f1 score)
#'
sim_metrics <- function(tp, fp, tn, fn){
  if (sum(tp, fp, na.rm = TRUE) > 0){
    precision <- tp/sum(tp, fp, na.rm = TRUE)
    fdr <- fp/sum(tp, fp, na.rm = TRUE)
  } else{
    precision <- NA
    fdr <- NA
  }
  if (sum(tp, fn, na.rm = TRUE) >0){
    recall <- tp/sum(tp, fn, na.rm = TRUE)
  } else{
    recall <- NA
  }
  if (sum(tn, fp, na.rm = TRUE) > 0){
    specificity <- tn/sum(tn, fp, na.rm = TRUE)
    fpr <- fp/sum(tn, fp, na.rm = TRUE)
  } else{
    specificity <- NA
    fpr <- NA
  }
  if (sum(tp, fp, tn, fn, na.rm = TRUE) > 0){
    accuracy <- sum(tp, tn, na.rm = TRUE)/sum(tp, tn, fp, fn, na.rm = TRUE)
  } else{
    accuracy <- NA
  }
  f1 <- prod(2,precision,recall, na.rm = TRUE)/sum(precision, recall, na.rm = TRUE)
  metric_list <- list(precision=precision,
                      recall=recall,
                      accuracy=accuracy,
                      specificity=specificity,
                      fdr = fdr,
                      fpr = fpr,
                      f1=f1)
  return (metric_list)
}