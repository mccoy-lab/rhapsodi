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
  if ((tp+fp) > 0){
    precision <- tp/(tp+fp)
    fdr <- fp/(tp+fp)
  } else{
    precision <- NA
    fdr <- NA
  }
  if ((tp+fn) >0){
    recall <- tp/(tp+fn)
  } else{
    recall <- NA
  }
  if ((tn+fp) > 0){
    specificity <- tn/(tn+fp)
    fpr <- fp/(tn+fp)
  } else{
    specificity <- NA
    fpr <- NA
  }
  if ((tp+fp+tn+fn) > 0){
    accuracy <- (tp + tn)/(tp + tn + fp + fn)
  } else{
    accuracy <- NA
  }
  f1 <- (2*precision*recall)/(precision+recall)
  metric_list <- list(precision=precision,
                      recall=recall,
                      accuracy=accuracy,
                      specificity=specificity,
                      fdr = fdr,
                      fpr = fpr,
                      f1=f1)
  return (metric_list)
}