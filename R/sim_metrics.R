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