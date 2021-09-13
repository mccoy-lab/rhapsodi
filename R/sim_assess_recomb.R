#' This function assesses the effectiveness of rhapsodi's meiotic recombination discovery
#' 
#' This function compares the simulated truth gamete specific recombination breakpoints to the rhapsodi predicted gamete specific breakpoints,
#' producing a named list with values for `precision`, `recall`, `accuracy`, `specificity`, `fdr` (false discovery rate), `fpr` (false positive rate), `f1` (f1 score),
#' `true_n` (number of true recombination breakpoints), `pred_n` (number of predicted recombination breakpoints), `tn` (# of true negatives), `fn` (# of false negatives), `tp` (# of true positives), `fp` (# of false positives)
#' 
#' @param true_recomb a data.table data table containing the true recombination breakpoints from the generative model with columns `gam`, `start`, `end`
#' @param pred_recomb a tibble containing the predicted recombination breakpoints from rhapsodi with columns `Ident` `Genomic_start`, `Genomic_end`
#' @param cons a bool; default = FALSE, If TRUE, compares recombination breakpoints in a conservative manner such that if two or more true breakpoints intersect a single predicted breakpoint, we only consider one intersection to be a tp and the rest to be fn. If FALSE, all are tp.
#'
#' @return metrics_out a named list with values for `precision`, `recall`, `accuracy`, `specificity`, `fdr`, `fpr`, `true_n`,`pred_n`, `tn`, `fn`, `tp`, `fp`
#'
#' @import data.table
#' @importFrom magrittr %>%
#'
#' @export
#'
sim_assess_recomb <- function(true_recomb, pred_recomb, cons=FALSE){
  ##Recombination Discovery assessment
  true_recomb <- data.table::data.table(true_recomb)
  true_recomb_nona <- true_recomb[!is.na(true_recomb$start),] %>% data.table::setkey()
  true_recomb_na <- true_recomb[is.na(true_recomb$start),]
  
  pred_recomb_dt <- data.table(gam=sapply(strsplit(pred_recomb$Ident, "_"), `[`, 3), start=pred_recomb$Genomic_start, end=pred_recomb$Genomic_end)
  pred_recomb_nona <- pred_recomb_dt[!is.na(pred_recomb_dt$start),] %>% data.table::setkey()
  pred_recomb_na <- pred_recomb_dt[is.na(pred_recomb_dt$start),]
  
  if (nrow(true_recomb_nona) > 0){
    no_truths = FALSE
  } else {no_truths = TRUE}
  if (nrow(pred_recomb_nona) > 0){
    no_preds = FALSE
  } else {no_preds = TRUE}
  
  if (!no_truths & !no_preds){
    truth_intersect <- data.table::foverlaps(true_recomb_nona, pred_recomb_nona) %>% `colnames<-`(c("gam", "Predicted_Start", "Predicted_End", "True_Start", "True_End"))
    pred_intersect <- data.table::foverlaps(pred_recomb_nona, true_recomb_nona) %>% `colnames<-`(c("gam", "True_Start", "True_End", "Predicted_Start", "Predicted_End"))
  } else{
    truth_intersect <- NULL
    pred_intersect <- NULL
  }
  
  tn <- sim_find_tn(true_recomb_na, pred_recomb_na)
  message(tn)
  fn <- sim_find_fn(truth_intersect, true_recomb_nona, cons=cons, no_truths=no_truths, no_preds=no_preds)
  tp <- sim_find_tp(truth_intersect, cons=cons, no_truths=no_truths, no_preds=no_preds)
  message(tp)
  fp <- sim_find_fp(pred_intersect, pred_recomb_nona, no_truths=no_truths, no_preds=no_preds)
  
  metrics_out <- sim_metrics(tp, fp, tn, fn)
  metrics_out$true_n <- nrow(true_recomb)
  metrics_out$pred_n <- nrow(pred_recomb_dt)
  metrics_out$tn <- tn
  metrics_out$fn <- fn
  metrics_out$tp <- tp
  metrics_out$fp <- fp
  
  return(metrics_out)
}