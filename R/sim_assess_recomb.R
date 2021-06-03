#'
#'
#'
sim_assess_recomb <- function(true_recomb, pred_recomb, cons){
  ##Recombination Discovery assessment
  true_recomb_nona <- true_recomb[!is.na(true_recomb$start),] %>% setkey()
  true_recomb_na <- true_recomb[is.na(true_recomb$start),]
  
  pred_recomb_dt <- data.table(gam=sapply(strsplit(pred_recomb$Ident, "_"), `[`, 3), start=pred_recomb$Genomic_start, end=pred_recomb$Genomic_end)
  pred_recomb_nona <- pred_recomb_dt[!is.na(pred_recomb_dt$start),] %>% setkey()
  pred_recomb_na <- pred_recomb_dt[is.na(pred_recomb_dt$start),]
  
  if (nrow(truth_recomb_nona) > 0){
    no_truths = FALSE
  } else {no_truths = TRUE}
  if (nrow(pred_recomb_nona) > 0){
    no_preds = FALSE
  } else {no_preds = TRUE}
  
  if (!no_truths & !no_preds){
    truth_intersect <- foverlaps(true_recomb_nona, pred_recomb_nona) %>% `colnames<-`(c("gam", "Predicted_Start", "Predicted_End", "True_Start", "True_End"))
    pred_intersect <- foverlaps(pred_recomb_nona, true_recomb_nona) %>% `colnames<-`(c("gam", "True_Start", "True_End", "Predicted_Start", "Predicted_End"))
  } else{
    truth_intersect <- NULL
    pred_intersect <- NULL
  }
  
  tn <- sim_find_tn(truth_recomb_na, pred_recomb_na)
  fn <- sim_find_fn(truth_intersect, true_recomb_nona, cons=cons, no_truths=no_truths, no_preds=no_preds)
  tp <- sim_find_tp(truth_intersect, cons=cons, no_truths=no_truths, no_preds=no_preds)
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