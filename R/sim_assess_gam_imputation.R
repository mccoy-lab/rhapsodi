#'
#'
#'
sim_assess_gam_imputation <- function(true_gam, pred_gam, num_snps, num_gametes){
  assess_gam <- list()
  
  ##Gamete Imputation assessment
  assess_gam$acc <- 100-sim_hamming_distance_ignoreNA(true_gam[,-1], pred_gam, num_snps)
  assess_gam$com <- sim_completeness(pred_gam, num_snps)
  assess_gam$lhs <- do.call(rbind, lapply(1:num_gametes,function(x) sim_lhs(true_gam[,-1][,x], pred_gam[,x])))
  assess_gam$ser <- do.call(rbind, lapply(1:num_gametes,function(x) sim_ser(true_gam[,-1][,x], pred_gam[,x], num_snps)))

  return(assess_gam)
}