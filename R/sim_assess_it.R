#'
#'
#'
sim_assess_it <- function(true_donor_haps, pred_donor_haps, true_recomb, pred_recomb, true_gam, pred_gam, num_snps, num_gametes, write_out_plot, cons=FALSE){
  assess_phasing_out <- sim_assess_phasing(true_donor_haps, pred_donor_haps, num_snps)
  
  assess_gam_imputation_out <- sim_assess_gam_imputation(true_gam, pred_gam, num_snps, num_gametes) 
  
  assess_recomb_out <- sim_assess_recomb(true_recomb, pred_recomb, cons)
  
  all_metrics <- list(phasing = assess_phasing_out, gam_imputation=assess_gam_imputation_out, recomb=assess_recomb_out)
  return (all_metrics)
}