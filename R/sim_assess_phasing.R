#'
#'
#'
sim_assess_phasing <- function(true_donor_haps, pred_donor_haps, num_snps){
  assess_phasing <- list()
  
  ##Phasing assessment
  assess_phasing$acc <- 100 - min(sim_hamming_distance_ignoreNA(true_donor_haps$donor1, pred_donor_haps$h1, num_snps),
                                          sim_hamming_distance_ignoreNA(true_donor_haps$donor2, pred_donor_haps$h1, num_snps))
  assess_phasing$com <- sim_completeness(pred_donor_haps$h1, num_snps)
  assess_phasing$lhs <- max(sim_lhs(true_donor_haps$donor1, pred_donor_haps$h1),
                                    sim_lhs(true_donor_haps$donor2, pred_donor_haps$h1))
  assess_phasing$ser <- min(sim_ser(true_donor_haps$donor1, pred_donor_haps$h1, num_snps),
                                    sim_ser(true_donor_haps$donor2, pred_donor_haps$h1, num_snps))
  return(assess_phasing)
}