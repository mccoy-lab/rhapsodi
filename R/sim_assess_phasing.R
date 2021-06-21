#' This function assesses the effectiveness of rhapsodi's diploid donor haplotype phasing
#' 
#' This function compares the simulated truth diploid donor haplotypes to the rhapsodi predicted phased data,
#' producing a named list of values for lhs (largest haplotype segment), ser (switch error rate), acc (accuracy), com (completeness)
#' 
#' @param true_donor_haps a data frame of phased donor haplotypes from the generative model with column names of `donor1` and `donor2`
#' @param pred_donor_haps a tibble of phased donor haplotypes from rhapsodi with column names `index`, `pos` (for SNP positions), `h1` (haplotype 1) & `h2` (haplotype 2)
#' @param num_snps the number of SNPs, should be equal to the number of rows of `true_donor_haps` or `pred_donor_haps` 
#' 
#' @return assess_phasing a named list with values for `acc`, `com`, `lhs`, & `ser`
#' 
#' @export
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