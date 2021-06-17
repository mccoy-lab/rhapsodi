#' This function drives the assessment of all 3 parts of rhapsodi: phasing, gamete imputation, recombination discovery
#' 
#' This function drives all 3 steps of assessment in comparing the simulated full truth data with the rhapsodi predicted data
#' First it assesses donor haplotypte phasing, producing a named list with single values for lhs (longest haplotype segment), ser (switch error rate), acc (accuracy), com (completeness)
#' Then it assesses gamete genotype imputation, producing a named list with vectors for lhs (longest haplotype segment), ser (switch error rate), acc (accuracy), com (completeness)
#' Then it assesses recombination discovery producing a named list with single values for precision, recall, accuracy, specificity, fdr (false discovery rate), fpr (false positive rate) f1 (f1 score), 
#' true_n (number of true recombination breakpoints), pred_n (number of predicted recombination breakpoints), tn (true negative), fn (false negative), tp (true positive), fp (false positive)
#' Finally, it returns a list of named lists where `phasing` contains the phasing assessment named list
#' `gam_imputation` contains the gamete genotype imputation assessment named list
#' and `recomb` contains the recombination breakpoint discovery assessment named list
#' 
#' @param true_donor_haps a data frame of phased donor haplotypeps from the generative model with column names of `donor1` and `donor2`
#' @param pred_donor_haps a tibble of phased donor haplotypes from rhapsodi with column names `index`, `pos` (for SNP positions) `h1` (haplotype 1), & `h2` (haplotype 2)
#' @param true_recomb a data.table data table containing the true recombination breakpoints from the generative model with columns `gam`, `start`, `end` 
#' @param pred_recomb a tibble containing the predicted recombination breakpoints from rhapsodi with columns `Ident`, `Genomic_start`, `Genomic_end`
#' @param true_gam a matrix, from the output of the generative model, the true/full gamete genotypes where the rows are the SNPs and the columns are the gametes
#' @param pred_gam a matrix, from the output of rhapsodi, the predicted/filled gamete genotypes where the rows are the SNPs and the columns are the gametes
#' @param write_out_plot a bool; default=FALSE, might not be applicable
#' @param cons a bool; default=FALSE, If TRUE, compares recombination breakpoints in a conservative manner such that if two or more true breakpoints intersect a single predicted breakpoint, we only consider one intersection to be a tp and the rest to be fn.  
#' 
#' @return all_metrics a named list of named lists with all the assessment metric values or vectors
#' 
#' @export
#'
sim_assess_it <- function(true_donor_haps, pred_donor_haps, true_recomb, pred_recomb, true_gam, pred_gam, write_out_plot=FALSE, cons=FALSE){
  num_snps <- nrow(true_gam)
  num_gametes <- ncol(true_gam[,-1])
  
  assess_phasing_out <- sim_assess_phasing(true_donor_haps, pred_donor_haps, num_snps)
  
  assess_gam_imputation_out <- sim_assess_gam_imputation(true_gam, pred_gam, num_snps, num_gametes) 
  
  assess_recomb_out <- sim_assess_recomb(true_recomb, pred_recomb, cons)
  
  all_metrics <- list(phasing = assess_phasing_out, gam_imputation=assess_gam_imputation_out, recomb=assess_recomb_out)
  return (all_metrics)
}