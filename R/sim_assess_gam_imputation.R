#' This function assesses the effectiveness of rhapsodi's gamete imputation 
#' 
#' This function compares the simulated full truth gamete data to the rhapsodi predicted filled gamete data,
#' producing a named list with vectors for lhs (largest haplotype segment), ser (switch error rate), acc (accuracy), com (completeness)
#' 
#' @param true_gam a matrix, from the output of the generative model, the true/full gamete genotypes where the rows are the SNPs and the columns are the gametes (except for the first column which is the SNP genomic positions)
#' @param pred_gam a matrix, from the output of rhapsodi, with the predicted/filled gamete genotypes where the rows are the SNPs and the columns are the gametes
#' @param num_snps the number of snps, should be equal to the number of rows of `true_gam` or `pred_gam`
#' @param num_gametes the number of gametes, should be equal to the number of cols - 1 of `true_gam` or just the number of cols of `pred_gam`
#' 
#' @return assess_gam a named list with vectors for `acc`, `com`, `lhs`, & `ser` 
#' 
#' @export
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