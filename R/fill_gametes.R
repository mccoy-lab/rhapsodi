#' A function to assign the haplotypes of each allele on every gamete
#' 
#' This function builds and applies a hidden Markov model to categorize each allele on each gamete. 
#' It then fills the positions missing data with the nearest haplotype assignment. It offers the option to 
#' superimpose the original haplotype. 
#' 
#' @param dt matrix of gametes
#' @param complete_haplotypes Inferred parental haplotypes 
#' @param sequencing_error User-input for expected error in sequencing (default = 0.005) 
#' @param avg_recomb User-input for average recombination rate that can be expected for a chromosome (default=1)
#' @param threads User-input value for calling `pbmclapply` and `pblapply` (default = 2)

fill_gametes <- function(dt, complete_haplotypes, sequencing_error=0.005, avg_recomb = 1, threads=2){
  dt_recoded <- recode_gametes(dt, complete_haplotypes)
  #build the HMM
  hmm <- build_hmm(dt, sequencing_error, avg_recomb)
  imputed_gametes <- as.tibble(do.call(cbind, pbmclapply(1:ncol(dt_recoded), 
                                                         function(x) run_hmm(dt_recoded, x, hmm),
                                                         mc.cores = getOption("mc.cores", threads))))
  
  filled_gametes <- as.tibble(do.call(cbind, pblapply(1:ncol(imputed_gametes),
                                                      function(x) fill_na(imputed_gametes, x),
                                                      mc.cores = getOption("mc.cores", threads))))
  colnames(filled_gametes) <- colnames(dt)
  return(filled_gametes)

}