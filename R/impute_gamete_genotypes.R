#' A function to assign the haplotypes of each allele on every gamete
#' 
#' This function builds and applies a hidden Markov model to categorize each allele on each gamete. 
#' It then fills the positions missing data with the nearest haplotype assignment.
#' 
#' @param dt matrix of gametes
#' @param complete_haplotypes Inferred parental haplotypes 
#' @param sequencing_error User-input for expected error in sequencing (default = 0.005) 
#' @param avg_recomb User-input for average recombination rate that can be expected for a chromosome (default=1)
#' @param threads User-input value for calling `pbmclapply` or `mclapply` (default = 2)
#' 
#' @return filled_gametes a tibble resulting from the HMM and fill_NAs function which has the imputed donor haplotypes for each gamete
#' 
#' @importFrom tibble as_tibble
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' 
#' @export

fill_gametes <- function(dt, complete_haplotypes, sequencing_error=0.005, avg_recomb = 1, threads=2){
  dt_recoded <- recode_gametes(dt, complete_haplotypes)
  #build the HMM
  hmm <- build_hmm(nrow(dt_recoded), sequencing_error, avg_recomb)
  if (requireNamespace("pbmcapply", quietly = TRUE)){
    
    imputed_gametes <- as_tibble(do.call(cbind, pbmcapply::pbmclapply(1:ncol(dt_recoded), 
                                                          function(x) run_hmm(dt_recoded, x, hmm),
                                                          mc.cores = getOption("mc.cores", threads))))
  } else {
    imputed_gametes <- as_tibble(do.call(cbind, mclapply(1:ncol(dt_recoded),
                                                        function(x) run_hmm(dt_recoded, x, hmm),
                                                        mc.cores = getOption("mc.cores", threads))))
  }
  
  if (requireNamespace("pbapply", quietly = TRUE)){
  
    filled_gametes <- as_tibble(do.call(cbind, pbapply::pblapply(1:ncol(imputed_gametes),
                                                        function(x) fill_na(imputed_gametes, x))))
  } else {
    filled_gametes <- as_tibble(do.call(cbind, lapply(1:ncol(imputed_gametes),
                                                      function(x) fill_na(imputed_gametes, x))))
  }
  colnames(filled_gametes) <- colnames(dt)
  return(filled_gametes)

}