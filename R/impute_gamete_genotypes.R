#' A function to assign the haplotypes of each allele on every gamete
#' 
#' This function builds and applies a hidden Markov model to categorize each allele on each gamete. 
#' It then fills the positions missing data with the nearest haplotype assignment.
#' If the user asks for unsmoothed genotypes (i.e. replacing original sequencing reads where HMM imputation disagrees with these original reads) by setting `smooth_imputed_genotypes` to FALSE
#' then the unsmooth function is called to replace imputed genotypes with original sequencing reads
#' 
#' @param original_gamete_data original matrix of gametes
#' @param complete_haplotypes Inferred parental haplotypes 
#' @param sequencing_error User-input for expected error in sequencing (default = 0.005) 
#' @param avg_recomb User-input for average recombination rate that can be expected for a chromosome (default=1)
#' @param smooth_imputed_genotypes a bool, default is FALSE, whether to use smoothed data for ending genotypes. If `TRUE`, doesn't replace with original reads, returning smoothed data only. If `FALSE`, will return both smoothed and unsmoothed
#' @param threads User-input value for calling `pbmclapply` or `mclapply` (default = 2)
#' 
#' @return gamete_data a named list with two data frames (names filled_gametes and unsmoothed_gametes) resulting from the HMM, fill_NAs, and potentially the unsmooth function, returning the imputed donor haplotypes for each gamete. In the filled_gametes output, the dataframe represents the direct output. In the unsmoothed_gametes output, if `smooth_imputed_genotypes` is TRUE, this is NULL; if FALSE, original sequencing reads replace imputed genotypes if they disagree
#' 
#' @importFrom tibble as_tibble
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' 
#' @export

impute_gamete_genotypes <- function(original_gamete_data, complete_haplotypes, sequencing_error=0.005, avg_recomb = 1, smooth_imputed_genotypes=FALSE, threads=2){
  complete_haplotypes <- as_tibble(complete_haplotypes)
  dt_recoded <- recode_gametes(original_gamete_data, complete_haplotypes)
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
  
    filled_gametes <- as.data.frame(do.call(cbind, pbapply::pblapply(1:ncol(imputed_gametes),
                                                        function(x) fill_na(imputed_gametes, x))))
  } else {
    filled_gametes <- as.data.frame(do.call(cbind, lapply(1:ncol(imputed_gametes),
                                                      function(x) fill_na(imputed_gametes, x))))
  }
  colnames(filled_gametes) <- colnames(original_gamete_data)
  filled_gametes_01 <- re_recode_gametes(filled_gametes, complete_haplotypes)
  if (!smooth_imputed_genotypes){
    filled_gametes_unsmooth <- unsmooth(recode_gametes(original_gamete_data, complete_haplotypes), as_tibble(filled_gametes)) %>% as.data.frame() %>% `colnames<-`(colnames(original_gamete_data))
    filled_gametes_unsmooth_01 <- re_recode_gametes(filled_gametes_unsmooth, complete_haplotypes)
  } else{ 
    filled_gametes_unsmooth <- NULL
    filled_gametes_unsmooth_01 <- NULL}
  gamete_data <- list(filled_gametes = filled_gametes_01,
                      filled_gametes_haps = filled_gametes,
                      unsmoothed_gametes = filled_gametes_unsmooth_01,
                      unsmoothed_gametes_haps = filled_gametes_unsmooth)
  return(gamete_data)

}