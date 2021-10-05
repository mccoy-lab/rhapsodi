#' A function to drive and report meiotic recombination breakpoint finding
#' 
#' This function takes as input two booleans controlling whether the gamete haplotypes and genotypes remain smoothed before crossover discovery (directly from HMM)
#' or unsmoothed (replacing inferred HMM state with the original reads if they disagree)
#' Then the function runs recombination finding
#' It offers the option to avoid oversmoothing by superimposing initial haplotype assignments over each gamete. For example, if an 
#' allele assigned to h1 was changed by the model to h2, the early functions can fill the NAs to h2, but then un-smoothing 
#' will replace the singular h1 at the correct allele. This singular h1 could be an example of gene conversion or non-crossover.
#' 
#' @param original_gamete_data original sparse gamete data matrix
#' @param complete_haplotypes dataframe of phased diploid donor genotypes in two columns, each column corresponding with a haplotype from the donor
#' @param filled_gamete_data_list the output list from `impute_gamete_genotypes` which contains each gamete data matrix with haplotype info from the HMM and `fill_NA` functions
#' @param positions the genomic positions corresponding to SNP indices
#' @param smooth_crossovers boolean, default is TRUE, whether to use smoothed data for recombination finding. If `TRUE`, doesn't replace with original reads
#' @param smooth_imputed_genotypes boolean, default is FALSE, whether to use smoothed data for ending genotypes. If `TRUE`, doesn't replace with original reads
#' @param sampleName sample name of sample given to rhapsodi, default is "sampleT"
#' @param chrom chromosome of sample given to rhapsodi, default is "chromT"
#' @param threads an integer, default = 2, the number of cores to use when we use mclapply or the like
#' 
#' @return out a dataframe recomb_breaks, specifying the predicted recombination breakpoints for each gamete
#'                                        
#' @export
#' 
#' @importFrom dplyr right_join
#' @importFrom tibble tibble as_tibble
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' 
discover_meiotic_recombination <- function(original_gamete_data, complete_haplotypes, filled_gamete_data_list, positions, smooth_crossovers = TRUE, smooth_imputed_genotypes = FALSE, sampleName = "sampleT", chrom = "chrT", threads=2){
  idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gamete_data_list$filled_gametes_haps))
  if (smooth_crossovers){
    filled_gamete_forrecomb <- as_tibble(filled_gamete_data_list$filled_gametes_haps)
  } else {
    if(smooth_imputed_genotypes){
      filled_gamete_forrecomb <- as_tibble(unsmooth(recode_gametes(original_gamete_data, complete_haplotypes), filled_gamete_data_list$filled_gametes_haps))
    } else{ filled_gamete_forrecomb <- as_tibble(filled_gamete_data_list$unsmoothed_gametes_haps)}
  }
  if (requireNamespace("pbmcapply", quietly = TRUE)){
    recomb_breaks <- do.call(rbind, pbmcapply::pbmclapply(1:ncol(filled_gamete_forrecomb),
                                                             function(x) find_recomb_spots(filled_gamete_forrecomb, x, idents_for_csv, positions),
                                                             mc.cores=getOption("mc.cores", threads))) %>%
      dplyr::right_join(tibble(Ident = idents_for_csv), by = "Ident") %>% as.data.frame()
  } else {
    recomb_breaks <- do.call(rbind, mclapply(1:ncol(filled_gamete_forrecomb),
                                                function(x) find_recomb_spots(filled_gamete_forrecomb, x, idents_for_csv, positions),
                                                mc.cores=getOption("mc.cores", threads))) %>%
      dplyr::right_join(tibble(Ident = idents_for_csv), by = "Ident")  %>% as.data.frame()
  }
  return(recomb_breaks)
}