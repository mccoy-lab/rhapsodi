#' A function to drive and report the unsmoothing (if desired) and recombination finding
#' 
#' This function takes as innput two booleans controlling wehther the reported genotypes are smoothed
#' or unsmoothed (replacing inferred HMM state with the original reads if they disagree)
#' and whether the data given to recombination finding is smoothed or unsmoothed. Then the function runs recombination finding
#' and eventually drives repporting/exporting the cool data
#' It offers the option to avoid oversmoothing by superimposing initial haplotype assignments over each gamete. For example, if an 
#' allele assigned to h1 was changed by the model to h2, the early functions can fill the NAs to h2, but then un-smoothing 
#' will replace the singular h1 at the correct allele. This singular h1 could be an example of gene conversion or non-crossover.
#' 
#' @param smooth_crossovers boolean whether to use smoothed data for recombination finding. If `TRUE`, doesn't replace with original reads
#' @param smooth_imputed_genotypes boolean whether to use smoothed data for ending genotypes. If `TRUE`, doesn't replace with original reads
#' @param complete_haplotypes dataframe of phased diploid donor genotypes in two columns, each column corresponding with a haplotype from the donor
#' @param original_gamete_data original sparse gamete data matrix
#' @param filled_gamete_data gamete data matrix from the HMM and `fill_NA` function
#' @param sampleName sample name of sample given to rhapsodi
#' @param chrom chromosome of sample given to rhapsodi
#' 
#' @export
#' 
#' @example
#' R code showing how my function works
#' 
report_gametes <- function(smooth_crossovers, smooth_imputed_genotypes, complete_haplotypes, original_gamete_data, filled_gamete_data, sampleName, chrom){
  if (!smooth_crossovers){
    filled_gamete_forrecomb <- unsmooth(original_gamete_data, filled_gamete_data) #filled_gamete_forrecomb is haplotypes
    idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gamete_forrecomb))
    recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gamete_forrecomb),
                                                  function(x) find_recomb_spots(filled_gamete_forrecomb, x, idents_for_csv, positions),
                                                  mc.cores=getOption("mc.cores", threads))) %>% 
      right_join(., tibble(Ident = idents_for_csv), by = "Ident")
    if (!smooth_imputed_genotypes){
      filled_gamete_recode <- re_recode_gametes(filled_gamete_forrecomb, complete_haplotypes) #filled_gamete_recode is 0's and 1's
      #want to report filled_gamete_forrecomb (the haplotypes)
    } else{ #smooth_imputed_genotypes is TRUE
      filled_gamete_recode <- re_recode_gametes(filled_gamete_data, complete_haplotypes)
      #want to report filled_gamete_data (the haplotypes)
    }
  } else { #smooth_crossovers is TRUE
    idents_for_csv <- paste0(paste0(samplenName, "_", chrom, "_"), colnames(filled_gamete_data))
    recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gamete_data),
                                                  function(x) find_recomb_spots(filled_gamete_data, x, idents_for_csv, positions),
                                                  mc.cores=getOption("mc.cores", threads))) %>%
      right_join(., tibble(Ident = idents_for_csv), by = "Ident")
    if (!smooth_imputed_genotypes){
      filled_gamete_data <- unsmooth(original_gamete_data, filled_gamete_data) #haplotypes
      filled_gamete_recode <- re_recode_gamtes(filled_gamete_data, complete_haplotypes) #0's and 1's
    } else { #smooth_imputed_genotypes is TRUE
      filled_gamete_recode <- re_recode_gamtes(filled_gamete_data, complete_haplotypes) #0's and 1's
    }
   #want to report filled_gamete_data (the haplotypes)
  }
  #want to report recomb_spots_all, filled_gamete_recode (the 0's and 1's)
}