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
#' @param positions the genomic positions corresponding to SNP indices
#' @param sampleName sample name of sample given to rhapsodi
#' @param chrom chromosome of sample given to rhapsodi
#' @param threads an integer, default = 2, the number of cores to use when we use mclapply or the like
#' 
#' @return out a named list which returns gamete_haps, or a data frame specifying from which donor haplotype each gamete position originates
#'                                        gamete_genotypes, or a matrix specifying the genotype in (0's and 1's) for each gamete position
#'                                        recomb_breaks, or a tibble specifying the predicted recombination breakpoints for each gamete
#'                                        donor_haps, or phased haplotypes as a tibble with column names: index, pos (for SNP positions), h1 (haplotype 1), & h2 (haplotype 2)
#'                                        
#' @export
#' 
#' @importFrom dplyr right_join
#' @importFrom tibble tibble
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' 
report_gametes <- function(smooth_crossovers, smooth_imputed_genotypes, complete_haplotypes, original_gamete_data, filled_gamete_data, positions, sampleName, chrom, threads=2){
  out <- list() #want to return out$gamete_haps filled gamete data haplotypes, out$recomb_breaks, out$
  if (!smooth_crossovers){
    filled_gamete_forrecomb <- unsmooth(original_gamete_data, filled_gamete_data) #filled_gamete_forrecomb is haplotypes
    idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gamete_forrecomb))
    if (requireNamespace("pbmcapply", quietly = TRUE)){
      recomb_spots_all <- do.call(rbind, pbmcapply::pbmclapply(1:ncol(filled_gamete_forrecomb),
                                                              function(x) find_recomb_spots(filled_gamete_forrecomb, x, idents_for_csv, positions),
                                                              mc.cores=getOption("mc.cores", threads))) %>% 
        dplyr::right_join(.data, tibble(Ident = idents_for_csv), by = "Ident")
    } else {
      recomb_spots_all <- do.call(rbind, mclapply(1:ncol(filled_gamete_forrecomb),
                                                  function(x) find_recomb_spots(filled_gamete_forrecomb, x, idents_for_csv, positions),
                                                  mc.cores=getOption("mc.cores", threads))) %>% 
        dplyr::right_join(.data, tibble(Ident = idents_for_csv), by = "Ident")
    }
    if (!smooth_imputed_genotypes){
      filled_gamete_recode <- re_recode_gametes(filled_gamete_forrecomb, complete_haplotypes) #filled_gamete_recode is 0's and 1's
      out$gamete_haps <- filled_gamete_forrecomb #want to report filled_gamete_forrecomb (the haplotypes)
      
    } else{ #smooth_imputed_genotypes is TRUE
      filled_gamete_recode <- re_recode_gametes(filled_gamete_data, complete_haplotypes)
      out$gamete_haps <- filled_gamete_data #want to report filled_gamete_data (the haplotypes)
    }
  } else { #smooth_crossovers is TRUE
    idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gamete_data))
    if (requireNamespace("pbmcapply", quietly = TRUE)){
      recomb_spots_all <- do.call(rbind, pbmcapply::pbmclapply(1:ncol(filled_gamete_data),
                                                              function(x) find_recomb_spots(filled_gamete_data, x, idents_for_csv, positions),
                                                              mc.cores=getOption("mc.cores", threads))) %>%
      dplyr::right_join(.data, tibble(Ident = idents_for_csv), by = "Ident")
    } else {
      recomb_spots_all <- do.call(rbind, mclapply(1:ncol(filled_gamete_data),
                                                               function(x) find_recomb_spots(filled_gamete_data, x, idents_for_csv, positions),
                                                               mc.cores=getOption("mc.cores", threads))) %>%
        dplyr::right_join(.data, tibble(Ident = idents_for_csv), by = "Ident")  
    }
    if (!smooth_imputed_genotypes){
      filled_gamete_data <- unsmooth(original_gamete_data, filled_gamete_data) #haplotypes
      filled_gamete_recode <- re_recode_gametes(filled_gamete_data, complete_haplotypes) #0's and 1's
    } else { #smooth_imputed_genotypes is TRUE
      filled_gamete_recode <- re_recode_gametes(filled_gamete_data, complete_haplotypes) #0's and 1's
    }
   out$gamete_haps <- filled_gamete_data #want to report filled_gamete_data (the haplotypes)
  }
  out$recomb_breaks <- recomb_spots_all #want to report recomb_spots_all
  out$donor_haps <- complete_haplotypes #want to report the diploid donor haplotypes (genotypes within each a single h1 or h2 column is a single haplotype)
  out$gamete_genotypes <- filled_gamete_recode #want to report filled_gamete_recode (the 0's and 1's)
  return(out)
}