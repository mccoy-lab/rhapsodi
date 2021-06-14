#' A function to impute haplotypes from diploid donor
#' 
#' This function infers the diploid donor haplotypes by first calling `split_with_overlap` to find the overlapping SNP segments, 
#' then `reconstruct_hap` for majority voting within each overlapping segment, finally stitching together the overlapping segments to phase both haplotypes
#' 
#' @param dt Matrix of gamete alleles 
#' @param positions vector of SNP positions 
#' @param window_length Size of window (default=3000)
#' @param overlap_denom User-input value for denominator in calculation of overlap (default = 2)
#' @param threads User-input value for calling `pbmclapply` (default = 2)
#' @param mcstop User-input value for continuing phasing even if mean_concordance isn't within advised bounds (default = FALSE)
#' @param stringent_stitch Boolean, whether or not to use the original bounds (0.1, 0.9 cutoff)
#' @param stitch_new_min Numeric/float, if `stringent_stitch` == FALSE, the new minimum value to use as a cutoff. If concordance is below this value, we assume the windows are opposite haplotypes, if greater than or equal, we assume the same haplotypes
#' 
#' @return complete_haplotypes phased haplotypes as a tibble with column names index, pos (for SNP positions), h1 (haplotype 1), & h2 (haplotype 2)
#' 
#' @export
#' 
#' @example 
#' R code here showing my function works 
#' 
impute_donor_haplotypes <- function(dt, positions, window_length=3000, overlap_denom=2, threads=2, mcstop=TRUE, stringent_stitch=TRUE, stitch_new_min=0.5){
  #Find overlapping windows
  windows <- split_with_overlap(rank(positions), window_length, overlap_denom)
  #Infer the haplotypes within the overlapping windows (a window at a time)
  if (requireNamespace("pbmcapply", quietly = TRUE)){
    inferred_haplotypes <- pbmcapply::pbmclapply(1:length(windows),
                                      function(x) reconstruct_hap(dt, positions, windows[[x]]),
                                      mc.cores = getOption("mc.cores", threads))
  } else {
    inferred_haplotypes <- mclapply(1:length(windows),
                                    function(x) reconstruct_hap(dt, positions, windows[[x]]),
                                    mc.cores = getOption("mc.cores", threads))
  }
  #Stitch together the haplotypes
  complete_haplotypes <- stitch_haplotypes(inferred_haplotypes, windows, mcstop=mcstop, stringent_stitch=stringent_stitch, stitch_new_min=stitch_new_min)
  return(complete_haplotypes)
}