#' A function to drive the phasing and reporting of the two haplotypes from a single diploid donor
#'
#' This function infers the diploid donor haplotypes by first calling `split_with_overlap` to construct the overlapping SNP segments,
#' then `reconstruct_hap` for majority voting within each overlapping segment, 
#' finally stitching together the overlapping segments to phase both haplotypes using `stitch_haplotypes`
#'
#' @param dt Matrix or dataframe of gamete alleles in 0/1/NA encoding
#' @param positions vector of SNP positions
#' @param window_length An integer; size of window for `split_with_overlap` (default=3000)
#' @param overlap_denom An integer; value for denominator in calculation of overlap for `split_with_overlap` (default = 2)
#' @param threads An integer; Number of threads to use when calling `pbmclapply` or the like (default = 2)
#' @param mcstop A bool; whether to stop phasing if mean_concordance isn't within advised bounds during stitching in `stitch_haplotypes` (default = TRUE ). TRUE stops phasing, FALSE allows it to continue, finding the minimum linear distance between the concordance and the thresholds.
#' @param stringent_stitch A bool; whether or not to use the original bounds (0.1, 0.9 cutoff) in `stitch_haplotypes`, or a user specified threshold value `stitch_new_min`.
#' @param stitch_new_min A numeric/float, if `stringent_stitch` == FALSE, the new minimum value to use as a cutoff in `stitch_haplotypes`. If concordance is below this value, we assume the windows are opposite haplotypes, if greater than or equal, we assume the same haplotypes
#' @param calculate_window_size_bool A bool; whether or not to calculate the window size based on characteristics of the input dataset; default = FALSE
#' @param cov A numeric; the expected sequencing depth coverage, necessary only if wanting rhapsodi to calculate the preferred window size; default = NULL
#' @param ger A numeric; the expected genotyping error rate, necessary only if wanting rhapsodi to calculate the preferred window size; default = NULL
#' @param avgr A numeric; the expected average recombination rate, necessary only if wanting rhapsodi to calculate the preferred window size; default = NULL
#'
#' @return complete_haplotypes phased haplotypes as a data frame with column names index, pos (for SNP positions), h1 (haplotype 1), & h2 (haplotype 2)
#'
#' @importFrom parallel mclapply
#'
#' @export
#'
phase_donor_haplotypes <- function(dt, positions, window_length=3000, overlap_denom=2, threads=2, mcstop=TRUE, stringent_stitch=TRUE, stitch_new_min=0.5, calculate_window_size_bool = FALSE, cov=NULL, ger = NULL, avgr = NULL){
  #Find overlapping windows
  if (calculate_window_size_bool){
    window_length <- calculate_window_size(ncol(dt), cov, nrow(dt), ger, avgr)
  }
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
