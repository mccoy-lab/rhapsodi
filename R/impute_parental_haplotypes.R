#' A function to impute haplotypes from diploid donor
#' 
#' This function infers the diploid donor haplotypes by first calling `split_with_overlap` to find the overlapping SNP segments, 
#' then `reconstruct_hap` for majority voting within each overlapping segment, finally stitching together the overlapping segments to phase both haplotypes
#' 
#' @export
#' 
#' @param dt Matrix of gamete alleles 
#' @param positions vector of SNP positions 
#' @param window_length Size of window (default=3000)
#' @param overlap_denom User-input value for denominator in calculation of overlap (default = 2)
#' @param threads User-input value for calling `pbmclapply` (default = 2)
#' @param mcstop User-input value for continuing phasing even if mean_concordance isn't within advised bounds (default = FALSE)
#' 
#' @return complete_haplotypes phased haplotypes as a tibble with column names index, pos (for SNP positions), h1 (haplotype 1), & h2 (haplotype 2)
#' 
#' @example 
#' R code here showing my function works 
#' 
impute_parental_haplotypes <- function(dt, positions, window_length=3000, overlap_denom=2, threads=2, mcstop=TRUE){
  #Find overlapping windows
  windows <- split_with_overlap(rank(positions), window_length, overlap_denom)
  #Infer the haplotypes within the overlapping windows (a window at a time)
  inferred_haplotypes <- pbmclapply(1:length(windows),
                                    function(x) reconstruct_hap(dt, positions, windows[[x]]),
                                    mc.cores = getOption("mc.cores", threads))
  #Stitch together the haplotypes
  initial_haplotype <- inferred_haplotypes[[1]]
  for (hap_window in 1:length(windows)){
    olap_haps <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index")
    olap_haps_complete <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index", all = TRUE)
    mean_concordance <- mean(olap_haps$h1.x == olap_haps$h1.y, na.rm=TRUE)
    if (mean_concordance < 0.1) {
      olap_haps_complete$h1.y <- invert_bits(olap_haps_complete$h1.y)
    } else if (mean_concordance < 0.9) {
      if (mcstop){
        stop(paste0("Haplotypes within overlapping windows are too discordant to merge with a mean concordance of ", mean_concordance, ". rhapsodi is exiting"))
      } else{
        message(paste0("Haplotypes within overlapping windows are too discordant for confident merging with a mean concordance of ", mean_concordance, ", but continuing."))
      }
    }
    initial_haplotype <- tibble(index = olap_haps_complete$index,
                                pos = c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$pos.x,
                                        olap_haps_complete[!is.na(olap_haps_complete$pos.x) &
                                                             !is.na(olap_haps_complete$pos.y),]$pos.x,
                                        olap_haps_complete[is.na(olap_haps_complete$pos.x),]$pos.y),
                                h1 = c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$h1.x,
                                       olap_haps_complete[!is.na(olap_haps_complete$pos.x) &
                                                            !is.na(olap_haps_complete$pos.y),]$h1.x,
                                       olap_haps_complete[is.na(olap_haps_complete$pos.x),]$h1.y))
  }
  complete_haplotypes <- initial_haplotype %>%
    mutate(h2 = invert_bits(h1))
  return(complete_haplotypes)
}