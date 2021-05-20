#' A function to impute parental haplotypes 
#' 
#' This function 
#' 
#' @export
#' 
#' @param dt Matrix of gamete alleles 
#' @param positions vector of SNP positions 
#' @param window_length Size of window (default=3000)
#' @param overlap_denom User-input value for denominator in calculation of overlap (default = 2)
#' @param threads User-input value for calling `pbmclapply` (default = 2)
#' 
#' @return complete_haplotypes
#' 
#' @example 
#' R code here showing my function works 
#' 
impute_parental_haplotypes <- function(dt, positions, window_length=3000, overlap_denom=2, threads=2){
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
      stop(paste0("Haplotypes within overlapping windows are too discordant to merge with a mean concordance of ", mean_concordance, ". rhapsodi is exiting"))
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