#'
#'
#'
stitch_haplotypes <- function(inferred_haplotypes, windows, mcstop=TRUE, stringent_stitch=TRUE, stitch_new_min=0.5){
  if (stringent_stitch){
    different_max = 0.1
    same_min = 0.9
  } else{
    different_max = stitch_new_min
    same_min = stitch_new_min
  }
  initial_haplotype <- inferred_haplotypes[[1]]
  for (hap_window in 1:length(windows)){
    olap_haps <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by="index")
    olap_haps_complete <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by="index", all=TRUE)
    mean_concordance <- mean(olap_haps$h1.x == olap_haps$h1.y, na.rm=TRUE)
    if (mean_concordance < different_max){
      olap_haps_complete$h1.y <- invert_bits(olap_haps_complete$h1.y)
    } else if (stringent_stitch & mean_concordance < same_min){
      if (mcstop){
        stop(paste0("Haplotypes within overlapping windows are too discordant to merge with a mean concordance of ", mean_concordance, ". rhapsodi is exiting."))
      } else{
        message(paste0("Haplotypes within overlapping windows are too discordant for confident merging with a mean concordance of ", mean_concordance, " , but continuing and merging as the same haplotype."))
      }
    }
    initial_haplotype <- tibble(index = olap_haps_complete$index,
                                pos=c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$pos.x,
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
  return (complete_haplotypes)
}