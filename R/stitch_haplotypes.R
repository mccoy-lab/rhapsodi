#' A function to stitch overlapping windows together into the two donor haplotypes
#' 
#' This function stitches together the haplotypes of overlapping windows by considering the amount of shared concordance between each window
#' Specific behavior depends on the parameters. There are two main modes by which this function operates. 
#' The first is stringent stitching, or the confident phasing mode, which uses two preset threshold values to determine whether or not the windows originate from the same donor
#' The second is non-stringent stitching which uses a single user-defined threshold value to determine whether or not the windows originate from the same donor
#' In both cases, consider there being a `different_max` threshold, where if the concordance between two windows is less than this value, we would say that the windows originate from different donors
#' and then there being a `same_min` threshold, where if the concordance between two windows is greater than this value, we would say that the windows originate from the same donor
#' In the first stringent stitching case, these are different values, 0.1 and 0.9 respectively. 
#' By default if mean concordance is between these two values, rhapsodi will exit
#' However, the user may set the `mcstop` argument to FALSE and rhapsodi will continue phasing, asking which threshold the mean concordance is closer to in the 1dimensional sense (the distance between two points on a line segment) and stitching continues, acting accordingly.
#' In the second non-stringent stitching case, these are the same value, recommended and default value of 0.5, but user may change.
#' Theoretically, stringent stitching with `mcstop` of FALSE is equivalent to non-stringent stitching with a stitch_new_min of 0.5
#' Mean concordance is found by finding the mean number of locations that have equal values, removing NAs from consideration since majority voting may return NA for some positions
#' After going through all the windows, the function finishes stitching together the whole first haplotype and finds the second by inverting the first haplotype.
#' 
#' @param inferred_haplotypes a list of tibbles from calling `reconstruct_hap` within `impute_donor_haplotypes` where each tibble is the inferred haplotype sorted by position for that window 
#' @param windows a list of lists from `split_with_overlap` within `impute_donor_haplotypes` where the number of lists equal to the number of windows and each inner list contains the SNP positions in that overlapping segment
#' @param stringent_stitch a bool, this parameter is used to determine the threshold values used in determining whether two windows originate from the same donor. If TRUE, the preset thresholds of 0.1 and 0.9 are used.
#' @param stitch_new_min a numeric >0, but <1; default is 0.5; this parameter is only evaluated if `stringent_stitch` is FALSE and is dually assigned as the `different_max` and `same_min` threshold values when considering the concordance between two windows and therefore which donors they originate from (same or different). 
#' @param mcstop a bool, only considered if `stringent_stitch` is TRUE; default is TRUE; this parameter is used to determine whether phasing continues or exits if the mean concordance between two windows is between 0.1 and 0.9. If TRUE, rhapsodi exits. If FALSE, rhapsodi and phasing continues, asking which threshold is closer and acting accordingly
#' 
#' @return complete_haplotypes a data frame with two columns, h1 and h2, containing the inferred haplotypes from each window stitched together in a single non-overlapping segment 
#' 
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' 
stitch_haplotypes <- function(inferred_haplotypes, windows, mcstop=TRUE, stringent_stitch=TRUE, stitch_new_min=0.5){
  if (stringent_stitch){
    different_max <- 0.1
    same_min <- 0.9
  } else{
    stopifnot(stitch_new_min < 1 & stitch_new_min > 0)
    different_max <- stitch_new_min
    same_min <- stitch_new_min
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
      } else{ #which threshold are we closer to?
        same_min_dif <- abs(mean_concordance - same_min)
        different_max_dif <- abs(mean_concordance - different_max)
        if (same_min_dif <= different_max_dif){
          message(paste0("Haplotypes within overlapping windows are too discordant for confident merging with a mean concordance of ", mean_concordance, " , but continuing and merging as the same haplotype."))
        } else {
          olap_haps_complete$h1.y <- invert_bits(olap_haps_complete$h1.y)
          message(paste0("Haplotypes within overlapping windows are too discordant for confident merging with a mean concordance of ", mean_concordance, " , but continuing and stitching windows as opposite haplotypes."))
        }
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
    dplyr::mutate(h2 = invert_bits(h1)) %>% as.data.frame()
  return (complete_haplotypes)
}