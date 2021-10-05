#' A function to turn a matrix of gamete haplotypes by position back to the 0/1 reads by
#' position, using the phased donor haplotypes as a guide
#' 
#' This function builds a matrix with position by row and gametes by column such that each cell
#' is a 0 or a 1 or an NA based on whether that cell in the input gamete matrix was from haplotype1, haplotype2, 
#' or was an NA. Then the 0 or 1 is found in the complete_haplotypes (phased donors) input at the corresponding
#' positions 
#'
#' @param dt input gamete haplotype data in tibble form
#' @param complete_haplotypes dataframe of phased donors with two columns (h1 and h2) and SNP positions as rows 
#'
#' @return to_return gamete data in dataframe form with read data (0's and 1's and NAs) instead of haplotype information
#' 
re_recode_gametes <- function(dt, complete_haplotypes) {
  to_return <- data.frame(matrix(NA_real_, nrow=nrow(dt), ncol=ncol(dt)))
  for (i in 1:ncol(dt)) {
    locs_h1 <- dt[,i] == "h1"
    locs_h1[which(is.na(locs_h1))] <- FALSE
    locs_h2 <- dt[,i] == "h2"
    locs_h2[which(is.na(locs_h2))] <- FALSE
    to_return[locs_h1, i] <- complete_haplotypes$h1[locs_h1]
    to_return[locs_h2, i] <- complete_haplotypes$h2[locs_h2]
    colnames(to_return) <- colnames(dt)
  }
  return (to_return)
}