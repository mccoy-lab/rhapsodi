#' A function to fill in missing data from each gamete 
#' 
#' This function fills in missing data (NAs) on each gamete. For each gamete, it fills the NA values with the nearest haplotype.
#' If the two adjacent haplotypes are not the same (i.e., at a recombination breakpoint), it leaves the values as NA. 
#' 
#' @param imputed_gametes Output of `run_hmm` which assigned a donor haplotype to each segment of each gamete
#' @param col_index Each column of `imputed_gametes`, pulled via `apply` function 
#' @param fill_ends a boolean; if TRUE, fills the NAs at the terminal edges of chromosomes with the last known or imputed SNP (for end of chromosome) and the first known or imputed SNP (for beginning of chromosome); if FALSE, leaves these genotypes as NA; default = TRUE 
#' 
#' @return gamete_sample_imputed Column with each gamete's imputed haplotypes 
#' 
#' @importFrom dplyr rename mutate 
#' @importFrom tidyr fill replace_na
#' @importFrom magrittr %>%
#' 
fill_na <- function(imputed_gametes, col_index, fill_ends = TRUE) {
  gamete_sample <- imputed_gametes[,col_index] %>%
    dplyr::rename(gamete = colnames(.)[1]) %>%
    dplyr::mutate(gamete_up = gamete) %>%
    dplyr::mutate(gamete_down = gamete) %>%
    tidyr::fill(gamete_up, .direction = "up") %>%
    tidyr::fill(gamete_down, .direction = "down") %>%
    dplyr::mutate(is_match = (gamete_up == gamete_down)) %>%
    tidyr::replace_na(list(is_match = FALSE))
  gamete_sample$gamete_imputed <- as.character(NA)
  gamete_sample[gamete_sample$is_match == TRUE,]$gamete_imputed <- gamete_sample[gamete_sample$is_match == TRUE,]$gamete_up
  if (fill_ends){
    #fill beginning of chromosome NAs
    first <- which(!is.na(gamete_sample$gamete_imputed))[1]
    if (!is.na(first)){
      gamete_sample$gamete_imputed[1:(first-1)] <- gamete_sample$gamete_imputed[first]
    }
    #fill end of chromosome NAs
    gamete_sample$gamete_imputed <- rev(gamete_sample$gamete_imputed)
    first <- which(!is.na(gamete_sample$gamete_imputed))[1]
    if (!is.na(first)){
      gamete_sample$gamete_imputed[1:(first-1)] <- gamete_sample$gamete_imputed[first]
    }
    #reverse chromosome imputation back so it faces the right way
    gamete_sample$gamete_imputed <- rev(gamete_sample$gamete_imputed)
  }
  gamete_sample_imputed <- gamete_sample$gamete_imputed
  return(gamete_sample_imputed)
}