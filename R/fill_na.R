#' A function to fill in missing data from each gamete 
#' 
#' This function fills in missing data (NAs) on each gamete. For each gamete, it fills the NA values with the nearest haplotype.
#' If the two adjacent haplotypes are not the same (i.e., at a recombination breakpoint), it leaves the values as NA. 
#' 
#' @param imputed_gametes Output of `run_hmm` which assigned a donor haplotype to each segment of each gamete
#' @param col_index Each column of `imputed_gametes`, pulled via `apply` function 
#' 
#' @return gamete_sample_imputed Column with each gamete's imputed haplotypes 
#' 
#' @importFrom dplyr rename mutate 
#' @importFrom tidyr fill replace_na
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' 
fill_na <- function(imputed_gametes, col_index) {
  gamete_sample <- imputed_gametes[,col_index] %>%
    dplyr::rename(gamete = colnames(.data)[1]) %>%
    dplyr::mutate(gamete_up = gamete_sample$gamete) %>%
    dplyr::mutate(gamete_down = gamete_sample$gamete) %>%
    tidyr::fill(gamete_sample$gamete_up, .direction = "up") %>%
    tidyr::fill(gamete_sample$gamete_down, .direction = "down") %>%
    dplyr::mutate(is_match = (gamete_sample$gamete_up == gamete_sample$gamete_down)) %>%
    tidyr::replace_na(list(is_match = FALSE))
  gamete_sample$gamete_imputed <- as.character(NA)
  gamete_sample[gamete_sample$is_match == TRUE,]$gamete_imputed <- gamete_sample[gamete_sample$is_match == TRUE,]$gamete_up
  #fill beginning of chromosome NAs
  first <- which(!is.na(gamete_sample$gamete_imputed))[1]
  gamete_sample$gamete_imputed[1:(first-1)] <- gamete_sample$gamete_imputed[first]
  #fill end of chromosome NAs
  gamete_sample$gamete_imputed <- rev(gamete_sample$gamete_imputed)
  first <- which(!is.na(gamete_sample$gamete_imputed))[1]
  gamete_sample$gamete_imputed[1:(first-1)] <- gamete_sample$gamete_imputed[first]
  #reverse chromosome imputation back so it faces the right way
  gamete_sample$gamete_imputed <- rev(gamete_sample$gamete_imputed)
  gamete_sample_imputed <- gamete_sample$gamete_imputed
  return(gamete_sample_imputed)
}