#' A function to fill in missing data from each gamete 
#' 
#' This function fills in missing data (NAs) on each gamete. For each gamete, it fills the NA values with the nearest haplotype.
#' If the two adjacent haplotypes are not the same (i.e., at a recombination breakpoint), it leaves the values as NA. 
#' It offers the option to avoid oversmoothing by superimposing initial haplotype assignments over each gamete. For example, if an 
#' allele assigned to h1 was changed by the model to h2, this function can fill the NAs to h2, but replace the singular h1
#' at the correct allele. This could be an example of gene conversion or non-crossover. 
#' 
#' @param imputed_gametes Output of `run_hmm` which assigned a parental haplotype to each segment of each gamete
#' @param col_index Each column of `imputed_gametes`, pulled via `apply` function 
#' 
#' @return gamete_sample_imputed Column with each gamete's imputed haplotypes 
#' 
#' @example 
#' R code here showing my function works 
#' 
fill_na <- function(imputed_gametes, col_index) {
  gamete_sample <- imputed_gametes[,col_index] %>%
    rename(gamete = colnames(.)[1]) %>%
    mutate(gamete_up = gamete) %>%
    mutate(gamete_down = gamete) %>%
    fill(gamete_up, .direction = "up") %>%
    fill(gamete_down, .direction = "down") %>%
    mutate(is_match = (gamete_up == gamete_down)) %>%
    replace_na(list(is_match = FALSE))
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