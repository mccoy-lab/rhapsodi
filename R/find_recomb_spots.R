#' A function to find recombination breakpoint regions for a single gamete
#' 
#' This function finds recombination breakpoint regions for a single gamete by finding 
#' where adjacent non-NA haplotypes switch. Most likely large regions will be returned 
#' rather than 2bp windows because the sparsity of the original data constrains how 
#' many NAs will remain after filling and therefore we can only say that the true
#' exchange point occurs somewhere between these two non-matching haplotypes
#' 
#' @param input_gamete_data (tibble) SNPs correspond to the rows and gametes correspond to the columns
#' @param x which gamete number to search
#' @param identities gamete identities vector
#' @param genomic_positions Genomic SNP positions 
#' 
#' @return recomb_spots tibble of the recombination spots for gamete x
#' 
#' @importFrom dplyr mutate row_number arrange
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats complete.cases
#' 
find_recomb_spots <- function(input_gamete_data, x, identities, genomic_positions) {
  ident <- identities[x]
  input_tibble <- input_gamete_data[, x] %>%
    dplyr::mutate(.data, index = dplyr::row_number()) %>%
    dplyr::mutate(.data, positions = genomic_positions)
  complete_cases_tibble <- input_tibble[complete.cases(input_tibble),]
  input_vec <- as.factor(complete_cases_tibble[[1]])
  switch_indices <- which(input_vec[-1] != input_vec[-length(input_vec)])
  switch_indices_input <- complete_cases_tibble[switch_indices,]$index
  crossover_start <- input_tibble[switch_indices_input,]$positions
  rev_input_tibble <- dplyr::arrange(input_tibble, -index) %>%
   dplyr:: mutate(.data, index = dplyr::row_number())
  complete_cases_rev_tibble <- rev_input_tibble[complete.cases(rev_input_tibble),]
  rev_input_vec <- as.factor(complete_cases_rev_tibble[[1]])
  rev_switch_indices <- which(rev_input_vec[-1] != rev_input_vec[-length(rev_input_vec)])
  rev_switch_indices_input <- complete_cases_rev_tibble[rev_switch_indices,]$index
  crossover_end <- rev(rev_input_tibble[rev_switch_indices_input,]$positions)
  recomb_spots <- tibble(Ident = ident, Genomic_start = crossover_start, Genomic_end = crossover_end)
  return (recomb_spots)
}