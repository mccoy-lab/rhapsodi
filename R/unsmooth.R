#' A function to unsmooth or replace original reads for the gametes
#' 
#' This function finds where the resulting haplotype assignments differ from the original reads and replaces the imputed data with
#' the originally observed data, hence unsmoothing the HMM signal.
#' 
#' @param original_gamete_df original gamete data with haplotype by position 
#' @param filled_gamete_data filled gamete data from `fill_gametes`
#'
#' @return filled_gamete_data filled gamete data in tibble form with replaced original reads 
#' 
#' @importFrom tibble as_tibble
#' 
unsmooth <- function(original_gamete_df, filled_gamete_data) {
  original_gamete_df[original_gamete_df == "hap1"] <- "h1"
  original_gamete_df[original_gamete_df == "hap2"] <- "h2"
  original_dt <- as.data.frame(original_gamete_df)
  filled_gamete_data <- as.data.frame(filled_gamete_data)
  filled_gamete_data[!is.na(original_dt)] <- original_dt[!is.na(original_dt)]
  filled_gamete_data <- as_tibble(filled_gamete_data)
  return (filled_gamete_data)
}