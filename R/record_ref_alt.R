#' This function can be used to append the reference and alternate alleles back to the output data 
#' 
#' This function can be called if the input data was originally nucleotide encoded, providing ref and alt alleles
#' Then this function will append a0 and a1 columns to 0/1 encoded outputs such that the nucleotide(s) associated with genotype 0 are recorded in a0
#' and the nucleotide(s) associated with genotype 1 are recorded in a1 
#' returning the original dataframes just with two new columns
#' 
#' @param to_output_df whatever dataframe you want to append the columns to
#' @param ref a vector of the reference nucleotides (those associated with the 0 genotypes)
#' @param alt a vector of the alternate nucleotides (those associated with the 1 genotypes)
#' 
#' @return to_return a dataframe with index, pos, ref, & alt columns before the remaining 0/1 encoded genotype columns. 
#' 
#' @importFrom magrittr %>%
#'
#' @export
#'

record_ref_alt <- function(to_output_df, ref, alt){
  if (!is.null(to_output_df)){
    to_return <- cbind(to_output_df$index, to_output_df$pos, ref, alt, to_output_df[,3:ncol(to_output_df)]) %>% `colnames<-`(c(c("index", "pos", "a0", "a1"), colnames(to_output_df[,3:ncol(to_output_df)])))
    return (to_return)}
  else{ return (to_output_df)}
}
 