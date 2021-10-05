#' This function adds an index and a pos column to dataframes
#'
#' @param input_df the dataframe which the user wants to append index and pos columns to
#' @param position_vec the vector of positions which will be appended to the dataframe
#' 
#' @return output_df the dataframe with two new columns at the beginning, index and pos, as well as the original columns
#' 
#' @export
#' 
#' @importFrom magrittr %>%
#' 
cbind_pos <- function(input_df, position_vec){
  if (!is.null(input_df)){
    orig_colnames <- colnames(input_df)
    output_df <- cbind(1:length(position_vec), position_vec, input_df) %>% `colnames<-`(c(c("index", "pos"), orig_colnames))
    return (output_df)
  } else {return (input_df)}
}