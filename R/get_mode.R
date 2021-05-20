#' A function to find the more common allele or NA at each SNP
#' 
#' This function gets the mode of a vector for majority voting or returns NA if there is no single mode
#' 
#' Adapted from  #from https://stackoverflow.com/questions/56552709/r-no-mode-and-exclude-na?noredirect=1#comment99692066_56552709
#' 
#' @param vector A subset of donor haplotypes
#' 
#' @return mode The most frequent value or NA at a position
#' 
#' @example 
#' R code here showing how my function works
#' 
get_mode <- function(vector){
  uniqv <- unique(na.omit(vector))
  tabv <- tabulate(match(vector, uniqv))
  if (length(uniqv) != 1 & sum(max(tabv) == tabv) > 1){
    if (is.character(uniqv)) return(NA_character_) else return(NA_real_)
  }
  max_tabv <- tabv == max(tabv)
  return(uniqv[max_tabv])
}