#' A function that assigns the position in each gamete to the correct donor haplotype
#' 
#' This function walks along each gamete and replaces the original observation with the inferred state.
#' This replaces any incorrect haplotype assignments (e.g., due to sequencing error).
#'
#' @param gamete_dt Matrix of gametes with 0 and 1
#' @param column_index `apply` function cycles through gamete_dt to act on each column (i.e., each gamete) 
#' @param hmm the HMM that viterbi will be applied on
#' 
#' @import HMM  
#'    
#' @return original_obs Replaced observed haplotype with that assigned by the model
#' 
#' @example 
#' R code here showing my function works 
#' 
run_hmm <- function(dt, column_index, hmm) {
  original_obs <- dt[,column_index]
  inferred_state <- HMM::viterbi(hmm, na.omit(dt[, column_index]))
  original_obs[!is.na(original_obs)] <- inferred_state
  return(original_obs)
}