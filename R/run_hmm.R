#' A function that assigns the position in each gamete to the correct donor haplotype
#' 
#' This function walks along each gamete and replaces the original observation with the inferred state.
#' If the number of original observations is 0 or 1, the viterbi algorithm is not called, as there is no "path" to evaluate.
#' If only one observation, we assign the one observation as the inferred state
#' If no observations, we assign an inferred state of NA, such that the inferred state (NA) simply matches the original observations (which were also NA in this situation)
#' This replaces any incorrect haplotype assignments (e.g., due to sequencing error).
#'
#' @param dt Matrix of gametes with 0 and 1
#' @param column_index `apply` function cycles through gamete_dt to act on each column (i.e., each gamete) 
#' @param hmm the HMM that viterbi will be applied on
#' 
#' @importFrom HMM viterbi
#' @importFrom stats na.omit
#'    
#' @return original_obs Replaced observed haplotype with that assigned by the model
#' 
run_hmm <- function(dt, column_index, hmm) {
  dict_list <- list("h1" = "haplotype1", "h2" = "haplotype2")
  original_obs <- dt[,column_index]
  na_omit_oo <- na.omit(original_obs)
  if (length(na_omit_oo == 1)){
    inferred_state <- dict_list[[na_omit_oo[1]]]
  } else if (length(na_omit_oo) == 0){
    inferred_state <- NA
  } else {
    inferred_state <- viterbi(hmm, na_omit_00)
  }
  original_obs[!is.na(original_obs)] <- inferred_state
  return(original_obs)
}