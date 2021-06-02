sim_ser <- function(truth, predicted, num_snps){
  match01_encoding <- predicted - truth #match is a 0, mismatch is a 1 | -1
  which_mismatch <- which(match01_encoding != 0)
  comp_with_before <- abs(diff(match01_encoding))[which_mismatch - 1]
  switch_errors <- sum(comp_with_before == 1, na.rm=TRUE) #if comp_with_before values are 0 | 2, then the value is a continuation following another error; if 1, it's a new switch error
  to_return <- switch_errors / num_snps
  return (to_return)
}