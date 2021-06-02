sim_hamming_distance_ignoreNA <- function(truth, predicted, num_snps){
  if (is.null(dim(truth))){
    num_mismatch <- sum((predicted - truth) != 0, na.rm=TRUE)
    to_return <- num_mismatch / num_snps * 100
  } else {
    num_mismatch_byCol <- colSums((predicted - truth) != 0, na.rm=TRUE)
    to_return <- num_mismatch_byCol / num_snps * 100
  }
  return (to_return)
}