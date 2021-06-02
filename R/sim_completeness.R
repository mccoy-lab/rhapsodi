sim_completeness <- function(predicted, num_snps){
  if (is.null(dim(predicted))){
    num_nas <- sum(is.na(predicted))
    to_return <- 1 - (num_nas/num_snps)
  } else {
    num_nas_byCol <- colSums(is.na(predicted))
    to_return <- 1 - (num_nas_byCol/num_snps)
  }
  return (to_return)
}