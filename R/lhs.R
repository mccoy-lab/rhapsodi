lhs <- function(truth, predicted){
  match01_encoding <- predicted - truth #match is a 0, mismatch is a 1 | -1
  rleres <- rle(match01_encoding) #any match location will be 0
  to_return <- max(rleres$lengths[which(rleres$values == 0)]) #find longest length of matches
  return (to_return)
}