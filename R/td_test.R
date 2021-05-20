td_test <- function(gam_matrix, row_index, hap_encoding = TRUE){
  test_row <- gam_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  if (hap_encoding){
    one_count <- sum(gt_vector == "haplotype1", na.rm=TRUE)
    two_count <- sum(gt_vector == "haplotype2", na.rm=TRUE)
  } else{
    one_count <- sum(gt_vector == 1, na.rm=TRUE)
    two_count <- sum(gt_vector == 2, na.rm=TRUE)
  }
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}