#'
#'
#'
sim_add_seq_error <- function(num_snps, num_gametes, seqError_add, gam_mat_with_na){
  num_genotypes <- num_snps * num_gametes
  num_nas <- sum(is.na(gam_mat_with_na))
  num_bits_to_flip <- as.integer(seqError_add * (num_genotypes - num_nas))
  switched_bit_mat <- c(abs(1-gam_mat_with_na)) #make a switched bit matrix compared to gam_mat_with_na
  where_locs <- sample(which(!is.na(switched_bit_mat)), size=num_bits_to_flip) #randomly pick num_bits_to_flip locations
  gam_mat_with_na <- c(gam_mat_with_na)
  gam_mat_with_na[where_locs] <- switched_bit_mat[where_locs] #take those locations from switched bit matrix and put them in place in the gam_mat_with_na matrix
  gam_mat_with_na <- matrix(gam_mat_with_na, nrow=num_snps, ncol=num_gametes)
  return(gam_mat_with_na)
}