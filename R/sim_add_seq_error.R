#' This function inverts a bit at random indices in effect adding sequencing error to data
#' 
#' This function takes as input a matrix and carries over its inverted data to a vector and its original data to a second vector
#' After subsetting for just non-NA locations, the number of available genotypes to change is multiplied by the sequencing error rate to find the number of bits to flip
#' Then these number of random (non NA) indices are sampled 
#' and the inverted bit data is used to replace the original input data in these locations
#' Finallly, this data with the introduced sequencing errors is reshaped to a matrix of equal dimension and size as the input and returned
#' 
#' @param num_snps an integer; the number of snps or the number of rows of the input and output
#' @param num_gametes an integer; the number of gametes or the number of columns of the input and output
#' @param seqError_add a numeric; the sequencing error rate
#' @param gam_mat_with_na a matrix; the input sparsified gamete data with nrow of `num_snps` and ncol of `num_gametes`
#'
#' @return gam_mat_with_na a matrix; reflecting the input with `num_bits_to_flip` inverted bits with nrow of `num_snps` and ncol of `num_gametes`
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