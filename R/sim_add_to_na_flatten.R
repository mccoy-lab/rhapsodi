#'
#'
#'
sim_add_to_na_flatten <- function(to_add_from, num_nas, num_gametes, num_snps){
  to_return <- rep(NA, (num_snps*num_gametes))
  coords_to_keep_genotype <- sample(1:(num_snps*num_gametes), size=((num_snps*num_gametes)-num_nas), replace= FALSE)
  to_return[coords_to_keep_genotype] <- as.vector(to_add_from)[coords_to_keep_genotype]
  to_return <- matrix(to_return, ncol=num_gametes, nrow=num_snps)
  return (to_return)
}