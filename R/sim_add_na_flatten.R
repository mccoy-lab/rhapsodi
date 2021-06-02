#'
#'
#'
sim_add_na_flatten <- function(to_change, num_nas, num_gametes, num_snps){
  to_change <- as.vector(to_change)
  coords_to_change <- sample(1:(num_snps*num_gametes), size=num_nas, replace = FALSE)
  to_change[coords_to_change] <- NA
  to_return <- matrix(to_change, ncol=num_gametes, nrow=num_snps)
  return (to_return)
}