#' This function simulates fully known gamete genotypes given a specified number of meiotic recombination events
#' 
#' This function runs portions of the generative model to simulate fully known gamete genotypes
#' for a specified number of gametes and number of SNPs
#' as well as a specified number of meiotic recombination events for each gamete
#'
#' @param num_gametes an integer, the number of gametes to simulate; also the number of columns -1 of the output
#' @param num_snps an integer, the number of SNPs to simulate; also the number of rows of the output
#' @param n_crossovers a vector of integers, whose length is `num_gametes`
#' @param random_seed an integer, the random seed which will be set for reproducibility; default is 42
#'
#' @return gam_full_df a dataframe of 0's and 1's with num_gametes + 1 columns, and num_snps rows. The first column is the SNP index and has a column name of "positions". The rest of the columns have a name corresponding to which gamete number it is (i.e. "gam101_")
#' 
#' @importFrom magrittr %>% 
#' 
#' @export

sim_generate_global <- function(num_gametes, num_snps, n_crossovers, random_seed=42){
  set.seed(random_seed)
  
  #n_crossovers <- sample(n_crossovers, size = length(n_crossovers))
  
  donor_haps <- sim_generate_donor_hap(num_snps)
  
  sim_gam <- lapply(1:num_gametes, function(x) sim_generate_gam(donor_haps, n_crossovers[x]))
  
  gam_mat <- sapply(sim_gam, "[[", 2)
  gam_full_df <- data.frame(pseudo_pos = 1:nrow(gam_mat), gam_mat) %>% `colnames<-`(c("positions", paste0("gam", 1:num_gametes, "_")))
  
  return(gam_full_df)
}