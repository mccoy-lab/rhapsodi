#'
#'
#'
sim_generate_gam <- function(donor_haplotypes, n_crossovers){
  init_hap_index <- sample(1:2, 1)
  init_hap <- donor_haplotypes[,init_hap_index]
  if (n_crossovers == 0){
    return(list(NA, init_hap, rep(init_hap_index, length(init_hap))))
  } else {
    n_snps <- length(init_hap)
    crossover_indices <- sample(2:(n_snps-1), n_crossovers) #limiting so that crossover locations can't be at the beginning or end to avoid an edge effect
    crossover_indices <- crossover_indices[order(crossover_indices)]
    recombined_hap <- init_hap
    recombined_hap_index <- rep(init_hap_index-1, n_snps)#recode 1's to 0's or 2's to 1's and store the initial haplotype at all SNP locations
    for (crossover in crossover_indices) { #switch haplotypes at some (unknown) location/exchange point between the indices crossover-1 and crossover
      recombined_hap <- c(recombined_hap[1:(crossover-1)], abs(1-recombined_hap[crossover:n_snps]))
      recombined_hap_index <- c(recombined_hap_index[1:(crossover-1)], abs(recombined_hap_index[crossover:n_snps]-1))
    }
    return (list(crossover_indices, recombined_hap, recombined_hap_index+1)) #re-recode recombined_hap_index so that 1's are 2's and 0' are 1's
  }
}