#' This function is used to generate the fully known gamete genotype data given the diploid donor haplotypes. It also records gamete specific crossover exchange points
#'
#' This function generates a single gamete given a diploid donor's two haplotypes and the number of crossovers that should occur for a given gamete
#' First the function picks which donor haplotype corresponds to the beginning of the gamete's chromosome, using a uniform distribution
#' Then the function assigns the genotypes from this corresponding donor hapltoype to the whole gamete
#' When tracking the originating haplotype, we assign the 1 or 2 to the whole gamete. However, we then re-code this to 0 or 1 (by subtracting 1)
#' and at the end re-recode these back to 1's and 2's by adding 1
#' For each crossover index, i, we assume between index i-1 and index i the recombination breakpoint occurs 
#' and switch the bits for the rest of the chromosome genotypes to the opposite donor haplotype or genotype compared to the one just previous by inverting the bits
#' We repeat this when tracking the originating haplotypes, which is why we recoded them.
#' 
#' @param donor_haplotypes a data frame with the diploid donor haplotypes in two columns
#' @param n_crossovers an integer for the number of crossovers that should occur for this gamete
#' 
#' @return list with generated crossover indices by gamete (a vector of integers such taht each integer without loss of generality, i, represents an index location of a crossover exchange point, such that the crossover occured somewhere between SNP index location i-1 and i), the 0/1 encoded genotypes for each SNP position in each genotype,the haplotypes from which each SNP in each gamete originates (encoded as 0s and 1s) 
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