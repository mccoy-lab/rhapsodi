#' This function generates the diploid donor haplotypes of length `num_snps` with only heterozygous sites
#'
#' This function generates the diploid donor haplotypes by first generating a single haplotype
#' using a uniform distribution to pick the reference (a 0) or the alternate (a 1) allele
#' for each position for `num_snps` positions.
#' Then for the second donor haplotype we assume that each of these SNP positions are heterozygous 
#' and just invert the bits of the first haplotypes to make this so. 
#'
#' @param num_snps an integer, the number of SNP positions for which we want to generate data
#'
#' @return donor_haps a data frame with the truth phased diploid donor haplotypes with two columns `donor1` and `donor2`
#'
#' @import tidyverse
#'
sim_generate_donor_hap <- function(num_snps){
  hap1 <- data.frame(V1 = sample(c(0,1), size=num_snps, replace=TRUE)) #simulate first donor chromosome
  hap2 <- abs(1-hap1) #switch bits to construct second donor
  donor_haps <- data.frame(cbind(hap1, hap2)) %>% `colnames<-`(c("donor1", "donor2"))
  return (donor_haps)
}