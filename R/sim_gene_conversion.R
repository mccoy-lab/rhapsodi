#' This function is used to run the generative model to simulate sparse gamete data with transmission distortion caused by gene conversion for rhapsodi
#' 
#
#' This function runs the generative model to simulate input gamete data for rhapsodi. To simulate transmission distortion, a SNP is chosen as the start of a gene conversion event, the length of which is sampled from the Poisson distribution. Haplotypes are bit inverted for a user-defined fraction of gametes containing this SNP/haplotype combination. 
#' In addition to returning gamete data, the function also returns the fully known generated gamete data, this function returns true donor haplotypes and the identity of the SNP and haplotype that are used for transmission distortion simulation
#'
#'
#'
#' @param num_gametes an integer, the number of gametes, or the number of columns for the sparse gamete data you want generated
#' @param num_snps an integer, the number of SNPs, or the number of rows for the sparse gamete data you want generated. Note: not all of these will be heterozygous due to the coverage and therefore this number won't necessarily equal the number of SNPs following filtering at the end of the generation 
#' @param p_convert a float, the probability that a gamete with the TD allele will undergo a gene conversion event
#' @param conversion_lambda a float, used as lambda in the Poisson distribution to determine the length of the gene conversion event
#' @param converted_snp an integer indicating the specific SNP which will be subject so TD. Randomly selected if not provided by the user. 
#' @param converted_haplotype an integer, 0 or 1, indicating which haplotype will be subject to transmission distortion. Randomly selected if not provided by the user.
#' @param recomb_lambda a numeric, the average recombination rate, used as the mean for the poisson distribution
#' 
#' @return a list of the following: sim_gam_filtered, the simulated gametes; converted_snp, the snp used to simulated gene conversion; converted_haplotype, the haplotype used to simulated gene conversion; donor_haps, the true donor haplotypes 
#'

sim_gene_conversion <- function(num_gametes = 500, num_snps = 5000, p_convert=0.5, conversion_lambda = 4, 
                                   converted_snp, converted_haplotype, recomb_lambda = 1){
  # From sim_run_generative_model - get number of crossovers; generate donor haps
  n_crossovers <- sim_find_num_crossovers(num_gametes, recomb_lambda)
  donor_haps <- sim_generate_donor_hap(num_snps)
  
  sim_gam <- lapply(1:num_gametes, function(x) sim_generate_gam(donor_haps, n_crossovers[x]))
  
  # Get length of gene conversion 
  conversion_length <- rpois(1, conversion_lambda)
  # Randomly select a SNP for TD if not provided by user
  if (is.null(converted_snp)){
    converted_snp <- sample(1:(length(sim_gam[[1]][[2]]) - conversion_length), 1)
  }
  
  # Randomly select a haplotype to be converted if not provided by user
  if (is.null(converted_haplotype)){
    converted_haplotype <- sample(1:2,1)
  }
  # Find which gametes have the converted haplotype at the correct snp
  gam_with_converted <- unlist(lapply(1:num_gametes, function(x) sim_gam[[x]][[3]][converted_snp] == converted_haplotype))
  # Get indexes of gametes with converted haplotype and sample which ones to convert
  gam_with_conveted_idx <- which(gam_with_converted == TRUE)
  num_converted = as.integer(p_convert * length(gam_with_conveted_idx))
  if (length(gam_with_conveted_idx) >= 1){
    TD_converted_gams <- sample(gam_with_conveted_idx, num_converted)
  }
  else {
    TD_converted_gams <- c()
  }
  
  for (i in TD_converted_gams){
    for (j in c(converted_snp:(converted_snp+conversion_length))){
      sim_gam[[i]][[2]][[j]] <- abs(1-sim_gam[[i]][[2]][[j]])
      sim_gam[[i]][[3]][[j]] <- 1 + abs(1 - (converted_haplotype - 1))
    }
  }
  
  return(list(sim_gam, converted_snp, converted_haplotype, donor_haps))
}