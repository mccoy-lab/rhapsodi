#' This function is used to run the generative model to simulate sparse gamete data with transmission distortion caused by gamete killing rhapsodi
#' 
#
#' This function runs the generative model to simulate input gamete data for rhapsodi. This function samples more gametes than necessary and then to simulate transmission distortion removes a fraction of gametes containing a specific haplotype at a specific SNP. This SNP/haplotype combination can be user-defined or randomly generated.
#' Remaining gametes are sampled such that the desired number remains. In addition to returning gamete data, the function also returns the fully known generated gamete data, this function returns true donor haplotypes and the identity of the SNP and haplotype that are used for transmission distortion simulation
#'
#'
#'
#' @param num_gametes an integer, the number of gametes, or the number of columns for the sparse gamete data you want generated
#' @param num_snps an integer, the number of SNPs, or the number of rows for the sparse gamete data you want generated. Note: not all of these will be heterozygous due to the coverage and therefore this number won't necessarily equal the number of SNPs following filtering at the end of the generation 
#' @param p_kill a float, the probability that a gamete containing the SNP subject to TD will be removed from the dataset through simulated gamete killing
#' @param killer_snp an integer indicating the specific SNP which will be subject to TD. Randomly selected if not provided by the user. 
#' @param killer_haplotype an integer, 0 or 1, indicating which haplotype will be subject to transmission distortion. Randomly selected if not provided by the user.
#' @param recomb_lambda a numeric, the average recombination rate, used as the mean for the poisson distribution 
#' 
#' @return a list of the following: sim_gam_filtered, the simulated gametes; killer_snp, the snp used to simulated gamete killing; killer_haplotype, the haplotype used to simulated gamete killing; donor_haps, the true donor haplotypes 
#'

sim_gamete_killing <- function(num_gametes = 500, num_snps = 5000, p_kill=0.5, killer_snp, killer_haplotype,
                                  recomb_lambda = 1){
  
  # Determine how many gametes must be simulated 
  num_gam_to_sim <- ceiling(1.5 * (num_gametes / (1 - p_kill/2)))
  
  # From sim_run_generative_model - get number of crossovers; generate donor haps
  n_crossovers <- sim_find_num_crossovers(num_gam_to_sim, recomb_lambda)
  donor_haps <- sim_generate_donor_hap(num_snps)
  
  sim_gam <- lapply(1:num_gam_to_sim, function(x) sim_generate_gam(donor_haps, n_crossovers[x]))
  
  # Randomly select a SNP for TD if not provided by user
  if (is.null(killer_snp)){
    killer_snp <- sample(1:length(sim_gam[[1]][[2]]), 1)
  }
  
  # Randomly select a haplotype to be killed if not provided by user
  if (is.null(killer_haplotype)){
    killer_haplotype <- sample(1:2,1)
  }
  
  # Find which gametes have the killer haplotype at the correct snp
  gam_with_killer <- unlist(lapply(1:num_gam_to_sim, function(x) sim_gam[[x]][[3]][killer_snp] == killer_haplotype))
  
  # Get indexes of gametes with killer haplotype and sample which ones to remove
  gam_with_killer_idx <- which(gam_with_killer == TRUE)
  num_killed = as.integer(p_kill * length(gam_with_killer_idx))
  if ((num_killed == 0) & (p_kill > 0)){
    num_killed = 1
  }  else if (num_killed == length(gam_with_killer_idx)){
    num_killed = num_killed - 1
  }
  if (length(gam_with_killer_idx) >= 1){
    TD_removed_gams <- sample(gam_with_killer_idx, num_killed)
  }
  else {
    TD_removed_gams <- c()
  }
  
  
  # If at least one gamete killed, remove killed gametes from dataset
  if ((p_kill * length(gam_with_killer_idx)) >= 1){
    sim_gam_filtered <- sim_gam[-c(TD_removed_gams)]
  } else{
    sim_gam_filtered <- sim_gam
  }
  
  # Sample dataset 
  num_to_keep <- sample(length(sim_gam_filtered), num_gametes)
  sim_gam_filtered <- sim_gam_filtered[c(sample(length(sim_gam_filtered), num_gametes))]
  
  return(list(sim_gam_filtered, killer_snp, killer_haplotype, donor_haps))
  
}