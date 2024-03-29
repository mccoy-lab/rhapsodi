#' This function runs the generative model to simulate input sparse gamete data for rhapsodi
#' 
#' This function runs the generative model to simulate input sparse gamete data for rhapsodi. In addition,
#' to returning the sparse gamete data, the function also returns the fully known generated gamete data,
#' the diploid donor phased haplotypes, and the true recombination break points for each gamete.
#' The following variables of the simulation can all be controlled:
#' the number of gametes, the number of SNPs, the sequencing coverage (or missing genotype rate), 
#' the average recombination rate, whether to simulate sequencing error, the sequencing error rate to use,
#' whether to add de novo mutations, values parameterizing how many de novo mutations there are and how many gametes are affected by the de novo mutations,
#' and the random seed for reproducibility
#'
#' @param num_gametes an integer, the number of gametes, or the number of columns for the sparse gamete data you want generated
#' @param num_snps an integer, the number of SNPs, or the number of rows for the sparse gamete data you want generated. Note: not all of these will be heterozygous due to the coverage and therefore this number won't necessarily equal the number of SNPs following filtering at the end of the generation 
#' @param coverage a numeric, input if input_cov is TRUE, suggested NULL otherwise
#' @param recomb_lambda a numeric, the average rate of recombination expected for the simulation
#' @param random_seed an integer, the random seed which will be set for the simulation, default=42
#' @param input_cov a logical, TRUE if coverage (i.e. like 0.01 (x)) will be input rather than missing genotype rate 
#' @param input_mgr a logical, TRUE if missing genotype rate (i.e. like 80 (%) or 0.8) will be inpupt rather than coverage, default = FALSE
#' @param missing_genotype_rate a numeric, input if input_mgr is TRUE and input_COV is FALSE, suggested NULL otherwise, default=NULL
#' @param add_seq_error a logical, TRUE if you want to add sequencing error to the generated data, default=TRUE
#' @param seqError_add a numeric, the sequencing error rate if adding sequencing error to the generated data, default=0.005
#' @param add_de_novo_mut a logical, TRUE if you want to add de novo mutations to the generated data, default=FALSE
#' @param de_novo_lambda an integer, default=5, parameterizes a poisson distribution to find the number of de novo mutations (DNM) total
#' @param de_novo_alpha a numeric, default=7.5, shape parameter for a gamma distribution to find the number of gametes affected per DNM
#' @param de_novo_beta a numeric, default=10, scale parameter for a gamma distribution to find the number of gametes affected per DNM
#'
#' @return generated_data a named list returning the generated input and full truth data, specifically `gam_na` for the sparse rhapsodi input, `gam_full` for the fully known gamete data input equivalent, `recomb_spots` for the true recombination spots for each gamete, and `donor_haps` for the diploid donor phased haplotypes
#'
#' @importFrom magrittr %>%
#' @import data.table
#'
#' @export
sim_run_generative_model <- function(num_gametes, num_snps, coverage, 
                                     recomb_lambda, random_seed=42, 
                                     input_cov=TRUE, input_mgr=FALSE, missing_genotype_rate=NULL,
                                     add_seq_error=TRUE, seqError_add=0.005,
                                     add_de_novo_mut=FALSE, de_novo_lambda=5, de_novo_alpha=7.5, de_novo_beta=10){
  ##is it necessary to verify argument input types?
 
  generated_data <- list()
  
  if (input_cov){
    missing_genotype_rate <- sim_find_mgr_from_cov(coverage)
  } else if (input_mgr){
    coverage <- sim_find_cov_from_mgr(missing_genotype_rate)
  }
  
  num_nas <- sim_find_num_nas(num_gametes, num_snps, missing_genotype_rate)
  
  set.seed(random_seed)
  n_crossovers <- sim_find_num_crossovers(num_gametes, recomb_lambda)
  
  donor_haps <- sim_generate_donor_hap(num_snps)
  
  sim_gam <- lapply(1:num_gametes, function(x) sim_generate_gam(donor_haps, n_crossovers[x]))
  
  gam_haps <- sapply(sim_gam, "[[", 3)
  
  crossover_indices <- sapply(sim_gam, "[[", 1)
  names(crossover_indices) <- paste0(rep("gam", num_gametes), 1:num_gametes, "_")
  unlist_ci <- unlist(crossover_indices, use.names=TRUE)
  tci_dt <- data.table::data.table(gam=sapply(strsplit(names(unlist_ci), "_"), `[`, 1), start =(unlist_ci-1), end=(unlist_ci))
  
  gam_mat <- sapply(sim_gam, "[[", 2)
  if (missing_genotype_rate > 0.5){
    gam_mat_with_na <- sim_add_to_na_flatten(gam_mat, num_nas, num_gametes, num_snps)
  } else if (missing_genotype_rate <= 0.5){
    gam_mat_with_na <- sim_add_na_flatten(gam_mat, num_nas, num_gametes, num_snps)
  }
  
  if (add_de_novo_mut){
    dnm_out <- sim_add_de_novo_mut(de_novo_lambda, de_novo_alpha, de_novo_beta, num_snps, num_gametes, gam_haps, gam_mat, gam_mat_with_na, donor_haps, unlist_ci, missing_genotype_rate)
    num_snps <- dnm_out$num_snps
    
    donor_haps <- dnm_out$donor_haps
    
    gam_mat_with_na <- dnm_out$gam_mat_with_na
    
    gam_mat <- dnm_out$gam_mat
    
    unlist_ci <- dnm_out$unlist_ci
    tci_dt <- data.table::data.table(gam=sapply(strsplit(names(unlist_ci), "_"), `[`, 1), start =(unlist_ci-1), end=(unlist_ci))
    
    new_dnm_rows <- dnm_out$new_rows
  } else {new_dnm_rows <- c() }
  if (add_seq_error){
    gam_mat_with_na <- sim_add_seq_error(num_snps, num_gametes, seqError_add, gam_mat_with_na)
  }
  gam_na_df <- data.frame(pseudo_pos = 1:nrow(gam_mat_with_na), gam_mat_with_na) %>% `colnames<-`(c("positions", paste0("gam", 1:num_gametes, "_")))
  gam_full_df <- data.frame(pseudo_pos = 1:nrow(gam_mat), gam_mat) %>% `colnames<-`(c("positions", paste0("gam", 1:num_gametes, "_")))
  
  #Filtering
  filtered_out <- sim_filter_generated_data(gam_na_df, gam_full_df, donor_haps, new_dnm_rows, add_de_novo_mut)
  gam_na_df <- filtered_out$gam_na_df
  gam_full_df <- filtered_out$gam_full_df
  donor_haps <- filtered_out$donor_haps
  num_snps <- filtered_out$num_snps
  
  generated_data$recomb_spots <- tci_dt
  generated_data$gam_full <- gam_full_df
  generated_data$gam_na <- gam_na_df
  generated_data$donor_haps <- donor_haps
  
  return(generated_data)
}