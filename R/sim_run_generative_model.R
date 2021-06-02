#'
#'
#'
sim_run_generative_model <- function(){
  #get argument inputs
  #input_cov a bool
  ##coverage a float
  #input_mgr a bool
  ##missing_genotype_rate a float
  #random_seed an int
  #num_gametes an int
  #recomb_lambda a float
  #num_snps an int
  #add_de_novo_mut a bool
  ##de_novo_lambda an int
  ##de_novo_alpha a float
  ##de_novo_beta a float
  #add_seq_error a bool
  ##seqError_add a float
  
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
  tci_dt <- data.table(gam=str_split(names(unlist_ci), "_", simplify=TRUE)[,1], start =(unlist_ci-1), end=(unlist_ci))
  
  gam_mat <- sapply(sim_gam, "[[", 2)
  if (missing_genotype_rate > 0.5){
    gam_mat_with_na <- sim_add_to_na_flatten(gam_mat, num_nas, num_gametes, num_snps)
  } else if (missing_genotype_rate <= 0.5){
    gam_mat_with_na <- sim_add_na_flatten(gam_mat, num_nas, num_gametes, num_snps)
  }
  
  if (add_de_novo_mut){
    dnm_out <- sim_add_de_novo_mut(de_novo_lambda, de_novo_alppha, de_novo_beta, num_snps, num_gametes, gam_haps, gam_mat, gam_mat_with_na, donor_haps, unlist_ci, missing_genotype_rate)
    #do I have to do anything fancy to get objects from this, or is it ok since it's a named list?
    
    num_snps <- dnm_out$num_snps
    
    gam_mat_with_na <- dnm_out$gam_mat_with_na
    
    gam_mat <- dnm_out$gam_mat
    
    unlist_ci <- dnm_out$unlist_ci
    tci_dt <- data.table(gam=str_split(names(unlist_ci), "_", simplify=TRUE)[,1], start =(unlist_ci-1), end=(unlist_ci))
  }
  if (add_seq_error){
    gam_mat_with_na <- sim_add_seq_error(num_snps, num_gametes, seqError_add, gam_mat_with_na)
  }
  gam_na_df <- data.frame(pseudo_pos = 1:nrow(gam_mat_with_na), gam_mat_with_na) %>% `colnames<-`(c("positions", paste0("gam", 1:num_gametes, "_")))
  gam_full_df <- data.frame(pseudo_pos = 1:nrow(gam_mat), gam_mat) %>% `colnames<-`(c("positions", paste0("gam", 1:num_gametes, "_")))
  
  #filtering?
}