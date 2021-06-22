#' This function filters the generated gamete data so that the sparse input data only reflects heterozygous SNPs
#' 
#' This function is used to filter the generated data such that the sparse gamete data only reflects heterozygous SNPs (one 0 and one 1)
#' Then we only keep these SNP indices for the full gamete data and the diploid donor haplotypes
#' Further, if de novo mutations were simulated, we filter the list of DNM locations to reflect those that weren't filtered out. 
#' However, we don't return the filtered dnm list, we just output it as a message
#' 
#' @param gam_na_df the generated, sparsified gamete data with SNP indices in the first column 
#' @param gam_full_df the generated full gamete data with SNP indices in the first column
#' @param donor_haps the generated diploid donor phased haplotypes
#' @param new_rows a vector of the integer SNP indices for where de novo mutation were added if DNMs were simulated, this should be NULL otherwise
#' @param add_de_novo_mut a bool, if TRUE de novo mutations were simulated in the generative model and then `new_rows` is also filtered
#'
#' @return out a named list with the filtered `gam_na_df`, filtered `gam_full_df`, filtered `donor_haps`, new number of snps after filtering `num_snps`
#' 
sim_filter_generated_data <- function(gam_na_df, gam_full_df, donor_haps, new_rows, add_de_novo_mut){
  out <- list()
  keep_bool <- unname((rowSums(gam_na_df[, 2:ncol(gam_na_df)] == 0, na.rm=TRUE) > 0) & (rowSums(gam_na_df[, 2:ncol(gam_na_df)] == 1, na.rm=TRUE) > 0))
  
  out$gam_na_df <- gam_na_df[keep_bool,]
  out$gam_full_df <- gam_full_df[keep_bool, ]
  out$donor_haps <- donor_haps[keep_bool, ]
  
  out$num_snps <- sum(keep_bool)
  message(paste0("new number of snps: ", out$num_snps))
  if (add_de_novo_mut){
    `%notin%` <- Negate(`%in%`)
    for (i in 1:length(new_rows)){
      message(paste0("dnm ", i, " is filtered out: ", new_rows[i] %notin% out$gam_na_df[,1]))
    }
  }
  
  return(out)
}