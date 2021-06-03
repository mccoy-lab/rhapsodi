#'
#'
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